/* hmmsearch: search profile HMM(s) against a sequence database.
 *
 * To do:
 *  - in MPI mode, add a check to make sure ncpus >= 2. If 1, then we
 *    only have a master, no workers. See Infernal commit r3972 on the
 *    same point; and same note in hmmscan.c's to do list.
 */
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#ifdef HMMER_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif 

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif 

#include "hmmer.h"










// SPLICING STUFF BEGINS!!!


static bool ALEX_DEBUG = true;


static char[21] AMINO_CHARS = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'];
static char[5]  DNA_CHARS   = ['A','C','G','T','-'];
static char[5]  RNA_CHARS   = ['A','C','G','U','-'];


// How many amino acids are we willing to extend
// to bridge two hits?
static int MAX_AMINO_EXT = 8;





/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DumpSeqData
 *
 *  Inputs:
 *
 *  Output:
 *
 */
void DumpSeqData
(ESL_DSQ * Aminos, int num_aminos, ESL_DSQ * Nucls, int num_nucls)
{

  printf("\nNum Aminos: %d\n",num_aminos);
  for (int i=0; i<num_aminos; i++) printf(" %c ",AMINO_CHARS[Aminos[i]]);
  printf("\nNum Nucls : %d\n",num_nucls);
  for (int i=0; i<num_aminos; i++) printf("%c",DNA_CHARS[Nucls[i]]);
  printf("\n");

}





/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DomainsAreSpliceComaptible
 *
 *  Inputs:
 *
 *  Output:
 *
 */
bool DomainsAreSpliceCompatible
(P7_ALIDISPLAY * upstream, P7_ALIDISPLAY * downstream)
{

  // Start by checking if we either have amino acid
  // overlap, or are close enough to consider extending
  int amino_start_1 = upstream->hmmfrom;
  int amino_end_1   = upstream->hmmto;

  int amino_start_2 = downstream->hmmfrom;
  int amino_end_2   = downstream->hmmto;

  // If the upstream ain't upstream, then obviously we can't treat
  // these as splice-compatible!
  if (!(amino_start_1 < amino_start_2 && amino_end_1 < amino_end_2))
    return false;

  // Do we have overlap OR sufficient proximity to consider
  // extending?
  if (!(amino_end_1 + MAX_AMINO_EXT >= amino_start_2))
    return false;


  // Fantastic!  The amino acid coordinates support splice
  // compatibility!  Now it's just time to confirm that
  // the nucleotides also look good.

  int  nucl_start_1 = upstream->sqfrom;
  int  nucl_end_1   = upstream->sqto;
  
  int revcomp1 = 0;
  if (nucl_start_1 > nucl_end_1)
    revcomp1 = 1;


  int nucl_start_2 = downstream->sqfrom;
  int nucl_end_2   = downstream->sqto;
 
  int revcomp2 = 0;
  if (nucl_start_2 > nucl_end_2)
    revcomp2 = 1;


  if (revcomp1 != revcomp2) 
    return false;


  // We want to make sure that these aren't unrealistically
  // close together on the genome...
  if (revcomp1) {

    if (nucl_start_2 + (3 * MAX_AMINO_EXT) >= nucl_end_1)
      return false;

  } else {

    if (nucl_start_2 - (3 * MAX_AMINO_EXT) <= nucl_end_1)
      return false;

  }


  // Looks like we've got a viable upstream / downstream pair!
  return true;


}




/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GatherViableDownstreamHits
 *
 *  Inputs:
 *
 *  Output:
 *
 */
int GatherViableDownstreamHits
(
  P7_TOPHITS * TopHits,
  uint64_t     upstream_hit_id,
  uint64_t  ** ValidCompsByHitID,
  int       ** ValidCompUpstreamDoms,
  int       ** ValidCompDownstreamDoms
)
{

  P7_HIT * UpstreamHit = TopHits->hit[upstream_hit_id];
  int num_upstream_domains = UpstreamHit->ndom;


  // We'll setup a temporary array to hold the indices
  // of hits that could serve as downstream exons for our
  // candidate (upstream) hit
  int vc_capacity  = 8;
  int vc_occupancy = 0;
  uint64_t * ViableComps        = malloc(vc_capacity * sizeof(uint64_t));
  int * ViableUpstreamDomains   = malloc(vc_capacity * sizeof(int));
  int * ViableDownstreamDomains = malloc(vc_capacity * sizeof(int));


  // For each hit, gather all of the indices of other hits that
  // could potentially be downstream exons.
  uint64_t num_hits = TopHits->N;
  for (uint64_t downstream_hit_id; downstream_hit_id < num_hits; downstream_hit_id++) {

    // No self-splicing, gentlemen!
    if (downstream_hit_id == upstream_hit_id)
      continue;

    // We only consider splicing when the two hits come 
    // from the same sequence.  Further, because of how
    // we determine splice viability (specifically, by
    // checking compatibility of amino acid start/end
    // coordinates), this wouldn't be a trivial check to
    // remove to allow for splicing across sequences.
    if (strcmp(UpstreamHit->name,TopHits->hit[downstream_hit_id]->name))
      continue;

    // Now that we know these are valid to check, check 'em!
    P7_HIT * DownstreamHit = TopHits->hit[downstream_hit_id];
    int num_downstream_domains = DownstreamHit->ndom;

    for (int upstream_domain_id = 0; upstream_domain_id < num_upstream_domains; upstream_domain_id++) {

      P7_DOMAIN * UpstreamDomain = &UpstreamHit->dcl[upstream_domain_id];

      for (int downstream_domain_id = 0; downstream_domain_id < num_downstream_domains; downstream_domain_id++) {

        P7_DOMAIN * DownstreamDomain = &DownstreamHit->dcl[downstream_domain_id];

        if (DomainsAreSpliceCompatible(UpstreamDomain,DownstreamDomain)) {

          // Do we need to resize before entering new info.?
          if (vc_occupancy == vc_capacity) {

            uint64_t * TmpComps  = malloc(vc_occupancy * sizeof(uint64_t));
            int * TmpViableUps   = malloc(vc_occupancy * sizeof(int));
            int * TmpViableDowns = malloc(vc_occupancy * sizeof(int));

            for (int i=0; i<vc_occupancy; i++) {
              TmpComps[i]       = ViableComps[i];
              TmpViableUps[i]   = ViableUpstreamDomains[i];
              TmpViableDowns[i] = ViableDownstreamDomains[i];
            }

            free(ViableComps);
            free(ViableUpstreamDomains);
            free(ViableDownstreamDomains);

            vc_capacity *= 2;
            ViableComps             = malloc(vc_capacity * sizeof(uint64_t));
            ViableUpstreamDomains   = malloc(vc_capacity * sizeof(int));
            ViableDownstreamDomains = malloc(vc_capacity * sizeof(int));

            for (int i=0; i<vc_occupancy; i++) {
              ViableComps[i]             = TmpComps[i];
              ViableUpstreamDomains[i]   = TmpViableUps[i];
              ViableDownstreamDomains[i] = TmpViableDowns[i];
            }

            free(TmpComps);
            free(TmpViableUps);
            free(TmpViableDowns);

          }

          // Record that splice compatibility!
          ViableComps[vc_occupancy]             = downstream_hit_id;
          ViableUpstreamDomains[vc_occupancy]   = upstream_domain_id;
          ViableDownstreamDomains[vc_occupancy] = downstream_domain_id;

          vc_occupancy++;

        }

      }

    }


  }


  // Did we find any viable downstream exons for this hit?
  // If so, copy 'em over!
  if (vc_occupancy > 0) {

    ValidCompsByHitID[upstream_hit_id]       = malloc(vdh_occupancy * sizeof(uint64_t));
    ValidCompUpstreamDoms[upstream_hit_id]   = malloc(vdh_occupancy * sizeof(int));
    ValidCompDownstreamDoms[upstream_hit_id] = malloc(vdh_occupancy * sizeof(int));

    for (int i=0; i<vc_occupancy; i++) {
      ValidCompsByHitID[upstream_hit_id][i]       = ViableComps[i];
      ValidCompUpstreamDoms[upstream_hit_id][i]   = ViableUpstreamDomains[i];
      ValidCompDownstreamDoms[upstream_hit_id][i] = ViableDownstreamDomains[i];
    }

  }


  free(ViableComps);
  free(ViableUpstreamDomains);
  free(ViableDownstreamDomains);


  return vc_occupancy;

}






/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: AttemptSpliceEdge
 *
 *  Inputs:  
 *
 *  Output:  (Eventually) splice graphs built around the original hits.
 *
 */
void AttemptSpliceEdge
(
  P7_HIT      * UpstreamHit,
  P7_DOMAIN   * UpstreamDomain,
  P7_HIT      * DownstreamHit,
  P7_DOMAIN   * DownstreamDomain,
  P7_OPROFILE * om, 
  ESL_GENCODE * gcode
)
{

  // We'll want to have the emission score array and transition
  // costs on-hand as we're scoring the net score change of splicing
  float * FwdEmissionScores  = malloc((om->abc->Kp * (om->M + 1)) * sizeof(float));
  float * FwdTransitionProbs = malloc((7           * (om->M + 2)) * sizeof(float));


  p7_oprofile_GetFwdEmissionScoreArray(om,FwdEmissionScores);


  int mm_transit = 0;
  int mi_transit =      om->M + 2 ;
  int md_transit = 2 * (om->M + 2);
  int ii_transit = 3 * (om->M + 2);
  int im_transit = 4 * (om->M + 2);
  int dd_transit = 5 * (om->M + 2);
  int dm_transit = 6 * (om->M + 2);

  p7_oprofile_GetFwdTransitionArray(om,1,&FwdTransitionProbs[mm_transit]);
  p7_oprofile_GetFwdTransitionArray(om,5,&FwdTransitionProbs[mi_transit]);
  p7_oprofile_GetFwdTransitionArray(om,4,&FwdTransitionProbs[md_transit]);
  p7_oprofile_GetFwdTransitionArray(om,6,&FwdTransitionProbs[ii_transit]);
  p7_oprofile_GetFwdTransitionArray(om,2,&FwdTransitionProbs[im_transit]);
  p7_oprofile_GetFwdTransitionArray(om,7,&FwdTransitionProbs[dd_transit]);
  p7_oprofile_GetFwdTransitionArray(om,3,&FwdTransitionProbs[dm_transit]);


  // Convert to log scores (natural)
  for (i=0; i<7*(om->M+2); i++)
    FwdTransitionProbs[i] = log(FwdTransitionProbs[i]);


  // We'll need a background model for... transition probabilities...?
  P7_BG * BackgroundModel = p7_bg_Create(om->abc);
  float null_transit_prob = bg->p1;



  // CHECKPOINT
  //
  // The above work sets up up to be able to interact with the model that
  // provided the hits we're going to attempt to splice into one another.
  //
  // Now we can access those hits themselves (as alignments) and do the work
  // of splicing them!



  // Use the alignment displays as we consider what these
  // domain alignments actually look like
  P7_ALIDISPLAY * UpstreamDisplay   = UpstreamDomain->ad;
  P7_ALIDISPLAY * DownstreamDisplay = DownstreamDomain->ad;


  // NOTE that these are (essentially) numerical codes for the aminos,
  // so we'll need to use AMINO_CHARS anytime we want the actual letter
  ESL_DSQ * UpstreamAminos;
  ESL_DSQ * DownstreamAminos;

  esl_abc_CreateDsq(om->abc,UpstreamDisplay->aseq  ,&UpstreamAminos  );
  esl_abc_CreateDsq(om->abc,DownstreamDisplay->aseq,&DownstreamAminos);



  // Get the length of the alignments (including indels, of course!)
  int upstream_amino_ali_len   = strlen(UpstreamDisplay->aseq);
  int downstream_amino_ali_len = strlen(DownstreamDisplay->aseq);

  int upstream_nucl_ali_len   = strlen(UpstreamDisplay->ntseq);
  int downstream_nucl_ali_len = strlen(DownstreamDisplay->ntseq);

  // We'll also need the model consensus sequence in order to check whether
  // or not an emission is an insertion (if the model consensus character is '.').
  ESL_DSQ * UpstreamConsensus;
  ESL_DSQ * DownstreamConsensus;

  esl_abc_CreateDsq(om->abc,UpstreamDisplay->model  ,&UpstreamConsensus  );
  esl_abc_CreateDsq(om->abc,DownstreamDisplay->model,&DownstreamConsensus);



  // Because we want to be able to accommodate nucleotide sequences that
  // are either DNA or RNA, we'll check the 
  int nucl_type = 0;


  for (int i=0; i<upstream_nucl_ali_len; i++) {
    
    if (UpstreamDisplay->ntseq[i] == 'T') {
      nucl_type = eslDNA;
      break;
    }

    if (UpstreamDisplay->ntseq[i] == 'U') {
      nucl_type = eslRNA;
      break;
    }

  }


  if (nucl_type == 0) {

    for (int i=0; i<downstream_nucl_ali_len; i++) {
    
      if (DownstreamDisplay->ntseq[i] == 'T') {
        nucl_type = eslDNA;
        break;
      }

      if (DownstreamDisplay->ntseq[i] == 'U') {
        nucl_type = eslRNA;
        break;
      }

    }

    if (nucl_type == 0) 
      nucl_type = eslDNA;

  }


  ESL_ALPHABET * nucl_alphabet = esl_alphabet_Create(nucl_type);



  // FINALLY we can grab the nucleotide sequences!
  ESL_DSQ * UpstreamNucls;
  ESL_DSQ * DownstreamNucls;

  esl_abc_CreateDsq(nucl_alphabet,UpstreamDisplay->ntseq  ,&UpstreamNucls  );
  esl_abc_CreateDsq(nucl_alphabet,DownstreamDisplay->ntseq,&DownstreamNucls);



  // If we're debugging, we might want to dump the sequence
  // data we've acquired to the ol' terminal...
  if (ALEX_DEBUG) {
    printf("\n\nUPSTREAM\n");
    DumpSeqData(UpstreamAminos  , upstream_amino_ali_len  , UpstreamNucls  , upstream_nucl_ali_len  );
    printf("\nDOWNSTREAM\n");
    DumpSeqData(DownstreamAminos, downstream_amino_ali_len, DownstreamNucls, downstream_nucl_ali_len);
    printf("\n\n");
  }



  // I think that's all we need to get to work!
  // So... GET TO WORK!


  // How many aminos do we want to extend each hit?
  // Our target is to have 2 aminos of overlap, just for
  // a bit o' fun
  int min_overlap_aminos  = 2;
  int init_overlap_aminos = (DownstreamDisplay->hmmfrom - UpstreamDisplay->hmmto) + 1;

  // If the two hits don't overlap (*require* extension)
  // then init_overlap_aminos will be negative.
  // If the two hits do    overlap, then we might dip below
  // zero (no extension, even to get to our min_overlap_aminos)
  int aminos_to_extend = min_overlap_aminos - init_overlap_aminos;

  if (aminos_to_extend > 0) {
  
    ExtendUpstreamHit();
    ExtendDownstreamHit();

  }



}






/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GenSpliceGraphs
 *
 *  Inputs:  
 *
 *  Output:  (Eventually) splice graphs built around the original hits.
 *
 */
void GenSpliceGraphs
(P7_TOPHITS * TopHits, P7_OPROFILE * om, ESL_GENCODE * gcode)
{


  uint64_t num_hits = TopHits->N;


  // For each hit, which other hits function as downstream
  // exons in the splice graph?
  //
  // The 'ValidCompsByHitID' gives us pairs of indices in TopHits,
  // and our indices *into* VCBHID correspond to the entries in the
  // ValidComp(Up/Down)streamDoms data.
  //
  uint64_t ** ValidCompsByHitID       = malloc(num_hits * sizeof(uint64_t *));
  int      ** ValidCompUpstreamDoms   = malloc(num_hits * sizeof(int *));
  int      ** ValidCompDownstreamDoms = malloc(num_hits * sizeof(int *));
  int       * NumCompsByHitID         = malloc(num_hits * sizeof(int));

  for (uint64_t hit_id = 0; hit_id < num_hits; hit_id++) {
    
    ValidCompsByHitID[hit_id]       = NULL;
    ValidCompUpstreamDoms[hit_id]   = NULL;
    ValidCompDownstreamDoms[hit_id] = NULL;

    NumDownstreamHits[hit_id] = 
      GatherViableDownstreamHits(TopHits,hit_id,ValidCompsByHitID,ValidCompUpstreamDoms,ValidCompDownstreamDoms);

  }


  // Now we can run through all of our paired domains and actually
  // splice 'em up (or at least try our best to)!
  for (uint64_t hit_id = 0; hit_id < num_hits; hit_id++) {

    P7_HIT * UpstreamHit = TopHits->hit[hit_id];

    for (int i=0; i<NumDownstreamHits[hit_id]; i++) {

      P7_DOMAIN * UpstreamDomain = &UpstreamHit->dcl[ValidCompUpstreamDoms[hit_id][i]];

      P7_HIT    * DownstreamHit    = TopHits->hit[ValidCompsByHitID[hit_id][i]];
      P7_DOMAIN * DownstreamDomain = &DownstreamHit->dcl[ValidCompDownstreamDoms[hit_id][i]];

      AttemptSpliceEdge(UpstreamHit,UpstreamDomain,DownstreamHit,DownstreamDomain,om,gcode);

    }

  }


}











typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif 
  P7_BG            *bg;	         /* null model                              */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
  P7_OPROFILE      *om;          /* optimized query profile                 */
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#if defined (HMMER_THREADS) && defined (HMMER_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <f>",              2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout", eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save table of hits and domains to file, in Pfam format <f>",   2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "report domains <= this E-value threshold in output",           4 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "report domains >= this score cutoff in output",                4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "consider domains <= this E-value threshold as significant",    5 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "consider domains >= this score threshold as significant",      5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",      7 },
  { "--F1",         eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F2",         eslARG_REAL,  "1e-3", NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F3",         eslARG_REAL,  "1e-5", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--nobias",     eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL, "--max",          "turn off composition bias filter",                             7 },

/* Other options */
  { "--nonull2",    eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "turn off biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert target <seqfile> is in format <s>: no autodetection",  12 },

#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT, p7_NCPU,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,      "number of parallel CPU workers to use for multithreads",      12 },
#endif
#ifdef HMMER_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--mpi", NULL,            "arrest after start: for debugging MPI under gdb",             12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  MPIOPTS,         "run as an MPI parallel program",                              12 },
#endif

  /* Restrict search to subset of database - hidden because these flags are
   *   (a) currently for internal use
   *   (b) probably going to change
   * Doesn't work with MPI
   */
  { "--restrictdb_stkey", eslARG_STRING, "0",  NULL, NULL,    NULL,  NULL,  NULL,       "Search starts at the sequence with name <s> (not with MPI)",     99 },
  { "--restrictdb_n",eslARG_INT,        "-1",  NULL, NULL,    NULL,  NULL,  NULL,       "Search <j> target sequences (starting at --restrictdb_stkey)",   99 },
  { "--ssifile",    eslARG_STRING,       NULL, NULL, NULL,    NULL,  NULL,  NULL,       "restrictdb_x values require ssi file. Override default to <s>",  99 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <hmmfile> <seqdb>";
static char banner[] = "search profile(s) against a sequence database";

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *dbfile;            /* target sequence database file                   */
  char            *hmmfile;           /* query HMM file                                  */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */

  char             *firstseq_key;     /* name of the first sequence in the restricted db range */
  int              n_targetseq;       /* number of sequences in the restricted range */
};

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp, int n_targetseqs);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, int n_targetseqs);
static void pipeline_thread(void *arg);
#endif 

#ifdef HMMER_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif 


static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      if (puts("\nOptions controlling reporting thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      if (puts("\nOptions controlling inclusion (significance) thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      if (puts("\nOptions controlling model-specific thresholding:")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 2)     { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0 && strcmp(*ret_seqfile, "-") == 0) 
    { if (puts("Either <hmmfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere most common options are:")                                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  p7_banner(ofp, go->argv[0], banner);
  
  if (fprintf(ofp, "# query HMM file:                  %s\n", hmmfile)                                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(ofp, "# target sequence database:        %s\n", seqfile)                                                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")           && fprintf(ofp, "# output directed to file:         %s\n",             esl_opt_GetString(go, "-o"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-A")           && fprintf(ofp, "# MSA of all hits saved to file:   %s\n",             esl_opt_GetString(go, "-A"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")     && fprintf(ofp, "# per-seq hits tabular output:     %s\n",             esl_opt_GetString(go, "--tblout"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domtblout")  && fprintf(ofp, "# per-dom hits tabular output:     %s\n",             esl_opt_GetString(go, "--domtblout"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pfamtblout") && fprintf(ofp, "# pfam-style tabular hit output:   %s\n",             esl_opt_GetString(go, "--pfamtblout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--acc")        && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")      && fprintf(ofp, "# show alignments in output:       no\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")    && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")      && fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")           && fprintf(ofp, "# sequence reporting threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "-E"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")           && fprintf(ofp, "# sequence reporting threshold:    score >= %g\n",    esl_opt_GetReal(go, "-T"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domE")       && fprintf(ofp, "# domain reporting threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--domE"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domT")       && fprintf(ofp, "# domain reporting threshold:      score >= %g\n",    esl_opt_GetReal(go, "--domT"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")       && fprintf(ofp, "# sequence inclusion threshold:    E-value <= %g\n",  esl_opt_GetReal(go, "--incE"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")       && fprintf(ofp, "# sequence inclusion threshold:    score >= %g\n",    esl_opt_GetReal(go, "--incT"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomE")    && fprintf(ofp, "# domain inclusion threshold:      E-value <= %g\n",  esl_opt_GetReal(go, "--incdomE"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incdomT")    && fprintf(ofp, "# domain inclusion threshold:      score >= %g\n",    esl_opt_GetReal(go, "--incdomT"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--cut_ga")     && fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--cut_nc")     && fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--cut_tc")     && fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")        && fprintf(ofp, "# Max sensitivity mode:            on [all heuristic filters off]\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")         && fprintf(ofp, "# MSV filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F1"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")         && fprintf(ofp, "# Vit filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F2"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")         && fprintf(ofp, "# Fwd filter P threshold:       <= %g\n",             esl_opt_GetReal(go, "--F3"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")     && fprintf(ofp, "# biased composition HMM filter:   off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") && fprintf(ofp, "# Restrict db to start at seq key: %s\n",            esl_opt_GetString(go, "--restrictdb_stkey"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_n")     && fprintf(ofp, "# Restrict db to # target seqs:    %d\n",            esl_opt_GetInteger(go, "--restrictdb_n")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--ssifile")          && fprintf(ofp, "# Override ssi file to:            %s\n",            esl_opt_GetString(go, "--ssifile"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--nonull2")    && fprintf(ofp, "# null2 bias corrections:          off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")           && fprintf(ofp, "# sequence search space set to:    %.0f\n",           esl_opt_GetReal(go, "-Z"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domZ")       && fprintf(ofp, "# domain search space set to:      %.0f\n",           esl_opt_GetReal(go, "--domZ"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                               fprintf(ofp, "# random number seed set to:       %d\n",             esl_opt_GetInteger(go, "--seed"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--tformat")    && fprintf(ofp, "# targ <seqfile> format asserted:  %s\n",             esl_opt_GetString(go, "--tformat"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")        && fprintf(ofp, "# number of worker threads:        %d\n",             esl_opt_GetInteger(go, "--cpu"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
#endif
#ifdef HMMER_MPI
  if (esl_opt_IsUsed(go, "--mpi")        && fprintf(ofp, "# MPI:                             on\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#endif
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;	
  struct cfg_s     cfg;        
  int              status   = eslOK;

  impl_Init();                  /* processor specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.hmmfile    = NULL;
  cfg.dbfile     = NULL;
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.dbfile);    

/* is the range restricted? */
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") )
    if ((cfg.firstseq_key = esl_opt_GetString(go, "--restrictdb_stkey")) == NULL)  p7_Fail("Failure capturing --restrictdb_stkey\n");

  if (esl_opt_IsUsed(go, "--restrictdb_n") )
    cfg.n_targetseq = esl_opt_GetInteger(go, "--restrictdb_n");

  if ( cfg.n_targetseq != -1 && cfg.n_targetseq < 1 )
    p7_Fail("--restrictdb_n must be >= 1\n");


  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HMMER_MPI
  /* pause the execution of the programs execution until the user has a
   * chance to attach with a debugger and send a signal to resume execution
   * i.e. (gdb) signal SIGCONT
   */
  if (esl_opt_GetBoolean(go, "--stall")) pause();

  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HMMER_MPI*/
    {
      status = serial_master(go, &cfg);
    }

  esl_getopts_Destroy(go);

  return status;
}


/* serial_master()
 * The serial version of hmmsearch.
 * For each query HMM in <hmmfile> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().  We also use the
 * ESL_EXCEPTION and ERROR: mechanisms, but only because we know we're
 * using a fatal exception handler.
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-dom (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam tabular output (--pfamtblout)    */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  


  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
    if (esl_opt_IsUsed(go, "--ssifile"))
      esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
    else
      esl_sqfile_OpenSSI(dbfp, NULL);
  }



  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout")); }
  if (esl_opt_IsOn(go, "--pfamtblout")){ if ((pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)  esl_fatal("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
  ncpus = ESL_MIN( esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, (ptrdiff_t) sizeof(*info) * infocnt);

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      esl_sqfile_SetDigital(dbfp, abc); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks

      for (i = 0; i < infocnt; ++i)
	{
	  info[i].bg    = p7_bg_Create(abc);
#ifdef HMMER_THREADS
	  info[i].queue = queue;
#endif
	}

#ifdef HMMER_THREADS
      for (i = 0; i < ncpus * 2; ++i)
	{
	  block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
	  if (block == NULL) 	      esl_fatal("Failed to allocate sequence block");

 	  status = esl_workqueue_Init(queue, block);
	  if (status != eslOK)	      esl_fatal("Failed to add block to work queue");
	}
#endif
    }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */

      nquery++;
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)
      {
        if (! esl_sqfile_IsRewindable(dbfp))
          esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);

        if (! esl_opt_IsUsed(go, "--restrictdb_stkey") )
          esl_sqfile_Position(dbfp, 0); //only re-set current position to 0 if we're not planning to set it in a moment
      }

      if ( cfg->firstseq_key != NULL ) { //it's tempting to want to do this once and capture the offset position for future passes, but ncbi files make this non-trivial, so this keeps it general
        sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
        if (sstatus != eslOK)
          p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
      }

      if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
      if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                  /* <om> is now p7_LOCAL, multihit */

      for (i = 0; i < infocnt; ++i)
      {
        /* Create processing pipeline and hit list */
        info[i].th  = p7_tophits_Create();
        info[i].om  = p7_oprofile_Clone(om);
        info[i].pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
        status = p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
        if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

#ifdef HMMER_THREADS
        if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

#ifdef HMMER_THREADS
      if (ncpus > 0)  sstatus = thread_loop(threadObj, queue, dbfp, cfg->n_targetseq);
      else            sstatus = serial_loop(info, dbfp, cfg->n_targetseq);
#else
      sstatus = serial_loop(info, dbfp, cfg->n_targetseq);
#endif
      switch(sstatus)
      {
      case eslEFORMAT:
        esl_fatal("Parse failed (sequence file %s):\n%s\n",
            dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        break;
      case eslEOF:
        /* do nothing */
        break;
      default:
        esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i)
      {
        p7_tophits_Merge(info[0].th, info[i].th);
        p7_pipeline_Merge(info[0].pli, info[i].pli);

        p7_pipeline_Destroy(info[i].pli);
        p7_tophits_Destroy(info[i].th);
        p7_oprofile_Destroy(info[i].om);
      }















      // NORD - START
      GenSpliceGraphs(info->th,om,gcode);
      // NORD - END

















      /* Print the results.  */
      p7_tophits_SortBySortkey(info->th);
      p7_tophits_Threshold(info->th, info->pli);
      p7_tophits_Targets(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, info->th, info->pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
      if (domtblfp)  p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, info->th, info->pli, (nquery == 1));
      if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp, hmm->name, hmm->acc, info->th, info->pli);
  
      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, info->pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if (p7_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
	  {
	    esl_msa_SetName     (msa, hmm->name, -1);
	    esl_msa_SetAccession(msa, hmm->acc,  -1);
	    esl_msa_SetDesc     (msa, hmm->desc, -1);
	    esl_msa_FormatAuthor(msa, "hmmsearch (HMMER %s)", HMMER_VERSION);

	    if (textw > 0) esl_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msafile_Write(afp, msa, eslMSAFILE_PFAM);
	  
	    if (fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
	  } 
	else { if (fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
	  
	esl_msa_Destroy(msa);
      }

      p7_pipeline_Destroy(info->pli);
      p7_tophits_Destroy(info->th);
      p7_oprofile_Destroy(info->om);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus) {
  case eslEOD:       p7_Fail("read failed, HMM file %s may be truncated?", cfg->hmmfile);      break;
  case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);      break;
  case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);      break;
  case eslEOF:       /* do nothing. EOF is what we want. */                                    break;
  default:           p7_Fail("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
  }


  /* Terminate outputs... any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp, "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (pfamtblfp) p7_tophits_TabularTail(pfamtblfp,"hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit
   */
  for (i = 0; i < infocnt; ++i)
    p7_bg_Destroy(info[i].bg);

#ifdef HMMER_THREADS
  if (ncpus > 0)
    {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
	esl_sq_DestroyBlock(block);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
#endif

  free(info);
  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  if (pfamtblfp)     fclose(pfamtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}

#ifdef HMMER_MPI

/* Define common tags used by the MPI master/slave processes */
#define HMMER_ERROR_TAG          1
#define HMMER_HMM_TAG            2
#define HMMER_SEQUENCE_TAG       3
#define HMMER_BLOCK_TAG          4
#define HMMER_PIPELINE_TAG       5
#define HMMER_TOPHITS_TAG        6
#define HMMER_HIT_TAG            7
#define HMMER_TERMINATING_TAG    8
#define HMMER_READY_TAG          9

/* mpi_failure()
 * Generate an error message.  If the clients rank is not 0, a
 * message is created with the error message and sent to the
 * master process for handling.
 */
static void
mpi_failure(char *format, ...)
{
  va_list  argp;
  int      status = eslFAIL;
  int      len;
  int      rank;
  char     str[512];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* format the error mesg */
  va_start(argp, format);
  len = vsnprintf(str, sizeof(str), format, argp);
  va_end(argp);

  /* make sure the error string is terminated */
  str[sizeof(str)-1] = '\0';

  /* if the caller is the master, print the results and abort */
  if (rank == 0)
    {
      if (fprintf(stderr, "\nError: ") < 0) exit(eslEWRITE);
      if (fprintf(stderr, "%s", str)   < 0) exit(eslEWRITE);
      if (fprintf(stderr, "\n")        < 0) exit(eslEWRITE);
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD, status);
      exit(1);
    }
  else
    {
      MPI_Send(str, len, MPI_CHAR, 0, HMMER_ERROR_TAG, MPI_COMM_WORLD);
      pause();
    }
}

#define MAX_BLOCK_SIZE (512*1024)

typedef struct {
  uint64_t  offset;
  uint64_t  length;
  uint64_t  count;
} SEQ_BLOCK;

typedef struct {
  int        complete;
  int        size;
  int        current;
  int        last;
  SEQ_BLOCK *blocks;
} BLOCK_LIST;

/* this routine parses the database keeping track of the blocks
 * offset within the file, number of sequences and the length
 * of the block.  These blocks are passed as work units to the
 * MPI workers.  If multiple hmm's are in the query file, the
 * blocks are reused without parsing the database a second time.
 */
int next_block(ESL_SQFILE *sqfp, ESL_SQ *sq, BLOCK_LIST *list, SEQ_BLOCK *block, int n_targetseqs)
{
  int      status   = eslOK;

  /* if the list has been calculated, use it instead of parsing the database */
  if (list->complete)
    {
      if (list->current == list->last)
      {
        block->offset = 0;
        block->length = 0;
        block->count  = 0;

        status = eslEOF;
      }
      else
      {
        int inx = list->current++;

        block->offset = list->blocks[inx].offset;
        block->length = list->blocks[inx].length;
        block->count  = list->blocks[inx].count;

        status = eslOK;
      }

      return status;
    }

  block->offset = 0;
  block->length = 0;
  block->count = 0;

  esl_sq_Reuse(sq);
  if (n_targetseqs == 0) status = eslEOF; //this is to handle the end-case of a restrictdb scenario, where no more targets are required, and we want to mark the list as complete
  while (block->length < MAX_BLOCK_SIZE && (n_targetseqs <0 || block->count < n_targetseqs) && (status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      if (block->count == 0) block->offset = sq->roff;
      block->length = sq->eoff - block->offset + 1;
      block->count++;
      esl_sq_Reuse(sq);
    }

  if (block->count > 0)
    if (status == eslEOF || block->count == n_targetseqs)
      status = eslOK;
  if (status == eslEOF) list->complete = 1;

  /* add the block to the list of known blocks */
  if (status == eslOK)
    {
      int inx;

      if (list->last >= list->size)
	{
	  void *tmp;
	  list->size += 500;
	  ESL_RALLOC(list->blocks, tmp, sizeof(SEQ_BLOCK) * list->size);
	}

      inx = list->last++;
      list->blocks[inx].offset = block->offset;
      list->blocks[inx].length = block->length;
      list->blocks[inx].count  = block->count;
    }

  return status;

 ERROR:
  return eslEMEM;
}

/* mpi_master()
 * The MPI version of hmmbuild.
 * Follows standard pattern for a master/worker load-balanced MPI program (J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * p7_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-dom (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam-style tabular output  (--pfamtblout) */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              dest;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  BLOCK_LIST      *list     = NULL;
  SEQ_BLOCK        block;

  int              i;
  int              size;
  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  int              n_targets;

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) mpi_failure("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
      if (esl_opt_IsUsed(go, "--ssifile"))
        esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
      else
        esl_sqfile_OpenSSI(dbfp, NULL);
  }


  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o") && (ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
    mpi_failure("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o"));

  if (esl_opt_IsOn(go, "-A") && (afp = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) 
    mpi_failure("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A"));

  if (esl_opt_IsOn(go, "--tblout") && (tblfp = fopen(esl_opt_GetString(go, "--tblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout"));

  if (esl_opt_IsOn(go, "--domtblout") && (domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout"));

  if (esl_opt_IsOn(go, "--pfamtblout") && (pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)
    mpi_failure("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout"));

  ESL_ALLOC(list, sizeof(BLOCK_LIST));
  list->complete = 0;
  list->size     = 0;
  list->current  = 0;
  list->last     = 0;
  list->blocks   = NULL;

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      dbsq = esl_sq_CreateDigital(abc);
      bg = p7_bg_Create(abc);
    }
  

  if ( cfg->firstseq_key != NULL ) { //it's tempting to want to do this once and capture the offset position for future passes, but ncbi files make this non-trivial, so this keeps it general
    sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
    if (sstatus != eslOK)
      p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
  }

  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
      P7_PIPELINE     *pli     = NULL;
      P7_TOPHITS      *th      = NULL;
      int              seq_cnt = 0;
      nquery++;
      esl_stopwatch_Start(w);

      n_targets = cfg->n_targetseq;

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)   list->current = 0;

      if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
      if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
      p7_oprofile_Convert(gm, om);

      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, hmm->M, 100, FALSE, p7_SEARCH_SEQS);
      p7_pli_NewModel(pli, om, bg);

      /* Main loop: */
      while ((n_targets==-1 || seq_cnt<=n_targets) && (sstatus = next_block(dbfp, dbsq, list, &block, n_targets-seq_cnt)) == eslOK )
      {
        seq_cnt += block.count;

        if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0)
          mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

        MPI_Get_count(&mpistatus, MPI_PACKED, &size);
        if (mpi_buf == NULL || size > mpi_size) {
          void *tmp;
          ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
          mpi_size = size;
        }

        dest = mpistatus.MPI_SOURCE;
        MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

        if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
          mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
        if (mpistatus.MPI_TAG != HMMER_READY_TAG)
          mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);

        MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
      }

      if (n_targets!=-1 && seq_cnt==n_targets)
        sstatus = eslEOF;

      switch(sstatus)
      {
      case eslEFORMAT:
        mpi_failure("Parse failed (sequence file %s):\n%s\n", dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        break;
      case eslEOF:
        break;
      default:
        mpi_failure("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      block.offset = 0;
      block.length = 0;
      block.count  = 0;

      /* wait for all workers to finish up their work blocks */
      for (i = 1; i < cfg->nproc; ++i)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	}

      /* merge the results of the search results */
      for (dest = 1; dest < cfg->nproc; ++dest)
	{
	  P7_PIPELINE     *mpi_pli   = NULL;
	  P7_TOPHITS      *mpi_th    = NULL;

	  /* send an empty block to signal the worker they are done */
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = p7_tophits_MPIRecv(dest, HMMER_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = p7_pipeline_MPIRecv(dest, HMMER_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  p7_tophits_Merge(th, mpi_th);
	  p7_pipeline_Merge(pli, mpi_pli);

	  p7_pipeline_Destroy(mpi_pli);
	  p7_tophits_Destroy(mpi_th);
	}

      /* Print the results.  */
      p7_tophits_SortBySortkey(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, th, pli, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, th, pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, th, pli, (nquery == 1));
      if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp, hmm->name, hmm->acc, th, pli);

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if (p7_tophits_Alignment(th, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
	  {
	    esl_msa_SetName     (msa, hmm->name, -1);
	    esl_msa_SetAccession(msa, hmm->acc,  -1);
	    esl_msa_SetDesc     (msa, hmm->desc, -1);
	    esl_msa_FormatAuthor(msa, "hmmsearch (HMMER %s)", HMMER_VERSION);

	    if (textw > 0) esl_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msafile_Write(afp, msa, eslMSAFILE_PFAM);
	  
	    if (fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
	  } 
	else { if (fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
	  
	esl_msa_Destroy(msa);
      }

      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus) {
  case eslEOD:       mpi_failure("read failed, HMM file %s may be truncated?", cfg->hmmfile);      break;
  case eslEFORMAT:   mpi_failure("bad file format in HMM file %s",             cfg->hmmfile);      break;
  case eslEINCOMPAT: mpi_failure("HMM file %s contains different alphabets",   cfg->hmmfile);      break;
  case eslEOF:       /* EOF is good, that's what we expect here */                                 break;
  default:           mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
  }

  /* monitor all the workers to make sure they have ended */
  for (i = 1; i < cfg->nproc; ++i)
    {
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

      MPI_Get_count(&mpistatus, MPI_PACKED, &size);
      if (mpi_buf == NULL || size > mpi_size) {
	void *tmp;
	ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	mpi_size = size; 
      }

      dest = mpistatus.MPI_SOURCE;
      MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

      if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
      if (mpistatus.MPI_TAG != HMMER_TERMINATING_TAG)
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    }

  /* Terminate outputs... any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,     "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp,  "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (pfamtblfp)p7_tophits_TabularTail(pfamtblfp, "hmmsearch", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (ofp)     { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit
   */
  free(list);
  if (mpi_buf != NULL) free(mpi_buf);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);

  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);
  if (pfamtblfp)     fclose(pfamtblfp);

  return eslOK;

 ERROR:
  return eslEMEM;
}


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures           */
  int              mpi_size = 0;                 /* size of the allocated buffer                    */

  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      dbsq = esl_sq_CreateDigital(abc);
      bg = p7_bg_Create(abc);
      esl_sqfile_SetDigital(dbfp, abc);
    }
  
  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */
      P7_PIPELINE     *pli     = NULL;
      P7_TOPHITS      *th      = NULL;

      SEQ_BLOCK        block;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL);
      p7_oprofile_Convert(gm, om);

      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      p7_pli_NewModel(pli, om, bg);

      /* receive a sequence block from the master */
      MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
      while (block.count > 0)
	{
	  uint64_t length = 0;
	  uint64_t count  = block.count;

	  status = esl_sqfile_Position(dbfp, block.offset);
	  if (status != eslOK) mpi_failure("Cannot position sequence database to %ld\n", block.offset);

	  while (count > 0 && (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
	    {
	      length = dbsq->eoff - block.offset + 1;

	      p7_pli_NewSeq(pli, dbsq);
	      p7_bg_SetLength(bg, dbsq->n);
	      p7_oprofile_ReconfigLength(om, dbsq->n);
      
	      p7_Pipeline(pli, om, bg, dbsq, NULL, th);

	      esl_sq_Reuse(dbsq);
	      p7_pipeline_Reuse(pli);

	      --count;
	    }

	  /* lets do a little bit of sanity checking here to make sure the blocks are the same */
	  if (count > 0)              mpi_failure("Block count mismatch - expected %ld found %ld at offset %ld\n",  block.count,  block.count - count, block.offset);
	  if (block.length != length) mpi_failure("Block length mismatch - expected %ld found %ld at offset %ld\n", block.length, length,              block.offset);

	  /* inform the master we need another block of sequences */
	  status = 0;
	  MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

	  /* wait for the next block of sequences */
	  MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
	}

      esl_stopwatch_Stop(w);

      /* Send the top hits back to the master. */
      p7_tophits_MPISend(th, 0, HMMER_TOPHITS_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);

      p7_pipeline_Destroy(pli);
      p7_tophits_Destroy(th);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);

      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus)
    {
    case eslEOF:
      /* do nothing */
      break;
    case eslEFORMAT:
      mpi_failure("bad file format in HMM file %s", cfg->hmmfile);
      break;
    case eslEINCOMPAT:
      mpi_failure("HMM file %s contains different alphabets", cfg->hmmfile);
      break;
    default:
      mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
    }

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, HMMER_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);

  p7_bg_Destroy(bg);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  return eslOK;
}
#endif /*HMMER_MPI*/

static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp, int n_targetseqs)
{
  int      sstatus;
  ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */
  int seq_cnt = 0;

  dbsq = esl_sq_CreateDigital(info->om->abc);

  /* Main loop: */
  while ( (n_targetseqs==-1 || seq_cnt<n_targetseqs) &&  (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
  {
      p7_pli_NewSeq(info->pli, dbsq);
      p7_bg_SetLength(info->bg, dbsq->n);
      p7_oprofile_ReconfigLength(info->om, dbsq->n);
      
      p7_Pipeline(info->pli, info->om, info->bg, dbsq, NULL, info->th);

      seq_cnt++;
      esl_sq_Reuse(dbsq);
      p7_pipeline_Reuse(info->pli);
  }

  if (n_targetseqs!=-1 && seq_cnt==n_targetseqs)
    sstatus = eslEOF;

  esl_sq_Destroy(dbsq);

  return sstatus;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, int n_targetseqs)
{
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK )
    {
      block = (ESL_SQ_BLOCK *) newBlock;

      if (n_targetseqs == 0)
      {
        block->count = 0;
        sstatus = eslEOF;
      } else {
        sstatus = esl_sqio_ReadBlock(dbfp, block, -1, n_targetseqs, /*max_init_window=*/FALSE, FALSE);
        n_targetseqs -= block->count;
      }

      if (sstatus == eslEOF)
      {
        if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
        ++eofCount;
      }

      if (sstatus == eslOK)
      {
        status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
        if (status != eslOK) esl_fatal("Work queue reader failed");
      }
    }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
    {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
    }

  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;

  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;
  
  impl_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (ESL_SQ_BLOCK *) newBlock;
  while (block->count > 0)
    {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
	{
	  ESL_SQ *dbsq = block->list + i;

	  p7_pli_NewSeq(info->pli, dbsq);
	  p7_bg_SetLength(info->bg, dbsq->n);
	  p7_oprofile_ReconfigLength(info->om, dbsq->n);
	  
	  p7_Pipeline(info->pli, info->om, info->bg, dbsq, NULL, info->th);
	  
	  esl_sq_Reuse(dbsq);
	  p7_pipeline_Reuse(info->pli);
	}

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (ESL_SQ_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* HMMER_THREADS */
 


