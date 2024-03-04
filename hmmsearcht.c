/* hmmsearcht: search profile HMM(s) against a sequence database.*/

#include "p7_config.h"

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

/* for hmmsearcht */
#include "esl_gencode.h"

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG                 *bg;	 /* null model                                                               */
  ESL_SQ                *ntqsq;  /* query or target sequence; this is a DNA sequence in the case of hmmscant */
  P7_PIPELINE           *pli;    /* work pipeline                                                            */
  P7_TOPHITS            *th;     /* top hit results                                                          */
  P7_OPROFILE           *om;     /* optimized query profile                                                  */
  ESL_GENCODE           *gcode;  /* used for translating ORFs                                                */
  ESL_GENCODE_WORKSTATE *wrk;    /* maintain state of nucleotide sequence in the midst of processing ORFs    */
} WORKER_INFO;

typedef struct {
  int    id;         /* internal sequence ID  */
  int    length;     /* length of sequence */
} ID_LENGTH;

typedef struct {
  ID_LENGTH  *id_lengths;
  int        count;
  int        size;
} ID_LENGTH_LIST;




















static int ALEX_DEBUG = 1;


static char AMINO_CHARS[21] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};
static char DNA_CHARS[5]    = {'A','C','G','T','-'};
static char RNA_CHARS[5]    = {'A','C','G','U','-'};

// How many amino acids are we willing to extend to bridge two hits?  
// How many overlapping aminos do we require to perform bridging?
static int MAX_AMINO_EXT     = 5;
static int MIN_AMINO_OVERLAP = 4;



//
//  TARGET_SEQ
//
typedef struct _target_seq {
  
  ESL_SQ  * esl_sq; // For freeing, we'll need this pointer
  ESL_DSQ * Seq;
  const ESL_ALPHABET * abc;

  int64_t start;
  int64_t end;

} TARGET_SEQ;



//
//  DOM_OVERLAP
//
typedef struct _dom_overlap {

  const ESL_ALPHABET * ntalpha;

  int amino_start;
  int amino_end;


  uint64_t upstream_hit_id;
  int      upstream_dom_id;
  int      upstream_nucl_start;
  int      upstream_nucl_end;

  P7_DOMAIN * UpstreamDomain;
  ESL_DSQ   * UpstreamNucls;
  

  uint64_t downstream_hit_id;
  int      downstream_dom_id;
  int      downstream_nucl_start;
  int      downstream_nucl_end;

  P7_DOMAIN * DownstreamDomain;
  ESL_DSQ   * DownstreamNucls;


  int optimal_splice_index;


} DOM_OVERLAP;





/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DumpDOM_OVERLAP
 *
 *  Inputs:
 *
 *  Output:
 *
 */
void DumpDOM_OVERLAP
(DOM_OVERLAP * D_O)
{

  printf("\n");
  printf("  Aminos  : %d..%d\n",D_O->amino_start,D_O->amino_end);

  // DOMAIN 1
  printf("\n");
  printf("  Domain 1: %d..%d\n",D_O->upstream_nucl_start,D_O->upstream_nucl_end);
  printf("            ");
  for (int i=1; i<=abs(D_O->upstream_nucl_start-D_O->upstream_nucl_end)+1; i++) {
    if (i % 60 == 0) printf("\n            ");
    printf("%c",DNA_CHARS[D_O->UpstreamNucls[i]]);
  }  
  printf("\n\n");

  // DOMAIN 2
  printf("  Domain 2: %d..%d\n",D_O->downstream_nucl_start,D_O->downstream_nucl_end);
  printf("            ");
  for (int i=1; i<=abs(D_O->downstream_nucl_start-D_O->downstream_nucl_end)+1; i++) {
    if (i % 60 == 0) printf("\n            ");
    printf("%c",DNA_CHARS[D_O->DownstreamNucls[i]]);
  }  
  printf("\n\n");

  printf("\n");

}






/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: InitEmptyDSQ
 *
 *  Inputs:
 *
 *  Output:
 *
 */
void InitEmptyDSQ
(const ESL_ALPHABET * alphabet, int len, ESL_DSQ * Target)
{

  char * DummyStr = malloc(len * sizeof(char));
  for (int i=0; i<len; i++) 
    DummyStr[i] = 'A';
  
  esl_abc_CreateDsq(alphabet,DummyStr,&Target);

  free(DummyStr);

}





/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GrabNuclRange
 *
 *  Inputs:
 *
 *  Output:
 *
 */
ESL_DSQ * GrabNuclRange
(TARGET_SEQ * TargetNuclSeq, int start, int end)
{

  int len = abs(end - start) + 1;

  char * Seq = malloc(len * sizeof(char));

  // Keep in mind that DSQs are [1..n]
  int read_index = (int)((uint64_t)start - TargetNuclSeq->start) + 1;

  if (start < end) {
  
    for (int i=0; i<len; i++)
      Seq[i] = DNA_CHARS[TargetNuclSeq->Seq[read_index++]];

  } else {

    for (int i=0; i<len; i++) {
      int fwd_nucl_code = TargetNuclSeq->Seq[read_index--];
      if (fwd_nucl_code < 4) Seq[i] = DNA_CHARS[3 - fwd_nucl_code];
      else                   Seq[i] ='N';
    }

  }


  ESL_DSQ * NuclSubseq;
  esl_abc_CreateDsq(TargetNuclSeq->abc,Seq,&NuclSubseq);

  free(Seq);

  return NuclSubseq;

}





/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DetermineNuclType
 *
 *  Inputs:
 *
 *  Output:
 *
 */
int DetermineNuclType 
(P7_ALIDISPLAY * AliDisplay)
{
  for (int i=0; i<strlen(AliDisplay->ntseq); i++) {
    if (AliDisplay->ntseq[i] == 'T') return eslDNA;
    if (AliDisplay->ntseq[i] == 'U') return eslRNA;
  }
  return 0;
}






/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindOptimalSpliceSite
 *
 *  Inputs:  
 *
 *  Output:  
 *
 */
void FindOptimalSpliceSite
(
  DOM_OVERLAP * SpliceEdge,
  int           upstream_hmm_to,
  int           downstream_hmm_to,
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
  int mi_transit =  om->M + 2;
  int md_transit = (om->M + 2) * 2;
  int ii_transit = (om->M + 2) * 3;
  int im_transit = (om->M + 2) * 4;
  int dd_transit = (om->M + 2) * 5;
  int dm_transit = (om->M + 2) * 6;

  p7_oprofile_GetFwdTransitionArray(om,1,&FwdTransitionProbs[mm_transit]);
  p7_oprofile_GetFwdTransitionArray(om,5,&FwdTransitionProbs[mi_transit]);
  p7_oprofile_GetFwdTransitionArray(om,4,&FwdTransitionProbs[md_transit]);
  p7_oprofile_GetFwdTransitionArray(om,6,&FwdTransitionProbs[ii_transit]);
  p7_oprofile_GetFwdTransitionArray(om,2,&FwdTransitionProbs[im_transit]);
  p7_oprofile_GetFwdTransitionArray(om,7,&FwdTransitionProbs[dd_transit]);
  p7_oprofile_GetFwdTransitionArray(om,3,&FwdTransitionProbs[dm_transit]);


  // Convert to log scores (natural)
  for (int i=0; i<7*(om->M+2); i++)
    FwdTransitionProbs[i] = log(FwdTransitionProbs[i]);
  

  // We'll also want to have some simple way of acknowledging whether a
  // splice site looks good (i.e., is canonical)
  float    canon_ss_score = log(0.98);
  float noncanon_ss_score = log(0.02);


  // We'll need a background model for... stuff?
  P7_BG * BackgroundModel = p7_bg_Create(om->abc);
  float null_transit_prob = BackgroundModel->p1;


  // The above work sets up up to be able to interact with the model that
  // provided the hits we're going to attempt to splice into one another.
  //
  // Now we can get to work with actually finding an optimal splice site!




}







/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: AttemptSpliceEdge
 *
 *  Inputs:  
 *
 *  Output:
 *
 */
void AttemptSpliceEdge
(
  P7_TOPHITS  * TopHits,
  DOM_OVERLAP * SpliceEdge,
  TARGET_SEQ  * TargetNuclSeq,
  P7_OPROFILE * om,
  ESL_GENCODE * gcode
)
{

  P7_DOMAIN *   UpDom = &(TopHits->hit[SpliceEdge->upstream_hit_id]->dcl[SpliceEdge->upstream_dom_id]);
  P7_DOMAIN * DownDom = &(TopHits->hit[SpliceEdge->downstream_hit_id]->dcl[SpliceEdge->downstream_dom_id]);


  int   upstream_hmm_from =   UpDom->ad->hmmfrom;
  int   upstream_hmm_to   =   UpDom->ad->hmmto;
  int downstream_hmm_from = DownDom->ad->hmmfrom;
  int downstream_hmm_to   = DownDom->ad->hmmto;
  
  int   upstream_nt_from =   UpDom->ad->sqfrom;
  int   upstream_nt_to   =   UpDom->ad->sqto;
  int downstream_nt_from = DownDom->ad->sqfrom;
  int downstream_nt_to   = DownDom->ad->sqto;
  

  int strand = 3;
  if (upstream_nt_from > upstream_nt_to) 
    strand = -3;


  SpliceEdge->upstream_nucl_end = upstream_nt_to;
  SpliceEdge->downstream_nucl_start = downstream_nt_from;


  int amino_overlap  = 1 + upstream_hmm_to - downstream_hmm_from;
  int codons_to_pull = MIN_AMINO_OVERLAP - amino_overlap;
  if (codons_to_pull > 0) {

    SpliceEdge->upstream_nucl_end     += strand * codons_to_pull;
    SpliceEdge->downstream_nucl_start -= strand * codons_to_pull;

    amino_overlap = MIN_AMINO_OVERLAP;

  } else {

    codons_to_pull = 0;

  }


  // Log the start and end amino acid coord.s
  SpliceEdge->amino_start = upstream_hmm_to - ((amino_overlap-1) - codons_to_pull);
  SpliceEdge->amino_end   = SpliceEdge->amino_start + (amino_overlap-1);


  // Prime these coordinates
  SpliceEdge->upstream_nucl_start = SpliceEdge->upstream_nucl_end     - (codons_to_pull * strand);
  SpliceEdge->downstream_nucl_end = SpliceEdge->downstream_nucl_start + (codons_to_pull * strand);


  // Find the upstream nucleotide start (we need to be mindful of gaps)
  P7_ALIDISPLAY * AD = UpDom->ad;
  int ad_pos = AD->N-1;
  int overlap_aminos_covered = 0;
  while (overlap_aminos_covered<amino_overlap) {
    
    if (AD->model[ad_pos] == '.') { // // // insertion relative to hmm

      SpliceEdge->upstream_nucl_start -= strand;

    } else if (AD->aseq[ad_pos] == '-') { // insertion relative to genome

      overlap_aminos_covered++;

    } else { // // // // // // // // // // // Match state! (easy!)

      SpliceEdge->upstream_nucl_start -= strand;
      overlap_aminos_covered++;

    }

    ad_pos--;
  
  }


  AD = DownDom->ad;
  ad_pos = 0;
  overlap_aminos_covered = 0;
  while (overlap_aminos_covered<amino_overlap) {

    if (AD->model[ad_pos] == '.') { // // // insertion relative to hmm

      SpliceEdge->downstream_nucl_end += strand;

    } else if (AD->aseq[ad_pos] == '-') { // insertion relative to genome

      overlap_aminos_covered++;

    } else { // // // // // // // // // // // Match state! (easy!)

      SpliceEdge->downstream_nucl_end += strand;
      overlap_aminos_covered++;

    }

    ad_pos++;

  }


  // We'll need to trim the last nucleotide to have inclusive bounds
  if (strand < 0) {
    SpliceEdge->upstream_nucl_start--;
    SpliceEdge->downstream_nucl_end++;
  } else {
    SpliceEdge->upstream_nucl_start++;
    SpliceEdge->downstream_nucl_end--;
  }


  // Grab them nucleotides!
  SpliceEdge->UpstreamNucls   = GrabNuclRange(TargetNuclSeq,SpliceEdge->upstream_nucl_start,SpliceEdge->upstream_nucl_end);
  SpliceEdge->DownstreamNucls = GrabNuclRange(TargetNuclSeq,SpliceEdge->downstream_nucl_start,SpliceEdge->downstream_nucl_end);


  // Finish off by adding these friendly little pointers
  SpliceEdge->ntalpha          = TargetNuclSeq->abc;
  SpliceEdge->UpstreamDomain   =   UpDom;
  SpliceEdge->DownstreamDomain = DownDom;

  // DEBUGGING
  DumpDOM_OVERLAP(SpliceEdge);

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
int DomainsAreSpliceCompatible
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
    return 0;

  // Do we have overlap OR sufficient proximity to consider
  // extending?
  if (!(amino_end_1 + MAX_AMINO_EXT >= amino_start_2))
    return 0;


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
    return 0;


  // We want to make sure that these aren't unrealistically
  // close together on the genome...
  if (revcomp1) {

    if (nucl_start_2 + (3 * MAX_AMINO_EXT) >= nucl_end_1)
      return 0;

  } else {

    if (nucl_start_2 - (3 * MAX_AMINO_EXT) <= nucl_end_1)
      return 0;

  }


  // Looks like we've got a viable upstream / downstream pair!
  return 1;


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
  uint64_t upstream_hit_id,
  DOM_OVERLAP *** SpliceEdges
)
{

  P7_HIT * UpstreamHit  = TopHits->hit[upstream_hit_id];
  int num_upstream_doms = UpstreamHit->ndom;


  // We'll setup a temporary array to hold the indices
  // of hits that could serve as downstream exons for our
  // candidate (upstream) hit
  int vc_capacity  = 8;
  int vc_occupancy = 0;
  int      * ViableUpstreamDoms   = malloc(vc_capacity * sizeof(int));
  uint64_t * ViableDownstreamHits = malloc(vc_capacity * sizeof(uint64_t));
  int      * ViableDownstreamDoms = malloc(vc_capacity * sizeof(int));


  // For each hit, gather all of the indices of other hits that
  // could potentially be downstream exons.
  uint64_t num_hits = TopHits->N;
  for (uint64_t downstream_hit_id=0; downstream_hit_id < num_hits; downstream_hit_id++) {


    // We only consider splicing when the two hits come 
    // from the same sequence.  Further, because of how
    // we determine splice viability (specifically, by
    // checking compatibility of amino acid start/end
    // coordinates), this wouldn't be a trivial check to
    // remove to allow for splicing across sequences.
    if (strcmp(UpstreamHit->name,TopHits->hit[downstream_hit_id]->name))
      continue;


    // Now that we know these are valid to check, check 'em!
    P7_HIT * DownstreamHit  = TopHits->hit[downstream_hit_id];
    int num_downstream_doms = DownstreamHit->ndom;


    for (int upstream_dom_id = 0; upstream_dom_id < num_upstream_doms; upstream_dom_id++) {

      P7_DOMAIN * UpstreamDomain = &UpstreamHit->dcl[upstream_dom_id];

      for (int downstream_dom_id = 0; downstream_dom_id < num_downstream_doms; downstream_dom_id++) {

        // NO SELF-SPLICING, YOU ABSOLUTE MANIAC!
        if (upstream_hit_id == downstream_hit_id && upstream_dom_id == downstream_dom_id)
          continue;


        P7_DOMAIN * DownstreamDomain = &DownstreamHit->dcl[downstream_dom_id];

        if (DomainsAreSpliceCompatible(UpstreamDomain->ad,DownstreamDomain->ad)) {

          // Do we need to resize before entering new info.?
          if (vc_occupancy == vc_capacity) {

            int      * TmpUpDoms   = malloc(vc_occupancy * sizeof(int));
            uint64_t * TmpDownHits = malloc(vc_occupancy * sizeof(uint64_t));
            int      * TmpDownDoms = malloc(vc_occupancy * sizeof(int));

            for (int i=0; i<vc_occupancy; i++) {
              TmpUpDoms[i]   = ViableUpstreamDoms[i];
              TmpDownHits[i] = ViableDownstreamHits[i];
              TmpDownDoms[i] = ViableDownstreamDoms[i];
            }

            free(ViableUpstreamDoms);
            free(ViableDownstreamHits);
            free(ViableDownstreamDoms);

            vc_capacity *= 2;
            ViableUpstreamDoms   = malloc(vc_capacity * sizeof(int));
            ViableDownstreamHits = malloc(vc_capacity * sizeof(uint64_t));
            ViableDownstreamDoms = malloc(vc_capacity * sizeof(int));

            for (int i=0; i<vc_occupancy; i++) {
              ViableUpstreamDoms[i]   = TmpUpDoms[i];
              ViableDownstreamHits[i] = TmpDownHits[i];
              ViableDownstreamDoms[i] = TmpDownDoms[i];
            }

            free(TmpUpDoms);
            free(TmpDownHits);
            free(TmpDownDoms);

          }

          // Record that splice compatibility!
          ViableUpstreamDoms[vc_occupancy]   = upstream_dom_id;
          ViableDownstreamHits[vc_occupancy] = downstream_hit_id;
          ViableDownstreamDoms[vc_occupancy] = downstream_dom_id;

          vc_occupancy++;

        }

      }

    }


  }


  // Did we find any viable downstream exons for this hit?
  // If so, copy 'em over!
  if (vc_occupancy > 0) {
    
    SpliceEdges[upstream_hit_id] = (DOM_OVERLAP **)malloc(vc_occupancy * sizeof(DOM_OVERLAP *));

    for (int i=0; i<vc_occupancy; i++) {
          
      SpliceEdges[upstream_hit_id][i] = (DOM_OVERLAP *)malloc(sizeof(DOM_OVERLAP));
      
      // Initialize this DOM_OVERLAP entry
      DOM_OVERLAP * OverlapIniter = SpliceEdges[upstream_hit_id][i];

      OverlapIniter->upstream_hit_id   = upstream_hit_id;
      OverlapIniter->upstream_dom_id   = ViableUpstreamDoms[i];
      OverlapIniter->downstream_hit_id = ViableDownstreamHits[i];
      OverlapIniter->downstream_dom_id = ViableDownstreamDoms[i];

    }

  }


  free(ViableUpstreamDoms);
  free(ViableDownstreamHits);
  free(ViableDownstreamDoms);


  return vc_occupancy;

}






/* * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Debugging Function: TestNuclGrab
 *
 */
void TestNuclGrab
(TARGET_SEQ * TargetNuclSeq)
{

  int start = (int)TargetNuclSeq->start + 100;
  int end   = (int)TargetNuclSeq->end   - 100;

  ESL_DSQ * FwdRead = GrabNuclRange(TargetNuclSeq,start,end);
  ESL_DSQ * RevRead = GrabNuclRange(TargetNuclSeq,end,start);

  printf("\n\n  %d..%d\n",start,end);

  int line_chars = 60;
  for (int i=1; i <= (end-start)+1; i += line_chars) {

    printf("\n  ");

    int read_index = i;
    while (read_index < i+line_chars && read_index <= (end-start)+1)
      printf("%c",DNA_CHARS[FwdRead[read_index++]]);
    printf("\n  ");

    read_index = i;
    while (read_index < i+line_chars && read_index <= (end-start)+1)
      printf("%c",DNA_CHARS[RevRead[read_index++]]);
    printf("\n");

  }
  printf("\n");

  exit(69);
}






/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetMinAndMaxCoords
 *
 *  Inputs:  
 *
 */
void GetMinAndMaxCoords
(P7_TOPHITS * TopHits, TARGET_SEQ * NuclTargetSeq)
{

  // First, let's figure out if we're revcomp
  int revcomp = 0;
  int64_t min,max;

  P7_HIT    * Hit;
  P7_DOMAIN * Dom;
  Dom = &(TopHits->hit[0]->dcl[0]);
  if (Dom->ad->sqfrom > Dom->ad->sqto) {

    min = Dom->ad->sqto;
    max = Dom->ad->sqfrom;
    revcomp = 1;

  } else {

    min = Dom->ad->sqfrom;
    max = Dom->ad->sqto;

  }


  if (revcomp) {

    for (uint64_t hit_id = 0; hit_id < TopHits->N; hit_id++) {

      Hit = TopHits->hit[hit_id];
      for (int dom_id = 0; dom_id < Hit->ndom; dom_id++) {

        Dom = &(Hit->dcl[dom_id]);

        if (Dom->ad->sqfrom > max)
          max = Dom->ad->sqfrom;

        if (Dom->ad->sqto < min)
          min = Dom->ad->sqto;

      }

    }

  } else {

    for (uint64_t hit_id = 0; hit_id < TopHits->N; hit_id++) {

      Hit = TopHits->hit[hit_id];
      for (int dom_id = 0; dom_id < Hit->ndom; dom_id++) {

        Dom = &(Hit->dcl[dom_id]);

        if (Dom->ad->sqto > max)
          max = Dom->ad->sqto;

        if (Dom->ad->sqfrom < min)
          min = Dom->ad->sqfrom;

      }

    }

  }

  NuclTargetSeq->start = min;
  NuclTargetSeq->end   = max;

}







/* * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetTargetNuclSeq
 *
 *  Inputs:  
 *
 */
TARGET_SEQ * GetTargetNuclSeq
(ESL_SQFILE * GenomicSeqFile, P7_TOPHITS * TopHits)
{

  TARGET_SEQ * TargetNuclSeq = (TARGET_SEQ *)malloc(sizeof(TARGET_SEQ));
  GetMinAndMaxCoords(TopHits,TargetNuclSeq);

  TargetNuclSeq->abc = GenomicSeqFile->abc;

  ESL_SQFILE * TmpSeqFile;
  esl_sqfile_Open(GenomicSeqFile->filename,GenomicSeqFile->format,NULL,&TmpSeqFile);
  esl_sqfile_OpenSSI(TmpSeqFile,NULL);

  TargetNuclSeq->esl_sq = esl_sq_CreateDigital(TargetNuclSeq->abc);
  int fetch_err_code    = esl_sqio_FetchSubseq(TmpSeqFile,TopHits->hit[0]->name,TargetNuclSeq->start,TargetNuclSeq->end,TargetNuclSeq->esl_sq);

  if (fetch_err_code != eslOK) {
    fprintf(stderr,"\n  ERROR fetching subsequence\n\n");
    exit(1);
  }

  TargetNuclSeq->Seq = TargetNuclSeq->esl_sq->dsq;

  return TargetNuclSeq;

  // DESTROY TmpSeqFile

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
(P7_TOPHITS * TopHits, ESL_SQFILE * GenomicSeqFile, P7_OPROFILE * om, ESL_GENCODE * gcode)
{

  // Very first thing we want to do is make sure that our hits are
  // organized by the 'Seqidx' (the target genomic sequence) and position
  // within that file.
  //
  // NOTE: I still need to take advantage of this sorting!
  //
  p7_tophits_SortBySeqidxAndAlipos(TopHits);

  uint64_t num_hits = TopHits->N;


  // For each hit, how many possible splice edges are there (all domains)?
  //
  int * ViableSpliceTargets = malloc(num_hits * sizeof(int));
  DOM_OVERLAP *** SpliceEdges = (DOM_OVERLAP ***)malloc(num_hits * sizeof(DOM_OVERLAP **));

  for (uint64_t hit_id = 0; hit_id < num_hits; hit_id++) {
    
    SpliceEdges[hit_id]         = NULL;
    ViableSpliceTargets[hit_id] = GatherViableDownstreamHits(TopHits,hit_id,SpliceEdges);

  }


  // Given that our hits are organized by target sequence, we can
  // be a bit more efficient in our file reading by only pulling
  // target sequences as they change (wrt the upstream hit)
  TARGET_SEQ * TargetNuclSeq = GetTargetNuclSeq(GenomicSeqFile,TopHits);


  // Now we can run through all of our paired domains and actually
  // splice 'em up (or at least try our best to)!
  for (uint64_t upstream_hit_id = 0; upstream_hit_id < num_hits; upstream_hit_id++) {

    for (int splice_edge_id = 0; splice_edge_id < ViableSpliceTargets[upstream_hit_id]; splice_edge_id++)
      AttemptSpliceEdge(TopHits,SpliceEdges[upstream_hit_id][splice_edge_id],TargetNuclSeq,om,gcode);

  }


}


























static ID_LENGTH_LIST* init_id_length( int size );
static void            destroy_id_length( ID_LENGTH_LIST *list );
static int             add_id_length(ID_LENGTH_LIST *list, int id, int L);
static int             assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list);

/* set the max residue count to 1/4 meg when reading a block */
#define HMMSEARCHT_MAX_RESIDUE_COUNT (1024 * 256)  /* 1/4 Mb */

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL

static ESL_OPTIONS options[] = {
  /* name             type            default env          range     toggles  reqs  incomp              help                                                      docgroup*/
  { "-h",             eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "show brief help on version and usage",                        1 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE, NULL,   NULL,        NULL,       NULL,  NULL,  NULL,             "direct output to file <f>, not stdout",                       2 },
  { "-A",             eslARG_OUTFILE, NULL,   NULL,        NULL,       NULL,  NULL,  NULL,             "save multiple alignment of all hits to file <f>",             2 },
  { "--tblout",       eslARG_OUTFILE, NULL,   NULL,        NULL,       NULL,  NULL,  NULL,             "save parseable table of per-sequence hits to file <f>",       2 },
  { "--domtblout",    eslARG_OUTFILE, NULL,   NULL,        NULL,       NULL,  NULL,  NULL,             "save parseable table of per-domain hits to file <f>",         2 },
  { "--acc",          eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "prefer accessions over names in output",                      2 },
  { "--noali",        eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "don't output alignments, so output is smaller",               2 },
  { "--notrans",      eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "don't show the translated DNA sequence in domain alignment",  2 }, /*for hmmsearcht */
  { "--vertcodon",    eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "show the DNA codon vertically in domain alignment",           2 }, /*for hmmsearcht */
  { "--notextw",      eslARG_NONE,    NULL,   NULL,        NULL,       NULL,  NULL,  "--textw",        "unlimit ASCII text output line width",                        2 },
  { "--textw",        eslARG_INT,     "120",  NULL,        "n>=120",   NULL,  NULL,  "--notextw",      "set max width of ASCII text output lines",                    2 },
  /* Tranlsation Options */
  { "-c",             eslARG_INT,     "1",    NULL,        NULL,       NULL,  NULL,  NULL,             "use alt genetic code of NCBI transl table <n>",               15 },
  { "-l",             eslARG_INT,     "20",   NULL,        NULL,       NULL,  NULL,  NULL,             "minimum ORF length",                                          15 },
  { "-m",             eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  "-M",             "ORFs must initiate with AUG only",                            15 },
  { "-M",             eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  "-m",             "ORFs must start with allowed initiation codon",               15 },
  { "--watson",       eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "only translate top strand",                                   15 },
  { "--crick",        eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  NULL,             "only translate bottom strand",                                15 },
  /* Control of reporting thresholds */
  { "-E",             eslARG_REAL,    "10.0", NULL,        "x>0",      NULL,  NULL,  REPOPTS,          "report sequences <= this E-value threshold in output",        4 },
  { "-T",             eslARG_REAL,    FALSE,  NULL,        NULL,       NULL,  NULL,  REPOPTS,          "report sequences >= this score threshold in output",          4 },
  { "--domE",         eslARG_REAL,    "10.0", NULL,        "x>0",      NULL,  NULL,  DOMREPOPTS,       "report domains <= this E-value threshold in output",          4 },
  { "--domT",         eslARG_REAL,    FALSE,  NULL,        NULL,       NULL,  NULL,  DOMREPOPTS,       "report domains >= this score cutoff in output",               4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",         eslARG_REAL,    "0.01", NULL,        "x>0",      NULL,  NULL,  INCOPTS,          "consider sequences <= this E-value threshold as significant", 5 },
  { "--incT",         eslARG_REAL,    FALSE,  NULL,        NULL,       NULL,  NULL,  INCOPTS,          "consider sequences >= this score threshold as significant",   5 },
  { "--incdomE",      eslARG_REAL,    "0.01", NULL,        "x>0",      NULL,  NULL,  INCDOMOPTS,       "consider domains <= this E-value threshold as significant",   5 },
  { "--incdomT",      eslARG_REAL,    FALSE,  NULL,        NULL,       NULL,  NULL,  INCDOMOPTS,       "consider domains >= this score threshold as significant",     5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",       eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  THRESHOPTS,       "use profile's GA gathering cutoffs to set all thresholding",  6 },
  { "--cut_nc",       eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  THRESHOPTS,       "use profile's NC noise cutoffs to set all thresholding",      6 },
  { "--cut_tc",       eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  THRESHOPTS,       "use profile's TC trusted cutoffs to set all thresholding",    6 },
  /* Control of acceleration pipeline */
  { "--max",          eslARG_NONE,    FALSE,  NULL,        NULL,       NULL,  NULL,  "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",     7 },
  { "--F1",           eslARG_REAL,    "0.02", NULL,        NULL,       NULL,  NULL,  "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",            7 },
  { "--F2",           eslARG_REAL,    "1e-3", NULL,        NULL,       NULL,  NULL,  "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",            7 },
  { "--F3",           eslARG_REAL,    "1e-5", NULL,        NULL,       NULL,  NULL,  "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",            7 },
  { "--nobias",       eslARG_NONE,    NULL,   NULL,        NULL,       NULL,  NULL,  "--max",          "turn off composition bias filter",                            7 },

  /* Other options */
  { "--nonull2",      eslARG_NONE,    NULL,   NULL,        NULL,       NULL,  NULL,  NULL,             "turn off biased composition score corrections",               12 },
  { "-Z",             eslARG_REAL,    FALSE,  NULL,        "x>0",      NULL,  NULL,  NULL,             "set # of comparisons done, for E-value calculation",          12 },
  { "--domZ",         eslARG_REAL,    FALSE,  NULL,        "x>0",      NULL,  NULL,  NULL,             "set # of significant seqs, for domain E-value calculation",   12 },
  { "--seed",         eslARG_INT,     "42",   NULL,        "n>=0",     NULL,  NULL,  NULL,             "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--tformat",      eslARG_STRING,  NULL,   NULL,        NULL,       NULL,  NULL,  NULL,             "assert target <seqfile> is in format <s>: no autodetection",  12 },

#ifdef HMMER_THREADS 
  { "--block_length", eslARG_INT,     NULL,   NULL,        "n>=50000", NULL,  NULL,  NULL,             "length of blocks read from target database (threaded) ",      12 },
  { "--cpu",          eslARG_INT,     NULL,   "HMMER_NCPU","n>=0",     NULL,  NULL,  CPUOPTS,          "number of parallel CPU workers to use for multithreads",      12 },
#endif

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[options] <hmmfile> <seqdb>";
static char banner[] = "search protein profile(s) against DNA sequence database";



/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads).
 */
struct cfg_s {
  char            *dbfile;            /* target sequence database file                   */
  char            *hmmfile;           /* query HMM file                                  */
};

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


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

      if (puts("\nTranslation options:")                                     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 15, 2, 80);
      
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
 
      if (puts("\nAvailable NCBI genetic code tables (for -c <id>):")        < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_gencode_DumpAltCodeTable(stdout); 

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
  if (esl_opt_IsUsed(go, "--acc")        && fprintf(ofp, "# prefer accessions over names:    yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")      && fprintf(ofp, "# show alignments in output:       no\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")    && fprintf(ofp, "# max ASCII text line length:      unlimited\n")                                             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")      && fprintf(ofp, "# max ASCII text line length:      %d\n",             esl_opt_GetInteger(go, "--textw"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notrans")    && fprintf(ofp, "# show translated DNA sequence:    no\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--vertcodon")  && fprintf(ofp, "# show DNA codon vertically:       no\n")                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
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

  if (esl_opt_IsUsed(go, "--nonull2")    && fprintf(ofp, "# null2 bias corrections:          off\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-Z")           && fprintf(ofp, "# sequence search space set to:    %.0f\n",           esl_opt_GetReal(go, "-Z"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--domZ")       && fprintf(ofp, "# domain search space set to:      %.0f\n",           esl_opt_GetReal(go, "--domZ"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--seed"))  
  {
    if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n")                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                                    fprintf(ofp, "# random number seed set to:       %d\n",        esl_opt_GetInteger(go, "--seed"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--tformat")    && fprintf(ofp, "# targ <seqfile> format asserted:  %s\n",             esl_opt_GetString(go, "--tformat"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")        && fprintf(ofp, "# number of worker threads:        %d\n",             esl_opt_GetInteger(go, "--cpu"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
#endif
  if (esl_opt_IsUsed(go, "-c")           && fprintf(ofp, "# use alt genetic code of NCBI transl table: %d\n",   esl_opt_GetInteger(go, "-c"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-l")           && fprintf(ofp, "# minimum ORF length: %d\n",                          esl_opt_GetInteger(go, "-l"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-m")           && fprintf(ofp, "# ORFs must initiate with AUG only:    yes\n")                                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-M")           && fprintf(ofp, "# ORFs must start with allowed initiation codon:    yes\n")                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--watson")     && fprintf(ofp, "# only translate top strand:    yes\n")                                                      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--crick")      && fprintf(ofp, "# only translate bottom strand:    yes\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (                                      fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
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

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.hmmfile    = NULL;
  cfg.dbfile     = NULL;

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.dbfile);    

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);

  return status;
}

/* translate_sequence()
 * For input DNA sequence, add all ORFs (6 frames) to wrk block
 */
static int
translate_sequence(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq)
{
      if (wrk->do_watson) {
        esl_gencode_ProcessStart(gcode, wrk, sq);
        esl_gencode_ProcessPiece(gcode, wrk, sq);
        esl_gencode_ProcessEnd  (wrk, sq);
      }

      if (wrk->do_crick) {
        esl_sq_ReverseComplement(sq);
        esl_gencode_ProcessStart(gcode, wrk, sq);
        esl_gencode_ProcessPiece(gcode, wrk, sq);
        esl_gencode_ProcessEnd  (wrk, sq);
        esl_sq_ReverseComplement(sq);
      }

  return eslOK;
}

/* serial_master()
 * The serial version of hmmsearcht.
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
  FILE            *ofp      = stdout;            /* results output file (-o)                               */
  FILE            *afp      = NULL;              /* alignment output file (-A)                             */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)           */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-dom (--domtblout)        */

  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                                    */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                               */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                          */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                       */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file                 */
  ESL_STOPWATCH   *w;

  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  /* used to keep track of the lengths of the sequences that are processed */
  ID_LENGTH_LIST  *id_length_list = NULL;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  /* hmmsearcht */
  P7_TOPHITS      *tophits_accumulator = NULL;      /* to hold the top hits information from all 6 frame translations     */
  P7_PIPELINE     *pipelinehits_accumulator = NULL; /* to hold the pipeline hit information from all 6 frame translations */
  ESL_ALPHABET    *abcDNA = NULL;                   /* DNA sequence alphabet                                              */
  ESL_ALPHABET    *abcAMINO = NULL;                 /* DNA sequence alphabet                                              */
  ESL_SQ          *qsqDNA = NULL;		    /* DNA query sequence                                                 */
  ESL_GENCODE     *gcode       = NULL;
  /* end hmmsearcht */

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

  /* Open the query profile HMM file */
  status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp       = fopen(esl_opt_GetString(go, "-o"), "w"))           == NULL) p7_Fail("Failed to open output file %s for writing\n",                      esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp       = fopen(esl_opt_GetString(go, "-A"), "w"))           == NULL) p7_Fail("Failed to open alignment file %s for writing\n",                   esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp     = fopen(esl_opt_GetString(go, "--tblout"),    "w"))  == NULL) esl_fatal("Failed to open tabular per-seq output file %s for writing\n",    esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp  = fopen(esl_opt_GetString(go, "--domtblout"), "w"))  == NULL) esl_fatal("Failed to open tabular per-dom output file %s for writing\n",    esl_opt_GetString(go, "--domtblout")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);

  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue     = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);


  /*the query sequence will be DNA but will be translated to amino acids */
  abcDNA   = esl_alphabet_Create(eslDNA); 
  abcAMINO = esl_alphabet_Create(eslAMINO); 
  qsqDNA   = esl_sq_CreateDigital(abcDNA);

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);

  if (hstatus == eslOK)
  {
      if (abc->type != eslAMINO) p7_Fail("hmmsearcht only supports amino acid HMMs; %s uses a different alphabet", cfg->hmmfile);

      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
      esl_sqfile_SetDigital(dbfp, abcDNA); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks

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
	  block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abcDNA);
	  if (block == NULL) 	      esl_fatal("Failed to allocate sequence block");

 	  status = esl_workqueue_Init(queue, block);
	  if (status != eslOK)	      esl_fatal("Failed to add block to work queue");
	}
#endif
  }

  /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any Amino Acid   */
  gcode = esl_gencode_Create(abcDNA, abcAMINO);
  esl_gencode_Set(gcode, esl_opt_GetInteger(go, "-c"));  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m or -M are set


  /* Outer loop: over each query HMM in <hmmfile>. */
  while (hstatus == eslOK) 
    {
      P7_PROFILE      *gm      = NULL;
      P7_OPROFILE     *om      = NULL;       /* optimized query profile                  */


      /* defining the maximum window overlap in case the target DNA sequence is very long and
       * multiple windows are taken for a single sequence
       */
      p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);
      hmm->max_length *=3;

      nquery++;
      esl_stopwatch_Start(w);

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1)
      {
        if (! esl_sqfile_IsRewindable(dbfp))
          esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);

        esl_sqfile_Position(dbfp, 0);
      }

      if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
      if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

      /* Convert to an optimized model */
      gm = p7_profile_Create (hmm->M, abc);
      om = p7_oprofile_Create(hmm->M, abc);
      p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gm, om);                        /* <om> is now p7_LOCAL, multihit */

      /* Create processing pipeline and hit list accumulators */
      tophits_accumulator      = p7_tophits_Create(); 
      pipelinehits_accumulator = p7_pipeline_Create(go, 100, 100, FALSE, p7_SEARCH_SEQS);

      pipelinehits_accumulator->nmodels       = 1;
      pipelinehits_accumulator->nnodes        = hmm->M;
      pipelinehits_accumulator->is_translated = TRUE;

      for (i = 0; i < infocnt; ++i)
      {
        /* Create processing pipeline and hit list */
        info[i].gcode = gcode;
        info[i].wrk                = esl_gencode_WorkstateCreate(go, gcode);
        info[i].th                 = p7_tophits_Create();
        info[i].om                 = p7_oprofile_Clone(om);
        info[i].pli                = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
        info[i].pli->is_translated = TRUE;
        status                     = p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
        if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

        if (  esl_opt_IsUsed(go, "--watson") )
          info[i].pli->strands = p7_STRAND_TOPONLY;
        else if (  esl_opt_IsUsed(go, "--crick") )
          info[i].pli->strands = p7_STRAND_BOTTOMONLY;
        else
          info[i].pli->strands = p7_STRAND_BOTH;

        if (  esl_opt_IsUsed(go, "--block_length") )
          info[i].pli->block_length = esl_opt_GetInteger(go, "--block_length");
        else
          info[i].pli->block_length = HMMSEARCHT_MAX_RESIDUE_COUNT;


#ifdef HMMER_THREADS
        if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

      /* establish the id_lengths data structutre */
      id_length_list = init_id_length(1000);

#ifdef HMMER_THREADS
      if (ncpus > 0)
          sstatus = thread_loop(info, id_length_list, threadObj, queue, dbfp);
      else
#endif
      sstatus = serial_loop(info, id_length_list, dbfp);

      switch(sstatus)
      {
          case eslOK:
           /* do nothing */
           break;
         case eslEOF:
           /* do nothing */
           break;
         default:
           esl_fatal("Unexpected error %d processing ORFs", sstatus);
      }

      /* merge the results of the search results */
      for (i = 0; i < infocnt; ++i)
      {
          p7_tophits_Merge(tophits_accumulator, info[i].th);
          p7_pipeline_Merge(pipelinehits_accumulator, info[i].pli);

          p7_pipeline_Destroy(info[i].pli);
          p7_tophits_Destroy(info[i].th);
          p7_oprofile_Destroy(info[i].om);
      }













      // NORD - START
      if (tophits_accumulator->N)
        GenSpliceGraphs(tophits_accumulator,dbfp,om,gcode);
      // NORD - END



















      esl_sq_Reuse(qsqDNA);
		
      /* Sort and remove duplicates */
      p7_tophits_SortBySeqidxAndAlipos(tophits_accumulator);
      assign_Lengths(tophits_accumulator, id_length_list);
      p7_tophits_RemoveDuplicates(tophits_accumulator, pipelinehits_accumulator->use_bit_cutoffs);


      /* Print the results.  */
      p7_tophits_SortBySortkey(tophits_accumulator);
      p7_tophits_Threshold(    tophits_accumulator, pipelinehits_accumulator);
      p7_tophits_Targets(ofp,  tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp,  tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));
      if (domtblfp)  p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pipelinehits_accumulator, w);
      if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      /* Output the results in an MSA (-A option) */
      if (afp) {
         ESL_MSA *msa = NULL;

         if (p7_tophits_Alignment(tophits_accumulator, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
         {
           if (textw > 0) esl_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
           else           esl_msafile_Write(afp, msa, eslMSAFILE_PFAM);
	  
           if (fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
         } 
         else { if (fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
	  
         esl_msa_Destroy(msa);
      }

	  
      p7_pipeline_Destroy(pipelinehits_accumulator);
      p7_tophits_Destroy(tophits_accumulator);
      p7_oprofile_Destroy(om);
      p7_profile_Destroy(gm);
      p7_hmm_Destroy(hmm);
      destroy_id_length(id_length_list);
      hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
    } /* end outer loop over query HMMs */

  switch(hstatus) {
  case eslEOD:       p7_Fail("read failed, HMM file %s may be truncated?", cfg->hmmfile);      break;
  case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->hmmfile);      break;
  case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->hmmfile);      break;
  case eslEOF:       /* do nothing. EOF is what we want. */                                    break;
  default:           p7_Fail("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->hmmfile);
  }


  /* Terminate outputs... any last words? */
  if (tblfp)    p7_tophits_TabularTail(tblfp,     "hmmsearcht", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp,  "hmmsearcht", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit */
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
  for (i = 0; i < infocnt; ++i) {
    p7_bg_Destroy(info[i].bg);
    esl_gencode_WorkstateDestroy(info[i].wrk);
  }

  free(info);
  
  esl_gencode_Destroy(gcode);

  esl_sq_Destroy(qsqDNA);  
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAMINO);
    
  p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}


static int
serial_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp)
{

  int      k;
  int      sstatus       = eslOK;
  int      seq_id        = 0;
  uint64_t prev_char_cnt = 0;
  ESL_ALPHABET *abc         = esl_alphabet_Create(eslDNA);
  ESL_SQ       *dbsq_dna    = esl_sq_CreateDigital(abc);   /* (digital) nucleotide sequence, to be translated into ORFs  */
  ESL_SQ_BLOCK *block       = NULL;                        /* for translated ORFs                                        */
  ESL_SQ       *dbsq_aa     = NULL;                        /* used to hold a current ORF                                 */
  ESL_SQ       *dbsq_dnatxt = esl_sq_Create();
  
  dbsq_dnatxt->abc = abc;

  sstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq_dna);

  info->wrk->orf_block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, info->om->abc);
  if (info->wrk->orf_block == NULL)          esl_fatal("Failed to allocate sequence block");
  
  while (sstatus == eslOK  ) {
      dbsq_dna->idx = seq_id;

      /* copy and convert the DNA sequence to text so we can print it in the domain alignment display */
      esl_sq_Copy(dbsq_dna, dbsq_dnatxt);
      info->ntqsq = dbsq_dnatxt; // for printing the DNA target sequence in the domain hits display
      
      /* translate DNA sequence to 6 frame ORFs */
      dbsq_dna->L = dbsq_dna->n; /* here, L is not the full length of the sequence in the db, just of the currently-active window;  required for esl_gencode machinations */
      esl_sq_ReuseBlock(info->wrk->orf_block);
      translate_sequence(info->gcode, info->wrk, dbsq_dna);

      block =  info->wrk->orf_block;

      /* Main loop: */
      for (k = 0; k < block->count; ++k)
      {
          dbsq_aa = &(block->list[k]);
      
          if (   (dbsq_aa->start < dbsq_aa->end    &&  dbsq_aa->end < dbsq_dna->C )  ||
                 (dbsq_aa->end < dbsq_aa->start    &&  dbsq_aa->start < dbsq_dna->C ) )
              continue; /* don't bother with an orf that showed up completely within a previous window */


          p7_pli_NewSeq(info->pli, dbsq_aa);
          /*we use overlapping windows to ensure that we don't miss something at a boundary, but now
           * we need to adjust for overcounting candidate ORFs
           */
          if (dbsq_aa->start < dbsq_aa->end) {
              if (dbsq_aa->start < dbsq_dna->C - 59 ){
                  info->pli->nseqs--;
                  info->pli->nres -= ESL_MIN(dbsq_aa->n,  (dbsq_dna->C - dbsq_aa->start + 1)/3);
              }
          } else {
              if (dbsq_aa->end < dbsq_dna->C - 58 ){
                  info->pli->nseqs--;
                  info->pli->nres -= ESL_MIN(dbsq_aa->n,  (dbsq_dna->C - dbsq_aa->end + 1)/3);
              }
          }

          /* Use the name, accession, and description from the DNA sequence and not 
           * from the ORF which is generated by gencode and only for internal use.
           * Set the orfid to a number that will be uniq and consistent across 
           * threading options */
          dbsq_aa->idx = prev_char_cnt+k;
          snprintf(dbsq_aa->orfid, dbsq_aa->orfalloc, "orf%" PRId64 "", dbsq_aa->idx);

          if ((sstatus = esl_sq_SetName     (dbsq_aa, info->ntqsq->name))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence name failed");
          if ((sstatus = esl_sq_SetAccession(dbsq_aa, info->ntqsq->acc))    != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence accession failed");
          if ((sstatus = esl_sq_SetDesc     (dbsq_aa, info->ntqsq->desc))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence description failed");
          
          dbsq_aa->idx = dbsq_dna->idx;

          p7_bg_SetLength(info->bg, dbsq_aa->n);
          p7_oprofile_ReconfigLength(info->om, dbsq_aa->n);

          p7_Pipeline(info->pli, info->om, info->bg, dbsq_aa, info->ntqsq, info->th, NULL);

          esl_sq_Reuse(dbsq_aa);
          p7_pipeline_Reuse(info->pli);
      }

      prev_char_cnt += dbsq_dna->n;

      /* "maxlength * 3" because we're looking for a protein of maxlength, and this is DNA */
      sstatus = esl_sqio_ReadWindow(dbfp, info->om->max_length * 3, info->pli->block_length, dbsq_dna);

      if (sstatus == eslEOD) { // no more left of this sequence ... move along to the next sequence.
          add_id_length(id_length_list, dbsq_dna->idx, dbsq_dna->L); 
          seq_id++;
          esl_sq_Reuse(dbsq_dna);
          sstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq_dna);
      }

  }
  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(dbsq_dna);
  esl_sq_Destroy(dbsq_dnatxt);

  return sstatus;
}

#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{

  int      i;
  int      status        = eslOK;
  int      sstatus       = eslOK;
  int      eofCount      = 0;
  int      seqid         = 0;
  uint64_t prev_char_cnt = 0;
  ESL_SQ_BLOCK *block; // block of 1 or more nucleotide sequences
  void         *newBlock;
  ESL_ALPHABET *abcdna    = esl_alphabet_Create(eslDNA);
  ESL_SQ       *tmpsq_dna = esl_sq_CreateDigital(abcdna ) ;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
  ((ESL_SQ_BLOCK *)newBlock)->complete = TRUE;

  /* Main loop: */
  while (sstatus == eslOK )
    {
      block = (ESL_SQ_BLOCK *) newBlock;
      sstatus = esl_sqio_ReadBlock(dbfp, block, info->pli->block_length, -1,  /*max_init_window=*/FALSE, TRUE);

      block->first_seqidx = info->pli->ndbseqs;
      seqid = block->first_seqidx;

      for (i=0; i<block->count; i++) {
        block->list[i].idx = seqid;
        add_id_length(id_length_list, seqid, block->list[i].L);
        seqid++;

        block->list[i].prev_n = prev_char_cnt;
        prev_char_cnt += block->list[i].n;
      }
       
      info->pli->ndbseqs += block->count  - (block->complete ? 0 : 1);
        
      if (sstatus == eslEOF) {
          if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
          ++eofCount;
      } else if (!block->complete ) {
          /* The final sequence on the block was an incomplete window of the active sequence,
           * so our next read will need a copy of it to correctly deal with overlapping
           * regions. We capture a copy of the sequence here before sending it off to the
           * pipeline to avoid odd race conditions that can occur otherwise.
           * Copying the entire sequence isn't really necessary, and is a bit heavy-
           * handed. Could accelerate if this proves to have any notable impact on speed. */
          esl_sq_Copy(block->list + (block->count - 1) , tmpsq_dna);
      }


      if (sstatus == eslOK) {

          status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
          if (status != eslOK) esl_fatal("Work queue reader failed");

          /* newBlock needs all this information so the next ReadBlock call will know what to do */
          ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
          if (!block->complete) {
              /* Push the captured copy of the previously-read sequence into the new block, 
               * in preparation for ReadWindow  (double copy ... slower than necessary) */
              esl_sq_Copy(tmpsq_dna, ((ESL_SQ_BLOCK *)newBlock)->list);
              esl_sq_Reuse(tmpsq_dna);

              if (  ((ESL_SQ_BLOCK *)newBlock)->list->n < info->om->max_length ) {
                /* no reason to search the final partial sequence on the block, as the next block will search this whole chunk */
                ((ESL_SQ_BLOCK *)newBlock)->list->C = ((ESL_SQ_BLOCK *)newBlock)->list->n;
                (((ESL_SQ_BLOCK *)newBlock)->count)--;
              } else {
                /* This sets the overlap that the next window will contain, relative to the current
                 * window.  "maxlength * 3" because we're looking for a protein of maxlength, and this is DNA */
                ((ESL_SQ_BLOCK *)newBlock)->list->C = info->om->max_length * 3;
              }
          }
      }
    }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
    {
      /* wait for all the threads to complete */

      esl_alphabet_Destroy(abcdna);
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
      esl_sq_Destroy(tmpsq_dna);
    }

  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int i, k;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  ESL_SQ_BLOCK  *orfblock    = NULL;
  void          *newBlock;
  ESL_SQ_BLOCK  *dnablock    = NULL;
  ESL_SQ        *dbsq_aa     = NULL;   /* used to hold a current ORF  */
  ESL_SQ        *dbsq_dna    = NULL;
  ESL_SQ        *dbsq_dnatxt = esl_sq_Create();

  impl_Init();

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  info->wrk->orf_block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, info->om->abc);
  orfblock = info->wrk->orf_block;
  if (orfblock == NULL)          esl_fatal("Failed to allocate sequence block");

  /* thread loops until all blocks have been processed */
  dnablock = (ESL_SQ_BLOCK *) newBlock; //block from threads

  while (dnablock->count > 0)
    {
      /* Main loop: */
      for (i = 0; i < dnablock->count; ++i)
      {

          dbsq_dna = dnablock->list + i;

          /* copy and convert the DNA sequence to text so we can print it in the domain alignment display */
          esl_sq_Reuse(dbsq_dnatxt);
          esl_sq_Copy(dbsq_dna, dbsq_dnatxt);
          dbsq_dnatxt->abc = dbsq_dna->abc;
          info->ntqsq = dbsq_dnatxt; // for printing the DNA target sequence in the domain hits display

          /* translate DNA sequence to 6 frame ORFs */
          dbsq_dna->L = dbsq_dna->n; /* here, L is not the full length of the sequence in the db, just of the currently-active window;  required for esl_gencode machinations */

          esl_sq_ReuseBlock(info->wrk->orf_block);
          translate_sequence(info->gcode, info->wrk, dbsq_dna);

          orfblock =  info->wrk->orf_block;
          
          /* Main loop: */
          for (k = 0; k < orfblock->count; ++k)
          {
              dbsq_aa = &(orfblock->list[k]);

              if (   (dbsq_aa->start < dbsq_aa->end    &&  dbsq_aa->end < dbsq_dna->C )  ||
                     (dbsq_aa->end < dbsq_aa->start    &&  dbsq_aa->start < dbsq_dna->C ) )
                  continue; /* don't bother with an orf that showed up completely within a previous window */


              /* Use the name, accession, and description from the DNA sequence and not 
               * from the ORF which is generated by gencode and only for internal use. 
               * Set the orfid to a number that will be uniq and consistent across 
               * threading options */
              dbsq_aa->idx = dbsq_dna->prev_n+k;
              snprintf(dbsq_aa->orfid, dbsq_aa->orfalloc, "orf%" PRId64 "", dbsq_aa->idx);
              if ((status = esl_sq_SetName     (dbsq_aa, info->ntqsq->name))   != eslOK)  esl_fatal("Set query sequence name failed");
              if ((status = esl_sq_SetAccession(dbsq_aa, info->ntqsq->acc))    != eslOK)  esl_fatal("Set query sequence accession failed");
              if ((status = esl_sq_SetDesc     (dbsq_aa, info->ntqsq->desc))   != eslOK)  esl_fatal("Set query sequence description failed");
              dbsq_aa->idx = dbsq_dna->idx;

              p7_pli_NewSeq(info->pli, dbsq_aa);
              /*we use overlapping windows to ensure that we don't miss something at a boundary, but now
               * we need to adjust for overcounting candidate ORFs
               */
              if (dbsq_aa->start < dbsq_aa->end) {
                  if (dbsq_aa->start < dbsq_dna->C - 59 ){
                      info->pli->nseqs--;
                      info->pli->nres -= ESL_MIN(dbsq_aa->n,  (dbsq_dna->C - dbsq_aa->start + 1)/3);
                  }
              } else {
                  if (dbsq_aa->end < dbsq_dna->C - 58 ){
                      info->pli->nseqs--;
                      info->pli->nres -= ESL_MIN(dbsq_aa->n,  (dbsq_dna->C - dbsq_aa->end + 1)/3);
                  }
              }

              p7_bg_SetLength(info->bg, dbsq_aa->n);
              p7_oprofile_ReconfigLength(info->om, dbsq_aa->n);

              p7_Pipeline(info->pli, info->om, info->bg, dbsq_aa, info->ntqsq, info->th, NULL);

              esl_sq_Reuse(dbsq_aa);
              p7_pipeline_Reuse(info->pli);
          }

      }

      status = esl_workqueue_WorkerUpdate(info->queue, dnablock, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      dnablock = (ESL_SQ_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, dnablock, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_sq_Destroy(dbsq_dnatxt);
  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* HMMER_THREADS */
 

/* helper functions for tracking id_lengths */

static ID_LENGTH_LIST *
init_id_length( int size )
{
  int status;
  ID_LENGTH_LIST *list;

  ESL_ALLOC (list, sizeof(ID_LENGTH_LIST));
  list->count = 0;
  list->size  = size;
  list->id_lengths = NULL;

  ESL_ALLOC (list->id_lengths, size * sizeof(ID_LENGTH));

  return list;

ERROR:
  return NULL;
}

static void
destroy_id_length( ID_LENGTH_LIST *list )
{

  if (list != NULL) {
    if (list->id_lengths != NULL) free (list->id_lengths);
    free (list);
  }

}

static int
add_id_length(ID_LENGTH_LIST *list, int id, int L)
{
   int status;

   if (list->count > 0 && list->id_lengths[list->count-1].id == id) {
     // the last time this gets updated, it'll have the sequence's actual length
     list->id_lengths[list->count-1].length = L;
   } else {

     if (list->count == list->size) {
       list->size *= 10;
       ESL_REALLOC(list->id_lengths, list->size * sizeof(ID_LENGTH));
     }

     list->id_lengths[list->count].id     = id;
     list->id_lengths[list->count].length = L;

     list->count++;
   }
   return eslOK;

ERROR:
   return status;
}

static int
assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list) {

  int i,d;
  int j = 0;
  for (i=0; i<th->N; i++) {
    while (th->hit[i]->seqidx != id_length_list->id_lengths[j].id) { j++; }
    for(d=0; d<th->hit[i]->ndom; d++) {
      th->hit[i]->dcl[d].ad->L = id_length_list->id_lengths[j].length;
    } 
  }
  return eslOK;
}


