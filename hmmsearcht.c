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















/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                           BEGIN SPLICING STUFF
 * 
 *  Catalogue of Ships
 *  ==================
 *
 *  Fleet 1: Pertaining to Data Preparation
 *
 *  + FloatHighLowSortIndex : 
 *  + FloatLowHighSortIndex :
 *  + GetMinAndMaxCoords    :
 *  + GetTargetNuclSeq      :
 *  + GrabNuclRange         :
 *  + DetermineNuclType     :
 *  
 *
 *
 *  Fleet 2: Pertaining to Splicing Hits
 *
 *  + GetSpliceOptions             :
 *  + FindOptimalSpliceSite        :
 *  + SpliceOverlappingDomains     :
 *  + GetNuclRangesFromAminoCoords :
 *  + SketchSpliceEdge             :
 *  + HitsAreSpliceComaptible      :
 *  + GatherViableSpliceEdges      :
 *
 *
 *
 *  Fleet 3: Pertaining to the Splice Graph
 *
 *  + InitSpliceNode        :
 *  + ConnectNodesByEdge    :
 *  + FindBestPathToNode    :
 *  + EvangelizePath        :
 *  + GatherNTermNodes      :
 *  + GatherCTermNodes      :
 *  + GenCumScoreSort       :
 *  + EvaluatePaths         :
 *  + FillOutGraphStructure :
 *  + FindBestFullPath      :
 *  + BuildSpliceGraph      :
 *
 *
 *
 *  Fleet 4: Pertaining to Filling in Gaps (Sub-Model Search)
 *
 *  + ExtractSubProfile       :
 *  + NodesAreDCCCompatible   :
 *  + GetBoundedSearchRegions :
 *  + SelectFinalSubHits      :
 *  + FindSubHits             :
 *  + IntegrateMissedHits     :
 *  + SeekMissingExons        :
 *  + AddMissingExonsToGraph  :
 *
 *
 *
 *  Fleet 5: Pertaining to the Final Alignment
 *
 *  + GetExonSetFromStartNode :
 *  + FindComponentBestStart  :
 *  + TranslateExonSetNucls   :
 *  + GrabExonCoordSetNucls   :
 *  + GetSplicedExonCoordSets :
 *  + ReportSplicedTopHits    :
 *  + RunModelOnExonSets      :
 *
 *
 *
 *  Flagship: SpliceHits
 *
 */



// Before we get to the fun stuff, let's just set up some
// bureaucratic stuff to make debugging relatively (hopefully)
// painless


static int DEBUGGING = 1; // Print debugging output?
int FUNCTION_DEPTH = 0;
void DEBUG_OUT (const char * message, const int func_depth_change) {

  if (func_depth_change > 0) 
    FUNCTION_DEPTH += func_depth_change;
  
  fprintf(stderr,"  SplDebug:");
  for (int i=0; i<FUNCTION_DEPTH; i++) 
    fprintf(stderr,"  ");
  fprintf(stderr,"%s\n",message);
  fflush(stderr);
  
  if (func_depth_change < 0) 
    FUNCTION_DEPTH += func_depth_change;

}


static float SSSCORE[2]      = {-0.7,0.7}; // Non-canon vs canon splice site
static float EDGE_FAIL_SCORE = -14773.0;   // Makes me thirsty for a latte!

static char AMINO_CHARS[21] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};
static char   DNA_CHARS[ 5] = {'A','C','G','T','-'};
static char   RNA_CHARS[ 5] = {'A','C','G','U','-'};

// How many amino acids are we willing to extend to bridge two hits?  
// How many overlapping aminos do we require to perform bridging?
static int MAX_AMINO_EXT     = 4;
static int MIN_AMINO_OVERLAP = 6;


int intMax (int a, int b) { if (a>b) return a; return b; }
int intMin (int a, int b) { if (a<b) return a; return b; }




//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: TARGET_SEQ
//
//  Desc. : 
//
typedef struct _target_seq {
  
  ESL_SQ  * esl_sq; // For freeing, we'll need this pointer
  ESL_DSQ * Seq;
  const ESL_ALPHABET * abc;

  int64_t start;
  int64_t end;

} TARGET_SEQ;




//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: DOMAIN_OVERLAP
//
//  Desc. :
//
typedef struct _domain_overlap {

  const ESL_ALPHABET * ntalpha;


  int amino_start;
  int amino_end;


  int upstream_hit_id;
  int upstream_dom_id;
  int upstream_nucl_start;
  int upstream_nucl_end;
  int upstream_ext_len;
  int upstream_disp_start;

  P7_TOPHITS    * UpstreamTopHits;
  P7_ALIDISPLAY * UpstreamDisplay;
  ESL_DSQ       * UpstreamNucls;
  

  int downstream_hit_id;
  int downstream_dom_id;
  int downstream_nucl_start;
  int downstream_nucl_end;
  int downstream_ext_len;
  int downstream_disp_end;

  P7_TOPHITS    * DownstreamTopHits;
  P7_ALIDISPLAY * DownstreamDisplay;
  ESL_DSQ       * DownstreamNucls;


  int   upstream_exon_terminus;
  int downstream_exon_terminus;

  int   upstream_spliced_nucl_end;
  int downstream_spliced_nucl_start;


  float score_density;
  float score;

} DOMAIN_OVERLAP;







//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: SPLICE_NODE
//
//  Desc. :
//
typedef struct _splice_node {


  int node_id;


  int was_missed;
  int hit_id;
  int dom_id;


  int in_edge_cap;
  int num_in_edges;
  int best_in_edge; // Relative to the following arrays
  DOMAIN_OVERLAP ** InEdges;
  struct _splice_node ** UpstreamNodes;


  int out_edge_cap;
  int num_out_edges;
  int best_out_edge; // Relative to the following arrays
  DOMAIN_OVERLAP ** OutEdges;
  struct _splice_node ** DownstreamNodes;


  int is_n_terminal;
  int is_c_terminal;


  float hit_score;
  float cumulative_score;
  float best_path_score;


} SPLICE_NODE;







//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: SPLICE_GRAPH
//
//  Desc. :
//
typedef struct _splice_graph {


  P7_TOPHITS  *    TopHits;
  P7_TOPHITS  * MissedHits; // This is a fake P7_TOPHITS!  Don't be fooled!

  P7_PROFILE  *  Model;
  P7_OPROFILE * OModel;

  int revcomp;


  // Because we build based on the DOMAIN_OVERLAP datastructure,
  // it's helpful to be able to find each node ID by way of
  // these lookups.
  int ** TH_HitDomToNodeID; //    TopHits
  int  * MH_HitToNodeID;    // MissedHits


  int num_nodes;
  int num_edges;
  SPLICE_NODE ** Nodes; // NOTE: these go from [1..num_nodes]

  int * CumScoreSort; // Yes, it says cum.  FOR CUMULATIVE!

  int   num_n_term;
  int * NTermNodeIDs; // Sorted by cumulative score

  int   num_c_term;
  int * CTermNodeIDs; // Also sorted by cum score. ~CUM~


  int   has_full_path;
  int   best_full_path_length;
  int   best_full_path_start;
  int   best_full_path_end;
  float best_full_path_score;


} SPLICE_GRAPH;












/* DEBUGGING FUNCTION: DumpNode  */
void DumpNode (SPLICE_NODE * Node)
{
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |       NODE %d\n",Node->node_id);
  fprintf(stderr,"   |     ,---------------------------------------------------\n");
  fprintf(stderr,"   |     |\n");
  fprintf(stderr,"   |     |  Source in P7_TOPHITS  :  Hit %d, Domain %d\n", Node->hit_id, Node->dom_id);
  fprintf(stderr,"   |     |\n");
  fprintf(stderr,"   |     |  Score of Source Hit      :  %f\n",Node->hit_score);
  fprintf(stderr,"   |     |  Score of Path Up To Node :  %f\n",Node->cumulative_score);
  fprintf(stderr,"   |     |  Total Score of Best Path\n");
  fprintf(stderr,"   |     |      that Uses this Node  :  %f\n",Node->best_path_score);
  fprintf(stderr,"   |     |\n");
  if (Node->num_in_edges > 0) {
    fprintf(stderr,"   |     |  Number of Incoming Edges :  %d\n",Node->num_in_edges);
    fprintf(stderr,"   |     |  Best Upstream Node       :  Node %d\n",Node->UpstreamNodes[Node->best_in_edge]->node_id);
  } else if (Node->is_n_terminal) {
    fprintf(stderr,"   |     |  * N-TERMINAL NODE\n");
  } else {
    fprintf(stderr,"   |     |  - No Incoming Edges\n");
  }
  fprintf(stderr,"   |     |\n");
  if (Node->num_out_edges > 0) {
    fprintf(stderr,"   |     |  Number of Outgoing Edges :  %d\n",Node->num_out_edges);
    fprintf(stderr,"   |     |  Best Downstream Node     :  Node %d\n",Node->DownstreamNodes[Node->best_out_edge]->node_id);
  } else if (Node->is_c_terminal) {
    fprintf(stderr,"   |     |  * C-TERMINAL NODE\n");
  } else {
    fprintf(stderr,"   |     |  - No Outgoing Edges\n");
  }
  fprintf(stderr,"   |     |\n");
  fprintf(stderr,"   |     '\n");
  fprintf(stderr,"   |\n");
  fflush(stderr);
}
/* DEBUGGING FUNCTION: DumpGraph */
void DumpGraph(SPLICE_GRAPH * Graph)
{
  fprintf(stderr,"\n\n");
  fprintf(stderr,"     SPLICE GRAPH\n");
  fprintf(stderr,"   +=========================================================+\n");
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |  Total Number of Nodes      : %d\n",Graph->num_nodes);
  fprintf(stderr,"   |  Total Number of Edges      : %d\n",Graph->num_edges);
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |  Number of N-terminal nodes : %d\n",Graph->num_n_term);
  fprintf(stderr,"   |  Number of C-terminal nodes : %d\n",Graph->num_c_term);
  fprintf(stderr,"   |\n");
  if (Graph->has_full_path)
    fprintf(stderr,"   |  * This graph has at least one full path through the HMM!\n");
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |\n");
  for (int i=1; i<=Graph->num_nodes; i++) DumpNode(Graph->Nodes[i]);
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   +=========================================================+\n");
  fprintf(stderr,"\n\n\n");
  fflush(stderr);
}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  TARGET_SEQ_Destroy
 *
 */
void TARGET_SEQ_Destroy
(TARGET_SEQ * TS)
{
  esl_sq_Destroy(TS->esl_sq); // This takes care of 'TS->Seq' too
  free(TS);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  DOMAIN_OVERLAP_Destroy
 *
 */
void DOMAIN_OVERLAP_Destroy
(DOMAIN_OVERLAP * DO)
{
  free(DO->UpstreamNucls);
  free(DO->DownstreamNucls);
  free(DO);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  SPLICE_NODE_Destroy
 *
 */
void SPLICE_NODE_Destroy
(SPLICE_NODE * Node)
{

  for (int in_edge_id = 0; in_edge_id < Node->num_in_edges; in_edge_id++) {
    if (Node->InEdges[in_edge_id]) {
      DOMAIN_OVERLAP_Destroy(Node->InEdges[in_edge_id]);
      Node->InEdges[in_edge_id] = NULL;
    }
  }

  for (int out_edge_id = 0; out_edge_id < Node->num_out_edges; out_edge_id++) {
    if (Node->OutEdges[out_edge_id]) {
      DOMAIN_OVERLAP_Destroy(Node->OutEdges[out_edge_id]);
      Node->OutEdges[out_edge_id] = NULL;
    }
  }

  free(Node->InEdges);
  free(Node->OutEdges);
  free(Node);
  Node = NULL;

}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  SPLICE_GRAPH_Destroy
 *
 */
void SPLICE_GRAPH_Destroy
(SPLICE_GRAPH * Graph)
{

  
  // Wipe them nodes *OUT*
  for (int node_id = 1; node_id <= Graph->num_nodes; node_id++)
    SPLICE_NODE_Destroy(Graph->Nodes[node_id]);


  free(Graph->NTermNodeIDs);
  free(Graph->CTermNodeIDs);

  Graph = NULL;

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FloatHighLowSortIndex
 *
 *  Desc. : Generate a sort index for an array of floats (in decreasing value)
 *          using a basic mergesort implementation.
 *
 *  Inputs:  1. Vals     : The floats that we want to produce a sorting for.
 *           2. num_vals : The number of values in the array.
 *
 *  Output:  An array of indices corresponding to a sorting of 'Vals'.
 *
 */
int * FloatHighLowSortIndex
(float * Vals, int num_vals)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'FloatHighLowSortIndex'",1);

  int * Write = malloc(num_vals * sizeof(int));
  int * Read  = malloc(num_vals * sizeof(int));
  int * Tmp;

  for (int i=0; i<num_vals; i++) {
    Write[i] = i;
    Read[i]  = i;
  }

  // Of course I'm doing a merge sort
  int writer,left_reader,left_break,right_reader,right_break;
  int ms_block_size = 1;
  while (ms_block_size < num_vals) {

    writer = 0;
    while (writer+ms_block_size < num_vals) {

      left_reader = writer;
      left_break  = writer + ms_block_size;

      right_reader = left_break;
      right_break  = right_reader + ms_block_size;
      if (right_break > num_vals)
        right_break = num_vals;

      while (left_reader < left_break && right_reader < right_break) {
        if (Vals[Read[left_reader]] > Vals[Read[right_reader]]) {
          Write[writer] = Read[left_reader];
          writer++;
          left_reader++;
        }
        else { 
          Write[writer] = Read[right_reader];
          writer++;
          right_reader++;
        }
      }
      while ( left_reader <  left_break) {
        Write[writer] = Read[left_reader];
        writer++;
        left_reader++;
      }
      while (right_reader < right_break) {
        Write[writer] = Read[right_reader];
        writer++;
        right_reader++;
      }

    }

    while (writer < num_vals) {
      Write[writer] = Read[writer];
      writer++;
    }

    // Flip 'em
    Tmp = Read;
    Read = Write;
    Write = Tmp;

    // Again! Again!
    ms_block_size *= 2;

  }

  if (DEBUGGING) DEBUG_OUT("'FloatHighLowSortIndex' Complete",-1);

  free(Write);
  return Read;

}
/*  Function: FloatLowHighSortIndex
 *
 *  Desc. : Invert the sort index produced by 'FloatHighLowSortIndex'
 *
 */
int * FloatLowHighSortIndex
(float * Vals, int num_vals)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'FloatLowHighSortIndex'",1);
  int * SortIndex = FloatHighLowSortIndex(Vals,num_vals);
  for (int i=0; i<num_vals/2; i++) {
    int tmp = SortIndex[i];
    SortIndex[i] = SortIndex[(num_vals-1)-i];
    SortIndex[(num_vals-1)-i] = tmp;
  }
  if (DEBUGGING) DEBUG_OUT("'FloatLowHighSortIndex' Complete",-1);
  return SortIndex;
}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetMinAndMaxCoords
 *
 *  Desc. :  Iterate over a collection of hits (from non-spliced hmmsearcht)
 *           and identify the minimum and maximum coordinates of hits to the
 *           genomic sequence.
 *
 *  Inputs:  1. TopHits       :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *
 *  Output:
 *
 */
void GetMinAndMaxCoords
(P7_TOPHITS * TopHits, TARGET_SEQ * TargetNuclSeq)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetMinAndMaxCoords'",1);

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

    for (int hit_id = 0; hit_id < (int)(TopHits->N); hit_id++) {

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

    for (int hit_id = 0; hit_id < (int)(TopHits->N); hit_id++) {

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

  TargetNuclSeq->start = min;
  TargetNuclSeq->end   = max;

  if (DEBUGGING) DEBUG_OUT("'GetMinAndMaxCoords' Complete",-1);

}













/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetTargetNuclSeq
 *
 *  Desc. :  Given a sequence file and a collection of translated search hits, extract a
 *           sub-region of the genome that contains the all of the nucleotides implicated in
 *           the set of hits.
 *
 *           This assumes (correctly!) that all of the hits are to the same nucleotide sequence.
 *
 *  Inputs:  1. GenomicSeqFile : An ESL_SQFILE struct representing the input genom(e/ic sequence) file.
 *           2.        TopHits : A collection of hits achieved by "standard" hmmsearcht (unspliced)
 *
 *  Output:  A TARGET_SEQ struct, containing the nucleotide subsequence within which all of the
 *           unspliced hmmsearcht hits reside.
 *
 */
TARGET_SEQ * GetTargetNuclSeq
(ESL_SQFILE * GenomicSeqFile, P7_TOPHITS * TopHits)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetTargetNuclSeq'",1);

  TARGET_SEQ * TargetNuclSeq = (TARGET_SEQ *)malloc(sizeof(TARGET_SEQ));
  GetMinAndMaxCoords(TopHits,TargetNuclSeq);

  TargetNuclSeq->abc = GenomicSeqFile->abc;

  ESL_SQFILE * TmpSeqFile;
  esl_sqfile_Open(GenomicSeqFile->filename,GenomicSeqFile->format,NULL,&TmpSeqFile);
  esl_sqfile_OpenSSI(TmpSeqFile,NULL);

  TargetNuclSeq->esl_sq = esl_sq_CreateDigital(TargetNuclSeq->abc);
  int fetch_err_code    = esl_sqio_FetchSubseq(TmpSeqFile,TopHits->hit[0]->name,TargetNuclSeq->start,TargetNuclSeq->end,TargetNuclSeq->esl_sq);

  esl_sqfile_Close(TmpSeqFile);

  if (fetch_err_code != eslOK) {
    fprintf(stderr,"\n  ERROR: Failed to fetch target subsequence (is there an .ssi index for the sequence file?)\n\n");
    exit(1);
  }

  TargetNuclSeq->Seq = TargetNuclSeq->esl_sq->dsq;

  if (DEBUGGING) DEBUG_OUT("'GetTargetNuclSeq' Complete",-1);

  return TargetNuclSeq;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GrabNuclRange
 *
 *  Desc. :  Extract a sub-sequence from our (already a sub-)sequence stored
 *           in a TARGET_SEQ datastructure.
 *
 *           Importantly, the 'start' and 'end' coordinate values are with respect
 *           to the full sequence of which TargetNuclSeq is a part.  This is because
 *           we're expecting to be working with the coordinates listed in an ALIDISPLAY,
 *           which are w.r.t. the full input sequence.
 *
 *  Inputs:  1. TargetNuclSeq : The nucleotide sub-sequence that we want to
 *                              extract a (further sub-)sequence from.
 *           2.         start : The (inclusive) start coordinate of the sub-sequence we
 *                              want to extract.
 *           3.           end : The (inclusive) end coordinate of the sub-sequence we
 *                              want to extract.
 *
 *  Output:  An ESL_DSQ containing the numeric nucleotide codes for the requested
 *           sub-sequence.
 *
 *           TODO: Add catches for out-of-bounds requests (these *shouldn't* happen)
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








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DetermineNuclType
 *
 *  Desc. :  Given a P7_ALIDISPLAY, determine the alphabet type of the nucleotide sequence.
 *
 *  Inputs:  1. AliDisplay : A P7_ALIDISPLAY object, assumed to be from hmmsearcht.
 *
 *  Output:  The easel code for the alphabet type of the nucleotide sequence (DNA or RNA).
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








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetSpliceOptions
 *
 *  Desc. :  Once we know the codon that would be a candidate for being split by a splice site,
 *           check each of the possible ways to split that codon (or assign it fully to one side)
 *           for (1.) the amino acid encoded by each candidate codon triple, and (2.) the
 *           strength of the splice signal for each candidate splicing (currently, just a boolean
 *           communicating whether the canonical AG/GT is being used).
 *
 *  Inputs:  1.        Overlap :
 *           2.  TargetNuclSeq :
 *           3.    upstream_ss :
 *           4.  downstream_ss :
 *           5.   SpliceCodons :
 *           6.    Canon5Prime :
 *           7.    Canon3Prime :
 *           8.          gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:  Nothing is returned, but the SpliceCodons, Canon5Prime, and Canon3Prime arrays are
 *           filled in with values representing each of the precise splice site options.
 *
 */
void GetSpliceOptions
(
  DOMAIN_OVERLAP * Overlap,
  TARGET_SEQ     * TargetNuclSeq,
  int   upstream_ss,
  int   downstream_ss,
  int * SpliceCodons,
  int * Canon5Prime,
  int * Canon3Prime,
  ESL_GENCODE * gcode
)
{
  
  if (DEBUGGING) DEBUG_OUT("Starting 'GetSpliceOptions'",1);


  // This is often too much unnecessary junk even for debugging...
  if (DEBUGGING && 0) {
    fprintf(stderr,"\n");
    fprintf(stderr,"  > GSO Dom1: %d / %d\n",(int)(Overlap->upstream_hit_id),(int)(Overlap->upstream_dom_id));
    fprintf(stderr,"            : %d..%d\n",(int)(Overlap->upstream_nucl_start),(int)(Overlap->upstream_nucl_end));
    fprintf(stderr,"  > GSO Dom2: %d / %d\n",(int)(Overlap->downstream_hit_id),(int)(Overlap->downstream_dom_id));
    fprintf(stderr,"            : %d..%d\n",(int)(Overlap->downstream_nucl_start),(int)(Overlap->downstream_nucl_end));
    fprintf(stderr,"    splice\n");
    fprintf(stderr,"    coord.s : %d , %d\n",upstream_ss,downstream_ss);
    fprintf(stderr,"\n");
    fflush(stderr);
  }



  // We need to convert the indices 'upstream_ss' and 'downstream_ss'
  // to actual nucleotide coordinates
  int us_option_range_start = Overlap->upstream_nucl_start;
  int us_option_range_end;
  int ds_option_range_start;
  int ds_option_range_end   = Overlap->downstream_nucl_start;
  if (Overlap->upstream_nucl_start < Overlap->upstream_nucl_end) {

    // Forward strand

    us_option_range_start += upstream_ss + 1; // Gets us into the "contested" zone
    us_option_range_end    = us_option_range_start + 4;

    ds_option_range_end   += downstream_ss - 1;
    ds_option_range_start  = ds_option_range_end - 4;

  } else {

    // Revcomp

    us_option_range_start -= upstream_ss + 1;
    us_option_range_end    = us_option_range_start - 4;

    ds_option_range_end   -= downstream_ss - 1;
    ds_option_range_start  = ds_option_range_end + 4;

  }



  // Grab them ranges!
  ESL_DSQ * UN = GrabNuclRange(TargetNuclSeq,us_option_range_start,us_option_range_end);
  ESL_DSQ * DN = GrabNuclRange(TargetNuclSeq,ds_option_range_start,ds_option_range_end);



  for (int i=0; i<4; i++) {
    Canon5Prime[i] = 0;
    Canon3Prime[i] = 0;
  }



  ESL_DSQ * Codon;
  esl_abc_CreateDsq(Overlap->ntalpha,"AAA",&Codon);


  // NOTE: Even though there's an intuitive loop structure to this,
  //       it doesn't kill us to unroll it, so I've unrolled it.


  // Option 1: |ABCxx|...yy|
  Codon[1] = UN[1];
  Codon[2] = UN[2];
  Codon[3] = UN[3];
  SpliceCodons[0] = esl_gencode_GetTranslation(gcode,&Codon[1]);

  if (SpliceCodons[0] >= 0 && SpliceCodons[0] <= 20) {
    if (UN[4] == 2 && UN[5] == 3) Canon5Prime[0] = 1;
    if (DN[4] == 0 && DN[5] == 2) Canon3Prime[0] = 1;
  } else {
    SpliceCodons[0] = -1; // Stop codon
  }


  // Option 2: |ABxx.|..yyD|
  Codon[3] = DN[5];  
  SpliceCodons[1] = esl_gencode_GetTranslation(gcode,&Codon[1]);

  if (SpliceCodons[1] >= 0 && SpliceCodons[1] <= 20) {
    if (UN[3] == 2 && UN[4] == 3) Canon5Prime[1] = 1;
    if (DN[3] == 0 && DN[4] == 2) Canon3Prime[1] = 1;
  } else {
    SpliceCodons[1] = -1; // Stop codon
  }


  // Option 3: |Axx..|.yyDE|
  Codon[2] = DN[4];
  SpliceCodons[2] = esl_gencode_GetTranslation(gcode,&Codon[1]);

  if (SpliceCodons[2] >= 0 && SpliceCodons[2] <= 20) {
    if (UN[2] == 2 && UN[3] == 3) Canon5Prime[2] = 1;
    if (DN[2] == 0 && DN[3] == 2) Canon3Prime[2] = 1;
  } else {
    SpliceCodons[2] = -1; // Stop codon
  }


  // Option 4: |xx...|yyDEF|
  Codon[1] = DN[3];
  SpliceCodons[3] = esl_gencode_GetTranslation(gcode,&Codon[1]);

  if (SpliceCodons[3] >= 0 && SpliceCodons[3] <= 20) {
    if (UN[1] == 2 && UN[2] == 3) Canon5Prime[3] = 1;
    if (DN[1] == 0 && DN[2] == 2) Canon3Prime[3] = 1;
  } else {
    SpliceCodons[3] = -1; // Stop codon
  }


  free(Codon);
  free(UN);
  free(DN);

  if (DEBUGGING) DEBUG_OUT("'GetSpliceOptions' Complete",-1);

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindOptimalSpliceSite
 *
 *  Desc. :
 *
 *  Inputs:  1.                 Overlap :
 *           2.                      gm : The straightforward profile for the protein / family.
 *           3.                   gcode : An ESL_GENCODE struct (mainly used for translation).
 *           4.   upstream_splice_index :
 *           5. downstream_splice_index :
 *           6. split_amino_model_index :
 *
 *  Output:  A float representing the total score of the optimal splicing of the two hits.
 *
 */
float FindOptimalSpliceSite
(
  DOMAIN_OVERLAP * Overlap,
  P7_PROFILE     * gm,
  ESL_GENCODE    * gcode,
  int *   upstream_splice_index,
  int * downstream_splice_index,
  int * split_amino_model_index
)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'FindOptimalSpliceSite'",1);



  //
  //  UPSTREAM 
  //


  int upstream_nucl_cnt = abs(Overlap->upstream_nucl_start - Overlap->upstream_nucl_end) + 1;

  // Oversize -- Deal with it!
  int   * USTrans    = malloc(upstream_nucl_cnt*sizeof(int));
  int   * USModelPos = malloc(upstream_nucl_cnt*sizeof(int));
  float * USScores   = malloc(upstream_nucl_cnt*sizeof(float));


  int us_trans_len   = 0; // Without messy indels, comes out to upstream_nucl_cnt/3
  int nucl_read_pos  = 1;
  int model_pos      = Overlap->amino_start;
  int display_pos    = Overlap->upstream_disp_start;
  while (model_pos <= Overlap->UpstreamDisplay->hmmto) {


    int amino_index = 27;
    if (Overlap->UpstreamDisplay->aseq[display_pos] != '-') {
      amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->UpstreamNucls[nucl_read_pos]));
      nucl_read_pos += 3;
    }


    USTrans[us_trans_len]    = amino_index;
    USModelPos[us_trans_len] = model_pos;
    if (amino_index == 27) 
      USScores[us_trans_len] = -0.7; // Let's try this...
    else 
      USScores[us_trans_len] = gm->rsc[amino_index][2*model_pos];


    if (Overlap->UpstreamDisplay->aseq[display_pos] != '.')
      model_pos++;
    display_pos++;
    us_trans_len++;


  }
  int us_pre_ext_end_pos = us_trans_len-1; // For scoring purposes


  // Incorporate the extension
  for (int i=1; i<=Overlap->upstream_ext_len; i++) {


    int amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->UpstreamNucls[nucl_read_pos]));
    nucl_read_pos += 3;


    USTrans[us_trans_len]    = amino_index;
    USModelPos[us_trans_len] = model_pos;
    USScores[us_trans_len]   = gm->rsc[amino_index][2*model_pos];


    model_pos++;
    us_trans_len++;


  }




  //
  //  DOWNSTREAM
  //


  int downstream_nucl_cnt = abs(Overlap->downstream_nucl_start - Overlap->downstream_nucl_end) + 1;

  // Oversize -- deal with it!
  int   * DSTrans    = malloc(downstream_nucl_cnt*sizeof(int));
  int   * DSModelPos = malloc(downstream_nucl_cnt*sizeof(int));
  float * DSScores   = malloc(downstream_nucl_cnt*sizeof(float));

  int ds_trans_len = 0;
  nucl_read_pos    = 1;
  model_pos        = Overlap->amino_start;
  while (ds_trans_len < Overlap->downstream_ext_len) {


    int amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->DownstreamNucls[nucl_read_pos]));
    nucl_read_pos += 3;


    DSTrans[ds_trans_len]    = amino_index;
    DSModelPos[ds_trans_len] = model_pos;
    DSScores[ds_trans_len]   = gm->rsc[amino_index][2*model_pos];


    model_pos++;
    ds_trans_len++;


  }



  display_pos = 0;
  while (model_pos <= Overlap->amino_end) {


    int amino_index = 27;
    if (Overlap->DownstreamDisplay->aseq[display_pos] != '-') {
      amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->DownstreamNucls[nucl_read_pos]));
      nucl_read_pos += 3;
    }


    DSTrans[ds_trans_len]    = amino_index;
    DSModelPos[ds_trans_len] = model_pos;
    if (amino_index == 27) 
      DSScores[ds_trans_len] = -0.7; // Let's try this...
    else 
      DSScores[ds_trans_len] = gm->rsc[amino_index][2*model_pos];


    if (Overlap->DownstreamDisplay->aseq[display_pos] != '.')
      model_pos++;
    display_pos++;
    ds_trans_len++;


  }



  // We're really just interested in the sum scores on each side,
  // so let's switch over to that
  for (int i=1             ; i<us_trans_len; i++) USScores[i] += USScores[i-1];
  for (int i=ds_trans_len-2; i>=0          ; i--) DSScores[i] += DSScores[i+1];


  // In order to determine the score difference created by the
  // splice, we'll find the sum score at each position right
  // before the extension, and those become the baseline to
  // remove for determining the contribution of our cut.
  float baseline_score = DSScores[Overlap->downstream_ext_len];
  baseline_score      += USScores[us_pre_ext_end_pos];


  // What position in the model are we splitting on?
  int   optimal_us_pos    = 0;
  int   optimal_ds_pos    = 0;
  int   optimal_model_pos = 0;
  float optimal_score     = 0.0;

  int us_start = 0;
  int ds_start = 0;
  for (model_pos = Overlap->amino_start; model_pos < Overlap->amino_end; model_pos++) {

    
    while (USModelPos[us_start] < model_pos) us_start++;
    while (DSModelPos[ds_start] < model_pos) ds_start++;


    int us_pos = us_start;
    int ds_pos;
    while (USModelPos[us_pos] == model_pos) {

      ds_pos = ds_start;

      while (DSModelPos[ds_pos] == model_pos) {

        float sum_score = USScores[us_pos] + DSScores[ds_pos];

        if (sum_score > optimal_score) {
          optimal_score     = sum_score;
          optimal_us_pos    = us_pos;
          optimal_ds_pos    = ds_pos;
          optimal_model_pos = model_pos;
        }

        ds_pos++;

      }

      us_pos++;

    }

    us_start = us_pos;
    ds_start = ds_pos;

  }
  free(USTrans);
  free(USModelPos);
  free(USScores);
  free(DSTrans);
  free(DSModelPos);
  free(DSScores);


  // Currently, the optimal positions share an amino,
  // which we're going to want to consider methods
  // for splitting.
  //
  // For that reason, what we return as the splice
  // indices (relative to the *nucleotide* sequences)
  // is the last / first "safe" nucelotide (inside the
  // exon).
  //
  *upstream_splice_index   = 3 * optimal_us_pos;
  *downstream_splice_index = 3 * optimal_ds_pos + 4;
  *split_amino_model_index = optimal_model_pos;


  if (DEBUGGING) DEBUG_OUT("'FindOptimalSpliceSite' Complete",-1);


  return optimal_score-baseline_score;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SpliceOverlappingDomains
 *
 *  Desc. :
 *
 *  Inputs:  1.       Overlap :
 *           2. TargetNuclSeq :
 *           3.            gm : The straightforward profile for the protein / family.
 *           4.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void SpliceOverlappingDomains
(
  DOMAIN_OVERLAP * Overlap,
  TARGET_SEQ     * TargetNuclSeq,
  P7_PROFILE     * gm, 
  ESL_GENCODE    * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SpliceOverlappingDomains'",1);


  // These splice indices are the terminal nucleotides *within* the
  // uncontested coding region, relative to the nucleotide sequences
  // in Overlap
  int   upstream_splice_index;
  int downstream_splice_index;
  int split_amino_model_index;
  float splice_score = FindOptimalSpliceSite(Overlap,gm,gcode,&upstream_splice_index,&downstream_splice_index,&split_amino_model_index);


  if (splice_score == EDGE_FAIL_SCORE) {
    Overlap->score         = EDGE_FAIL_SCORE;
    Overlap->score_density = EDGE_FAIL_SCORE;
    if (DEBUGGING) DEBUG_OUT("'SpliceOverlappingDomains' Complete (splicing score was -inf...?)",-1);
    return;
  }


  // At this point, the optimal amino acid is the one that's
  // contested between the two hits, so we need to decide how
  // best to splice it.
  int * SpliceCodons = malloc(4*sizeof(int));
  int * Canon5Prime  = malloc(4*sizeof(int)); // GT
  int * Canon3Prime  = malloc(4*sizeof(int)); // AG

  GetSpliceOptions(Overlap,TargetNuclSeq,upstream_splice_index,downstream_splice_index,SpliceCodons,Canon5Prime,Canon3Prime,gcode);



  // Let's see which SpliceCodon has the best match to the model at
  // the "model_ss," while also factoring in splice site signals.
  int   best_split_opt   = -1;
  float best_split_score = EDGE_FAIL_SCORE;
  for (int i=0; i<4; i++) {
    
    // Did we end up with a stop codon?
    if (SpliceCodons[i] == -1) continue;

    float split_score = gm->rsc[SpliceCodons[i]][split_amino_model_index * 2] + SSSCORE[Canon5Prime[i]] + SSSCORE[Canon3Prime[i]];
    if (split_score > best_split_score) {
      best_split_score = split_score;
      best_split_opt   = i;
    }

  }
  free(SpliceCodons);
  free(Canon5Prime);
  free(Canon3Prime);
  


  // This would be really bizarre (and worth checking to see if it
  // ever actually happens...)
  if (best_split_opt == -1) {
    Overlap->score         = EDGE_FAIL_SCORE;
    Overlap->score_density = EDGE_FAIL_SCORE;
    if (DEBUGGING) DEBUG_OUT("'SpliceOverlappingDomains' Complete (BUT WITH TERRIBLE OPTIONS?!)",-1);
    return;
  }



  if (Overlap->upstream_nucl_start < Overlap->upstream_nucl_end) {
    Overlap->upstream_spliced_nucl_end = Overlap->upstream_nucl_start + (upstream_splice_index + (3 - best_split_opt) - 1);
    Overlap->downstream_spliced_nucl_start = Overlap->downstream_nucl_start + (downstream_splice_index - best_split_opt - 1);
  } else {
    Overlap->upstream_spliced_nucl_end = Overlap->upstream_nucl_start - (upstream_splice_index + (3 - best_split_opt) - 1);
    Overlap->downstream_spliced_nucl_start = Overlap->downstream_nucl_start - (downstream_splice_index - best_split_opt - 1);
  }


  if (best_split_opt < 2) {
    Overlap->upstream_exon_terminus   = split_amino_model_index;
    Overlap->downstream_exon_terminus = split_amino_model_index + 1;
  } else {
    Overlap->upstream_exon_terminus   = split_amino_model_index - 1;
    Overlap->downstream_exon_terminus = split_amino_model_index;    
  }


  Overlap->score = splice_score + best_split_score;
  Overlap->score_density = Overlap->score / (1 + Overlap->amino_end - Overlap->amino_start);


  if (DEBUGGING) DEBUG_OUT("'SpliceOverlappingDomains' Complete",-1);

}






/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetNuclRangesFromAminoCoords
 *
 *  Desc. :
 *
 *  Inputs:  1. Edge :
 *
 *  Output:
 *
 */
void GetNuclRangesFromAminoCoords 
(DOMAIN_OVERLAP * Edge)
{

  P7_ALIDISPLAY *   UpDisp = Edge->UpstreamDisplay;
  P7_ALIDISPLAY * DownDisp = Edge->DownstreamDisplay;


  int strand = 1;
  if (UpDisp->sqfrom > UpDisp->sqto)
    strand = -1;


  //
  // UPSTREAM
  //

  Edge->upstream_disp_start = UpDisp->N - 1;
  int disp_amino = UpDisp->hmmto;

  Edge->upstream_nucl_end   = UpDisp->sqto + (3 * strand * (Edge->amino_end - disp_amino));
  Edge->upstream_nucl_start = UpDisp->sqto + strand;

  while (Edge->upstream_disp_start >= 0 && disp_amino >= Edge->amino_start) {

    Edge->upstream_nucl_start -= 3 * strand;
    disp_amino--;

    // If we were presumptuous, undo the previous work
    if (UpDisp->model[Edge->upstream_disp_start] == '.') // Insertion relative to pHMM
      disp_amino++;
    if (UpDisp->aseq[Edge->upstream_disp_start] == '-') // Insertion relative to genome
      Edge->upstream_nucl_start += 3 * strand;

    Edge->upstream_disp_start -= 1;

  }
  disp_amino++;
  Edge->upstream_ext_len     = disp_amino - Edge->amino_start;
  Edge->upstream_nucl_start -= 3 * strand * Edge->upstream_ext_len;
  
  // We'll have overstepped by one
  Edge->upstream_disp_start += 1;



  //
  // DOWNSTREAM
  //

  Edge->downstream_disp_end = 0;
  disp_amino = DownDisp->hmmfrom;

  Edge->downstream_nucl_start = DownDisp->sqfrom - (3 * strand * (disp_amino - Edge->amino_start));
  Edge->downstream_nucl_end   = DownDisp->sqfrom - strand;

  while (Edge->downstream_disp_end < DownDisp->N && disp_amino <= Edge->amino_end) {

    Edge->downstream_nucl_end += 3 * strand;
    disp_amino++;

    // Presumptuous? I hardly knumptuous!
    if (DownDisp->model[Edge->downstream_disp_end] == '.') // Insertion relative to pHMM
      disp_amino--;
    if (DownDisp->aseq[Edge->downstream_disp_end] == '-') // Insertion relative to genome
      Edge->downstream_nucl_end -= 3 * strand;

    Edge->downstream_disp_end += 1;

  }
  disp_amino--;
  Edge->downstream_ext_len   = Edge->amino_end - disp_amino;
  Edge->downstream_nucl_end += 3 * strand * Edge->downstream_ext_len;

  // We'll have overstepped by one
  Edge->upstream_disp_start -= 1;


}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SketchSpliceEdge
 *
 *  Desc. :
 *
 *  Inputs:  1.          Edge :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.            gm : The straightforward profile for the protein / family.
 *           4.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void SketchSpliceEdge
(
  DOMAIN_OVERLAP * Edge,
  TARGET_SEQ     * TargetNuclSeq,
  P7_PROFILE     * gm,
  ESL_GENCODE    * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SketchSpliceEdge'",1);


  P7_ALIDISPLAY *   UpDisp = Edge->UpstreamDisplay;
  P7_ALIDISPLAY * DownDisp = Edge->DownstreamDisplay;


  Edge->amino_start = DownDisp->hmmfrom;
  Edge->amino_end   =   UpDisp->hmmto;


  // Do we need to extend beyond the bounds of these
  // hits to have the required number of overlapping
  // amino acids?
  int num_ext_aminos = MIN_AMINO_OVERLAP - (1 + Edge->amino_end - Edge->amino_start);
  if (num_ext_aminos > 0) {

    // Extending the overlap region means going in
    // both directions, so we'll extend each side by
    // half of what's required
    num_ext_aminos = (num_ext_aminos+1)/2;

    Edge->amino_start -= num_ext_aminos;
    Edge->amino_end   += num_ext_aminos;

  }


  // Now we can do the work of finding the (indel-aware)
  // upstream_start and downstream_end coordinates.
  GetNuclRangesFromAminoCoords(Edge);


  // Grab them nucleotides!
  Edge->UpstreamNucls   = GrabNuclRange(TargetNuclSeq,Edge->upstream_nucl_start,Edge->upstream_nucl_end);
  Edge->DownstreamNucls = GrabNuclRange(TargetNuclSeq,Edge->downstream_nucl_start,Edge->downstream_nucl_end);


  // Finish off by adding this friendly little pointer
  Edge->ntalpha = TargetNuclSeq->abc;


  if (DEBUGGING && 0) {
    fprintf(stderr,"\n");
    fprintf(stderr,"  Overlap  Nucl. Range:   Upstream : %d ... %d\n",  Edge->upstream_nucl_start,  Edge->upstream_nucl_end);
    fprintf(stderr,"                                   : ");
    for (int i=1; i<=abs(Edge->upstream_nucl_start-Edge->upstream_nucl_end)+1; i++)
      fprintf(stderr,"%c",DNA_CHARS[Edge->UpstreamNucls[i]]);
    fprintf(stderr,"\n");
    fprintf(stderr,"                      : Downstream : %d ... %d\n",Edge->downstream_nucl_start,Edge->downstream_nucl_end);
    fprintf(stderr,"                                   : ");
    for (int i=1; i<=abs(Edge->downstream_nucl_start-Edge->downstream_nucl_end)+1; i++)
      fprintf(stderr,"%c",DNA_CHARS[Edge->DownstreamNucls[i]]);
    fprintf(stderr,"\n\n");
  }



  SpliceOverlappingDomains(Edge,TargetNuclSeq,gm,gcode);


  if (DEBUGGING) DEBUG_OUT("'SketchSpliceEdge' Complete",-1);

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: HitsAreSpliceComaptible
 *
 *  Desc. :
 *
 *  Inputs:  1.   Upstream :
 *           2. Downstream :
 *
 *  Output:
 *
 */
int HitsAreSpliceCompatible
(P7_ALIDISPLAY * Upstream, P7_ALIDISPLAY * Downstream)
{

   
  if (DEBUGGING) DEBUG_OUT("Starting 'HitsAreSpliceCompatible'",1);

  
  // Start by checking if we either have amino acid
  // overlap, or are close enough to consider extending
  int amino_start_1 = Upstream->hmmfrom;
  int amino_end_1   = Upstream->hmmto;

  int amino_start_2 = Downstream->hmmfrom;
  int amino_end_2   = Downstream->hmmto;

  // If the upstream ain't upstream, then obviously we can't treat
  // these as splice-compatible!
  if (!(amino_start_1 < amino_start_2 && amino_end_1 < amino_end_2)) {
    if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }

  // Do we have overlap OR sufficient proximity to consider
  // extending?
  if (!(amino_end_1 + MAX_AMINO_EXT >= amino_start_2)) {
    if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }


  // Fantastic!  The amino acid coordinates support splice
  // compatibility!  Now it's just time to confirm that
  // the nucleotides also look good.

  int  nucl_start_1 = Upstream->sqfrom;
  int  nucl_end_1   = Upstream->sqto;
  
  int revcomp1 = 0;
  if (nucl_start_1 > nucl_end_1)
    revcomp1 = 1;


  int nucl_start_2 = Downstream->sqfrom;
  int nucl_end_2   = Downstream->sqto;
 
  int revcomp2 = 0;
  if (nucl_start_2 > nucl_end_2)
    revcomp2 = 1;


  if (revcomp1 != revcomp2) {
    if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }


  // We want to make sure that these aren't unrealistically
  // close together on the genome...
  if (revcomp1) {

    if (nucl_start_2 + (3 * MAX_AMINO_EXT) >= nucl_end_1) {
      if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
      return 0;
    }

  } else {

    if (nucl_start_2 - (3 * MAX_AMINO_EXT) <= nucl_end_1) {
      if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
      return 0;
    }

  }


  if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);


  // Looks like we've got a viable upstream / downstream pair!
  return 1;

}






/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GatherViableSpliceEdges
 *
 *  Desc. :
 *
 *  Inputs:  1.          TopHits :
 *           2.    TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.               gm : The straightforward profile for the protein / family.
 *           4.            gcode : An ESL_GENCODE struct (mainly used for translation).
 *           5. num_splice_edges :
 *
 *  Output:
 *
 */
DOMAIN_OVERLAP ** GatherViableSpliceEdges
(
  P7_TOPHITS  * TopHits,
  TARGET_SEQ  * TargetNuclSeq,
  P7_PROFILE  * gm,
  ESL_GENCODE * gcode,
  int * num_splice_edges
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GatherViableSpliceEdges'",1);


  int num_hits = (int)(TopHits->N);

  int edge_capacity = 2 * num_hits;
  DOMAIN_OVERLAP ** SpliceEdges = (DOMAIN_OVERLAP **)malloc(edge_capacity * sizeof(DOMAIN_OVERLAP *));

  int num_edges = 0;
  for (int upstream_hit_id = 0; upstream_hit_id < num_hits; upstream_hit_id++) {


    P7_HIT * UpstreamHit  = TopHits->hit[upstream_hit_id];
    int num_upstream_doms = UpstreamHit->ndom;


    // For each hit, gather all of the indices of other hits that
    // could potentially be downstream exons.
    for (int downstream_hit_id=0; downstream_hit_id < (int)(TopHits->N); downstream_hit_id++) {


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

        P7_ALIDISPLAY * UpstreamDisplay = (&UpstreamHit->dcl[upstream_dom_id])->ad;

        for (int downstream_dom_id = 0; downstream_dom_id < num_downstream_doms; downstream_dom_id++) {

          // NO SELF-SPLICING, YOU ABSOLUTE MANIAC!
          if (upstream_hit_id == downstream_hit_id && upstream_dom_id == downstream_dom_id)
            continue;


          P7_ALIDISPLAY * DownstreamDisplay = (&DownstreamHit->dcl[downstream_dom_id])->ad;

          if (HitsAreSpliceCompatible(UpstreamDisplay,DownstreamDisplay)) {

            // MUST WE RESIZE?!
            if (num_edges == edge_capacity) {

              edge_capacity *= 2;

              DOMAIN_OVERLAP ** MoreSpliceEdges = (DOMAIN_OVERLAP **)malloc(edge_capacity * sizeof(DOMAIN_OVERLAP *));
              for (int j=0; j<num_edges; j++)
                MoreSpliceEdges[j] = SpliceEdges[j];

              free(SpliceEdges);
              SpliceEdges = MoreSpliceEdges;

            }


            // Record that splice compatibility!
            SpliceEdges[num_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
            DOMAIN_OVERLAP * Edge  = SpliceEdges[num_edges];

            Edge->upstream_hit_id   = upstream_hit_id;
            Edge->upstream_dom_id   = upstream_dom_id;
            Edge->downstream_hit_id = downstream_hit_id;
            Edge->downstream_dom_id = downstream_dom_id;

            Edge->UpstreamTopHits   = TopHits;
            Edge->UpstreamDisplay   = UpstreamDisplay;
            Edge->DownstreamTopHits = TopHits;
            Edge->DownstreamDisplay = DownstreamDisplay;

            num_edges++;

          }

        }

      }

    }

  }


  //
  //  Now that we have our splice edges, we can more fully
  //  sketch out how they connect!
  //
  //
  //  This *could* be part of the above loop, but what's the
  //  rush?
  //


  // We'll run through all of our paired domains and actually
  // splice 'em up (or at least try our best to)!
  for (int splice_edge_id = 0; splice_edge_id < num_edges; splice_edge_id++) {

    SketchSpliceEdge(SpliceEdges[splice_edge_id],TargetNuclSeq,gm,gcode);

    // If we failed to find a reasonable splice site, then we'll
    // just rip this edge outta consideration.
    if (SpliceEdges[splice_edge_id]->score == EDGE_FAIL_SCORE) {
      free(SpliceEdges[splice_edge_id]);
      SpliceEdges[splice_edge_id] = NULL;
    }

  }



  if (DEBUGGING) DEBUG_OUT("'GatherViableSpliceEdges' Complete",-1);


  *num_splice_edges = num_edges;
  return SpliceEdges;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: InitSpliceNode
 *
 *  Desc. :
 *
 *  Inputs:  1.      Graph :
 *           2.    node_id :
 *           3.     hit_id :
 *           4.     dom_id :
 *           5. was_missed :
 *
 *  Output:
 *
 */
SPLICE_NODE * InitSpliceNode
(
  SPLICE_GRAPH * Graph,
  int node_id, 
  int hit_id, 
  int dom_id,
  int was_missed
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'InitSpliceNode'",1);
  

  SPLICE_NODE * NewNode = (SPLICE_NODE *)malloc(sizeof(SPLICE_NODE));

  NewNode->node_id = node_id;
  NewNode->hit_id  = hit_id;
  NewNode->dom_id  = dom_id;


  NewNode->out_edge_cap    = 10;
  NewNode->num_out_edges   =  0;
  NewNode->best_out_edge   = -1;
  NewNode->OutEdges        = (DOMAIN_OVERLAP **)malloc(NewNode->out_edge_cap*sizeof(DOMAIN_OVERLAP *));
  NewNode->DownstreamNodes = (SPLICE_NODE    **)malloc(NewNode->out_edge_cap*sizeof(SPLICE_NODE    *));


  NewNode->in_edge_cap    = 10;
  NewNode->num_in_edges   =  0;
  NewNode->best_in_edge   = -1;
  NewNode->InEdges        = (DOMAIN_OVERLAP **)malloc(NewNode->in_edge_cap*sizeof(DOMAIN_OVERLAP *));
  NewNode->UpstreamNodes  = (SPLICE_NODE    **)malloc(NewNode->in_edge_cap*sizeof(SPLICE_NODE    *));


  NewNode->was_missed = was_missed;


  if (was_missed) {
  
    NewNode->is_n_terminal = 0;
    if (Graph->MissedHits->hit[hit_id]->dcl->ad->hmmfrom == 1)
      NewNode->is_n_terminal = 1;

    NewNode->is_c_terminal = 0;
    if (Graph->MissedHits->hit[hit_id]->dcl->ad->hmmto == Graph->Model->M)
      NewNode->is_c_terminal = 1;

    NewNode->hit_score = Graph->MissedHits->hit[hit_id]->dcl->bitscore;
  
  } else {
  
    NewNode->is_n_terminal = 0;
    if ((&(Graph->TopHits->hit[hit_id]->dcl[dom_id]))->ad->hmmfrom == 1)
      NewNode->is_n_terminal = 1;

    NewNode->is_c_terminal = 0;
    if ((&(Graph->TopHits->hit[hit_id]->dcl[dom_id]))->ad->hmmto == Graph->Model->M)
      NewNode->is_c_terminal = 1;

    NewNode->hit_score = Graph->TopHits->hit[hit_id]->dcl[dom_id].bitscore;

  }


  // Initialize these scores to recognizable low values
  NewNode->cumulative_score = EDGE_FAIL_SCORE; // Best score up to and including this node
  NewNode->best_path_score  = EDGE_FAIL_SCORE; // Best full path score using this node


  if (DEBUGGING) DEBUG_OUT("'InitSpliceNode' Complete",-1);

  return NewNode;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ConnectNodesByEdge
 *
 *  Desc. :
 *
 *  Inputs:  1.  Edge :
 *           2. Graph :
 *
 *  Output:
 *
 */
void ConnectNodesByEdge
(DOMAIN_OVERLAP * Edge, SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'ConnectNodesByEdge'",1);


  int upstream_node_id;
  if (Edge->UpstreamTopHits == Graph->TopHits) {
    upstream_node_id = Graph->TH_HitDomToNodeID[Edge->upstream_hit_id][Edge->upstream_dom_id];
  } else {
    upstream_node_id = Graph->MH_HitToNodeID[Edge->upstream_hit_id];
  }
  SPLICE_NODE * UpstreamNode = Graph->Nodes[upstream_node_id];


  int downstream_node_id;
  if (Edge->DownstreamTopHits == Graph->TopHits) {
    downstream_node_id = Graph->TH_HitDomToNodeID[Edge->downstream_hit_id][Edge->downstream_dom_id];
  } else {
    downstream_node_id = Graph->MH_HitToNodeID[Edge->downstream_hit_id];
  }
  SPLICE_NODE * DownstreamNode = Graph->Nodes[downstream_node_id];



  UpstreamNode->OutEdges[UpstreamNode->num_out_edges] = Edge;
  UpstreamNode->DownstreamNodes[UpstreamNode->num_out_edges] = DownstreamNode;
  UpstreamNode->num_out_edges += 1;


  DownstreamNode->InEdges[DownstreamNode->num_in_edges] = Edge;
  DownstreamNode->UpstreamNodes[DownstreamNode->num_in_edges]  = UpstreamNode;
  DownstreamNode->num_in_edges += 1;



  // Resize?
  if (UpstreamNode->num_out_edges == UpstreamNode->out_edge_cap) {

    UpstreamNode->out_edge_cap *= 2;

    DOMAIN_OVERLAP ** NewOutEdges = (DOMAIN_OVERLAP **)malloc(UpstreamNode->out_edge_cap*sizeof(DOMAIN_OVERLAP *)); 
    SPLICE_NODE    ** NewDSNodes  = (SPLICE_NODE    **)malloc(UpstreamNode->out_edge_cap*sizeof(SPLICE_NODE    *));

    for (int i = 0; i < UpstreamNode->num_out_edges; i++) {
      NewOutEdges[i] = UpstreamNode->OutEdges[i];
      NewDSNodes[i]  = UpstreamNode->DownstreamNodes[i];
    }

    free(UpstreamNode->OutEdges);
    free(UpstreamNode->DownstreamNodes);

    UpstreamNode->OutEdges = NewOutEdges;
    UpstreamNode->DownstreamNodes = NewDSNodes;

  }
  // Resize?
  if (DownstreamNode->num_in_edges == DownstreamNode->in_edge_cap) {

    DownstreamNode->in_edge_cap *= 2;

    DOMAIN_OVERLAP ** NewInEdges = (DOMAIN_OVERLAP **)malloc(DownstreamNode->in_edge_cap*sizeof(DOMAIN_OVERLAP *)); 
    SPLICE_NODE    ** NewUSNodes = (SPLICE_NODE    **)malloc(DownstreamNode->in_edge_cap*sizeof(SPLICE_NODE    *));

    for (int i = 0; i < DownstreamNode->num_in_edges; i++) {
      NewInEdges[i] = DownstreamNode->InEdges[i];
      NewUSNodes[i] = DownstreamNode->UpstreamNodes[i];
    }

    free(DownstreamNode->InEdges);
    free(DownstreamNode->UpstreamNodes);

    DownstreamNode->InEdges = NewInEdges;
    DownstreamNode->UpstreamNodes = NewUSNodes;

  }



  Graph->num_edges += 1;


  if (DEBUGGING) DEBUG_OUT("'ConnectNodesByEdge' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindBestPathToNode
 *
 *  Desc. :
 *
 *  Inputs:  1. Node :  
 *
 *  Output:
 *
 */
void FindBestPathToNode
(SPLICE_NODE * Node)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'FindBestPathToNode'",1);

  
  // Have we already examined this node?
  if (Node->cumulative_score != EDGE_FAIL_SCORE) {
    if (DEBUGGING) DEBUG_OUT("'FindBestPathToNode' Complete",-1);
    return;
  }


  // For each of the nodes that feed into this node,
  // recursively find *their* best path, then determine
  // across those options which is our favorite.
  for (int in_edge_id = 0; in_edge_id < Node->num_in_edges; in_edge_id++) {
  

    FindBestPathToNode(Node->UpstreamNodes[in_edge_id]);


    // What's the score of the path up to the current node using
    // this as our input edge?
    float edge_sum_score = Node->UpstreamNodes[in_edge_id]->cumulative_score;
    edge_sum_score += Node->InEdges[in_edge_id]->score;


    if (edge_sum_score > Node->cumulative_score) {
      Node->cumulative_score = edge_sum_score;
      Node->best_in_edge     = in_edge_id;
    }


  }


  // If this node doesn't have any inputs, its cumulative score
  // should be initialized to 'EDGE_FAIL_SCORE'
  if (Node->cumulative_score == EDGE_FAIL_SCORE) 
    Node->cumulative_score  = Node->hit_score;
  else
    Node->cumulative_score += Node->hit_score;


  // DEBUGGING
  fprintf(stderr,"  :: Node ID:%d --> [ %f / %f / %f ]\n",Node->node_id,Node->hit_score,Node->cumulative_score,Node->best_path_score);


  if (DEBUGGING) DEBUG_OUT("'FindBestPathToNode' Complete",-1);

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: EvangelizePath
 *
 *  Desc. :
 *
 *  Inputs:  1. Node :
 *
 *  Output:
 *
 */
void EvangelizePath
(SPLICE_NODE * Node)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'EvangelizePath'",1);

  if (Node->best_path_score == EDGE_FAIL_SCORE) 
    Node->best_path_score = Node->cumulative_score;

  for (int i=0; i<Node->num_in_edges; i++) {

    SPLICE_NODE * UpstreamNode = Node->UpstreamNodes[i];

    if (Node->best_path_score > UpstreamNode->best_path_score) {

      // NOTE: This isn't necessarily the most rigorous
      //       definition of the best outgoing edge, but
      //       for our purposes it's fine.
      for (int j=0; j<UpstreamNode->num_out_edges; j++) {
        if (UpstreamNode->DownstreamNodes[j] == Node) {
          UpstreamNode->best_out_edge = j;
          break;
        }
      }

      // We'll only set the upstream node's best_path_score
      // to this node's bps *if* the upstream node is how we
      // achieved our best_path_score!
      if (i == Node->best_in_edge) {
        UpstreamNode->best_path_score = Node->best_path_score;    
        EvangelizePath(UpstreamNode);
      }
    
    }
  }

  if (DEBUGGING) DEBUG_OUT("'EvangelizePath' Complete",-1);

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GatherNTermNodes
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void GatherNTermNodes
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GatherNTermNodes'",1);

  // We'll re-count N terminal hits, just to be certain
  int num_n_term = 0;

  int   * NTermIDs    = malloc(Graph->num_nodes * sizeof(int  ));
  float * NTermScores = malloc(Graph->num_nodes * sizeof(float));

  for (int i=1; i<=Graph->num_nodes; i++) {
    if (Graph->Nodes[i]->is_n_terminal) {
      NTermIDs[num_n_term] = i;
      NTermScores[num_n_term] = Graph->Nodes[i]->best_path_score;
      num_n_term++;
    }
  }

  int * NTermScoreSort = FloatHighLowSortIndex(NTermScores,num_n_term);
  
  if (Graph->NTermNodeIDs)
    free(Graph->NTermNodeIDs);

  Graph->NTermNodeIDs = malloc(num_n_term * sizeof(int));
  Graph->num_n_term   = num_n_term;
  for (int i=0; i<num_n_term; i++)
    Graph->NTermNodeIDs[i] = NTermIDs[NTermScoreSort[i]];

  free(NTermIDs);
  free(NTermScores);
  free(NTermScoreSort);

  if (DEBUGGING) DEBUG_OUT("'GatherNTermNodes' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GatherCTermNodes
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void GatherCTermNodes
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GatherCTermNodes'",1);

  // We'll re-count C terminal hits, just to be certain
  int num_c_term = 0;

  int   * CTermIDs    = malloc(Graph->num_nodes * sizeof(int  ));
  float * CTermScores = malloc(Graph->num_nodes * sizeof(float));

  for (int i=1; i<=Graph->num_nodes; i++) {
    if (Graph->Nodes[i]->is_c_terminal) {
      CTermIDs[num_c_term] = i;
      CTermScores[num_c_term] = Graph->Nodes[i]->best_path_score;
      num_c_term++;
    }
  }

  int * CTermScoreSort = FloatHighLowSortIndex(CTermScores,num_c_term);
  
  if (Graph->CTermNodeIDs)
    free(Graph->CTermNodeIDs);

  Graph->CTermNodeIDs = malloc(num_c_term * sizeof(int));
  Graph->num_c_term   = num_c_term;
  for (int i=0; i<num_c_term; i++)
    Graph->CTermNodeIDs[i] = CTermIDs[CTermScoreSort[i]];

  free(CTermIDs);
  free(CTermScores);
  free(CTermScoreSort);

  if (DEBUGGING) DEBUG_OUT("'GatherNTermNodes' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GenCumScoreSort
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void GenCumScoreSort
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GenCumScoreSort'",1);

  float * Scores = malloc(Graph->num_nodes * sizeof(float));
  for (int i=0; i<Graph->num_nodes; i++)
    Scores[i] = Graph->Nodes[i+1]->cumulative_score;

  int * SortIndex = FloatHighLowSortIndex(Scores,Graph->num_nodes);
  for (int i=0; i<Graph->num_nodes; i++)
    SortIndex[i]++;

  if (Graph->CumScoreSort != NULL)
    free(Graph->CumScoreSort);
  Graph->CumScoreSort = SortIndex;

  free(Scores);

  if (DEBUGGING) DEBUG_OUT("'GenCumScoreSort' Complete",-1);

}






/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: EvaluatePaths
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void EvaluatePaths
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'EvaluatePaths'",1);


  // We'll reset all of our pathfinding, in case we're
  // running this after adding missing exons to the graph.
  Graph->has_full_path         = 0;
  Graph->best_full_path_start  = 0;
  Graph->best_full_path_end    = 0;
  Graph->best_full_path_length = 0;
  Graph->best_full_path_score  = EDGE_FAIL_SCORE;
  for (int node_id=1; node_id<=Graph->num_nodes; node_id++) {
    Graph->Nodes[node_id]->best_in_edge     = -1;
    Graph->Nodes[node_id]->best_out_edge    = -1;
    Graph->Nodes[node_id]->cumulative_score =  EDGE_FAIL_SCORE;
    Graph->Nodes[node_id]->best_path_score  =  EDGE_FAIL_SCORE;
  }


  // Find the best path to each node
  for (int node_id=1; node_id<=Graph->num_nodes; node_id++)
    FindBestPathToNode(Graph->Nodes[node_id]);

  GenCumScoreSort(Graph);

  for (int sort_id=0; sort_id<Graph->num_nodes; sort_id++) {
    
    int node_id = Graph->CumScoreSort[sort_id];
    SPLICE_NODE * Node = Graph->Nodes[node_id];
    EvangelizePath(Node);

  }


  // We gather the N- and C-terminal nodes at the end of this
  // function because these are fundamentally an ordering of
  // paths through the graph
  GatherNTermNodes(Graph);
  GatherCTermNodes(Graph);


  if (DEBUGGING) DEBUG_OUT("'EvaluatePaths' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FillOutGraphStructure
 *
 *  Desc. :
 *
 *  Inputs:  1.            Graph :
 *           2.      SpliceEdges :
 *           3. num_splice_edges :
 *
 *  Output:
 *
 */
void FillOutGraphStructure
(SPLICE_GRAPH * Graph, DOMAIN_OVERLAP ** SpliceEdges, int num_splice_edges)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'FillOutGraphStructure'",1);


  // We'll want to be a little careful, just because this is our flagship
  // datastructure...
  if (Graph->Nodes != NULL) {
    for (int i=1; i<=Graph->num_nodes; i++)
      if (Graph->Nodes[i] != NULL)
        free(Graph->Nodes[i]);
    free(Graph->Nodes);
  }


  // We'll want this lookup table to be able to go from
  // DOMAIN_OVERLAP content
  int num_hits = (int)(Graph->TopHits->N);
  int max_doms = 0;
  int sum_doms = 0;
  for (int hit_id=0; hit_id<num_hits; hit_id++) {
    int hit_doms = Graph->TopHits->hit[hit_id]->ndom;
    if (hit_doms > max_doms)
      max_doms = hit_doms;
    sum_doms += hit_doms;
  }

  // Allocate space for the maximum number of domains
  Graph->Nodes = (SPLICE_NODE **)malloc((sum_doms+1)*sizeof(SPLICE_NODE *));


  Graph->TH_HitDomToNodeID = malloc(num_hits * sizeof(int *));

  int node_id = 0;
  for (int hit_id=0; hit_id<num_hits; hit_id++) {

    Graph->TH_HitDomToNodeID[hit_id] = malloc(max_doms * sizeof(int));
    
    for (int dom_id=0; dom_id<max_doms; dom_id++) {

      Graph->TH_HitDomToNodeID[hit_id][dom_id] = 0;

      if (dom_id < Graph->TopHits->hit[hit_id]->ndom) {

        Graph->TH_HitDomToNodeID[hit_id][dom_id] = ++node_id;
        
        Graph->Nodes[node_id] = InitSpliceNode(Graph,node_id,hit_id,dom_id,0);

      }

    }

  }
  Graph->num_nodes = node_id;


  // Are your hits to the reverse complement of the query nucleotide seq?
  Graph->revcomp = 0;
  if (Graph->TopHits->N) {
    P7_DOMAIN * Dom = &(Graph->TopHits->hit[Graph->Nodes[1]->hit_id]->dcl[Graph->Nodes[1]->dom_id]);
    if (Dom->ad->sqfrom > Dom->ad->sqto)
      Graph->revcomp = 1;
  }


  for (int edge_id=0; edge_id<num_splice_edges; edge_id++) {
    if (SpliceEdges[edge_id] != NULL)
      ConnectNodesByEdge(SpliceEdges[edge_id],Graph);
  }


  EvaluatePaths(Graph);


  if (DEBUGGING) DEBUG_OUT("'FillOutGraphStructure' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindBestFullPath
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void FindBestFullPath
(SPLICE_GRAPH * Graph)
{
  
  if (DEBUGGING) DEBUG_OUT("Starting 'FindBestFullPath'",1);


  // NOTE: Because NTermNodeIDs is sorted by best_path_score,
  //       we can break the first time we find a full path.
  SPLICE_NODE * Walker;
  for (int i=0; i<Graph->num_n_term; i++) {

    Walker = Graph->Nodes[Graph->NTermNodeIDs[i]];

    Graph->best_full_path_length = 1;
    while (Walker->num_out_edges) {
      Walker = Walker->DownstreamNodes[Walker->best_out_edge];
      Graph->best_full_path_length += 1;
    }

    if (Walker->is_c_terminal) {
      Graph->has_full_path        = 1;
      Graph->best_full_path_start = Graph->NTermNodeIDs[i];
      Graph->best_full_path_end   = Walker->node_id;
      Graph->best_full_path_score = Walker->best_path_score;
      if (DEBUGGING) DEBUG_OUT("'FindBestFullPath' Complete",-1);
      return;
    } else {
      // RESET!
      Graph->best_full_path_length = 0;
    }

  }

  if (DEBUGGING) DEBUG_OUT("'FindBestFullPath' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: BuildSpliceGraph
 *
 *  Desc. :
 *
 *  Inputs:  1.       TopHits :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.            gm : The straightforward profile for the protein / family.
 *           4.            om : The optimized profile (assumed to be built on 'gm').
 *           5.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
SPLICE_GRAPH * BuildSpliceGraph
(
  P7_TOPHITS  * TopHits, 
  TARGET_SEQ  * TargetNuclSeq,
  P7_PROFILE  * gm,
  P7_OPROFILE * om,
  ESL_GENCODE * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'BuildSpliceGraph'",1);

  // We'll just make an unordered list of our splice edges for now
  int num_splice_edges = 0;
  DOMAIN_OVERLAP ** SpliceEdges = GatherViableSpliceEdges(TopHits,TargetNuclSeq,gm,gcode,&num_splice_edges);


  SPLICE_GRAPH * Graph = (SPLICE_GRAPH *)malloc(sizeof(SPLICE_GRAPH));


  Graph->TopHits    = TopHits;
  Graph->MissedHits = NULL;
  
  Graph->Model  = gm;
  Graph->OModel = om;

  Graph->TH_HitDomToNodeID = NULL;
  Graph->MH_HitToNodeID    = NULL;


  Graph->Nodes        = NULL;
  Graph->CumScoreSort = NULL;
  Graph->NTermNodeIDs = NULL;
  Graph->CTermNodeIDs = NULL;


  // Initialize our basic metadata
  Graph->num_nodes  = 0;
  Graph->num_edges  = 0;
  Graph->num_n_term = 0;
  Graph->num_c_term = 0;


  // Eventually, we'll need to know if we have a full
  // path through this graph.
  Graph->has_full_path        = 0;
  Graph->best_full_path_start = 0;
  Graph->best_full_path_end   = 0;
  Graph->best_full_path_score = 0.0;


  // Build that stinky graph!
  FillOutGraphStructure(Graph,SpliceEdges,num_splice_edges);
  FindBestFullPath(Graph);


  free(SpliceEdges);

  if (DEBUGGING) DEBUG_OUT("'BuildSpliceGraph' Complete",-1);

  return Graph;

}













/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ExtractSubProfile
 *
 *  Desc. :
 *
 *  Inputs:  1.     FullModel :
 *           2. hmm_start_pos :
 *           3.   hmm_end_pos :
 *
 *  Output:
 *
 */
P7_PROFILE * ExtractSubProfile
(
  P7_PROFILE * FullModel, 
  int hmm_start_pos, 
  int hmm_end_pos
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'ExtractSubProfile'",1);


  int fullM = FullModel->M;
  int  subM = 1 + hmm_end_pos - hmm_start_pos;
  P7_PROFILE * SubModel = p7_profile_Create(subM,FullModel->abc);


  // 1. TRANSITION SCORES
  //
  for (int trans_type_id = 0; trans_type_id < p7P_NTRANS; trans_type_id++) {

    int  sub_trans_base =  subM * trans_type_id;
    int full_trans_base = fullM * trans_type_id;

    // Position 0, always playing pranks
    SubModel->tsc[sub_trans_base] = FullModel->tsc[full_trans_base];

    full_trans_base += hmm_start_pos-1; // Minus 1 because we're doing 1-indexing in the loop

    for (int sub_model_pos = 1; sub_model_pos <= subM; sub_model_pos++)
      SubModel->tsc[sub_trans_base+sub_model_pos] = FullModel->tsc[full_trans_base+sub_model_pos];

  }


  // 2. EMISSION SCORES
  //
  for (int nr = 0; nr < p7P_NR; nr++) {
    for (int residue_id = 0; residue_id < FullModel->abc->Kp; residue_id++) {

      // Position 0 is a special little baby
      SubModel->rsc[residue_id][0+nr] = FullModel->rsc[residue_id][0+nr];
    
      for (int sub_model_pos = 1; sub_model_pos <= subM; sub_model_pos++)
        SubModel->rsc[residue_id][p7P_NR * sub_model_pos + nr] = FullModel->rsc[residue_id][p7P_NR * (sub_model_pos + hmm_start_pos-1) + nr];
    
    }
  }


  // 3. SPECIAL STATES
  //
  for (int i=0; i<p7P_NXSTATES; i++) {
    for (int j=0; j<p7P_NXTRANS; j++)
      SubModel->xsc[i][j] = FullModel->xsc[i][j];
  }


  // 4. CONSENSUS SEQUENCE
  //
  SubModel->consensus[0] = FullModel->consensus[0];
  for (int sub_model_pos = 1; sub_model_pos <= subM; sub_model_pos++)
    SubModel->consensus[sub_model_pos] = FullModel->consensus[(sub_model_pos-1)+hmm_start_pos];


  // 5. The REST!
  //
  SubModel->mode       = FullModel->mode;
  SubModel->M          = subM;
  SubModel->L          = FullModel->L;
  SubModel->max_length = FullModel->max_length;
  SubModel->nj         = FullModel->nj;
  SubModel->roff       = FullModel->roff;
  SubModel->eoff       = FullModel->eoff;
  for (int i=0; i< p7_NOFFSETS; i++) SubModel->offs[i]    = FullModel->offs[i];
  for (int i=0; i< p7_NEVPARAM; i++) SubModel->evparam[i] = FullModel->evparam[i];
  for (int i=0; i< p7_NCUTOFFS; i++) SubModel->cutoff[i]  = FullModel->cutoff[i];
  for (int i=0; i< p7_MAXABET;  i++) SubModel->compo[i]   = FullModel->compo[i];



  if (DEBUGGING && hmm_start_pos == 1 && hmm_end_pos == FullModel->M) {
    fprintf(stderr,"\n  Full Model Copy Validation:  ");
    if (p7_profile_Compare(SubModel,FullModel,0.0) != eslOK) { fprintf(stderr,"Failed\n\n");   }
    else                                                     { fprintf(stderr,"Success!\n\n"); }
    fflush(stderr);
  }


  if (DEBUGGING) DEBUG_OUT("'ExtractSubProfile' Complete",-1);


  return SubModel;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: NodesAreDCCCompatible
 *
 *  Desc. :
 *
 *  Inputs:  1.      upstream_hmm_to :
 *           2.     upstream_nucl_to :
 *           3.  downstream_hmm_from :
 *           4. downstream_nucl_from :
 *
 *  Output:
 *
 */
int NodesAreDCCCompatible
(
  SPLICE_GRAPH * Graph,
  SPLICE_NODE  * UpstreamNode,
  SPLICE_NODE  * DownstreamNode
)
{


  // How large of an area are we willing to consider searching?
  int max_hmm_dist  =    50;
  int max_nucl_dist = 50000;


  // Pull the relevant data from the upstream and downstream nodes
  P7_DOMAIN * USDom     = &(Graph->TopHits->hit[UpstreamNode->hit_id]->dcl[UpstreamNode->dom_id]);
  int upstream_hmm_to   = USDom->ad->hmmto;
  int upstream_nucl_to  = USDom->ad->sqto;


  P7_DOMAIN * DSDom        = &(Graph->TopHits->hit[DownstreamNode->hit_id]->dcl[DownstreamNode->dom_id]);
  int downstream_hmm_from  = DSDom->ad->hmmfrom;
  int downstream_nucl_from = DSDom->ad->sqfrom;


  // Is the "upstream" hit even really upstream?
  int hmm_dist = 1 + downstream_hmm_from - upstream_hmm_to;
  if (hmm_dist <= 0)
    return 0;


  // I mean *REALLY* upstream?
  int nucl_dist = 1;
  if (Graph->revcomp) nucl_dist += upstream_nucl_to - downstream_nucl_from;
  else                nucl_dist += downstream_nucl_from - upstream_nucl_to;
  if (nucl_dist <= 0)
    return 0;


  // Hmmm, I suppose you're oriented correctly...
  // But are you close enough?!
  if (hmm_dist <= max_hmm_dist && nucl_dist <= max_nucl_dist)
    return 1;


  // So close, but not quite ready to work together!
  return 0;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetBoundedSearchRegions
 *
 *  Desc. :
 *
 *  Inputs:
 *
 *  Output:
 *
 */
int * GetBoundedSearchRegions
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetBoundedSearchRegions'",1);


  // Who all doesn't have an incoming / outgoing edge?
  // (not including "genuine" terminal nodes, of course!)
  int * NoOutEdgeNodes = malloc(Graph->num_nodes * sizeof(int));
  int * NoInEdgeNodes  = malloc(Graph->num_nodes * sizeof(int));
  int num_no_out_edge  = 0;
  int num_no_in_edge   = 0;


  for (int node_id = 1; node_id < Graph->num_nodes; node_id++) {

    SPLICE_NODE * Node = Graph->Nodes[node_id];

    if (Node->num_out_edges == 0 && Node->is_c_terminal == 0)
      NoOutEdgeNodes[num_no_out_edge++] = node_id;

    if (Node->num_in_edges == 0 && Node->is_n_terminal == 0)
      NoInEdgeNodes[num_no_in_edge++] = node_id;

  }


  // This could *surely* be optimized, but I think it's going to
  // be reasonably fast to do an all-versus-all sorta thang...
  int * DisConnCompOuts = malloc(num_no_out_edge * sizeof(int));
  int * DisConnCompIns  = malloc(num_no_in_edge  * sizeof(int));
  for (int i=0; i<num_no_out_edge; i++) DisConnCompOuts[i] = 0;
  for (int i=0; i<num_no_in_edge ; i++) DisConnCompIns[i]  = 0;


  int dcc_ids_issued   = 0; // Note that this is more like "num_dcc_ids issued"
  int num_live_dcc_ids = 0;
  int * LiveDCCIDs = malloc(num_no_out_edge * sizeof(int));
  for (int i=0; i<num_no_out_edge; i++) {


    // Grab the next node without any outgoing edges
    SPLICE_NODE * UpstreamNode = Graph->Nodes[NoOutEdgeNodes[i]];

    
    int dcc_id    = dcc_ids_issued+1; // What's our component ID?
    int connected = 0;                // Did we connect to anyone?

    
    // Scan through the nodes without incoming edges and
    // see who's available for friendship
    for (int j=0; j<num_no_in_edge; j++) {


      // Has this node already heard the good news?
      if (DisConnCompIns[j] == dcc_id)
        continue;


      SPLICE_NODE * DownstreamNode = Graph->Nodes[NoInEdgeNodes[j]];


      if (!NodesAreDCCCompatible(Graph,UpstreamNode,DownstreamNode));
        continue;


      // Compatibility achieved!
      // These might be redundantly set multiple times, but that's life!
      DisConnCompOuts[i] = dcc_id;
      connected = 1;


      // This DCC ID is live, baby!
      LiveDCCIDs[num_live_dcc_ids++] = dcc_id;


      // Has this downstream node already made friends?
      // If so, they're our friends now!
      //
      // This is the lazy way to do this comparison and adjustment,
      // but it's 
      //
      if (DisConnCompIns[j]) {

        int dcc_id_to_replace = DisConnCompIns[j];

        for (int x=0; x<num_no_out_edge; x++) {
          if (DisConnCompOuts[x] == dcc_id_to_replace)
            DisConnCompOuts[x] = dcc_id;
        }

        for (int x=0; x<num_no_in_edge; x++) {
          if (DisConnCompIns[x] == dcc_id_to_replace) {
            DisConnCompIns[x] = dcc_id;
          }
        }

        // The replaced ID is no longer alive :'(
        num_live_dcc_ids--;
        for (int x=0; x<num_live_dcc_ids; x++) {
          if (LiveDCCIDs[x] == dcc_id_to_replace) {
            LiveDCCIDs[x] = dcc_id;
          }
        }


      } else {

        // Just the one friend... for now!
        DisConnCompIns[j] = dcc_id;
      
      }

    }


    // Did we register a new dcc_id?
    //
    // NOTE that it's theoretically possible DCCs could be merged,
    //   so we could end up with an overcount of the number of DCCs.
    //
    if (connected)
      dcc_ids_issued++;

  }


  // Now we can go through all of the hits in each of our
  // "disconnected components" and determine a maximal 
  // search area
  int * SearchRegionAggregate = malloc((1 + 4 * num_live_dcc_ids) * sizeof(int));
  SearchRegionAggregate[0] = num_live_dcc_ids;


  // Once again, there is plenty of room for optimization,
  // but I really doubt this is going to be a bottleneck
  // in Splash...
  P7_DOMAIN * DomPtr;
  for (int meta_dcc_id = 0; meta_dcc_id < num_live_dcc_ids; meta_dcc_id++) {


    int dcc_id = LiveDCCIDs[meta_dcc_id];


    // 
    // 1. What should we set as the bounds for our search region, based
    //    on the nodes without outgoing edges?
    //

    int dcc_hmm_start  = -1;
    int dcc_nucl_start = -1;

    for (int i=0; i<num_no_out_edge; i++) {

      if (DisConnCompOuts[i] != dcc_id)
        continue;

      int node_id = NoOutEdgeNodes[i];
      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[Graph->Nodes[node_id]->dom_id]);


      // NOTE that the START of the search region will be defined by the
      // minimal END position in the model.
      // It looks like these variable names are at odds, but they aren't,
      // I swear!
      if (dcc_hmm_start == -1 || dcc_hmm_start > DomPtr->ad->hmmto)
        dcc_hmm_start = DomPtr->ad->hmmto;


      if (Graph->revcomp) {
        if (dcc_nucl_start == -1 || dcc_nucl_start < DomPtr->ad->sqto)
          dcc_nucl_start = DomPtr->ad->sqto;
      } else {
        if (dcc_nucl_start == -1 || dcc_nucl_start > DomPtr->ad->sqto)
          dcc_nucl_start = DomPtr->ad->sqto;
      }

    }



    // 
    // 2. What should we set as the bounds for our search region, based
    //    on the nodes without incoming edges?
    //

    int dcc_hmm_end  = -1;
    int dcc_nucl_end = -1;

    for (int i=0; i<num_no_in_edge; i++) {

      if (DisConnCompIns[i] != dcc_id)
        continue;

      int node_id = NoInEdgeNodes[i];
      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[Graph->Nodes[node_id]->dom_id]);


      if (dcc_hmm_end == -1 || dcc_hmm_end < DomPtr->ad->hmmfrom)
        dcc_hmm_end = DomPtr->ad->hmmfrom;


      if (Graph->revcomp) {
        if (dcc_nucl_end == -1 || dcc_nucl_end > DomPtr->ad->sqfrom)
          dcc_nucl_end = DomPtr->ad->sqfrom;
      } else {
        if (dcc_nucl_end == -1 || dcc_nucl_end < DomPtr->ad->sqfrom)
          dcc_nucl_end = DomPtr->ad->sqfrom;
      }

    }



    //
    // 3. Log it!
    //
    SearchRegionAggregate[4*meta_dcc_id + 1] = dcc_hmm_start;
    SearchRegionAggregate[4*meta_dcc_id + 2] = dcc_hmm_end;
    SearchRegionAggregate[4*meta_dcc_id + 3] = dcc_nucl_start;
    SearchRegionAggregate[4*meta_dcc_id + 4] = dcc_nucl_end;

  }



  free(NoOutEdgeNodes);
  free(NoInEdgeNodes);
  free(DisConnCompOuts);
  free(DisConnCompIns);
  free(LiveDCCIDs);


  if (DEBUGGING) DEBUG_OUT("'GetBoundedSearchRegions' Complete",-1);


  return SearchRegionAggregate;

}














/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SelectFinalSubHits
 *
 *  Desc. :
 *
 *  Inputs:  1.          SubHitADs :
 *           2.       SubHitScores :
 *           3.       num_sub_hits :
 *           4.        sub_hmm_len :
 *           5.          hmm_start :
 *           6. final_num_sub_hits :
 *
 *  Output:
 *
 */
P7_DOMAIN ** SelectFinalSubHits
(
  P7_ALIDISPLAY ** SubHitADs,
  float * SubHitScores,
  int     num_sub_hits,
  int     sub_hmm_len,
  int     hmm_start,
  int   * final_num_sub_hits
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SelectFinalSubHits'",1);


  // Get a sorting of the hits by their scores
  // NOTE that these are the uncorrected scores
  //  so they aren't useful for external purposes,
  //  but can still rank our hits *internally*
  int * ADSort = FloatHighLowSortIndex(SubHitScores,num_sub_hits);


  // Now that we have all of our hits, let's find which
  // one(s) will give us coverage of the part of the pHMM
  // we need covered.
  int * Coverage = malloc(sub_hmm_len*sizeof(int));
  for (int i=0; i<sub_hmm_len; i++) 
    Coverage[i] = 0;


  *final_num_sub_hits = 0;
  for (int sort_id=0; sort_id<num_sub_hits; sort_id++) {

        
    int sub_hit_id = ADSort[sort_id];
    P7_ALIDISPLAY * AD = SubHitADs[sub_hit_id];


    // Do we get any fresh coverage out of this hit?
    int added_coverage = 0;
    for (int i = AD->hmmfrom - hmm_start; i <= AD->hmmto - hmm_start; i++) {
      if (Coverage[i] == 0) {
        Coverage[i] = 1;
        added_coverage++;
      }
    }


    // No coverage?! Away with you, filth!
    if (!added_coverage) {
      p7_alidisplay_Destroy(AD);
      SubHitADs[sub_hit_id] = NULL;
      continue;
    }


    // THERE WE GO!
    *final_num_sub_hits += 1;


  }
  free(Coverage);
  free(ADSort);



  P7_DOMAIN ** FinalSubHits = malloc(*final_num_sub_hits * sizeof(P7_DOMAIN *));
  *final_num_sub_hits = 0;
  for (int sub_hit_id = 0; sub_hit_id < num_sub_hits; sub_hit_id++) {

    if (SubHitADs[sub_hit_id] == NULL)
      continue;

    FinalSubHits[*final_num_sub_hits]           = p7_domain_Create_empty();
    FinalSubHits[*final_num_sub_hits]->ad       = SubHitADs[sub_hit_id];
    FinalSubHits[*final_num_sub_hits]->bitscore = SubHitScores[sub_hit_id];

    *final_num_sub_hits += 1;

  }



  if (DEBUGGING) DEBUG_OUT("'SelectFinalSubHits' Complete",-1);


  return FinalSubHits;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindSubHits
 *
 *  Desc. :
 *
 *  Inputs:  1.              Graph :
 *           2.      TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.       SearchRegion :
 *           4.              gcode : An ESL_GENCODE struct (mainly used for translation).
 *           5. final_num_sub_hits :
 *
 *  Output:
 *
 *  NOTE: This could easily be improved by someone with a keener
 *        understanding of the HMMER internals, but for the purposes
 *        of plugging holes in a splice graph I think this is
 *        adequate... maybe...
 *
 */
P7_DOMAIN ** FindSubHits
(
  SPLICE_GRAPH * Graph,
  TARGET_SEQ   * TargetNuclSeq,
  int          * SearchRegion,
  ESL_GENCODE  * gcode,
  int          * final_num_sub_hits
)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'FindSubHits'",1);


  int hmm_start  = SearchRegion[0];
  int hmm_end    = SearchRegion[1];
  int nucl_start = SearchRegion[2];
  int nucl_end   = SearchRegion[3];


  int sub_hmm_len  = 1 + hmm_end - hmm_start;
  int min_exon_len = 5;


  // Generate a sub-model and an optimized sub-model
  // for the part of the pHMM we need to fill in
  P7_PROFILE  *  SubModel = ExtractSubProfile(Graph->Model,hmm_start,hmm_end);
  P7_OPROFILE * OSubModel = p7_oprofile_Create(SubModel->M,SubModel->abc);
  int submodel_create_err = p7_oprofile_Convert(SubModel,OSubModel);

  // Required for alidisplay generation
  OSubModel->name = malloc(9*sizeof(char));
  strcpy(OSubModel->name,"SubModel");


  // Grab the nucleotides we're searching our sub-model against
  ESL_DSQ * SubNucls = GrabNuclRange(TargetNuclSeq,nucl_start,nucl_end);
  int nucl_seq_len   = abs(nucl_end-nucl_start)+1;


  // Create a whole mess of objects that we'll need to get our
  // our sub-model alignments to fit into an alidisplay...
  ESL_SQ   * ORFAminos      = esl_sq_Create();
  ESL_SQ   * ORFNucls       = esl_sq_Create();
  P7_GMX   * ViterbiMatrix  = p7_gmx_Create(SubModel->M,1024);
  P7_TRACE * Trace          = p7_trace_Create();
  float viterbi_score;


  // Finally, the array where we'll record all of our incredible triumphs!
  int num_sub_hits =  0;
  int max_sub_hits = 10;
  P7_ALIDISPLAY ** SubHitADs = malloc(max_sub_hits * sizeof(P7_ALIDISPLAY *));
  float * SubHitScores = malloc(max_sub_hits * sizeof(float));


  // Loop over each full reading frame
  for (int frame=0; frame<3; frame++) {

    int orf_len = 0;
    int frame_start = 1 + frame;


    // As we walk through the reading frame,
    // we'll build up ORFs and search them
    // against our sub-model
    for (int frame_end = frame_start; frame_end+2 <= nucl_seq_len; frame_end += 3) {


      // Translate and add to the current frame
      int next_amino_index = esl_gencode_GetTranslation(gcode,&(SubNucls[frame_end]));
      esl_sq_CAddResidue(ORFAminos,AMINO_CHARS[next_amino_index]);

      esl_sq_CAddResidue(ORFNucls,DNA_CHARS[SubNucls[frame_end  ]]);
      esl_sq_CAddResidue(ORFNucls,DNA_CHARS[SubNucls[frame_end+1]]);
      esl_sq_CAddResidue(ORFNucls,DNA_CHARS[SubNucls[frame_end+2]]);

      orf_len++;


      // Have we hit a stop codon? Are we at the end of the road?
      if (next_amino_index > 21 || frame_end+5 >= nucl_seq_len) {

        
        // No point throwing a stop codon into the mix...
        if (next_amino_index > 21) 
          orf_len--;



        // Is this ORF long enough to be worth considering as an exon?
        if (orf_len >= min_exon_len) {


          int dsq_err_code  = esl_sq_Digitize(SubModel->abc,ORFAminos);
          if (dsq_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while digitizing ORF amino sequence\n\n");
          }



          int vit_err_code  = p7_GViterbi(ORFAminos->dsq,orf_len,SubModel,ViterbiMatrix,&viterbi_score);
          if (vit_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while running Viterbi\n\n");
          }
          


          p7_trace_Reuse(Trace);
          int gtrace_err_code  = p7_GTrace(ORFAminos->dsq,orf_len,SubModel,ViterbiMatrix,Trace);
          if (gtrace_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while generating a generic P7_TRACE for the Viterbi matrix\n\n");
          }


          p7_trace_Index(Trace);
          for (int trace_dom = 0; trace_dom < Trace->ndom; trace_dom++) {


            P7_ALIDISPLAY * AD = p7_alidisplay_Create(Trace,0,OSubModel,ORFAminos,NULL); // Currently the translated version isn't playing nicely...
            if (AD == NULL) {
              fprintf(stderr,"\n  ERROR (FindSubHits): Failed while generating P7_ALIDISPLAY\n\n");
            }


            if (AD->hmmto - AD->hmmfrom < min_exon_len)
              continue;


            // Oh, boy! Let's add this alidisplay to our array!


            // (but first... resize?)
            if (num_sub_hits == max_sub_hits) {
              max_sub_hits *= 2;
              P7_ALIDISPLAY ** NewSubHitADs = malloc(max_sub_hits*sizeof(P7_ALIDISPLAY));
              float * NewSubHitScores = malloc(max_sub_hits*sizeof(float));
              for (int sub_hit_id=0; sub_hit_id<num_sub_hits; sub_hit_id++) {
                NewSubHitADs[sub_hit_id] = SubHitADs[sub_hit_id];
                NewSubHitScores[sub_hit_id] = SubHitScores[sub_hit_id];
              }
              free(SubHitADs);
              free(SubHitScores);
              SubHitADs = NewSubHitADs;
              SubHitScores = NewSubHitScores;
            }


            AD->hmmfrom += hmm_start - 1;
            AD->hmmto   += hmm_start - 1;


            // Before we adjust the 'sqfrom' and 'sqto' values to
            // represent genomic coordinates, let's make sure we
            // have the nucleotide sequence on-hand
            int ntseq_len = 3 * (1 + AD->sqto - AD->sqfrom);
            AD->ntseq  = malloc(ntseq_len * sizeof(char));
            for (int i=0; i<ntseq_len; i++)
              AD->ntseq[i] = ORFNucls->seq[3 * (AD->sqfrom - 1) + i];


            // We'll need to do some quick math to figure out exactly
            // where this hit is in terms of nucleotide coordinates
            if (Graph->revcomp) {
              AD->sqfrom = (nucl_start + 1) - frame_start - 3*(AD->sqfrom - 1);
              AD->sqto   = (nucl_start + 1) - frame_start - 3*(AD->sqto   - 1) - 2;
            } else {
              AD->sqfrom = (nucl_start - 1) + frame_start + 3*(AD->sqfrom - 1);
              AD->sqto   = (nucl_start - 1) + frame_start + 3*(AD->sqto   - 1) + 2;
            }


            // As one last thing, we'll set the score to more closely
            // correspond to other scores.  NOTE that this is not the
            // true score (we'd need to have a null model score)!
            // Instead (like I've done elsewhere...) we're approximating
            // to get the methods fully sketched, after which we can refine.
            viterbi_score = (1 + AD->hmmto - AD->hmmfrom) * (float)exp(viterbi_score);


            // Welcome to the list, fella!
            SubHitADs[num_sub_hits]    = AD;
            SubHitScores[num_sub_hits] = viterbi_score;
            num_sub_hits++;


          }

        }


        // Reset and advance!
        esl_sq_Reuse(ORFAminos);
        esl_sq_Textize(ORFAminos);
        esl_sq_Reuse(ORFNucls);
        frame_start = frame_end + 3;
        orf_len = 0;


      }


    }

  }


  // It takes a lot of fun to need this much cleanup ;)
  esl_sq_Destroy(ORFAminos);
  esl_sq_Destroy(ORFNucls);
  free(SubNucls);
  p7_profile_Destroy(SubModel);
  p7_oprofile_Destroy(OSubModel);
  p7_gmx_Destroy(ViterbiMatrix);
  p7_trace_Destroy(Trace);


  // Quick check: Did we find *anything* worth considering?
  if (num_sub_hits == 0) {
    free(SubHitADs);
    free(SubHitScores);
    if (DEBUGGING) DEBUG_OUT("'FindSubHits' Complete (albeit, no hits)",-1);
    return NULL;
  }


  // Reduce down to just the hits we're really excited about
  P7_DOMAIN ** FinalSubHits = SelectFinalSubHits(SubHitADs,SubHitScores,num_sub_hits,sub_hmm_len,hmm_start,final_num_sub_hits);


  free(SubHitADs);
  free(SubHitScores);


  if (DEBUGGING) DEBUG_OUT("'FindSubHits' Complete",-1);


  // Happy days!
  return FinalSubHits;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: IntegrateMissedHits
 *
 *  Desc. :
 *
 *  Inputs:  1.          Graph :
 *           2. NewSpliceEdges :
 *           3.  num_new_edges :
 *
 *  Output:
 *
 */
void IntegrateMissedHits
(SPLICE_GRAPH * Graph, DOMAIN_OVERLAP ** NewSpliceEdges, int num_new_edges)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'IntegrateMissedHits'",1);


  // First off, we need to copy over all of the existing nodes
  int new_num_nodes = Graph->num_nodes + Graph->MissedHits->N;

  SPLICE_NODE ** NewNodeArray = (SPLICE_NODE **)malloc((new_num_nodes+1)*sizeof(SPLICE_NODE *));

  int node_id;  
  for (node_id = 1; node_id <= Graph->num_nodes; node_id++)
    NewNodeArray[node_id] = Graph->Nodes[node_id];
  for (node_id = Graph->num_nodes+1; node_id <= new_num_nodes; node_id++)
    NewNodeArray[node_id] = NULL;

  free(Graph->Nodes);
  Graph->Nodes = NewNodeArray;


  // Now we can actually integrate the new hits!
  node_id = Graph->num_nodes;
  Graph->num_nodes = new_num_nodes;

  Graph->MH_HitToNodeID = malloc(Graph->MissedHits->N * sizeof(int));
  for (int missed_hit_id = 0; missed_hit_id < Graph->MissedHits->N; missed_hit_id++) {

    Graph->MH_HitToNodeID[missed_hit_id] = ++node_id;

    Graph->Nodes[node_id] = InitSpliceNode(Graph,node_id,missed_hit_id,0,1);

  }


  for (int new_edge_id=0; new_edge_id<num_new_edges; new_edge_id++) {
    if (NewSpliceEdges[new_edge_id] != NULL)
      ConnectNodesByEdge(NewSpliceEdges[new_edge_id],Graph);
  }


  if (DEBUGGING) DEBUG_OUT("'IntegrateMissedHits' Complete",-1);

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SeekMissingExons
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
P7_TOPHITS * SeekMissingExons
(
  SPLICE_GRAPH * Graph, 
  TARGET_SEQ   * TargetNuclSeq, 
  ESL_GENCODE  * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SeekMissingExons'",1);


  // Grab the sub-regions of our conceptual DP zone that we want
  // to search for missing exons.
  //
  // This 'SearchRegionAggregate' is indexed in groups consisting
  // of: the starting position in the model, the ending position in
  // the model, the starting coordinate on the genome, and the ending
  // coordinate on the genome.
  //
  int * SearchRegionAggregate = GetBoundedSearchRegions(Graph);
  int   num_search_regions    = SearchRegionAggregate[0];

  if (SearchRegionAggregate == NULL) {
    if (DEBUGGING) DEBUG_OUT("'SeekMissingExons' Complete (none found)",-1);
    return NULL;
  }


  if (DEBUGGING) {
    fprintf(stderr,"\n  Num Search Regions: %d\n",num_search_regions);
    for (int i=0; i<num_search_regions; i++) {
      fprintf(stderr,"    - Search Region : %d\n",i+1);
      fprintf(stderr,"      HMM Positions : %d..%d\n",SearchRegionAggregate[i*4 + 1],SearchRegionAggregate[i*4 + 2]);
      fprintf(stderr,"      Nucl. Coord.s : %d..%d\n",SearchRegionAggregate[i*4 + 3],SearchRegionAggregate[i*4 + 4]);
    }
    fprintf(stderr,"\n");
    DEBUG_OUT("'SeekMissingExons' Complete (early kill)",-1);
    free(SearchRegionAggregate);
    return NULL;
  }


  // Now we can iterate over our list of search regions and,
  // for each:
  //
  //   1. Extract the appropriate sub-region of the model
  //   2. Pull the appropriate portion of the target sequence
  //   3. Search the sub-model against the sub-target,
  //        using the Viterbi algorithm
  //
  //   and...
  //
  //   4. Integrate any new hits into the graph!
  //
  int num_sub_hits = 0;
  int sub_hits_capacity = 20;
  P7_DOMAIN ** SubHits = malloc(sub_hits_capacity*sizeof(P7_DOMAIN *));
  for (int search_region_id = 0; search_region_id < num_search_regions; search_region_id++) {


    int num_new_sub_hits;
    P7_DOMAIN ** NewSubHits = FindSubHits(Graph,TargetNuclSeq,&SearchRegionAggregate[(4*search_region_id)+1],gcode,&num_new_sub_hits);


    if (num_new_sub_hits == 0)
      continue;


    // New sub-hit(s) alert!


    // Resize?
    if (num_sub_hits + num_new_sub_hits > sub_hits_capacity) {
      sub_hits_capacity = intMax(sub_hits_capacity*2,num_sub_hits+num_new_sub_hits);
      P7_DOMAIN ** MoreSubHits = malloc(sub_hits_capacity * sizeof(P7_DOMAIN *));
      for (int sub_hit_id = 0; sub_hit_id < num_sub_hits; sub_hit_id++)
        MoreSubHits[sub_hit_id] = SubHits[sub_hit_id];
      free(SubHits);
      SubHits = MoreSubHits;
    }


    // Pop those shrimps on the barbie, as they say up there in Argentina
    for (int new_sub_hit_id = 0; new_sub_hit_id < num_new_sub_hits; new_sub_hit_id++)
      SubHits[num_sub_hits++] = NewSubHits[new_sub_hit_id];


    free(NewSubHits); // clear the pointer

  }
  free(SearchRegionAggregate);


  if (num_sub_hits == 0) {
    free(SubHits);
    return NULL;
  }


  //
  //  NOTE:  From here, what we do is aggregate our new 'SubHits' into
  //         a 'P7_TOPHITS' datastructure.
  //
  //         IMPORTANTLY, we're faking these in order to take advantage
  //         of some of the functions that are already available.
  //         DO NOT ASSUME HMMER FUNCTIONS WILL WORK WITH THIS DATASTRUCTURE!
  //
  P7_TOPHITS * MissingHits = p7_tophits_Create();

  for (int sub_hit_id=0; sub_hit_id<num_sub_hits; sub_hit_id++) {

    P7_HIT * NewHit;
    int hit_create_err  = p7_tophits_CreateNextHit(MissingHits,&NewHit);
    if (hit_create_err != eslOK) {
      fprintf(stderr,"\n  ERROR (SeekMissingExons): Failed while attempting P7_HIT allocation\n\n");
    }

    NewHit->ndom  = 1;
    NewHit->dcl   = SubHits[sub_hit_id];
    NewHit->score = SubHits[sub_hit_id]->bitscore;

  }
  free(SubHits);

 
  if (DEBUGGING) DEBUG_OUT("'SeekMissingExons' Complete",-1);


  return MissingHits;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: AddMissingExonsToGraph
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void AddMissingExonsToGraph
(
  SPLICE_GRAPH * Graph,
  TARGET_SEQ   * TargetNuclSeq,
  ESL_GENCODE  * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'AddMissingExonsToGraph'",1);

  

  // As I note in 'FindMissingExons' (but is worth emphasizing),
  // this datastructure is basically a cheap way to pass around
  // the new hits that we find.
  //
  // MUCH OF THE STANDARD P7_HIT, P7_DOMAIN, and P7_ALIDISPLAY DATA
  // IS NOT CONTAINED IN THESE!  BE CAREFUL BEFORE USING HMMER INTERNAL
  // FUNCTIONS TO INTERROGATE THEM!
  //
  Graph->MissedHits = SeekMissingExons(Graph,TargetNuclSeq,gcode);

  
  if (Graph->MissedHits == NULL) {
    if (DEBUGGING) DEBUG_OUT("'AddMissingExonsToGraph' Complete",-1);
    return;
  }


  // This is basically the same code from 'GatherViableSpliceEdges,'
  // but reworked to integrate new exons into the splice graph.
  //
  // It probably wouldn't be hard to convert the two to a single function,
  // but the annoying thing is that there's the whole dereferencing thing
  // to get *real* P7_DOMAINs from *actual* P7_HITs (but not my goofy ones)
  //
  int new_splice_edge_cap = Graph->MissedHits->N * 2;
  DOMAIN_OVERLAP ** NewSpliceEdges = (DOMAIN_OVERLAP **)malloc(new_splice_edge_cap*sizeof(DOMAIN_OVERLAP *));

  int num_missing_hits = (int)(Graph->MissedHits->N);
  int num_new_edges    = 0;
  for (int missed_hit_id = 0; missed_hit_id < num_missing_hits; missed_hit_id++) {


    // Lucky us!  Cheap and easy access!
    P7_ALIDISPLAY * MissedAD = Graph->MissedHits->hit[missed_hit_id]->dcl->ad;


    for (int node_id = 1; node_id <= Graph->num_nodes; node_id++) {


      int node_hit_id = Graph->Nodes[node_id]->hit_id;
      int node_dom_id = Graph->Nodes[node_id]->dom_id;

      P7_ALIDISPLAY * NodeAD = (&Graph->TopHits->hit[node_hit_id]->dcl[node_dom_id])->ad;


      // Because order matters for 'HitsAreSpliceCompatible' we
      // need to have a catch for either possibility.
      if (HitsAreSpliceCompatible(NodeAD,MissedAD)) {

        NewSpliceEdges[num_new_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
        DOMAIN_OVERLAP * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = node_hit_id;
        Edge->upstream_dom_id   = node_dom_id;
        Edge->downstream_hit_id = missed_hit_id;
        Edge->downstream_dom_id = 0;

        Edge->UpstreamTopHits   = Graph->TopHits;
        Edge->UpstreamDisplay   = NodeAD;
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MissedAD;

        num_new_edges++;

      } else if (HitsAreSpliceCompatible(MissedAD,NodeAD)) {

        NewSpliceEdges[num_new_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
        DOMAIN_OVERLAP * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = missed_hit_id;
        Edge->upstream_dom_id   = 0;
        Edge->downstream_hit_id = node_hit_id;
        Edge->downstream_dom_id = node_dom_id;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MissedAD;
        Edge->DownstreamTopHits = Graph->TopHits;
        Edge->DownstreamDisplay = NodeAD;

        num_new_edges++;

      }


      // A little wasteful, but better than writing this twice, am I right?!
      if (num_new_edges == new_splice_edge_cap) {

        new_splice_edge_cap *= 2;

        DOMAIN_OVERLAP ** MoreNewEdges = (DOMAIN_OVERLAP **)malloc(new_splice_edge_cap*sizeof(DOMAIN_OVERLAP *));
        for (int i=0; i<num_new_edges; i++)
          MoreNewEdges[i] = NewSpliceEdges[i];

        free(NewSpliceEdges);
        NewSpliceEdges = MoreNewEdges;

      }

    }

  }



  // Back at it again!
  for (int splice_edge_id = 0; splice_edge_id < num_new_edges; splice_edge_id++) {
  
    SketchSpliceEdge(NewSpliceEdges[splice_edge_id],TargetNuclSeq,Graph->Model,gcode);
  
    if (NewSpliceEdges[splice_edge_id]->score == EDGE_FAIL_SCORE) {
      free(NewSpliceEdges[splice_edge_id]);
      NewSpliceEdges[splice_edge_id] = NULL;
    }

  }



  //
  //  If I'm not mistaken.... IT'S PARTY TIME!!!!
  //
  IntegrateMissedHits(Graph,NewSpliceEdges,num_new_edges);


  if (DEBUGGING) DEBUG_OUT("'AddMissingExonsToGraph' Complete",-1);

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetExonSetFromStartNode
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. start_node_id :
 *
 *  Output:
 *
 */
int * GetExonSetFromStartNode
(SPLICE_GRAPH * Graph, int start_node_id)
{
  
  if (DEBUGGING) DEBUG_OUT("Starting 'GetExonSetFromStartNode'",1);


  // In case we use *every* node
  // BestPathCoords[0] is the number of nodes along the
  //  path (Min:1,Max:num_nodes)
  int * BestPathCoords  = malloc((4*Graph->num_nodes + 1)*sizeof(int));
  int num_path_elements = 0;


  // N-Terminal is special
  SPLICE_NODE * Node = Graph->Nodes[start_node_id];
  if (Node->was_missed) {
    BestPathCoords[++num_path_elements] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->sqfrom;
    BestPathCoords[++num_path_elements] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->hmmfrom;
  } else {
    BestPathCoords[++num_path_elements] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->sqfrom;
    BestPathCoords[++num_path_elements] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->hmmfrom;
  }


  while (Node->num_out_edges) {

    BestPathCoords[++num_path_elements] = Node->OutEdges[Node->best_out_edge]->upstream_spliced_nucl_end;
    BestPathCoords[++num_path_elements] = Node->OutEdges[Node->best_out_edge]->upstream_exon_terminus;

    Node = Node->DownstreamNodes[Node->best_out_edge];

    BestPathCoords[++num_path_elements] = Node->InEdges[Node->best_in_edge]->downstream_spliced_nucl_start;
    BestPathCoords[++num_path_elements] = Node->InEdges[Node->best_in_edge]->downstream_exon_terminus;

  }


  // C-Terminal is also special
  if (Node->was_missed) {
    BestPathCoords[++num_path_elements] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->sqto;
    BestPathCoords[++num_path_elements] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->hmmto;
  } else {
    BestPathCoords[++num_path_elements] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->sqto;
    BestPathCoords[++num_path_elements] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->hmmto;
  }


  BestPathCoords[0] = num_path_elements/4;


  if (DEBUGGING) DEBUG_OUT("'GetExonSetFromStartNode' Complete",-1);

  return BestPathCoords;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindComponentBestStart
 *
 *  Desc. :
 *
 *  Inputs:  1.            Node :
 *           2.    ComponentIDs :
 *           3.    component_id :
 *           4. best_comp_start :
 *           5. best_comp_score :
 *
 *  Output:
 *
 */
void FindComponentBestStart
(
  SPLICE_NODE * Node,
  int   * ComponentIDs,
  int     component_id,
  int   * best_comp_start,
  float * best_comp_score
)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'FindComponentBestStart'",1);


  if (ComponentIDs[Node->node_id]) {
    if (DEBUGGING) DEBUG_OUT("'FindComponentBestStart' Complete",-1);
    return;
  }

  ComponentIDs[Node->node_id] = component_id;

  if (Node->num_in_edges == 0 && Node->best_path_score > *best_comp_score) {
    *best_comp_start = Node->node_id;
    *best_comp_score = Node->best_path_score;
  }

  for (int i=0; i<Node->num_in_edges; i++)
    FindComponentBestStart(Node->UpstreamNodes[i],ComponentIDs,component_id,best_comp_start,best_comp_score);

  for (int i=0; i<Node->num_out_edges; i++)
    FindComponentBestStart(Node->DownstreamNodes[i],ComponentIDs,component_id,best_comp_start,best_comp_score);


  if (DEBUGGING) DEBUG_OUT("'FindComponentBestStart' Complete",-1);

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: TranslateExonSetNucls
 *
 *  Desc. :
 *
 *  Inputs:  1.      ExonSetNucls :
 *           2. coding_region_len :
 *           3.             gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
ESL_DSQ * TranslateExonSetNucls
(ESL_DSQ * ExonSetNucls, int coding_region_len, ESL_GENCODE * gcode)
{

  ESL_DSQ * ExonSetTrans = malloc((1 + coding_region_len/3) * sizeof(ESL_DSQ));

  ESL_DSQ * Codon = malloc(3*sizeof(ESL_DSQ));
  int trans_index = 1;

  for (int nucl_index=1; nucl_index<coding_region_len; nucl_index += 3) {
    Codon[0] = ExonSetNucls[nucl_index    ];
    Codon[1] = ExonSetNucls[nucl_index + 1];
    Codon[2] = ExonSetNucls[nucl_index + 2];
    ExonSetTrans[trans_index++] = esl_gencode_GetTranslation(gcode,&(Codon[0]));
  }

  free(Codon);

  return ExonSetTrans;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GrabExonCoordSetNucls
 *
 *  Desc. :
 *
 *  Inputs:  1.      ExonCoordSet :
 *           2.     TargetNuclSeq :
 *           3. coding_region_len :
 *
 *  Output:
 *
 */
ESL_DSQ * GrabExonCoordSetNucls
(int * ExonCoordSet, TARGET_SEQ * TargetNuclSeq, int * coding_region_len)
{

  int num_exons = ExonCoordSet[0];
  int num_nucls = 0;


  for (int exon_id = 0; exon_id < num_exons; exon_id++)
    num_nucls += abs(ExonCoordSet[exon_id*4 + 1] - ExonCoordSet[exon_id*4 + 3]) + 1;

  ESL_DSQ * ExonSetNucls = malloc(num_nucls * sizeof(ESL_DSQ));
  *coding_region_len = num_nucls;

  int nucl_placer = 1;
  for (int exon_id = 0; exon_id < num_exons; exon_id++) {

    int range_start = ExonCoordSet[exon_id*4 + 1];
    int range_end   = ExonCoordSet[exon_id*4 + 3];
    ESL_DSQ * ExonNucls = GrabNuclRange(TargetNuclSeq,range_start,range_end);


    for (int nucl_id = 1; nucl_id <= abs(range_end-range_start)+1; nucl_id++)
      ExonSetNucls[nucl_placer++] = ExonNucls[nucl_id];


    free(ExonNucls);

  }


  return ExonSetNucls;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  DEBUGGING Function: DumpExonSetSequence
 *
 *  Desc. :
 *
 *  Inputs:  1. ExonCoordSets :
 *           2. num_exon_sets :
 *           3. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           4.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void DumpExonSets
(int ** ExonCoordSets, int num_exon_sets, TARGET_SEQ * TargetNuclSeq, ESL_GENCODE * gcode)
{
  

  fprintf(stderr,"\n\n+");
  for (int i=0; i<60; i++)
    fprintf(stderr,"=");
  fprintf(stderr,"+\n\n");


  
  for (int exon_set_id=0; exon_set_id<num_exon_sets; exon_set_id++) {


    int * ExonCoords = ExonCoordSets[exon_set_id];
    int   num_exons  = ExonCoords[0];

    fprintf(stderr,">ExonSet__%d/%d:Nucls__",exon_set_id+1,num_exon_sets);
    for (int i=0; i<num_exons; i++) {
      if (i) fprintf(stderr,",");
      fprintf(stderr,"%d-%d",ExonCoords[i*4 + 1],ExonCoords[i*4 + 3]);
    }


    fprintf(stderr,":Aminos__");
    for (int i=0; i<num_exons; i++) {
      if (i) fprintf(stderr,",");
      fprintf(stderr,"%d-%d",ExonCoords[i*4 + 2],ExonCoords[i*4 + 4]);
    }


    int num_nucls;
    ESL_DSQ * NuclSeq  = GrabExonCoordSetNucls(ExonCoordSets[exon_set_id],TargetNuclSeq,&num_nucls);
    ESL_DSQ * TransSeq = TranslateExonSetNucls(NuclSeq,num_nucls,gcode);
    

    int line_length = 60;
    for (int i=1; i<=num_nucls; i++) {
      if (i % line_length == 1)
        fprintf(stderr,"\n");
      fprintf(stderr,"%c",DNA_CHARS[NuclSeq[i]]);
    }
    fprintf(stderr,"\n");


    fprintf(stderr,">Translation");
    for (int i=1; i<=num_nucls/3; i++) {
      if (i % line_length == 1)
        fprintf(stderr,"\n");
      fprintf(stderr,"%c",AMINO_CHARS[TransSeq[i]]);
    }
    fprintf(stderr,"\n");


    free(NuclSeq);
    free(TransSeq);

    fprintf(stderr,"\n");

  }
  fprintf(stderr,"\n");


  fprintf(stderr,"+");
  for (int i=0; i<60; i++)
    fprintf(stderr,"=");
  fprintf(stderr,"+\n\n\n");

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetSplicedExonCoordSets
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. num_exon_sets :
 *
 *  Output:
 *
 */
int ** GetSplicedExonCoordSets
(SPLICE_GRAPH * Graph, int * num_exon_sets)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetSplicedExonCoordSets'",1);


  // How many connected components do we have?
  // (In the case that we have a full path through
  //  the graph, we just report that one component.)
  int num_conn_comps = 0;
  int * StartNodes = malloc(Graph->num_nodes*sizeof(int));

  if (Graph->has_full_path) {

    StartNodes[0]  = Graph->best_full_path_start;
    num_conn_comps = 1;

  } else {

    int * ComponentIDs = malloc((Graph->num_nodes+1)*sizeof(int));
    for (int node_id=1; node_id<=Graph->num_nodes; node_id++)
      ComponentIDs[node_id] = 0;
    
    for (int node_id=1; node_id<=Graph->num_nodes; node_id++) {

      if (!ComponentIDs[node_id]) {

        int   best_comp_start = 0;
        float best_comp_score = 0.0;
        FindComponentBestStart(Graph->Nodes[node_id],ComponentIDs,num_conn_comps+1,&best_comp_start,&best_comp_score);

        StartNodes[num_conn_comps++] = best_comp_start;

      }

    }

    free(ComponentIDs);

  }


  // Cool!  Now let's grab the spliced coordinates
  int ** ExonCoordSets = malloc(num_conn_comps * sizeof(int *));
  for (int conn_comp_id = 0; conn_comp_id < num_conn_comps; conn_comp_id++)
    ExonCoordSets[conn_comp_id] = GetExonSetFromStartNode(Graph,StartNodes[conn_comp_id]);
  free(StartNodes);


  if (DEBUGGING) DEBUG_OUT("'GetSplicedExonCoordSets' Complete",-1);


  // Easy!
  *num_exon_sets = num_conn_comps;
  return ExonCoordSets;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ReportSplicedTopHits
 *
 *  Desc. :
 *
 *  Inputs:  1.  ExonSetTopHits :
 *           2. ExonSetPipeline :
 *           3.    ExonCoordSet :
 *           4.             ofp : An open file pointer (specificially, the target for output).
 *           5.           textw : The desired line width for output.
 *
 *  Output:
 *
 */
int ReportSplicedTopHits
(
  P7_TOPHITS  * ExonSetTopHits, 
  P7_PIPELINE * ExonSetPipeline, 
  int         * ExonCoordSet,
  int           exon_set_name_id, // just exon_set_id+1, for output
  FILE        * ofp,
  int           textw
)
{

  //
  // FOR NOW: I'm just going to use this printing style until
  //          we're getting to this point with *every* test case,
  //          and then I'll get the 'ExonCoordSet' data integrated. 
  //
  p7_tophits_SortBySeqidxAndAlipos(ExonSetTopHits);
  p7_tophits_RemoveDuplicates(ExonSetTopHits,ExonSetPipeline->use_bit_cutoffs);
  p7_tophits_SortBySortkey(ExonSetTopHits);
  p7_tophits_Threshold(ExonSetTopHits,ExonSetPipeline);

  fprintf(ofp,"\n\n");
  fprintf(ofp,"Exon Set %d (%d exons)\n",exon_set_name_id,ExonCoordSet[0]);
  for (int exon_id=0; exon_id<ExonCoordSet[0]; exon_id++) {
    fprintf(ofp,"  + Exon %d\n",exon_id+1);
    fprintf(ofp,"    Model Range:  %d..%d\n",ExonCoordSet[exon_id*4+2],ExonCoordSet[exon_id*4+4]);
    fprintf(ofp,"    Nucl. Range:  %d..%d\n",ExonCoordSet[exon_id*4+1],ExonCoordSet[exon_id*4+3]);
  }
  fprintf(ofp,"\n\n");

  p7_tophits_Targets(ofp, ExonSetTopHits, ExonSetPipeline, textw); 
  fprintf(ofp,"\n\n");
  
  p7_tophits_Domains(ofp, ExonSetTopHits, ExonSetPipeline, textw);
  fprintf(ofp,"\n\n");


  return 0;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: RunModelOnExonSets
 *
 *  Desc. :  This is the culmination of all of our work!
 *           Given a splice graph, use the spliced nucleotide coordinates for each
 *           connected component (ideally there will only be one...) to pull a final
 *           nucleotide sequence and re-run the model on that (spliced) sequence.
 *
 *  Inputs:  1.         Graph : The final SPLICE_GRAPH struct built on the set of unspliced hits.
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *           4.            go : An ESL_GETOPTS struct (in case of options... duh).
 *           5.           ofp : An open file pointer (specificially, the target for output).
 *           6.         textw : The desired line width for output.
 *
 *  Output:  Nothing is returned.  The results of running the model against the spliced
 *           nucleotide sequence are printed to ofp.
 *
 */
void RunModelOnExonSets
(
  SPLICE_GRAPH * Graph, 
  TARGET_SEQ   * TargetNuclSeq, 
  ESL_GENCODE  * gcode, 
  ESL_GETOPTS  * go,
  FILE         * ofp,
  int            textw
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'RunModelOnExonSets'",1);


  int num_exon_sets;
  int ** ExonCoordSets = GetSplicedExonCoordSets(Graph,&num_exon_sets);


  if (DEBUGGING) DumpExonSets(ExonCoordSets,num_exon_sets,TargetNuclSeq,gcode);


  // It's better to be re-using these than destroying
  // and re-allocating every time
  ESL_SQ      * NuclSeq           = esl_sq_Create();
  ESL_SQ      * AminoSeq          = esl_sq_Create();
  P7_PIPELINE * ExonSetPipeline   = NULL;
  P7_TOPHITS  * ExonSetTopHits    = NULL;
  P7_OPROFILE * ExonSetOModel     = p7_oprofile_Clone(Graph->OModel);
  P7_BG       * ExonSetBackground = p7_bg_Create(Graph->OModel->abc);


  for (int exon_set_id = 0; exon_set_id < num_exon_sets; exon_set_id++) {

    
    // If there's only one exon in this set of exons, we'll
    // skip reporting it (maybe have this be a user option?)
    //
    if (ExonCoordSets[exon_set_id][0] == 1)
      continue;


    // Grab the nucleotides for this set of exons and
    // convert them to a textized ESL_SQ
    //
    int coding_region_len;
    ESL_DSQ * ExonSetNucls = GrabExonCoordSetNucls(ExonCoordSets[exon_set_id],TargetNuclSeq,&coding_region_len);
    ESL_SQ  * NuclSeq      = esl_sq_CreateDigitalFrom(TargetNuclSeq->abc,"Exon Set",ExonSetNucls,(int64_t)coding_region_len,NULL,NULL,NULL);
    NuclSeq->idx = exon_set_id+1;
    esl_sq_Textize(NuclSeq);


    // Translate the nucleotide sequence and convert to
    // a digitized ESL_SQ
    //
    int trans_len = coding_region_len / 3;
    ESL_DSQ * ExonSetTrans = TranslateExonSetNucls(ExonSetNucls,coding_region_len,gcode);
    ESL_SQ  * AminoSeq     = esl_sq_CreateDigitalFrom(Graph->OModel->abc,"Trans",ExonSetTrans,(int64_t)trans_len,NULL,NULL,NULL);
    AminoSeq->idx = exon_set_id+1;
    strcpy(AminoSeq->orfid,"orf");



    // Prep a P7_PIPELINE datastructure and all of the other
    // friends that we need in order to produce our full evaluation
    // of the set of exons as a coding region for this protein model
    //
    P7_PIPELINE * ExonSetPipeline  = p7_pipeline_Create(go,Graph->OModel->M,coding_region_len,FALSE,p7_SEARCH_SEQS);
    ExonSetPipeline->is_translated = TRUE;
    ExonSetPipeline->strands       = p7_STRAND_TOPONLY;
    ExonSetPipeline->block_length  = coding_region_len;


    int pipeline_create_err = p7_pli_NewModel(ExonSetPipeline,ExonSetOModel,ExonSetBackground);
    if (pipeline_create_err == eslEINVAL) 
      p7_Fail(ExonSetPipeline->errbuf);

    p7_pli_NewSeq(ExonSetPipeline,AminoSeq);
    p7_bg_SetLength(ExonSetBackground,AminoSeq->n);
    p7_oprofile_ReconfigLength(ExonSetOModel,AminoSeq->n);



    // Create a P7_TOPHITS datastructure to capture the results
    // for this set of exons and run that rowdy ol' p7_Pipeline!
    //
    P7_TOPHITS * ExonSetTopHits = p7_tophits_Create();
    int pipeline_execute_err  = p7_Pipeline(ExonSetPipeline,ExonSetOModel,ExonSetBackground,AminoSeq,NuclSeq,ExonSetTopHits,NULL);
    if (pipeline_execute_err != eslOK) {
      fprintf(stderr,"\n  * PIPELINE FAILURE DURING EXON SET RE-ALIGNMENT *\n\n");
      exit(101);
    }



    // If we were successful in our search, report the hit(s)
    // that were produced!
    // This *should* always just be a single hit, but we'll
    // allow for the weird possibility that multiple hits were
    // acquired...
    //
    if (ExonSetTopHits->N)
      ReportSplicedTopHits(ExonSetTopHits,ExonSetPipeline,ExonCoordSets[exon_set_id],exon_set_id+1,ofp,textw);



    // DESTRUCTION AND REBIRTH!
    free(ExonSetNucls);
    free(ExonSetTrans);
    esl_sq_Reuse(NuclSeq);  // Takes care of 'ExonSetNucls'
    esl_sq_Reuse(AminoSeq); // Takes care of 'ExonSetTrans'
    p7_tophits_Reuse(ExonSetTopHits);
    p7_pipeline_Reuse(ExonSetPipeline);
    free(ExonCoordSets[exon_set_id]);

  }


  // That's all we needed!  Good work, team!
  if (ExonSetTopHits) p7_tophits_Destroy(ExonSetTopHits);
  if (ExonSetPipeline) p7_pipeline_Destroy(ExonSetPipeline);
  p7_oprofile_Destroy(ExonSetOModel);
  p7_bg_Destroy(ExonSetBackground);
  free(ExonCoordSets);


  if (DEBUGGING) DEBUG_OUT("'RunModelOnExonSets' Complete",-1);

}













/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SpliceHits
 *
 *  Desc. :  This is the top-level function for spliced hmmsearcht
 *
 *  Inputs:  1.        TopHits : A set of hits generated by hmmsearcht.
 *                               These serve as the initial collection of likely exons that
 *                               we try to find ways to stitch together into full-model hits.
 *           2. GenomicSeqFile : An ESL_SQFILE struct holding the genomic file that contains
 *                               all of the target sequences for search.
 *           3.             gm : The straightforward profile for the protein / family.
 *           4.             om : The optimized profile (assumed to be built on 'gm').
 *           5.          gcode : An ESL_GENCODE struct (mainly used for translation).
 *           6.             go : An ESL_GETOPTS struct (in case of options... duh).
 *           7.            ofp : An open file pointer (specificially, the target for output).
 *           8.          textw : The desired line width for output.
 *
 *  Output:  Nothing is returned, but if we find a way to splice together some of the
 *           input hits they're written out to the ofp.
 *
 */
void SpliceHits
(
  P7_TOPHITS  * TopHits,
  ESL_SQFILE  * GenomicSeqFile,
  P7_PROFILE  * gm,
  P7_OPROFILE * om,
  ESL_GENCODE * gcode,
  ESL_GETOPTS * go,
  FILE        * ofp,
  int           textw
)
{


  if (DEBUGGING) DEBUG_OUT("Starting 'SpliceHits'",1);


  // Very first thing we want to do is make sure that our hits are
  // organized by the 'Seqidx' (the target genomic sequence) and position
  // within that file.
  //
  // NOTE: I still need to ensure I'm taking advantage of this sorting!
  //
  p7_tophits_SortBySeqidxAndAlipos(TopHits);



  // Given that our hits are organized by target sequence, we can
  // be a bit more efficient in our file reading by only pulling
  // target sequences as they change (wrt the upstream hit)
  //
  TARGET_SEQ * TargetNuclSeq = GetTargetNuclSeq(GenomicSeqFile,TopHits);



  // This function encapsulates a *ton* of the work we do.
  // In short, take the collection of unspliced hits and build
  // a splice graph representing all (reasonable) ways of splicing
  // them.
  //
  SPLICE_GRAPH * Graph = BuildSpliceGraph(TopHits,TargetNuclSeq,gm,om,gcode);



  // If we're debugging, it might be useful to get a (quite rough)
  // picture of what our splice graph looks like.
  //
  if (DEBUGGING) DumpGraph(Graph);



  // If the graph doesn't currently have a complete path,
  // we'll see if we can find some holes to plug with small
  // (or otherwise missed) exons
  //
  if (Graph->has_full_path == 0)
    AddMissingExonsToGraph(Graph,TargetNuclSeq,gcode);



  // Re-run the model on the extracted nucleotide sequence(s)
  // for each connected component of the graph.
  //
  RunModelOnExonSets(Graph,TargetNuclSeq,gcode,go,ofp,textw);



  // CLEANUP
  //SPLICE_GRAPH_Destroy(Graph);
  TARGET_SEQ_Destroy(TargetNuclSeq);


  if (DEBUGGING) {
    DEBUG_OUT("'SpliceHits' Complete",-1);
    fprintf(stderr,"\n\n");
  }

}
/*
 *                                                                                              END SPLICING STUFF
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */






















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


                                    ////////////////////////////////////////////////////////////////////
                                    //                                                                //
                                    //                          * SPLICING *                          //
                                    //                                                                //
      if (tophits_accumulator->N > 1) SpliceHits(tophits_accumulator,dbfp,gm,om,gcode,go,ofp,textw);  //
                                    //                                                                //
                                    //                                                                //
                                    ////////////////////////////////////////////////////////////////////


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


