 /* hmmsearcht: search profile HMM(s) against a sequence database.
 *
 */
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







/* NORD'S SUPER COOL CLUBHOUSE STARTS HERE!
 */


int GetCodonFirstNucl (int aapos) { return ((3*aapos)-2); }
int GetCodonSecondNucl(int aapos) { return ((3*aapos)-1); }
int GetCodonThirdNucl (int aapos) { return  (3*aapos);    }


void FillCodonFromAAPos
(ESL_DSQ * ntsq, int aa_pos, ESL_DSQ * codon)
{
  codon[0] = ntsq[3 * aa_pos - 2];
  codon[1] = ntsq[3 * aa_pos - 1];
  codon[2] = ntsq[3 * aa_pos];
}


/*  FUNCTION:  DomainsAreSpliceCompatible
 *
 */
int DomainsAreSpliceCompatible
(P7_ALIDISPLAY * ad1, P7_ALIDISPLAY * ad2)
{
  int max_amino_overlap = 10;

  // NOTE: We need at least one overlapping amino acid so that each sequence
  //       can "sacrifice" an amino's encoding nucleotides to provide us with
  //       splice sites.

  // We'll begin with the most basic checks:
  if (ad1->hmmto >= ad2->hmmto
      || ad1->hmmfrom >= ad2->hmmfrom
      || ad1->hmmto-1 < ad2->hmmfrom) // Do we have an overlap?
    return 0;
  
  // If we have an unacceptable amount of amino overlap we'll
  // skip out on this lucrative splicing opportunity
  if (ad1->hmmto - max_amino_overlap >= ad2->hmmfrom)
    return 0;

  // Great! Next, let's check strand direction
  int revcomp = 0;
  if (ad1->sqfrom > ad1->sqto) revcomp = 1;

  // If they have different strand directions, then they definitely
  // aren't gonna splice nicely.
  if ((revcomp && ad2->sqfrom < ad2->sqto) || (!revcomp && ad2->sqfrom > ad2->sqto))
    return 0;

  // Oh, delightful!  Now we just need to make sure our nucleotide positions
  // are all in order
  if (revcomp) {
    if (ad1->sqto < ad2->sqfrom) return 0;
  } else {
    if (ad1->sqto > ad2->sqfrom) return 0;
  }

  printf("\n  + Such handsome gentlecoordinates: %d->%d // %d->%d\n\n",
	 ad1->hmmfrom,ad1->hmmto,ad2->hmmfrom,ad2->hmmto);
  
  // You two domains are very compatible ;)
  return 1;

}



/*  FUNCTION:  SpliceDomains
 *
 */
int SpliceDomains
(P7_DOMAIN * dom1, P7_DOMAIN * dom2, P7_OPROFILE *om, ESL_GENCODE * gcode)
{

  int i,j,k;

  // NOTE that these go from 1..M, (not from 0..M-1)
  // To get character c at position i, search for fwd_emits[om->abc->Kp * i + c]
  float * fwd_emits = malloc((om->abc->Kp * (om->M + 1))*sizeof(float));
  p7_oprofile_GetFwdEmissionScoreArray(om,fwd_emits);

  // Allocating for our transition probabilities
  float * fwd_trans = malloc(7*(om->M + 2)*sizeof(float));

  // NOTE: Here we need '+2' because there will be a buffering '0' on both the front
  //       and end of the array.
  //
  // Types for MY access     // Type enumeration ('enum p7o_tsc_e' in impl_sse.h)
  int MM = 0;                // 1
  int MI = om->M + 2;        // 5
  int MD = 2 * (om->M + 2);  // 4
  int II = 3 * (om->M + 2);  // 6
  int IM = 4 * (om->M + 2);  // 2
  int DD = 5 * (om->M + 2);  // 7
  int DM = 6 * (om->M + 2);  // 3

  // Pull in the transition probabilities
  p7_oprofile_GetFwdTransitionArray(om,1,&fwd_trans[MM]);
  p7_oprofile_GetFwdTransitionArray(om,5,&fwd_trans[MI]);
  p7_oprofile_GetFwdTransitionArray(om,4,&fwd_trans[MD]);
  p7_oprofile_GetFwdTransitionArray(om,6,&fwd_trans[II]);
  p7_oprofile_GetFwdTransitionArray(om,2,&fwd_trans[IM]);
  p7_oprofile_GetFwdTransitionArray(om,7,&fwd_trans[DD]);
  p7_oprofile_GetFwdTransitionArray(om,3,&fwd_trans[DM]);

  // We want these to be logs (base e)
  for (i=0; i<7*(om->M+2); i++)
    fwd_trans[i] = log(fwd_trans[i]);

  
  // Cool!  Now we can access our transition probabilities as follows:
  //             fwd_trans[<type> + <position>]
  // Where <type> is one of the seven variables above and <position> is
  // the 1..M position in the matrix.
  
  /* In case it's ever valuable to go back to alphas from numbers
   */
  char * aa_chars = malloc(21*sizeof(char));
  aa_chars[0]  = 'A';
  aa_chars[1]  = 'C';
  aa_chars[2]  = 'D';
  aa_chars[3]  = 'E';
  aa_chars[4]  = 'F';
  aa_chars[5]  = 'G';
  aa_chars[6]  = 'H';
  aa_chars[7]  = 'I';
  aa_chars[8]  = 'K';
  aa_chars[9]  = 'L';
  aa_chars[10] = 'M';
  aa_chars[11] = 'N';
  aa_chars[12] = 'P';
  aa_chars[13] = 'Q';
  aa_chars[14] = 'R';
  aa_chars[15] = 'S';
  aa_chars[16] = 'T';
  aa_chars[17] = 'V';
  aa_chars[18] = 'W';
  aa_chars[19] = 'Y';
  aa_chars[20] = '-';
  char * dna_chars = malloc(5*sizeof(char));
  dna_chars[0] = 'A';
  dna_chars[1] = 'C';
  dna_chars[2] = 'G';
  dna_chars[3] = 'T';
  dna_chars[4] = '-';
  /*
  */

  // Pull the alignment displays from the domains
  P7_ALIDISPLAY * ad1 = dom1->ad;
  P7_ALIDISPLAY * ad2 = dom2->ad;

  // Turn our translated (aligned) characters into [0..19] values, so we can look them
  // up in score matrices (and stuff).
  //
  // NOTE: Anything that isn't a canonical amino residue character will be a '20'
  //
  ESL_DSQ * dsqt1;
  ESL_DSQ * dsqt2;
  esl_abc_CreateDsq(om->abc,ad1->aseq,&dsqt1);
  esl_abc_CreateDsq(om->abc,ad2->aseq,&dsqt2);

  // How long are each of these alignments (in aminos, including indels)?
  int ali1_len = strlen(ad1->aseq);
  int ali2_len = strlen(ad2->aseq);

  // We'll also need the model consensus sequence in order to check whether
  // or not an emission is an insertion (if the model consensus character is '.').
  ESL_DSQ * dsqm1;
  ESL_DSQ * dsqm2;
  esl_abc_CreateDsq(om->abc,ad1->model,&dsqm1);
  esl_abc_CreateDsq(om->abc,ad2->model,&dsqm2);
  
  // As a final piece of the sequence puzzle, we'll need to pull in the translated
  // nucleotide sequence (so that we can look for splice signals).

  // NOTE: We'll need to figure out whether this is RNA or DNA. To do this
  //       we just do a quick scan of the sequence looking for a 'U' or 'T'
  //
  int ntype = 0;
  for (i=0; i<strlen(ad1->ntseq); i++) {
    if      (ad1->ntseq[i] == 'T') { ntype = eslDNA; break; }
    else if (ad1->ntseq[i] == 'U') { ntype = eslRNA; break; }    
  }

  // Still don't know? Check the second Alignment
  if (!ntype) {
    for (i=0; i<strlen(ad2->ntseq); i++) {
      if      (ad2->ntseq[i] == 'T') { ntype = eslDNA; break; }
      else if (ad2->ntseq[i] == 'U') { ntype = eslRNA; break; }    
    }
    // Okay, if we still don't know then let's call it DNA
    if (!ntype) ntype = eslDNA;
  }

  // Is it called an alphabet because you're starting to list off the letters?
  // Like, "Alpha, beta, etc."? That's pretty rad, if it is the truth.
  ESL_ALPHABET * ntalpha = esl_alphabet_Create(ntype);
  
  // Swag! Now we can finally make those digital boys!
  ESL_DSQ * ntdsq1;
  ESL_DSQ * ntdsq2;
  esl_abc_CreateDsq(ntalpha,ad1->ntseq,&ntdsq1);
  esl_abc_CreateDsq(ntalpha,ad2->ntseq,&ntdsq2);

  /*
   *  CHECKPOINT!
   *
   *  We now have all of the sequences and scoring information on-hand!
   *
   *  Recall: Because each alignment will need to "sacrifice" at least one
   *  amino acid's worth of nucleotides to give us room to search for a splice
   *  signal we've required that there is at least one amino of overlap.
   *
   *  Note that we'll need to start at least 2 nucleotides into the overlap
   *  and stop 2 nucleotides from the end of the overlap:
   *
   *             E  G  G  T  |x|  L  S  D  P   : Translated
   *            gaaggaggtacc |x| ctttcagaccct  : Nucleotides
   *             E  G  D  P  |x|  E  G  D  P   : Model (consensus)
   *                         |x|             
   *            -|ag         |x| ct|-S--D--P-
   *            -E|gg        |x|  tt|S--D--P-
   *             ....        |x|        .... 
   *            -E--G-|gt    |x|      ag|--P-   <==[this is the one we like!]
   *             ....        |x|        .... 
   *            -E--G--G-|cc |x|         cc|-
   *
   */

  // Ha!  You THOUGHT I would forget the background model!  Well, joke's on you,
  // because I did, but then I REMEMBERED IT!  (The main reason why we want to
  // have a background model built is just to get the transition probabilities).
  P7_BG * bg = p7_bg_Create(om->abc);
  float nulltransprob = bg->p1;

  // Next up, figure out the region of overlap so that we can unravel
  // the scores and look for splicables (spliceables? probably spliceables).
  //
  // DANGER: We use these values in ways that suggest lengths of translated
  //         nucleotides to traverse, but these numbers ignore indels!
  //
  int hmm_overlap_start = ad2->hmmfrom;
  int hmm_overlap_end   = ad1->hmmto;
  int hmm_overlap_len   = hmm_overlap_end - hmm_overlap_start;

  // Compute the edges of 'safe' parts of the alignments (the parts that aren't
  // at risk of being eaten into to find splice sites -- BUT THEN EAT ONE POS. IN!
  // (that is to say, these are INCLUSIVE edges)
  //
  // NOTE: We can compute the position of the codon center for a given amino using
  //       the function 'nucl = (amino * 3)-1' (-2 for the codon's first nucleotide)
  //

  

  //int ali1_left_edge_aa  = ali1_len - hmm_overlap_len; // << assuming all matches
  int ali1_right_edge_aa = ali1_len;
  int ali1_left_edge_aa  = ali1_len;
  i = 0;
  while (i<hmm_overlap_len) {
    if (dsqm1[ali1_left_edge_aa] != 20) {
      i++;
    }
    ali1_left_edge_aa--;
  }
  int ali1_left_edge_dna  = 3 * ali1_left_edge_aa - 2;
  int ali1_right_edge_dna = 3 * ali1_right_edge_aa;

  
  int ali2_left_edge_aa  = 1;
  //int ali2_right_edge_aa = 1 + hmm_overlap_len; // << assuming all matches
  int ali2_right_edge_aa = 1;
  i=0;
  while (i<hmm_overlap_len) {
    if (dsqm2[ali2_left_edge_aa] != 20) {
      i++;
    }
    ali2_right_edge_aa++;
  }
  int ali2_left_edge_dna  = 1; // 1 * 3 - 2 (hardcore math)
  int ali2_right_edge_dna = 3 * ali2_right_edge_aa;

  // Now we know just how many overlapping nucleotides we have (creation of variable
  // mainly just for readability).
  int dna_overlap_len = ali2_right_edge_dna;

  // Time to memoize the per-position score contributions of the two alignments!
  int ali1_true_aa_len = ali1_right_edge_aa - ali1_left_edge_aa + 1;
  int ali2_true_aa_len = ali2_right_edge_aa - ali2_left_edge_aa + 1;
  float * ali1_pos_scores = malloc(ali1_true_aa_len * sizeof(float));
  float * ali2_pos_scores = malloc(ali2_true_aa_len * sizeof(float));

  // We'll also store whether a given position was a [M]atch, [I]nsert, or [B]elete
  // (just kidding -- [D])
  char * ali1_states = malloc(ali1_true_aa_len * sizeof(char));
  char * ali2_states = malloc(ali2_true_aa_len * sizeof(char));
  
  // We'll also pre-compute the scores contributed at each position of the left-hand
  // side of the overlap (and have the right-hand side on-hand for fun).
  float ali1_sumscore = dom1->envsc;
  float ali2_sumscore = dom2->envsc;

  // Note: (Luckily) we know where in the HMM we start
  int hmm_pos = hmm_overlap_start;

  // Left side first, kids
  char last_state = 0;
  for (i=ali1_left_edge_aa; i<=ali1_right_edge_aa; i++) {

    // What state are we in?
    char state;
    if      (dsqt1[i] < 20 && dsqm1[i] < 20) { state = 'M'; } // A happy little match!
    else if (dsqt1[i] == 20)                 { state = 'I'; } // Not translating DNA
    else                                     { state = 'D'; } // Model isn't emitting

    float trans_score = 0.0;
    if (last_state == 'M' && state == 'M') { trans_score = fwd_trans[MM + hmm_pos]; }
    if (last_state == 'M' && state == 'I') { trans_score = fwd_trans[MI + hmm_pos]; }
    if (last_state == 'M' && state == 'D') { trans_score = fwd_trans[MD + hmm_pos]; }
    if (last_state == 'I' && state == 'I') { trans_score = fwd_trans[II + hmm_pos]; }
    if (last_state == 'I' && state == 'M') { trans_score = fwd_trans[IM + hmm_pos]; }
    if (last_state == 'D' && state == 'D') { trans_score = fwd_trans[DD + hmm_pos]; }
    if (last_state == 'D' && state == 'M') { trans_score = fwd_trans[DM + hmm_pos]; }

    float emit_score = 0.0;
    if (state != 'D') { emit_score = fwd_emits[om->abc->Kp * hmm_pos + dsqt1[i]]; }

    ali1_pos_scores[i-ali1_left_edge_aa] = trans_score + emit_score;
    ali1_states[i-ali1_left_edge_aa]     = state;

    ali1_sumscore -= trans_score + emit_score;

    // Increment model position, unless we're cycling in an insert state
    if (state != 'I' || dsqt1[i] != 20) hmm_pos++;
    last_state = state;
    
  }

  // Now for the right side
  hmm_pos = hmm_overlap_start;
  last_state = 0;
  for (i=ali2_left_edge_aa; i<=ali2_right_edge_aa; i++) {

    // What state are we in?
    char state;
    if      (dsqt2[i] < 20 && dsqm2[i] < 20) { state = 'M'; } // A happy little match!
    else if (dsqt2[i] == 20)                 { state = 'I'; } // Not translating DNA
    else                                     { state = 'D'; } // Model isn't emitting

    float trans_score = 0.0;
    if (last_state == 'M' && state == 'M') { trans_score = fwd_trans[MM + hmm_pos]; }
    if (last_state == 'M' && state == 'I') { trans_score = fwd_trans[MI + hmm_pos]; }
    if (last_state == 'M' && state == 'D') { trans_score = fwd_trans[MD + hmm_pos]; }
    if (last_state == 'I' && state == 'I') { trans_score = fwd_trans[II + hmm_pos]; }
    if (last_state == 'I' && state == 'M') { trans_score = fwd_trans[IM + hmm_pos]; }
    if (last_state == 'D' && state == 'D') { trans_score = fwd_trans[DD + hmm_pos]; }
    if (last_state == 'D' && state == 'M') { trans_score = fwd_trans[DM + hmm_pos]; }

    float emit_score = 0.0;
    if (state != 'D') { emit_score = fwd_emits[om->abc->Kp * hmm_pos + dsqt2[i]]; }

    ali2_pos_scores[i-ali2_left_edge_aa] = trans_score + emit_score;
    ali2_states[i-ali2_left_edge_aa]     = state;

    // Increment model position, unless we're cycling in an insert state
    if (state != 'I' || dsqt2[i] != 20) hmm_pos++;
    last_state = state;
    
  }

  // Just to make some of the lookup stuff easier
  // * NOTE: We need these to be ESL_DSQs because that lets us use the quick translate
  //         function 'esl_gencode_GetTranslation(gcode,codon)' which returns the
  //         index of the encoded amino acid (as would be used in aa_chars).
  ESL_DSQ * codon;
  esl_abc_CreateDsq(ntalpha,"AAA",&codon);


  // DEBUGGING OUTPUT
  //
  /*  This will print out the region where we're looking for spliceables
   */
  printf("  Model     :   "); // MODEL AMINOS -- these are the same on both sides
  for (i=ali1_left_edge_aa; i<=ali1_right_edge_aa; i++)
    printf(" %c ",aa_chars[dsqm1[i]]);
  printf("  |x|  ");
  for (i=ali2_left_edge_aa; i<=ali2_right_edge_aa; i++)
    printf(" %c ",aa_chars[dsqm2[i]]);
  printf("   : Model\n");
  //
  printf("  AL1 Nucls :   "); // TARGET NUCLEOTIDES
  for (i=ali1_left_edge_dna; i<ali1_right_edge_dna; i+=3) {
    codon[0] = ntdsq1[i];
    codon[1] = ntdsq1[i + 1];
    codon[2] = ntdsq1[i + 2];
    printf("%c%c%c",dna_chars[codon[0]],dna_chars[codon[1]],dna_chars[codon[2]]);
  }
  printf("  |x|  ");
  for (i=ali2_left_edge_dna; i<ali2_right_edge_dna; i+=3) {
    codon[0] = ntdsq2[i];
    codon[1] = ntdsq2[i + 1];
    codon[2] = ntdsq2[i + 2];
    printf("%c%c%c",dna_chars[codon[0]],dna_chars[codon[1]],dna_chars[codon[2]]);
  }
  printf("   : AL2 Nucls\n");
  //
  printf("  AL1 Trans :   "); // TARGET TRANSLATION
  for (i=ali1_left_edge_aa; i<=ali1_right_edge_aa; i++)
    printf(" %c ",aa_chars[dsqt1[i]]);
  printf("  |x|  ");
  for (i=ali2_left_edge_aa; i<=ali2_right_edge_aa; i++)
    printf(" %c ",aa_chars[dsqt2[i]]);
  printf("   : AL2 Trans\n");
  printf("\n");
  /*
   */

  
  // The log probability of beginning an intron
  float intron_score = log(0.05);

  // The log probabilities of canonical / noncanonical splice sites:
  float canon_ss_score    = log(0.95);
  float noncanon_ss_score = log(0.05 / 15.0);


  // NOTE: Currently, this is written to concurrently advance through both
  //       nucleotide sequences, but that only works if we're in a match or delete
  //       state (where the nucleotide sequence is being translated) --
  //       I'll need to figure out exactly how to augment this so that it handle
  //       weirder cases (I'm guessing this will just involve having a way to
  //       allow stepping backwards in one of the nucleotide sequences).

  // Alright, friends -- here's where the party really begins.  We walk
  // along the DNA, translating it and seeing if we get anything better
  // than the original alignment score by introducing a splice site.
  int model_pos = hmm_overlap_start;
  int nucls_from_ali1 = 2;

  // Indices into the position scores / state data for each alignment
  int ali1_scoring_index = 0;
  int ali2_scoring_index = 0;

  // Variables to track where we want to stick the splice site
  float best_spliced_score = log(0); // -inf
  int   best_splicing_pos  = -1;
  int   best_nucls_from_a1 = -1;
  for (i=3; i<dna_overlap_len; i++) { // how far we are along the overlapping dna

    // Pull in the characters contributed from alignment 1
    j = ali1_left_edge_dna + i - nucls_from_ali1 - 1;
    k = 0;
    while (k<nucls_from_ali1) { codon[k++] = ntdsq1[j++]; }

    // What's the 5' splice signal like? (looking for GT)
    float donor_score = noncanon_ss_score;
    if (ntdsq1[j]==2 && ntdsq1[j+1]==3)
      donor_score = canon_ss_score;

    // Now the characters contributed from alignment 2
    j = i;
    while (k<3) { codon[k++] = ntdsq2[j++]; }

    // Codon codoon!
    int encoded_aa = esl_gencode_GetTranslation(gcode,codon);

    // What's the 3' splice signal like? (looking for AG)
    float acceptor_score = noncanon_ss_score;
    if (ntdsq2[i-2]==0 && ntdsq2[i-1]==2)
      acceptor_score = canon_ss_score;

    // DOUBLE-CHECK THAT THIS IS IN THE RIGHT PLACE!!!
    // I think it's time to re-attribute that codon using the original score:
    if (nucls_from_ali1 == 3) {
      ali1_scoring_index++;
      ali2_scoring_index++;
      ali1_sumscore += ali1_pos_scores[ali1_scoring_index];
      ali2_sumscore -= ali2_pos_scores[ali2_scoring_index];
    }

    // If we're getting 3 nucleotides from alignment 1, then that means
    float spliced_score = ali1_sumscore + ali2_sumscore;
    spliced_score += donor_score + acceptor_score;
    spliced_score += fwd_emits[om->abc->Kp * model_pos + encoded_aa];

    if (spliced_score > best_spliced_score) {
      best_spliced_score = spliced_score;
      best_splicing_pos  = i;
      best_nucls_from_a1 = nucls_from_ali1;
      printf("          <<:>> "); // DEBUGGING
    } else {                      // DEBUGGING
      printf("            :   "); // DEBUGGING
    }

    // Increment the number of nucleotides from alignment 1 (or cycle)
    // Should this be 3->0, or 4->1? Does it matter?
    nucls_from_ali1++;
    if (nucls_from_ali1 == 4) {
      nucls_from_ali1 = 1;
      model_pos++;
    }

    // DEBUGGING OUTPUT
    //
    /* Print out the codon we currently have, what it translates into, and the amino 
     * we want it to be.
     */
    printf("%c%c%c  ==  ",dna_chars[codon[0]],dna_chars[codon[1]],dna_chars[codon[2]]);
    printf("%c (want %c)",aa_chars[encoded_aa],aa_chars[dsqm2[1+model_pos-hmm_overlap_start]]);
    //
    printf("     5'SS [");
    if (donor_score == canon_ss_score) { printf("x"); }
    else                               { printf(" "); }
    printf("]   3'SS [");
    if (acceptor_score == canon_ss_score) { printf("x"); }
    else                                  { printf(" "); }
    printf("]\n");
    /*
    */
    
  }


  // WE SHOULD PRINT OUT WHAT OUR SPLICED ALIGNMENT LOOKS LIKE FOR VERIFICATION!!!!
    
    
  /* Alright, we've had our fun here.  Let's say we clean up and head off, yeah?
   */
  free(fwd_emits);
  free(fwd_trans);
  free(ali1_pos_scores);
  free(ali2_pos_scores);
  free(ali1_states);
  free(ali2_states);

  // Catalog of stuff that I'll need to look up the custom destructors for
  //dsq1   (ESL_DSQ *)
  //dsq2   (ESL_DSQ *)
  //ntdsq1 (ESL_DSQ *)
  //ntdsq2 (ESL_DSQ *)
  //codon  (ESL_DSQ *)
  //bg     (P7_BG *  )

  // Final bit of debugging fun (just for formatting)
  printf("\n");
  
  return eslOK;

}


/*  FUNCTION:  CheckSpliceViability
 *
 */
int CheckSpliceViability
(P7_TOPHITS *hitlist, P7_OPROFILE *om, ESL_GENCODE * gcode)
{
  
  printf("\n\n");
  printf("------------------- VIABILITY CHECK BEGIN -----------------------\n");

  // NOTE: Because we've sorted based on 'SeqidxAndAlipos' we know that
  //       hits to the same sequence (the only ones we'll consider for
  //       splicing) will be grouped together and in order.
  //
  //       HOWEVER, we're going to be a little neurotic about fully checking
  //       hits for compatible domains.
  //
  uint64_t upstream_hit = 0;
  while (upstream_hit < hitlist->N) {

    P7_HIT * hit1 = hitlist->hit[upstream_hit];
    
    int upstream_dom = 0;
    while (upstream_dom < hit1->ndom) {

      P7_DOMAIN * dom1 = &hit1->dcl[upstream_dom];
    
      uint64_t downstream_hit = upstream_hit;
      while (downstream_hit < hitlist->N
	     && !strcmp(hit1->name,hitlist->hit[downstream_hit]->name)) {

	P7_HIT * hit2 = hitlist->hit[downstream_hit];

	int downstream_dom = 0;
	if (downstream_hit == upstream_hit && upstream_dom == 0)
	  downstream_dom = 1;

	while (downstream_dom < hit2->ndom) {
	  
	  P7_DOMAIN * dom2 = &hit2->dcl[downstream_dom];

	  // If the stars align (i.e., if it's possible that these
	  // domains might contain exons we can splice together),
	  // then splice away!
	  if (DomainsAreSpliceCompatible(dom1->ad,dom2->ad)) {
	    SpliceDomains(dom1,dom2,om,gcode);
	  }

	  // Behold: a horrible cascade of incrementing values
	  downstream_dom++;
	  if (downstream_hit == upstream_hit && downstream_dom == upstream_dom)
	    downstream_dom++;
	}
	downstream_hit++;
      }
      upstream_dom++;
    }
    upstream_hit++;
  }
  
  printf("------------------- VIABILITY CHECK  END  -----------------------\n");
  printf("\n\n");

  return eslOK;
  
}



/* NORD'S SUPER COOL CLUBHOUSE ENDS HERE :(
 */







typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG            *bg;	         /* null model                              */
  ESL_SQ           *ntqsq;      /* query or target sequence; this is a DNA sequence in the case of nhmmscant */
  P7_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
  P7_OPROFILE      *om;          /* optimized query profile                 */
  P7_SCOREDATA     *scoredata;   /* hmm-specific data used by nhmmer */
} WORKER_INFO;



#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <f>",              2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <f>",        2 },
  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <f>",          2 },
  { "--pfamtblout", eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save table of hits and domains to file, in Pfam format <f>",  99 }, /* not for hmmsearcht */
  //{ "--aliscoresout", eslARG_OUTFILE, NULL,NULL,NULL,    NULL,  NULL,  NULL,              "save scores for each position in each alignment to <f>",       2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notrans",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't show the translated DNA sequence in domain alignment",   2 }, /*for hmmsearcht */
  { "--vertcodon",  eslARG_NONE,   FALSE,  NULL, NULL,   NULL,  NULL,  NULL,            "show the DNA codon vertically in domain alignment",            2 }, /*for hmmsearcht */
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
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",      12 },
#endif

  /* name           type        default  env  range toggles reqs incomp  help                                          docgroup*/
  { "-c",         eslARG_INT,       "1", NULL, NULL, NULL,  NULL, NULL,  "use alt genetic code of NCBI transl table <n>", 15 },
  { "-l",         eslARG_INT,      "20", NULL, NULL, NULL,  NULL, NULL,  "minimum ORF length",                            15 },
  { "-m",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, "-M",  "ORFs must initiate with AUG only",              15 },
  { "-M",         eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, "-m",  "ORFs must start with allowed initiation codon", 15 },
  { "--informat", eslARG_STRING,  FALSE, NULL, NULL, NULL,  NULL, NULL,  "specify that input file is in format <s>",      15 },
  { "--watson",   eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "only translate top strand",                     15 },
  { "--crick",    eslARG_NONE,    FALSE, NULL, NULL, NULL,  NULL, NULL,  "only translate bottom strand",                  15 },

  /* Restrict search to subset of database - hidden because these flags are
   *   (a) currently for internal use
   *   (b) probably going to change
   */
  { "--restrictdb_stkey", eslARG_STRING, "0",  NULL, NULL,    NULL,  NULL,  NULL,       "Search starts at the sequence with name <s> ",     99 },
  { "--restrictdb_n",eslARG_INT,        "-1",  NULL, NULL,    NULL,  NULL,  NULL,       "Search <j> target sequences (starting at --restrictdb_stkey)",   99 },
  { "--ssifile",    eslARG_STRING,       NULL, NULL, NULL,    NULL,  NULL,  NULL,       "restrictdb_x values require ssi file. Override default to <s>",  99 },

  
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

  char             *firstseq_key;     /* name of the first sequence in the restricted db range */
  int              n_targetseq;       /* number of sequences in the restricted range */
};

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp, ESL_SQ_BLOCK  *orf_block);
#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, ESL_SQ_BLOCK  *orf_block);
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
      
      if (puts("\nTranslation options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 15, 2, 80); 

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
//  if (esl_opt_IsUsed(go, "--aliscoresout")  && fprintf(ofp, "# alignment scores output:         %s\n",          esl_opt_GetString(go, "--aliscoresout")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

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
  if (esl_opt_IsUsed(go, "-c")           && fprintf(ofp, "# use alt genetic code of NCBI transl table: %d\n",             esl_opt_GetInteger(go, "-c")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-l")           && fprintf(ofp, "# minimum ORF length: %d\n",                           esl_opt_GetInteger(go, "-l"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-m")           && fprintf(ofp, "# ORFs must initiate with AUG only:    yes\n")                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-M")           && fprintf(ofp, "# ORFs must start with allowed initiation codon:    yes\n")                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--informat")   && fprintf(ofp, "# specify that input file is in format %s\n",         esl_opt_GetString(go, "--informat"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--watson")           && fprintf(ofp, "# only translate top strand:    yes\n")                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--crick")           && fprintf(ofp, "# only translate bottom strand:    yes\n")                                               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");



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
  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;

  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.dbfile);    

/* is the range restricted? */

#ifndef eslAUGMENT_SSI
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")  || esl_opt_IsUsed(go, "--ssifile")  )
    p7_Fail("Unable to use range-control options unless an SSI index file is available. See 'esl_sfetch --index'\n");
#else
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") )
    if ((cfg.firstseq_key = esl_opt_GetString(go, "--restrictdb_stkey")) == NULL)  p7_Fail("Failure capturing --restrictdb_stkey\n");

  if (esl_opt_IsUsed(go, "--restrictdb_n") )
    cfg.n_targetseq = esl_opt_GetInteger(go, "--restrictdb_n");

  if ( cfg.n_targetseq != -1 && cfg.n_targetseq < 1 )
    p7_Fail("--restrictdb_n must be >= 1\n");

#endif


  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);

  return status;
}

static int
do_sq_by_sequences(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq)
{
      if (wrk->do_watson) {
	esl_gencode_ProcessStart(gcode, wrk, sq);
	esl_gencode_ProcessPiece(gcode, wrk, sq);
	esl_gencode_ProcessEnd(wrk, sq);
      }

      if (wrk->do_crick) {
	esl_sq_ReverseComplement(sq);
	esl_gencode_ProcessStart(gcode, wrk, sq);
	esl_gencode_ProcessPiece(gcode, wrk, sq);
	esl_gencode_ProcessEnd(wrk, sq);
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
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *domtblfp = NULL;              /* output stream for tabular per-dom (--domtblout) */
  FILE            *pfamtblfp= NULL;              /* output stream for pfam tabular output (--pfamtblout)    */
//  FILE            *aliscoresfp  = NULL;            /* output stream for alignment scores (--aliscoresout)   */

  P7_HMMFILE      *hfp      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  P7_HMM          *hmm      = NULL;              /* one HMM query                                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  int              dbfmt    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */
  ESL_STOPWATCH   *w;
  P7_SCOREDATA    *scoredata = NULL;

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

  /* hmmsearcht */
  P7_TOPHITS       *tophits_accumulator = NULL; /* to hold the top hits information from all 6 frame translations */
  P7_PIPELINE      *pipelinehits_accumulator = NULL; /* to hold the pipeline hit information from all 6 frame translations */
  ESL_ALPHABET    *abcDNA = NULL;       /* DNA sequence alphabet                               */
  ESL_ALPHABET    *abcAMINO = NULL;       /* DNA sequence alphabet                               */
  ESL_SQ          *qsqDNA = NULL;		 /* DNA query sequence                                  */
  ESL_SQ          *qsqDNATxt = NULL;    /* DNA query sequence that will be in text mode for printing */
  int             n_targetseqs = 0;
  ESL_GENCODE     *gcode       = NULL;
  ESL_GENCODE_WORKSTATE *wrk    = NULL;
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


  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
    if (esl_opt_IsUsed(go, "--ssifile"))
      esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
    else
      esl_sqfile_OpenSSI(dbfp, NULL);
  }



  /* Open the query profile HMM file */
  status = p7_hmmfile_OpenE(cfg->hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(esl_opt_GetString(go, "--domtblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-dom output file %s for writing\n", esl_opt_GetString(go, "--domtblout")); }
  if (esl_opt_IsOn(go, "--pfamtblout")){ if ((pfamtblfp = fopen(esl_opt_GetString(go, "--pfamtblout"), "w")) == NULL)  esl_fatal("Failed to open pfam-style tabular output file %s for writing\n", esl_opt_GetString(go, "--pfamtblout")); }
//  if (esl_opt_IsOn(go, "--aliscoresout"))  { if ((aliscoresfp  = fopen(esl_opt_GetString(go, "--aliscoresout"),"w")) == NULL)  esl_fatal("Failed to open alignment scores output file %s for writing\n", esl_opt_GetString(go, "--aliscoresout")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);

  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  
   /*the query sequence will be DNA but will be translated to amino acids */
  /* TODO can we detect the type???? */
  abcDNA = esl_alphabet_Create(eslDNA); 
  abcAMINO = esl_alphabet_Create(eslAMINO); 
  qsqDNA = esl_sq_CreateDigital(abcDNA);
  qsqDNATxt = esl_sq_Create();

  /* <abc> is not known 'til first HMM is read. */
  hstatus = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (abc->type != eslAMINO) p7_Fail("hmmsearcht only supports amino acid HMMs; %s uses a different alphabet", cfg->hmmfile);

  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->hmmfile, cfg->dbfile);
//      esl_sqfile_SetDigital(dbfp, abc); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks
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
	  block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
	  if (block == NULL) 	      esl_fatal("Failed to allocate sequence block");

 	  status = esl_workqueue_Init(queue, block);
	  if (status != eslOK)	      esl_fatal("Failed to add block to work queue");
	}
#endif
    }

  /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any aa   */
  gcode = esl_gencode_Create(abcDNA, abcAMINO);
  esl_gencode_Set(gcode, esl_opt_GetInteger(go, "-c"));  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m or -M are set


  /* Set up the workstate structure, which contains both stateful 
   * info about our position in <sqfp> and the DNA <sq>, as well as
   * one-time config info from options
   */
  wrk = esl_gencode_WorkstateCreate(go, gcode);

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
      
      /*
      if (esl_opt_IsOn(go, "--aliscoresout") ) {
        scoredata = p7_hmm_ScoreDataCreate(om, NULL);
        p7_hmm_ScoreDataComputeRest(om, scoredata);
      }
      */

      /* Create processing pipeline and hit list accumulators */
      tophits_accumulator  = p7_tophits_Create(); 
      pipelinehits_accumulator = p7_pipeline_Create(go, 100, 100, FALSE, p7_SEARCH_SEQS);

      /* Outside loop: over each query sequence in <seqfile>. */
      n_targetseqs = 0;
      while ((cfg->n_targetseq < 0 || (cfg->n_targetseq > 0 && n_targetseqs < cfg->n_targetseq)) && (sstatus = esl_sqio_Read(dbfp, qsqDNA)) == eslOK )
      {
	n_targetseqs++;
        if (qsqDNA->n < 3) continue; /* do not process sequence of less than 1 codon */

	/* copy and convert the DNA sequence to text so we can print it in the domain alignment display */
         esl_sq_Copy(qsqDNA, qsqDNATxt);
	 
	 //printf("Creating 6 frame translations\n");
         /* create sequence block to hold translated ORFs */
         wrk->orf_block = esl_sq_CreateDigitalBlock(3, abcAMINO);
	 
         /* translate DNA sequence to 6 frame ORFs */
         do_sq_by_sequences(gcode, wrk, qsqDNA);
	 
	 
         for (i = 0; i < infocnt; ++i)
         {
           /* Create processing pipeline and hit list */
           info[i].th  = p7_tophits_Create();
           info[i].om  = p7_oprofile_Clone(om);
           info[i].ntqsq = qsqDNATxt; /* for printing the DNA target sequence in the domain hits display */
           info[i].pli = p7_pipeline_Create(go, om->M, 100, FALSE, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
           status = p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
           if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

//           info[i].pli->do_alignment_score_calc = esl_opt_IsOn(go, "--aliscoresout") ;
	   
//           if (esl_opt_IsOn(go, "--aliscoresout") )
//             info[i].scoredata = p7_hmm_ScoreDataClone(scoredata, om->abc->Kp);
	   
	   
#ifdef HMMER_THREADS
           if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
         }

#ifdef HMMER_THREADS
         if (ncpus > 0)  sstatus = thread_loop(threadObj, queue, dbfp, wrk->orf_block);
         else            sstatus = serial_loop(info, dbfp, wrk->orf_block);
#else
         sstatus = serial_loop(info, dbfp, wrk->orf_block);
#endif
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

         if(wrk->orf_block != NULL)
         {
            esl_sq_DestroyBlock(wrk->orf_block);
            wrk->orf_block = NULL;
         }
		 
         esl_sq_Reuse(qsqDNATxt);
         esl_sq_Reuse(qsqDNA);
		 
      } /* while ((cfg->n_targetseq < 0 || (cfg->n_targetseq > 0 &&... loop */



      
      // START NORD ////////////////////////////////////////////////////////////////

      // Okay, back to using tophits_accumulator to track down domains with overlaps
      // and then reverse-engineering the scores of individual alignment decisions...
      p7_tophits_SortBySeqidxAndAlipos(tophits_accumulator);
      CheckSpliceViability(tophits_accumulator,om,gcode);
      
      //  END  NORD ////////////////////////////////////////////////////////////////


      
      
      /* Print the results.  */
      p7_tophits_SortBySortkey(tophits_accumulator);
      p7_tophits_Threshold(tophits_accumulator, pipelinehits_accumulator);
      p7_tophits_Targets(ofp, tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
      p7_tophits_Domains(ofp, tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

      if (tblfp)     p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));
      if (domtblfp)  p7_tophits_TabularDomains(domtblfp, hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));
      if (pfamtblfp) p7_tophits_TabularXfam(pfamtblfp, hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator);
//      if (aliscoresfp) p7_tophits_AliScores(aliscoresfp, hmm->name, tophits_accumulator);

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


      if (scoredata) {
        for (i = 0; i < infocnt; ++i)
          p7_hmm_ScoreDataDestroy(info[i].scoredata);
        p7_hmm_ScoreDataDestroy(scoredata);
      }
	  
      p7_pipeline_Destroy(pipelinehits_accumulator);
      p7_tophits_Destroy(tophits_accumulator);
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
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmsearcht", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp, "hmmsearcht", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
  if (pfamtblfp) p7_tophits_TabularTail(pfamtblfp,"hmmsearcht", p7_SEARCH_SEQS, cfg->hmmfile, cfg->dbfile, go);
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
  
  esl_gencode_WorkstateDestroy(wrk);
  esl_gencode_Destroy(gcode);

  esl_sq_Destroy(qsqDNA);  
  esl_sq_Destroy(qsqDNATxt);
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
  if (pfamtblfp)     fclose(pfamtblfp);
//  if (aliscoresfp)   fclose(aliscoresfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}


static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp, ESL_SQ_BLOCK  *orf_block)
{
  int  sstatus = eslOK;
  ESL_SQ   *dbsq     = NULL;   /* one target sequence (digital)  */
  int      k;

  /* Main loop: */
  for (k = 0; k < orf_block->count; ++k)
  {
      dbsq = &(orf_block->list[k]);
      /* 
      use the name, accession, and description from the DNA sequence and
      not from the ORF which is generated by gencode and only for internal use
      */
      if ((sstatus = esl_sq_SetName     (dbsq, info->ntqsq->name))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence name failed");
      if ((sstatus = esl_sq_SetAccession(dbsq, info->ntqsq->acc))    != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence accession failed");
      if ((sstatus = esl_sq_SetDesc     (dbsq, info->ntqsq->desc))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence description failed");
	 	  
      p7_pli_NewSeq(info->pli, dbsq);
      p7_bg_SetLength(info->bg, dbsq->n);
      p7_oprofile_ReconfigLength(info->om, dbsq->n);

      p7_Pipeline(info->pli, info->om, info->bg, dbsq, info->ntqsq, info->th, info->scoredata);

      esl_sq_Reuse(dbsq);
      p7_pipeline_Reuse(info->pli);
  }

  return sstatus;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, ESL_SQ_BLOCK  *orf_block)
{
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;
  ESL_SQ       *orf_seq;
  int          num_processed_seq = 0;
  int          k;
  
  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK )
    {
      block = (ESL_SQ_BLOCK *) newBlock;
      block->count = 0;
		  
      /* fill up the new block with ORF sequences from the ORF sequence block */
      for (k = num_processed_seq; (k < orf_block->count) && (block->count < block->listSize); k++, num_processed_seq++)
      {
         orf_seq = &(orf_block->list[k]);
         esl_sq_Copy(orf_seq, &(block->list[block->count]));
         block->count++;
      }

      if(num_processed_seq == orf_block->count)	
      {
         sstatus = eslEOF;
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

          /* 
          use the name, accession, and description from the DNA sequence and
          not from the ORF which is generated by gencode and only for internal use
          */
          if ((status = esl_sq_SetName     (dbsq, info->ntqsq->name))   != eslOK)  esl_fatal("Set query sequence name failed");
          if ((status = esl_sq_SetAccession(dbsq, info->ntqsq->acc))    != eslOK)  esl_fatal("Set query sequence accession failed");
          if ((status = esl_sq_SetDesc     (dbsq, info->ntqsq->desc))   != eslOK)  esl_fatal("Set query sequence description failed");
	  
	  p7_pli_NewSeq(info->pli, dbsq);
	  p7_bg_SetLength(info->bg, dbsq->n);
	  p7_oprofile_ReconfigLength(info->om, dbsq->n);

	  p7_Pipeline(info->pli, info->om, info->bg, dbsq, info->ntqsq, info->th, info->scoredata);

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
 

