/**
 *
 * museq.c: Main program control
 
 * Author: A.R.Subramanian
 *         
 */


/**
 * TODO LIST:
 *  - option: multi diags as optional argument
 */


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "parameters.h"
#include "struct.h"
#include "translate.h"
#include "orf.h"
#include "io.h"

/**
 * external functions definitions
 */

// diag.c
//extern struct seq_part* create_seq_part(int num, struct seq* aSeq, 
//					unsigned int startpos);
//extern struct diag* create_diag(struct seq_part* part1, struct seq_part* part2,
//				int dlength);
//extern void calc_weight(struct diag* dg, struct scr_matrix* smatrix, 
//		 struct prob_dist *pdist);
//extern struct diag_col *create_diag_col(int seq_amount);
extern void free_diag(struct diag* dg);
extern void free_diag_col(struct diag_col* dcol);
extern struct diag_col *find_all_diags(struct scr_matrix *smatrix, 
				struct prob_dist *pdist, 
				struct seq_col *in_seq_col, struct alignment *algn, int round);


// prob.c
//extern struct seq* create_random_seq(struct scr_matrix *smatrix, int length);
//extern double* approx_prob(struct scr_matrix *smatrix, struct prob_dist *sdist,
//                    int diaglen, int seqlen);
extern struct prob_dist* calc_score_dist(struct scr_matrix *smatrix, int mxdlen);

// alig.c
extern void free_alignment(struct alignment *algn);
extern struct alignment* create_empty_alignment(struct seq_col *scol);
//extern char adapt_diag(struct alignment *algn, struct scr_matrix *smatrix, struct diag* dg);
extern int  simple_aligner(struct seq_col *scol, struct diag_col *dcol, 
			    struct scr_matrix* smatrix, 
			    struct prob_dist *pdist, 
		     struct alignment *algn, int round);



extern struct simple_diag_col* read_anchors(char *filename, struct seq_col* scol);

extern struct alignment* guided_aligner(struct alignment *palgn,
                                 struct seq_col *scol, struct diag_col *dcol,
                                 struct scr_matrix* smatrix,
                                 struct prob_dist *pdist, struct gt_node *gtn,
					int round);


/**
 * main program routine
 */
int main(int argc, char **argv) 
{
  parameters(argc, argv);
  version();

  // read similarity matrix
  char *smatrixfile = (char *)build_pathname(para->conf_dir,para->SCR_MATRIX_FILE_NAME);
  struct scr_matrix *smatrix = read_scr_matrix(smatrixfile);

  // print the score matrix
  if( para->DEBUG >5) {
    print_scr_matrix(smatrix);
	
  }

  // read the probability distribution for diagonals
  char *pdistfilename = (char *)build_pathname(para->conf_dir,para->DIAG_PROB_FILE_NAME);
 
  int i;
  struct seq_col *in_seq_col; 
  struct prob_dist *pdist;
  if(!para->COMPUTE_PROB){
    pdist = read_diag_prob_dist(smatrix, pdistfilename);
    
    
    in_seq_col = read_fasta(para->in_file);
    
    if(para->DNA_TRANSLATION){	
      if(para->FIND_ORF){
	if(para->ORF_FRAME){ 
	  translate_sequence_collection_orf_frame(in_seq_col);
	}
	else {
	  in_seq_col = set_longest_orf(in_seq_col);
	}
      }
      else{
	translate_sequence_collection_default(in_seq_col);
      }
    }
  } 
  // print read input sequences
  if(para->DEBUG>5) {
    int sc,scl = in_seq_col->length;
    for(sc=0; sc < scl;sc++) {
      print_seq(&(in_seq_col->seqs[sc]));
      printf("\n");
    }
  }

  // ----------------------------------------------------------------------
  // probability table generation area
  // ----------------------------------------------------------------------
  if(para->COMPUTE_PROB)
    {
      char *prob_table_output_file = (char *)build_pathname(para->conf_dir,"prob_table");
      struct prob_dist *sdist = calc_score_dist(smatrix, 100);
      print_pdist_matrix(sdist, prob_table_output_file);	
      exit(0);
    }
  
  
  double tim = clock(),tim2;

  // fast mode has higher threshold weights
  if(para->FAST_MODE) {
    para->PROT_SIM_SCORE_THRESHOLD += 0.25; 
  }

  // ----------------------------------------------------------------------
  // Consider Anchors
  // ----------------------------------------------------------------------
  struct simple_diag_col *anchors = NULL;
  struct diag_col adcol;
  struct alignment *algn= NULL;
  if(! para->FAST_MODE) {
    algn = create_empty_alignment(in_seq_col);
  }
  struct alignment *salgn = create_empty_alignment(in_seq_col);

  if(para->DO_ANCHOR>0) {
    anchors = read_anchors(para->ANCHOR_FILE_NAME, in_seq_col);

    adcol.diags = anchors->data;
    adcol.diag_amount = anchors->length;


    simple_aligner(in_seq_col, &adcol,smatrix,pdist,salgn,1);
    if(! para->FAST_MODE) simple_aligner(in_seq_col, &adcol,smatrix,pdist,algn,1);


    if(anchors!=NULL) {
      for(i=0;i<adcol.diag_amount;i++) {
	free_diag(adcol.diags[i]);
      }
      free(adcol.diags);
      free(anchors);
    }
  }

  // ----------------------------------------------------------------------
  // Compute pairwise diagonals
  // ----------------------------------------------------------------------
  //if(para->DNA_PARAMETERS && (para->SENS_MODE>0)) {
  //  para->DIAG_THRESHOLD_WEIGHT = -log(0.875);
  //}

  struct diag_col *all_diags = find_all_diags(smatrix, pdist, in_seq_col,salgn,1);
  double duration = (clock()-tim)/CLOCKS_PER_SEC;
  if(para->DEBUG >1) printf("Found %i diags in %f secs\n", all_diags->diag_amount, duration);
  int diag_amount = all_diags->diag_amount;
  
  // ----------------------------------------------------------------------
  // Compute alignment
  // ----------------------------------------------------------------------
  tim2=clock();

  if(! para->FAST_MODE) {
    struct diag *cp_diags[all_diags->diag_amount];
    for(i=0;i<diag_amount;i++) {
      cp_diags[i] = malloc(sizeof (struct diag));
      *(cp_diags[i]) = *(all_diags->diags[i]);
    }
    guided_aligner(algn, in_seq_col, all_diags,smatrix,pdist,all_diags->gt_root, 1);
    
    
    for(i=0;i<diag_amount;i++) {
      all_diags->diags[i] = cp_diags[i];
    }
    all_diags->diag_amount = diag_amount;
  }
  //struct alignment *algn = salgn;
  simple_aligner(in_seq_col, all_diags,smatrix,pdist,salgn,1);
  duration = (clock()-tim2)/CLOCKS_PER_SEC;

  if(! para->FAST_MODE) {
    if(para->DEBUG >1) printf("First alignment after %f secs. simple: %f guided: %f\n", duration, salgn->total_weight, algn->total_weight);
  } else {
    if(para->DEBUG >1) printf("First alignment after %f secs. simple: %f \n", duration, salgn->total_weight);
  }
  //if(para->DEBUG >1) printf("First alignment after %f secs. guided: %f\n", duration, algn->total_weight);
  
  free_diag_col(all_diags);

  para->DO_ANCHOR = 0; // anchors done

  // round 2+
  int round;
  char newFound = 0;
  int type;
  //printf(" SENS %i\n", para->SENS_MODE);

  // consider sensitivity level
  if(! para->FAST_MODE) {
    if(para->SENS_MODE==0) {
      para->DIAG_THRESHOLD_WEIGHT = 0.0;
    } else if(para->SENS_MODE==1) {
      if(para->DNA_PARAMETERS)
	para->DIAG_THRESHOLD_WEIGHT = -log(0.75);//-log(.875+0.125/2.0);
      else
	para->DIAG_THRESHOLD_WEIGHT = -log(0.75);
    }else if(para->SENS_MODE==2) {
      if(para->DNA_PARAMETERS)
	para->DIAG_THRESHOLD_WEIGHT = -log(0.5);//-log(0.875);
      else
	para->DIAG_THRESHOLD_WEIGHT = -log(0.5);
    }
  }

  int stype = (para->FAST_MODE ? 1 : 0);
  for(type=stype;type<2;type++) {
    for(round=2;round<=20;round++) {
      //for(round=2;round<=1;round++) {
      tim2=clock();
      all_diags = find_all_diags(smatrix, pdist, in_seq_col,(type ? salgn : algn), round);
      //all_diags = find_all_diags(smatrix, pdist, in_seq_col, algn);
      duration = (clock()-tim2)/CLOCKS_PER_SEC;
      if(para->DEBUG >1) printf("Found %i diags after %f secs\n", all_diags->diag_amount, duration);
      if(all_diags->diag_amount ==0) {
	free_diag_col(all_diags);
	break;
      } else {
	// round 2 and further we use the simple aligner
	newFound = simple_aligner(in_seq_col, all_diags,smatrix,pdist,(type ? salgn : algn),round);
	//newFound = simple_aligner(in_seq_col, all_diags,smatrix,pdist,algn,round);
	//      newFound = complex_aligner(in_seq_col, all_diags,smatrix,pdist,algn,round);
	free_diag_col(all_diags);
	if(!newFound) break;
      }
    }
  }
  if(para->DEBUG >1) 
    printf("Alignment ready!\n");

  if(! para->FAST_MODE) {
    if(para->DEBUG >1) printf("Final alignment simple: %f guided: %f\n", salgn->total_weight, algn->total_weight);
  } else {
    if(para->DEBUG >1) printf("Final alignment simple: %f \n", salgn->total_weight);
  }
  //if(para->DEBUG >1) printf("Final alignment guided: %f\n", algn->total_weight);

  
  if( ( para->FAST_MODE) || (salgn->total_weight > algn->total_weight)) {
    if(! para->FAST_MODE) free_alignment(algn);
    algn = salgn;
  }
  
  //algn = salgn;
  

  if(para->out_file==NULL) {
    if(para->OUTPUT){ 
      if(para->DNA_TRANSLATION) 
	simple_print_alignment_dna_retranslate(algn);
      else 
	simple_print_alignment_default(algn);
    }
    else{
      simple_print_alignment_default(algn);
    }
  }else {
    if(para->OUTPUT){ 
      if(para->DNA_TRANSLATION) {
	fasta_print_alignment_dna_retranslate(algn, para->out_file);
      }
      else { 
	fasta_print_alignment_default(algn, para->out_file);
      }
    }
    else{
      fasta_print_alignment_default(algn, para->out_file);
    }
  }
  duration = (clock()-tim)/CLOCKS_PER_SEC;
  printf("Total time:   %f secs\n", duration);
  printf("Total weight: %f \n", algn->total_weight);
  exit(EXIT_SUCCESS);
}



