#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "parameters.h"
#include "struct.h"

extern void error(char *message);
extern void merror(char *msg1, char *msg2);
extern inline void calc_weight(struct diag* dg, struct scr_matrix* smatrix, 
			struct prob_dist *pdist);
extern inline void calc_ov_weight(struct diag* dg, struct diag_col *dcol, struct scr_matrix* smatrix, 
		    struct prob_dist *pdist);
//extern struct seq_part* create_seq_part(int num, struct seq* aSeq, unsigned int startpos);
extern long double** create_tmp_pdist(struct prob_dist *pdist);
extern void free_tmp_pdist(long double **dist, int length);

extern struct diag* create_diag(int n1, struct seq* sq1, unsigned int sp1,
                         int n2, struct seq* sq2, unsigned int sp2,
                         int dlength);
extern void free_diag(struct diag* dg);
extern inline struct simple_diag_col* find_diags_guided(struct scr_matrix *smatrix,  
							struct prob_dist *pdist, 
							struct gt_node* n1,  
							struct gt_node* n2, 
							struct alignment *algn, 
							double thres_weight,
							int *** diag_info);



extern struct alignment* create_empty_alignment(struct seq_col *scol);
extern void free_alignment(struct alignment *algn);
extern inline struct algn_pos *find_eqc(struct algn_pos **ap, int seqnum, int pos);
extern struct alignment* copy_alignment( struct alignment *o_algn, struct alignment *algn, char doDgc);
//extern char adapt_diag(struct alignment *algn, struct scr_matrix *smatrix, struct diag* dg);
extern inline char align_diag(struct alignment *algn, struct scr_matrix *smatrix, struct diag* dg);
//extern inline struct diag_cont* enter_sorted(struct diag_cont* backlog_diags, struct diag_cont *cand);
//extern inline char fit_fpos_diag(struct alignment *algn, struct diag* dg);

//extern unsigned long allocss;
//extern unsigned long freess;

/**
 *
 * assemble.c: assembles the alignment from the diagonals
 *
 */


//double tim = 0;


/**
 * enters the candidate in a sorted way to the list and returns pointer to the first element
 */
inline struct diag_cont* enter_sorted(struct diag_cont* backlog_diags, struct diag_cont *cand) {
  if(backlog_diags==NULL) {
    cand->next=NULL;
    return cand;
  }

  struct diag_cont *prev = backlog_diags;
  struct diag_cont *actual = backlog_diags;

  //printf(" before %f %i\n", cand->dg->total_weight, backlog_diags->dg);
  //  if( (cand==backlog_diags) || (backlog_diags->next==cand)) error(" enter_sorted(): critical error");

  if((cand->dg->total_weight >= backlog_diags->dg->total_weight)) {
    cand->next = backlog_diags;
    return cand;
  }

  while( (actual!=NULL) && (cand->dg->total_weight < actual->dg->total_weight)) {
    prev = actual;
    actual = actual->next;
  }
  prev->next = cand;
  cand->next = actual;

  //printf(" after\n");

  return backlog_diags;
}




/**
 * returns whether the given first position of the diag fits intho the given alignment
 */
inline char fit_fpos_diag(struct alignment *algn, struct diag* dg) {
  unsigned int s1 = dg->seq_p1.num;
  unsigned int s2 = dg->seq_p2.num;

  char o1 = algn->seq_is_orphane[s1];
  char o2 = algn->seq_is_orphane[s2];
  // if any sequence is orphane the diag is always consistent
  if( o1 || o2) return 1; 

  int sp1 = dg->seq_p1.startpos;
  int sp2 = dg->seq_p2.startpos;
  int ep1 = sp1+dg->length-1;
  int ep2 = sp2+dg->length-1;

  struct algn_pos **ap=algn->algn;
  int predF, succF;
  
  
  struct algn_pos *tap;

  tap = find_eqc(ap, s1, sp1);
  if(tap->predF!=NULL) { //&& tap->succF!=NULL) {
    predF = tap->predF[s2];
  } else {
    predF=-1; 
  }
  if(tap->succF!=NULL) {
    succF = tap->succF[s2];
  } else {
    succF = ep2+1;
  }
  if( ( (predF>=sp2)|| (succF<=sp2) ) && !(predF==sp2 && succF==sp2)) {
    //printf(" leave fit_fpos 0 \n");
    return 0;
  } else {
    //    printf(" leave fit_fpos 1 \n");
    return 1;
  }
  
}



/**
 * changes the startpos's and length of the diag such that
 * it becomes consistent with the given alignment. If the diag
 * has more than one part that is consistent the leftmost will be chosen
 * The return value is 1 if any changes were made otherwise 0. 
 * (The weight of the diag is not recalculated !).
 *
 */
char adapt_diag(struct alignment *algn, struct scr_matrix *smatrix, struct diag* dg) {
//char adapt_diag(struct alignment *algn, struct diag* dg) {
  //printf(" ENTER adapt\n");

  if(dg->multi_dg) {
    char adapted = 0;
    int i;
    for(i=0;i<dg->multi_length;i++) {
      if(dg->multi_cont[i]!=NULL) {
	if(dg->multi_cont[i]->length > 0) {
	  //printf("  adapted before %i %i\n", adapted, dg->multi_cont[i]->length);
	  adapted = adapted || adapt_diag(algn, smatrix, dg->multi_cont[i]);
	  /*
	  if(dg->multi_cont[i]->length == 0) {
	    dg->multi_cont[i]->meetsThreshold = 0;
	    }*/
	  //printf("   adapted after %i %i\n", adapted, dg->multi_cont[i]->length);
	}
	if(dg->multi_cont[i]->length==0) {
	  free(dg->multi_cont[i]);
	  dg->multi_cont[i]=NULL;
	}
      }
    }
    return adapted;
  }

  unsigned int s1 = dg->seq_p1.num;
  unsigned int s2 = dg->seq_p2.num;

  char o1 = algn->seq_is_orphane[s1];
  char o2 = algn->seq_is_orphane[s2];
  // if any sequence is orphane the diag is always consistent
  if( o1 || o2) return 0; 

  int sp1 = dg->seq_p1.startpos;
  int sp2 = dg->seq_p2.startpos;
  int ep1 = sp1+dg->length-1;
  int ep2 = sp2+dg->length-1;

  struct algn_pos **ap=algn->algn;
  int predF, succF;
  int rlen;
  
  
  // cut off the beginning of the diag
  struct algn_pos *tap;// = find_eqc(ap, s1, sp1);

  sp1--;
  sp2--;
  // jump to a consistent position
  //char included = 1;
  do {
    sp1++; sp2++;
    if(sp1<=ep1) {
      tap = find_eqc(ap, s1, sp1);
      if(tap->predF!=NULL) { //&& tap->succF!=NULL) {
	predF = tap->predF[s2];
      } else {
	predF=-1; 
      }
      if(tap->succF!=NULL) {
	succF = tap->succF[s2];
      } else {
	succF = ep2+1;
      }
    }
    //if(predF!=sp2 || succF!=sp2) included = 0;
  } while( ( (predF>=sp2)|| (succF<=sp2) ) && !(predF==sp2 && succF==sp2) && sp1<=ep1);

  // cutoff low scoring positions at the beginning
  /*
  while( (sp1<=ep1)) {
      a1 = c2n[data1[sp1]];
      a2 = c2n[data2[sp2]];
      score1 = sdata[smatrixlen*a1+a2];
      
      if(score1<PROT_SIM_SCORE_THRESHOLD) {
	//printf(" SHORT1 ALARM %i %i %i %i %i %i\n",sp1,ep1, sp2,ep2,dg->length, score1);
	sp1++;
	sp2++;
      } else {
	break;
      }
  }
  */
  
  // check whether the diagonal has been cut off to zero length
  if(sp1>ep1) {
    //    printf(" OUT OF RANGE %i %i %i \n",predF, succF, sp2);
    //if(included) 
    //  dg->weight = 0.0;
    //else
    //  dg->weight = -1.0;
		  
    dg->length=0;
    return 1;
  }

  // cut off the end of the diag
  rlen=0;
  do {
    rlen++;
    //printf("   rlen: %i %i %i %i %i\n", ep1, sp1, sp2, dg->length, sp1+rlen-1);
    if((sp1+rlen-1)>ep1) {
      break;
    }else {
      tap = find_eqc(ap, s1, sp1+rlen-1);
      //printf("   after rlen: %i\n", rlen);
      if(tap->predF!=NULL) { 
	predF = tap->predF[s2];
      } else {
	predF=-1; 
      }
      if(tap->succF!=NULL) {
	succF = tap->succF[s2];
      } else {
	succF = sp2+rlen;
      }
    }
  } while( ( ((succF>=(sp2+rlen)) && (predF<(sp2+rlen-1))) || ((succF==(sp2+rlen-1)) && (predF==(sp2+rlen-1))))
	   && ((sp1+rlen-1)<=ep1));
  rlen--;
  
  // cutoff low scoring positions at the end
  /*
  while(rlen>0) {
    a1 = c2n[data1[sp1+rlen-1]];
    a2 = c2n[data2[sp2+rlen-1]];
    //printf(" %i %i %i\n", a1, a2,smatrixlen);
    score1 = sdata[smatrixlen*a1+a2];
    break;
    
    if(score1<PROT_SIM_SCORE_THRESHOLD) {
      //printf(" SHORTENING ALARM %i %i %i %i !\n", sp1,sp2,rlen,score1);
      rlen--;
    } else {
      break;
    }
  }
  */
  
  if(rlen<=0) {
    dg->length=0;
    //printf("sp1: %i\n", sp1);
    return 1;
  }
  int oldlen = dg->length;
  //printf("sp1: %i\n", sp1);
  dg->length = rlen;
  dg->seq_p1.startpos=sp1;
  dg->seq_p2.startpos=sp2;
  
  if(oldlen==rlen) {
    return 0;
  }

  return 1;
}



/**
 * heapify 
void heapify_diag_array(struct diag **diags, int pos, int length, int up) {
  struct diag *dg = diags[pos];
  if(up) {
    if(pos<=3) return; 
    int parent = (pos-1)/2;
    struct diag *pdg = diags[parent];
    if(pdg->total_weight<dg->total_weight) {
      diags[parent]=dg;
      diags[pos] = pdg;
      //      printf("heapify: %i %i %i %i %i\n", pos,parent, length, dg, pdg);
            heapify_diag_array(diags, parent,length,up);
    }

  } else {
    int lchild = 2*pos+1;
    if( (lchild)>=length) return;
    int rchild = lchild+1;
    //struct diag *dg = diags[pos];
    struct diag *ldg = diags[lchild];
    struct diag *rdg = (rchild>=length ? ldg : diags[rchild]);
    int greatest = pos;
    if(ldg->total_weight > diags[greatest]->total_weight) greatest = lchild;
    if( (rchild<length) && (rdg->total_weight > diags[greatest]->total_weight)) greatest = rchild;
    if(greatest != pos) {
      diags[pos] = diags[greatest];
      diags[greatest] = dg;
      
      heapify_diag_array(diags, greatest,length,up);
      
    }
  }
}
 */
 

/**
 * heapify 
 */
void heapify_diag_array(struct diag **diags, int pos, int length, int up) {
  struct diag *dg = diags[pos];
  if(up) {
    if(pos<=3) return; 
    int parent = (pos-1)/2;
    struct diag *pdg = diags[parent];
    if( (pdg->total_weight*pdg->weight_fac)< (dg->total_weight*dg->weight_fac)) {
      //if( (pow(pdg->total_weight,2)*pdg->weight_fac)< (pow(dg->total_weight,2)*dg->weight_fac)) {
      diags[parent]=dg;
      diags[pos] = pdg;
      //      printf("heapify: %i %i %i %i %i\n", pos,parent, length, dg, pdg);
            heapify_diag_array(diags, parent,length,up);
    }

  } else {
    int lchild = 2*pos+1;
    if( (lchild)>=length) return;
    int rchild = lchild+1;
    //struct diag *dg = diags[pos];
    struct diag *ldg = diags[lchild];
    struct diag *rdg = (rchild>=length ? ldg : diags[rchild]);
    int greatest = pos;
    if( (ldg->total_weight*ldg->weight_fac) > (diags[greatest]->total_weight * diags[greatest]->weight_fac)) greatest = lchild;
    //if( (pow(ldg->total_weight,2)*ldg->weight_fac) > (pow(diags[greatest]->total_weight,2) * diags[greatest]->weight_fac)) greatest = lchild;
    if( (rchild<length) && ( (rdg->total_weight*rdg->weight_fac) > (diags[greatest]->total_weight*diags[greatest]->weight_fac))) greatest = rchild;
    //if( (rchild<length) && ( (pow(rdg->total_weight,2)*rdg->weight_fac) > (pow(diags[greatest]->total_weight,2) *diags[greatest]->weight_fac))) greatest = rchild;
    if(greatest != pos) {
      diags[pos] = diags[greatest];
      diags[greatest] = dg;
      
      heapify_diag_array(diags, greatest,length,up);
      
    }
  }
}
 


/**
 *------------------------------------------------------------------------------------------
 *                                SIMPLE ALIGNER SECTION
 *------------------------------------------------------------------------------------------
 */

/**
 * test function that constructs an arbitrary consistent alignment.this function
 * is used to check the the correctness of adapt_diag() and align_diag()
 *
 * Returns whether something new could be aligned
 */
char simple_aligner(struct seq_col *scol, struct diag_col *dcol, 
			    struct scr_matrix* smatrix, 
			    struct prob_dist *pdist, 
			    struct alignment *algn, int round) {
  int dlen = dcol->diag_amount;
  int i;
  
  struct diag * dg,*tdg;
  char changed;
  char alignedSomething = 0;
  // compute counter weights

  // heapify
  struct diag **diags = dcol->diags; //calloc(dlen, sizeof(struct diag *));
  int alloc_dlen = dlen;
  for(i= (dlen+1)/2-1;i>=0;i--) {
    heapify_diag_array(diags, i, dlen,0);
  }
  //memset(algn->redo_seqs, 0, sizeof(char)*slen*slen);

  double oldweight=0.0, prevweight;
  
  i=0;
  double total_weight = 0.0;
  double alig_time = 0.0;
  double tclock;
  struct diag tmp_diag;
  int tmp_end;
  int offset;
  struct diag *hookdg;
  while(dlen>0) {
    //    printf(" dlen %i\n", dlen);
    dg = diags[0];
    //if(dg->score==217) print_diag(dg);
    changed = 0;
    //    print_diag(dg);
    //    hookdg = NULL;
    if(dg!=NULL) {
      if((dg->length>0) && (dg->meetsThreshold || dg->anchor) ) {
	/*
	if(oldweight > 0.0)
	  if(dg->weight > oldweight) {
	    //printf(" ALARM %.20f %.20f\n", oldweight, dg->weight);
	    //print_diag(&odg);
	    //print_diag(dg);
	    //intf("       %i %i      \n", odg, dg);
	  }
	*/
	//odg = *dg;
	//printf(" pre changed\n");
	prevweight = dg->weight;
	tmp_diag = *dg;
	changed= adapt_diag(algn,smatrix, dg);
	//changed= adapt_diag(algn,NULL, dg);
	//print_diag(dg);
	//printf(" after changed %i\n", dg->length);
	//if(dg!=NULL) printf(" diag %i %i %i %i %.20f %i %i\n", dg, dg->multi_dg, dg->length, dg->meetsThreshold,dg->weight, dg->multi_length,changed);
	if(changed) {
	  //printf("\nCHANGED\n");
	  //print_diag(dg);
	  //printf("   pre recalc\n");
	  calc_weight(dg, smatrix, pdist);

	  if(dg->anchor) {
	    *dg = tmp_diag;
	  }

	  if( (dg->length > 0) && !(dg->anchor)) {
	    tmp_end = tmp_diag.seq_p1.startpos+tmp_diag.length-1;
	    if( ((dg->seq_p1.startpos+dg->length-1)< tmp_end) && 
		!(dg->multi_dg)) {
	      offset = dg->seq_p1.startpos+dg->length-tmp_diag.seq_p1.startpos;
	      tmp_diag.seq_p1.startpos += offset;
	      tmp_diag.seq_p2.startpos += offset;
	      tmp_diag.length -= offset;
	      
	      adapt_diag(algn,smatrix, &tmp_diag);
	      tmp_diag.length = tmp_end - tmp_diag.seq_p1.startpos+1;
	      calc_weight(&tmp_diag, smatrix, pdist);
	      
	      if((tmp_diag.length>0)&& tmp_diag.meetsThreshold) {
		
		hookdg = malloc(sizeof(struct diag));
		*hookdg = tmp_diag;
		hookdg->marked = 1;

		dcol->diag_amount++;
		if(dcol->diag_amount>alloc_dlen) {
		  //printf("\n\n\nresize %i %i\n\n\n", dlen, alloc_dlen);
		  alloc_dlen += 8;
		  dcol->diags = (diags = realloc(diags, sizeof(struct diag*)*alloc_dlen));
		  if(diags==NULL) error("Error increasing diag heap durign aligning.");
		}
		//print_diag(hookdg);
		//printf("dlen %i damount %i %i\n", dlen, dcol->diag_amount, hookdg);
		//free(diags[dcol->diag_amount-1]);
		dlen++;
		diags[dcol->diag_amount-1] = diags[dlen-1];
		diags[dlen-1]=hookdg;
		
		if(dlen>=2) {
		  //printf("heapify: %i\n", dlen-1);
		  heapify_diag_array(diags,dlen-1,dlen,1);
		}
		
		//print_diag(hookdg);
	      }
	    }
	  } else {
	    dg->meetsThreshold=0;
	  }
	  
	//printf("   \nafter recalc %i %.20f %.20f\n",dg->length,  oldweight, dg->weight);
	  
	//printf("%.20f %.20f\n", oldweight, dg->weight);
	  
	//if(dg->weight<prevweight*0.5) dg->meetsThreshold = 0;
	  
	// TODO: reactivate !!!
	//(dg->weight<oldweight*0.5) {
	//if(dg->weight>oldweight) printf(" WEIGHT %e %e\n", dg->weight, oldweight);
	//algn->redo_seqs[dg->seq_p1.num*slen+dg->seq_p2.num] = 1;
	//algn->redo_seqs[dg->seq_p2.num*slen+dg->seq_p1.num] = 1;
	  
	// DELETE THIS:
	//dg->meetsThreshold = 0;
	  
	  //if(para->DO_OVERLAP) calc_ov_weight(dg, dcol, smatrix, pdist); 
	} else {
	  //printf("  Pre align\n");
	  //print_diag(dg);
	  if(para->DEBUG >1) tclock = clock();
	  
	  //if(dg->anchor)  printf(" ANCHOR %.20f %.20f\n", dg->weight, dg->total_weight);
	  alignedSomething =  align_diag(algn, smatrix, dg) || alignedSomething;
	  if(para->DEBUG >1) alig_time += clock()-tclock;
	  
	  if(para->DEBUG >1) total_weight += dg->total_weight;
	  //printf("  After align\n");
	  //if(para->DEBUG >2) printf("  aligned diag %i %e\n", i, dg->weight);
	  //dg->length = 0;
	  //dg->weight = 0.0;
	  dg->meetsThreshold = 0;
	}
      } else {
	//      printf("ALARM %i %i %Le\n", i, dg->length, dg->weight);
	oldweight = dg->weight;
      }
    }

    
    if((dg==NULL) || (!dg->meetsThreshold)) {
      tdg = diags[dlen-1];
      diags[dlen-1]=dg;
      diags[0] = tdg;
      dlen--;
    } 
    heapify_diag_array(diags, 0, dlen,0);
  }
  //  if(para->DEBUG >1)   printf("  Total Weight: %.20f (total %f) with pure alignment time %f \n", total_weight, algn->total_weight, alig_time/CLOCKS_PER_SEC);
  if(para->DEBUG >1)   printf("  Total Weight: %.20f with pure alignment time %f \n", total_weight, alig_time/CLOCKS_PER_SEC);
  //  return algn;
  return alignedSomething;
}


/**
 *------------------------------------------------------------------------------------------
 *                                COMPLEX ALIGNER SECTION
 *------------------------------------------------------------------------------------------
 */



/**
 * returns a value >0 if the given diags are in conflict within the given alignment 
 * returns a value <0 if there is an non-conflicting overlap
 * returns 0 in all other non-conflicting cases
 */
static char confl_diag(struct alignment *algn, char *layer, struct diag *dg1, struct diag *dg2) {
  //  if(dg1->multi_dg || dg2->multi_dg) error(" confl_diag(): cannot accept multi dgs!");
  int s1_1 = dg1->seq_p1.num;
  int s1_2 = dg1->seq_p2.num;
  int s2_1 = dg2->seq_p1.num;
  int s2_2 = dg2->seq_p2.num;
  int ts;
  
  int sp1_1 = dg1->seq_p1.startpos;
  int sp1_2 = dg1->seq_p2.startpos;
  int sp2_1 = dg2->seq_p1.startpos;
  int sp2_2 = dg2->seq_p2.startpos;
  int tsp;

  int sl1_1 = dg1->seq_p1.sq->length;
  int sl1_2 = dg1->seq_p2.sq->length;
  int sl2_1 = dg2->seq_p1.sq->length;
  int sl2_2 = dg2->seq_p2.sq->length;
  int tsl;

  struct algn_pos* ap1_1; 
  struct algn_pos* ap1_2; 
  struct algn_pos* ap2_1; 
  struct algn_pos* ap2_2; 
  int p1_1, p1_2, p2_1, p2_2;
  int off1, off2;
  int l1 = dg1->length;
  int l2 = dg2->length;

  int pF1, pF2, sF1, sF2;
  int opF1, opF2, osF1, osF2;
  int pos2[4];
  int lpos2 = 0;

  int ret = 0;
  signed char ucmp,lcmp;


  //l1=1;l2=1;
  /*
  int step1 = (l1>10) ? 3 : 1;
  int step2 = (l2>10) ? 3 : 1;;
  if(step1==0) step1 = 1;
  if(step2==0) step2 = 1;
  */
  if(layer[dg1->seq_p1.num]!=1) {
    ts = s1_2;
    s1_2 = s1_1;
    s1_1 = ts;

    tsl = sl1_2;
    sl1_2 = sl1_1;
    sl1_1 = tsl;

    tsp = sp1_2;
    sp1_2 = sp1_1;
    sp1_1 = tsp;
  }
  if(layer[dg2->seq_p1.num]!=1) {
    ts = s2_2;
    s2_2 = s2_1;
    s2_1 = ts;

    tsl = sl2_2;
    sl2_2 = sl2_1;
    sl2_1 = tsl;

    tsp = sp2_2;
    sp2_2 = sp2_1;
    sp2_1 = tsp;
  }

  // if one is included in the other we define it as a conflict
  //if( (s1_1==s2_1) && (s1_2==s2_2)) {
  //  
  //}

  opF1 = -1;
  osF1 = sl2_1;
  opF2 = -1;
  osF2 = sl2_2;
  

  for(off1=0;off1<l1;off1++) {

    p1_1 = sp1_1; // dg1->seq_p1.startpos+off1;
    p1_2 = sp1_2; // dg1->seq_p2.startpos+off1;
    ap1_1 = find_eqc(algn->algn, s1_1, p1_1);
    ap1_2 = find_eqc(algn->algn, s1_2, p1_2);
    
    // calculate the positions in dg2 that have to be considered for conflicts
    if(ap1_1->predF!=NULL) {
      pF1 = ap1_1->predF[s2_1];
    } else {
      pF1 = opF1;
    }
    if(ap1_1->succF!=NULL) {
      sF1 = ap1_1->succF[s2_1];
    } else {
      sF1 = osF1;
    }
    if(ap1_2->predF!=NULL) {
      pF2 = ap1_2->predF[s2_2];
    } else {
      pF2 = opF2;
    }
    if(ap1_2->succF!=NULL) {
      sF1 = ap1_2->succF[s2_2];
    } else {
      sF2 = osF2;
    }

    /*
    //if(ret==0)
      if( (pF1>=sp2_1) && (pF1!=opF1) &&  (pF1<sp2_1+l2)) 
	ret--;//ret = -1;
      else
	if( (sF1>=sp2_1) && (sF1!=osF1)&& (sF1<sp2_1+l2)) 
	  ret--;//ret = -1;
	else
	  if( (pF2>=sp2_2) && (pF2!=opF2) && (pF2<sp2_2+l2)) 
	    ret--;//ret = -1;
	  else
	    if( (sF2>=sp2_2) && (sF2!=osF2)&& (sF2<sp2_2+l2)) 
	      ret--;//ret = -1;
    */

    if(pF1 < sp2_1) {
      pF1 = sp2_1;
    }
    if(pF1 >= sp2_1+l2) { 
      pF1 = sp2_1+l2-1;
      off1 = l1;
    }
    if(sF1 < sp2_1) {
      sF1 = sp2_1;
      off1 = l1;
    }
    if(sF1 >= sp2_1+l2) {
      sF1 = sp2_1+l2-1;
    }
    if(pF2 < sp2_2) {
      pF2 = sp2_2;
    }
    if(pF2 >= sp2_2+l2) {
      pF2 = sp2_2+l2-1;
      off1 = l1;
    }
    if(sF2 < sp2_2) {
      sF2 = sp2_2;
      off1 = l1;
    }
    if(sF2 >= sp2_2+l2) {
      sF2 = sp2_2+l2-1;
    }

    lpos2 = 0;
    if((pF1!=opF1)) {
      pos2[lpos2++] = pF1-sp2_1;
    }
    if((pF2!=opF2)) {
      pos2[lpos2++] = pF2-sp2_2;
    }
    if((sF1!=osF1)) {
      pos2[lpos2++] = sF1-sp2_1;
    }
    if((sF2!=opF2)) {
      pos2[lpos2++] = sF2-sp2_2;
    }

    opF1 = pF1;
    opF2 = pF2;
    osF1 = sF1;
    osF2 = sF2;
    //for(off2=0;off2<l2;off2++) {
    for(off2=0;off2<lpos2;off2++) {

      //p2_1 = sp2_1+off2;
      //p2_2 = sp2_2+off2;
      p2_1 = sp2_1+pos2[off2];
      p2_2 = sp2_2+pos2[off2];
      ap2_1 = find_eqc(algn->algn, s2_1, p2_1);
      ap2_2 = find_eqc(algn->algn, s2_2, p2_2);
      
      ucmp = -1;
      lcmp = -1;

      // upper compare
      if(s1_1==s2_1) {
	if(p1_1 == p2_1) ucmp = 0;
	if(p1_1 <  p2_1) ucmp = 1;
	if(p1_1 >  p2_1) ucmp = 2;
      } else if( (ap1_1->succF!=NULL) && (ap1_1->predF!=NULL) && 
	  (ap1_1->succF[s2_1]==ap1_1->predF[s2_1]) &&
	  (ap1_1->predF[s2_1]==p2_1)) {
	ucmp = 0;
      }  else if( ( (ap1_1->succF!=NULL) && ( ap1_1->succF[s2_1]<=p2_1))) {
	ucmp = 1;
      } else if( ( (ap1_1->predF!=NULL) && ( ap1_1->predF[s2_1]>=p2_1))) {
	ucmp = 2;
      } 

      // lower compare
      if(s1_2==s2_2) {
	if(p1_2 == p2_2) lcmp = 0;
	if(p1_2 <  p2_2) lcmp = 1;
	if(p1_2 >  p2_2) lcmp = 2;
      } else if( (ap1_2->succF!=NULL) && (ap1_2->predF!=NULL) && 
	  (ap1_2->succF[s2_2]==ap1_2->predF[s2_2]) &&
	  (ap1_2->predF[s2_2]==p2_2)) {
	lcmp = 0;
      }  else if( ( (ap1_2->succF!=NULL) && ( ap1_2->succF[s2_2]<=p2_2))) {
	lcmp = 1;
      } else if( ( (ap1_2->predF!=NULL) && ( ap1_2->predF[s2_2]>=p2_2))) {
	lcmp = 2;
      } 

      if( (ucmp>=0) && (lcmp>=0)) {
	if(ucmp!=lcmp) return 1;
      } 
    }
  }
  /*
  ret = -ret;
  if((ret<=l1*0.1)&&(ret<=l2*0.1)) {
    return 0;
  } else {
    return -1; 
  }
  */
  return ret;
}


/**
 * determine the best diags that fit into the alignment
struct simple_diag_col *determine_diags((struct alignment *algn, struct seq_col *scol, 
					 struct scr_matrix *smatrix,
					 struct prob_dist *pdist, char *layer,
					 struct diag **diags, 
					 int diag_amount) {
  struct simple_diag_col scol;
  // TODO
  return scol;
} 
 */
 


/**
 * guided aligner recursion
 *
 * doHigher: >0 - only higher and equal than threshold weight
 * doHigher: 0  - only lower than threshold weight
 */
struct alignment* guided_aligner_rec(struct alignment *palgn,
				     struct seq_col *scol, struct diag_col *dcol, 
				     struct scr_matrix* smatrix, 
				     struct prob_dist *pdist, 
				     struct gt_node *gtn, 
				     double thres_weight, char doHigher,
				 int round) {
  if(palgn==NULL) {
    palgn = create_empty_alignment(scol);
  }

  // recursion as per the guide tree
  if(gtn->isLeaf) return palgn;
  palgn = guided_aligner_rec(palgn,scol,dcol,smatrix,pdist, gtn->succ1, thres_weight,doHigher,round);

  palgn = guided_aligner_rec(palgn,scol,dcol,smatrix,pdist, gtn->succ2, thres_weight,doHigher,round);



  // now align the fragments that are between succ1 and succ2 of gtn
  int scol_len = scol->length;
  struct diag *dg, *tdg, *sdg, *stdg;
  int i,j,si,sj,k,l, l1,l2;
  struct gt_node *n1, *n2;
  n1 = gtn->succ1;
  n2 = gtn->succ2;
  int diag_amount = dcol->diag_amount;//n1->seq_num_length + n2->seq_num_length;  
  double multi_weight_fac = 0.0;
  int diag_p=0;
  struct diag **all_diags = malloc(sizeof( struct diag*)*diag_amount);
  //struct diag *looser_diags[diag_amount];
  int diag_l=0;
  struct simple_diag_col *sdcol;
  //char crossing[scol_len*scol_len];
  char changed;
  char layer[scol_len];
  //  memset(crossing, 0, sizeof(char)*scol_len*scol_len);
  memset(layer, 0, sizeof(char)*scol_len);


  for(i=0;i<n1->seq_num_length;i++) {
    si = n1->seq_num[i];
    layer[si]=1;
    for(j=0;j<n2->seq_num_length;j++) {
      sj = n2->seq_num[j];

      if(i==0) layer[sj]=2;
      //printf("   %i %i %i %i\n", i,j,si,sj);
      //crossing[si*scol_len+sj] = 1;
      //crossing[sj*scol_len+si] = 1;
      if(sj>si) {
	sdcol = dcol->diag_matrix[scol_len*si+sj];
      } else {
	sdcol = dcol->diag_matrix[scol_len*sj+si];
      }
      multi_weight_fac += pow(sdcol->weight_fac,0.5);
      for(k=0;k<sdcol->length;k++) {
	dg = (sdcol->data[k]);
	//changed = adapt_diag(palgn, smatrix, dg);
      
	//if(changed) {
	//calc_weight(dg, smatrix, pdist);
	//} 
	if(dg->meetsThreshold) {
	  if(doHigher>0) {
	    if((dg->total_weight) >= thres_weight) {
	      all_diags[diag_p++] = dg;
	    }
	  } else {
	    if((dg->total_weight) < thres_weight) {
	      all_diags[diag_p++] = dg;
	    }
	  }
	}
      }
    }
  }

  multi_weight_fac = pow(multi_weight_fac/((double)n1->seq_num_length*n2->seq_num_length),2.0);
  struct diag_col tdcol;
  int survive;


  /*
  tdcol.diags = all_diags;
  tdcol.diag_amount = diag_p;
  
  simple_aligner(scol, &tdcol, smatrix, pdist, palgn, 0);
  return palgn;
  */
  /*
  struct simple_diag_col ssdcol = split_diags(palgn, scol, 
					     smatrix,
					     pdist, all_diags, 
					     diag_p);
  struct diag **sall_diags = ssdcol.data;
  int sdiag_p = ssdcol.length;
  

  printf(" diags split before:%i after:%i\n", diag_p, sdiag_p);
  */
  
  int sdiag_p = diag_p;
  struct diag **sall_diags = all_diags;

  // if both successor nodes are leaf align directly
  if(n1->isLeaf && n2->isLeaf) {
    tdcol.diags = all_diags;
    tdcol.diag_amount = diag_p;

    //printf("  before inner rec %i %i %i\n",all_diags,diag_p,diag_amount);
    simple_aligner(scol, &tdcol, smatrix, pdist, palgn, 0);
    //printf("  after inner rec %i\n",all_diags);
    free(all_diags);
    //printf(" after inner rec\n");
    return palgn;
  }

  //if( (n1->seq_num_length > 1) && (n2->seq_num_length > 1)) {
  // OMIT for the time being
  /*
  if( ((n1->seq_num_length>1) || (n2->seq_num_length>1)) && 
      ((n1->seq_num_length*n2->seq_num_length)>=4.0) &&
     ((n1->seq_num_length > sqrt(scol->length)) || (n2->seq_num_length > sqrt(scol->length)))) {

    int ***diag_info;
    diag_info = malloc(sizeof(int **)*n1->seq_num_length);
    for(i=0;i<n1->seq_num_length;i++) {
      diag_info[i] = malloc(sizeof(int *)*n2->seq_num_length);
    }
    struct seq *seq1, *seq2;
    int maxk;

    for(i=0;i<n1->seq_num_length;i++) {
      si = n1->seq_num[i];
      for(j=0;j<n2->seq_num_length;j++) {
	sj = n2->seq_num[j];
	
	seq1 = &(scol->seqs[si]);
	seq2 = &(scol->seqs[sj]);
	// prepare diag info
	diag_info[i][j] = malloc(sizeof(int)*seq1->length);


	for(k=0;k<seq1->length;k++) {
	  diag_info[i][j][k] = -1;
	}
	if(sj>si) {
	  sdcol = dcol->diag_matrix[scol_len*si+sj];
	} else {
	  sdcol = dcol->diag_matrix[scol_len*sj+si];
	}
	for(k=0;k<sdcol->length;k++) {
	  dg = (sdcol->data[k]);
	  if(dg->meetsThreshold) {
	    for(l=0;l<dg->length;l++) {
	      if( (dg->weight >= thres_weight)) {
		if(dg->seq_p1.num==si) {
		  //printf(" 1) %i %i %i\n",dg->seq_p1.startpos+l, dg->length, seq1->length);
		  diag_info[i][j][dg->seq_p1.startpos+l] = dg->seq_p2.startpos+l;
		} else {
		  //printf(" 2) %i %i %i \n",dg->seq_p2.startpos+l,dg->length, seq1->length);
		  diag_info[i][j][dg->seq_p2.startpos+l] = dg->seq_p1.startpos+l;
		}
	      }
	    }
	  }
	}
      }
    }
    
    //printf(" before guided diags\n");
    struct simple_diag_col *mscol = find_diags_guided(smatrix, pdist,n1, n2, palgn, thres_weight, diag_info);
    //printf(" after guided diags\n");

    for(i=0;i<n1->seq_num_length;i++) {
      si = n1->seq_num[i];
      for(j=0;j<n2->seq_num_length;j++) {
	sj = n2->seq_num[j];
	free(diag_info[i][j]);
	//      free_tmp_pdist(tmp_dist[sj][si], pdist->max_dlen);
      }
      free(diag_info[i]);
    }
    free(diag_info);
    
    diag_amount += mscol->length;
    all_diags = realloc(all_diags, sizeof( struct diag*)*diag_amount);
    sall_diags = all_diags;
    for(i=0;i<mscol->length;i++) {
      //printf(" multi dg %i %i %.20f\n", mscol->length,mscol->data[i], mscol->data[i]->weight);
      //survive = 1;
      dg = mscol->data[i];
      
      survive = 0;
      dg->weight = 0.0;
      for(j=0;j<dg->multi_length;j++) {
	if( (dg->multi_cont[j]!=NULL) && (dg->multi_cont[j]->weight >= thres_weight)) {
	  survive++;
	  dg->weight += dg->multi_cont[j]->weight;
	} else {
	  free_diag(dg->multi_cont[j]);
	  dg->multi_cont[j]=NULL;
	}
      }
      
      
      dg->total_weight = dg->weight;
      if(survive>1) {
	sall_diags[sdiag_p++] = mscol->data[i];
	mscol->data[i]->weight_fac = multi_weight_fac;
	mscol->data[i]->pred_diag = NULL;
      } else {
	free_diag(mscol->data[i]);
      }
    }

    // free data 
    free(mscol->data);
    free(mscol);

  }
  if(para->DEBUG>1) printf(" Found %i relevant multi diags\n", sdiag_p - diag_p);
*/

  diag_p = sdiag_p;
  struct diag **looser_diags = malloc(sizeof( struct diag*)*diag_amount);


  /*
  // sort diags
  for(i=0;i<sdiag_p;i++) {
    //printf("   damn %i %i\n", sall_diags[i], i);
    for(j=i+1;j<sdiag_p;j++) {
      dg = sall_diags[i];
      tdg = sall_diags[j];
      //if( (tdg->total_weight*dg->weight_fac) > (dg->total_weight*dg->weight_fac)) {
      if( (tdg->total_weight) > (dg->total_weight)) {
	sall_diags[i] = tdg;
	sall_diags[j]= dg;
      }
    }
  }

  // align the good diags first
  diag_l=0;
  for(i=0;i<sdiag_p*0.75;i++) {
    looser_diags[diag_l++] = sall_diags[i];
  }

  tdcol.diags = looser_diags;
  tdcol.diag_amount = diag_l;
  simple_aligner(scol, &tdcol, smatrix, pdist, palgn, 0);

  for(i=diag_l;i<sdiag_p;i++) {
    sall_diags[i-diag_l] = sall_diags[i];
  }
  sdiag_p -= diag_l;
  diag_l = 0;
  */
  //if(para->DEBUG>2) printf("   BEGIN: calc conflict of #diags: %i\n", sdiag_p);

  // build vertex cover graph
  char confl;
  int multilen1; int multipos1;
  int multilen2; int multipos2;
  char multi1, multi2;
  for(i=0;i<sdiag_p;i++) {
    dg = sall_diags[i];
    //dg->weight_fac = 1.0;
    dg->weight_sum = dg->total_weight;//*dg->weight_fac;
    multilen1 = 1;
    multipos1 = 0;
    multi1 = dg->multi_dg;
    if(multi1) multilen1 = dg->multi_length;

    //dg->weight_sum = dg->total_weight;
    //printf(" degree %i\n", dg->degree);
    sdg = dg;
    while(multipos1<multilen1) {
      if(multi1) dg = sdg->multi_cont[multipos1];

      if(dg!=NULL) {
	for(j=i+1;j<sdiag_p;j++) {
	  tdg = sall_diags[j];
	  
	  multilen2 = 1;
	  multipos2 = 0;
	  multi2 = tdg->multi_dg;
	  if(multi2) multilen2 = tdg->multi_length;
	  stdg = tdg;
	  
	  while(multipos2<multilen2) {
	    if(multi2) tdg = stdg->multi_cont[multipos2];
	    if(tdg!=NULL) {

	      confl = confl_diag(palgn, layer, dg, tdg) ;

	      if( (confl>0) ) {
		//printf(" conflict %i %i !\n",i,j);
		sdg->degree++;
		if(sdg->degree > sdg->max_degree) {
		  sdg->max_degree += 64;
		  sdg->neighbours = realloc(sdg->neighbours, sizeof(struct diag *)*sdg->max_degree);
		}
		sdg->neighbours[sdg->degree - 1] = stdg;
		
		stdg->degree++;
		if(stdg->degree > stdg->max_degree) {
		  stdg->max_degree += 64;
		  stdg->neighbours = realloc(stdg->neighbours, sizeof(struct diag *)*stdg->max_degree);
		}
		stdg->neighbours[stdg->degree - 1] = sdg;
		//printf(" CONFLICT FOUND %i %i!\n", i,j);
	      } 
	    }
	    multipos2++;
	  }
	}
      }
      multipos1++;
    }
  }
  // if(para->DEBUG>2) printf("   END: calc conflict of #diags: %i\n", sdiag_p);

  // perform clarkson vertex cover
  int max_n;
  double max_d;
  double t_d;
  while(1) {
    max_n = -1;
    max_d = -1.0;
    for(i=0;i<sdiag_p;i++) {
      dg = sall_diags[i];
      if(dg != NULL) {
	//printf(" degree %i\n", dg->degree);
	//printf(" before %i %i %i\n",i,sdiag_p,dg);
	//printf("   degree %i\n", dg->degree);
	//printf(" after\n");
	if(dg->degree > 0) {
	  t_d = dg->degree / dg->weight_sum;
	  if((max_d<t_d) || (max_n<0)){
	    max_n = i;
	    max_d = t_d;
	  }
	}
      }
    }

    //printf(" max_n: %i\n", max_n);
    if(max_n < 0) break;
    dg = sall_diags[max_n];
    sall_diags[max_n] = NULL;
    looser_diags[diag_l++] = dg;
    //    printf("  wlooser %i %i %i %f\n", max_n, dg, dg->multi_dg, dg->total_weight);
    //printf("          0sall[0] %i \n",sall_diags[0]);
    //i=11;
    //printf("          00sall[0] %i \n",sall_diags[0]);
    for(i=0; i < dg->degree; i++) {
      //printf("            1sall[0] %i i=%i\n",sall_diags[0],i);
      //printf(" %i %i\n",i,dg->degree);
      dg->neighbours[i]->weight_sum -= max_d;
      tdg = dg->neighbours[i];

      //printf(" before loop %i %i\n",i,dg->degree);
      for(j=0;j<tdg->degree;j++) {
	if(tdg->neighbours[j]==dg) {
	  //printf(" MATCH \n");
	  tdg->neighbours[j]=tdg->neighbours[tdg->degree-1];
	  tdg->degree--;
	  j=tdg->degree;//break;
	}
      }
      //printf("            2sall[0] %i\n",sall_diags[0]);
    }
    free(dg->neighbours);
    dg->neighbours=NULL;
    dg->degree = 0;
    dg->max_degree=0;
  }

  /*
  struct algn_pos *ap1 = find_eqc(palgn->algn,4,19);
  if(ap1->predF!=NULL)printf(" 4-19 predF-0 %i\n", ap1->predF[0]);
  ap1 = find_eqc(palgn->algn,0,4);
  if(ap1->succF!=NULL)printf(" 0-4  succF-4 %i\n", ap1->succF[4]);
  */

  //if( (sdiag_p>60) && (sall_diags[60]!=NULL)) sall_diags[0]=sall_diags[60];
  /*
  if((sdiag_p>60) && (sall_diags[60]!=NULL) && (sall_diags[52]!=NULL)) {
    printf(" conflict pre %i\n", confl_diag(palgn, layer, sall_diags[52], sall_diags[60]));
  }
  */
  /*
  for(i=0;i<sdiag_p;i++) {
    dg = sall_diags[i];
    if(dg!=NULL) {
      //printf(" num: %i\n",i);
      //printf(" a) %i %i %i\n", dg->seq_p1.num, dg->seq_p1.startpos, dg->length);
      //printf(" b) %i %i %i\n\n", dg->seq_p2.num, dg->seq_p2.startpos, dg->length);
      //rdg = *dg;
      //printf(" before adapt %i %i\n", i, sdiag_p);
      changed = adapt_diag(palgn, smatrix, dg);
      //printf(" after adapt %i %i %i %i %i\n", changed, i, sdiag_p, rdg.length,dg->length);

      if(changed) {
	calc_weight(dg, smatrix, pdist);
	looser_diags[diag_l++] = dg;
      } else {
	//if(dg->meetsThreshold) {
	align_diag(palgn, smatrix, dg);
	//}
	dg->meetsThreshold = 0;
      }

      if(dg->neighbours!=NULL) {
	free(dg->neighbours);
	dg->neighbours = NULL;
      }
    } else {
    }
  }
  */
 

  // align the winners first
  tdcol.diags = all_diags;
  tdcol.diag_amount = 0;
  for(i=0;i<sdiag_p;i++) {
    if(sall_diags[i]!=NULL) {
      //printf("  winner %i %i %f\n", sall_diags[i], sall_diags[i]->multi_dg, sall_diags[i]->total_weight);
      all_diags[tdcol.diag_amount++] = sall_diags[i];
    }
  }

  simple_aligner(scol, &tdcol, smatrix, pdist, palgn, 0);
  all_diags = tdcol.diags;
  for(i=0;i<tdcol.diag_amount;i++) {
    if(all_diags[i]->marked ) {
      //printf(" free %i %i %i\n", all_diags[i], all_diags[i]->multi_dg, looser_diags[0]);
      free_diag(all_diags[i]);
      //printf(" free end%i\n", all_diags[i]);
    } else {
      if(all_diags[i]->neighbours!=NULL) {
	free(all_diags[i]->neighbours);
	all_diags[i]->neighbours=NULL;
      }
      all_diags[i]->meetsThreshold = 0;
    }
  }

  // align the loosers afterwards
  tdcol.diags = looser_diags;
  tdcol.diag_amount = diag_l;
  simple_aligner(scol, &tdcol, smatrix, pdist, palgn, 0);
  looser_diags = tdcol.diags;
  for(i=0;i<tdcol.diag_amount;i++) {
    if(looser_diags[i]->marked ) {
      //printf(" lfree %i %i %i %i %i %i\n", i, diag_l, tdcol.diag_amount, looser_diags[i], all_diags[i]->multi_dg, looser_diags[0]);
      free_diag(looser_diags[i]);
      //printf(" lfree end%i\n", all_diags[i]);
    } else {
      if(looser_diags[i]->neighbours!=NULL) {
	free(looser_diags[i]->neighbours);
	looser_diags[i]->neighbours=NULL;
      }
      looser_diags[i]->meetsThreshold = 0;
    }
  }


  //printf("-----------------END GUIDED ALIGNER RECURSION STEP-------------------\n");

  // TODO: free/remove used diags
  //free(sall_diags);
  free(all_diags);
  free(looser_diags);

  return palgn;
}


/**
 * guided aligner method
 * NOTE: all diags are freeded within this routine
 */
struct alignment* guided_aligner(struct alignment *palgn,
				 struct seq_col *scol, struct diag_col *dcol, 
				 struct scr_matrix* smatrix, 
				 struct prob_dist *pdist, 
				 struct gt_node *gtn,
				 int round) {
  struct diag **all_diags = dcol->diags;
  int diag_p = dcol->diag_amount;
  struct diag *dg, *tdg;
  int i,j,k;


  /*
  // sort diags
  for(i=0;i<diag_p;i++) {
    for(j=i+1;j<diag_p;j++) {
      dg = all_diags[i];
      tdg = all_diags[j];
      //if( (tdg->total_weight*dg->weight_fac) > (dg->total_weight*dg->weight_fac)) {
      if( (tdg->total_weight) > (dg->total_weight)) {
	all_diags[i] = tdg;
	all_diags[j]= dg;
      }
    }
  }
  */


  int slen = scol->length;

  double tavg, ttavg;
  double avg_weight = 0.0;
  struct simple_diag_col *sdcol;
  int tsdlen=0;
  ttavg = 0;
  
  for(i=0;i<slen;i++) {
    for(j=i+1;j<slen;j++) {

      tavg=0.0;
      sdcol = dcol->diag_matrix[slen*i+j];
      for(k=0; k<sdcol->length; k++) {
	tavg +=sdcol->data[k]->total_weight;
	ttavg += sdcol->data[k]->total_weight;
	tsdlen++;
      }

      if(sdcol->length >0) {
	tavg = (tavg/sdcol->length);
      }
      /*
      tavg = 0.0;
      tsdlen = 0;
      for(k=0; k<sdcol->length; k++) {
	if(sdcol->data[k]->total_weight >= ttavg) {
	  tavg +=sdcol->data[k]->total_weight;
	  tsdlen++;
	}
      }

      tavg = tavg / tsdlen;
      */
      avg_weight += tavg;
    }
  }

  avg_weight = avg_weight/ (slen*(slen-1)/2.0);
  double thres_weight = ttavg/tsdlen;//1.0*avg_weight; //all_diags[ 0]->total_weight*0.05;

  palgn = guided_aligner_rec(palgn,scol,dcol,smatrix,pdist,gtn,thres_weight,1,round);
  
  //palgn = guided_aligner_rec(palgn,scol,dcol,smatrix,pdist,gtn,thres_weight,0,round);
  //printf("   intermediate alignment weight %f %f\n",palgn->total_weight, thres_weight);
  
  
  struct diag_col tdcol;
  struct diag **sall_diags = malloc(sizeof(struct diag*)*diag_p);
  int sdiag_p=0;
  
  for(i=0;i<diag_p;i++) {
    if(all_diags[i]->meetsThreshold && (all_diags[i]->total_weight < thres_weight)) {
      sall_diags[sdiag_p++] = all_diags[i];
    } else {
      free(all_diags[i]);
    }
  }

  tdcol.diags = sall_diags;
  tdcol.diag_amount = sdiag_p;
  simple_aligner(scol, &tdcol, smatrix, pdist, palgn, 0);
  sall_diags = tdcol.diags;

  for(i=0;i<tdcol.diag_amount;i++) {
    free(sall_diags[i]);
  }
  free(sall_diags);
  //diag_p -= thres_pos;
  //  dcol->diag_amount = diag_p;
  return palgn;
}
