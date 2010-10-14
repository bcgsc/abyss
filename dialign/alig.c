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
extern struct diag* create_diag(struct seq_part* part1, struct seq_part* part2, 
			 int dlength);
extern void free_diag(struct diag* dg);

//  long balance = 0;
//unsigned long allocss = 0;
//unsigned long freess = 0;

long sslen;


/**
 *
 * alig.c: Takes care of the alignment data structure
 *
 * 2003-10-31  A.R.Subramanian
 *             (Initial)
 */

/**
 * creates initial empty alignment data structure
 *
 * The pointer returned (and the ones included in the struct) 
 * has to be deallocted explicitely from memory.
 */
struct alignment* create_empty_alignment(struct seq_col *scol) {
  /*
  printf("before\n");
  sleep(5);
  */
  struct alignment* algn = malloc(sizeof(struct alignment));
  sslen = scol->length;
  //allocss += sizeof(struct alignment);
  if(algn==NULL) error("create_empty_alignment(): (1) Out of memory !");

  
  //long xsize = sizeof(struct alignment);

  algn->next = NULL;
  //  algn->prev = NULL;
  algn->total_weight = 0.0;
  //algn->pos = 0;
  algn->max_pos = -1;
  algn->scol = scol;
  //algn->aligned_diags_amount=0;
  //algn->aligned_diags = malloc(sizeof(struct diag*)*diag_amount);
  //algn->max_aligned_diags_amount = diag_amount;
  //algn->orig_max_aligned_diags_amount = diag_amount;
  //algn->backlog_diags = NULL;
  ////allocss += (sizeof(struct diag*)*diag_amount);

  //if(algn->aligned_diags==NULL)
  //  error("create_empty_alignment(): (1.5) Out of memory !");


  unsigned int slen = scol->length;
  algn->seq_is_orphane = calloc(slen, sizeof(char));
  //allocss += (sizeof(char)*slen);

  //xsize += slen*sizeof(char);

  if(algn->seq_is_orphane==NULL) error("create_empty_alignment(): (2) Out of memory !");
  //  memset(algn->seq_is_orphane, 1, slen*sizeof(char));
  
  algn->algn = malloc(sizeof(struct algn_pos *)*slen);
  //allocss += sizeof(struct algn_pos *)*slen;

  //xsize += slen*sizeof(struct algn_pos *);

  if(algn->algn==NULL) error("create_empty_alignment(): (3) Out of memory !");
  //algn->redo_seqs = calloc(slen*slen, sizeof(char));
  int i,j;
  struct seq* sq;
  for(i=0;i<slen;i++) {
    sq = &(scol->seqs[i]);
    algn->seq_is_orphane[i]=1;
    algn->algn[i] = malloc(sizeof(struct algn_pos)*sq->length );
    //allocss += sizeof(struct algn_pos )*sq->length;
    //xsize += sq->length*sizeof(struct algn_pos *);

    if(algn->algn[i]==NULL) error("create_empty_alignment(): (4) Out of memory !");

    for(j=0;j<sq->length;j++) {
      algn->algn[i][j].state = para->STATE_ORPHANE;
      //      algn->algn[i][j].isInherited = 0;
      algn->algn[i][j].predFPos = -1;
      algn->algn[i][j].succFPos = -1;

      algn->algn[i][j].eqcParent= &(algn->algn[i][j]);
      //if(j==442) printf(" parent: %i\n", algn->algn[i][j].eqcParent);
      algn->algn[i][j].eqcRank= 0;
      algn->algn[i][j].eqcAlgnPos=calloc(1, sizeof(int));;
      //allocss += sizeof(int);
      *algn->algn[i][j].eqcAlgnPos=j;
      algn->algn[i][j].proceed=calloc(1, sizeof(char));;
      //allocss += sizeof(char);
      *algn->algn[i][j].proceed = 0;

      algn->algn[i][j].predF = NULL;
      algn->algn[i][j].succF = NULL;

      algn->algn[i][j].row = i;
      algn->algn[i][j].col = j;
      algn->algn[i][j].dg_cont = NULL;
      //xsize += sizeof(char) + sizeof(int);
    }
  }
  /*
  printf("after\n");
  sleep(5);
  printf("gone\n");
  */
  //  printf(" algnsize=%i\n",xsize);
  return algn;
}


/**
 * free alignment
 *
 */
void free_alignment(struct alignment *algn) {
  struct seq_col *scol = algn->scol;
  int slen = scol->length;
  int i,j;
  struct algn_pos *apos;
  struct algn_pos *o_apos;
  struct algn_pos *tpos;
  struct seq *sq;
  struct diag_cont *dgc,*ndgc;

  if(algn->seq_is_orphane!=NULL) {
    free(algn->seq_is_orphane);
    //freess += sizeof(char)*slen;
  }
  for(i=0;i<slen;i++) {
    sq = &(scol->seqs[i]);
    for(j=0;j<sq->length;j++) {
      apos = &(algn->algn[i][j]);
      if(! (apos->state & para->STATE_INHERITED) && ! (apos->state & para->STATE_ORPHANE)) {
	//if(! (apos->state & para->STATE_INHERITED) && ! (apos->state & para->STATE_ORPHANE)) {
	if(apos->predF!=NULL) {
	  free(apos->predF);
	  //freess += sizeof(int)*slen;
	}
	if(apos->succF!=NULL) {
	  free(apos->succF);
	  //freess += sizeof(int)*slen;
	}
	
      }
      if(! (apos->state & para->STATE_INHERITED)  || (apos->state & para->STATE_ORPHANE) ) {
	free(apos->eqcAlgnPos);
	free(apos->proceed);
	//freess += sizeof(int)+sizeof(char);
      }
      dgc = apos->dg_cont;
      while(dgc!=NULL) {
	ndgc = dgc->next;
	//printf(" free %i %i %i %i\n", i,j,dgc, dgc->next);
	free(dgc);
	//freess += sizeof(struct diag_cont);
	dgc = ndgc;
      }
      
    }
    free(algn->algn[i]);
    //freess += sizeof(struct algn_pos)*sq->length;
  }
  //dgc = algn->backlog_diags;
  //while(dgc!=NULL) {
  //  ndgc = dgc->next;
  //  free(dgc);
  //  //freess += sizeof(struct diag_cont);
  //  dgc = ndgc;
  //}



  //free(algn->aligned_diags);
  ////freess += sizeof(struct diag *)*algn->;
  free(algn->algn);
  free(algn);
  ////freess += sizeof(struct algn_pos *)*slen + sizeof(struct alignment);
}


/**
 * returnes the representative of the equivalence class
 */
struct algn_pos *_find_eqc(struct algn_pos *ap) {
  if(ap!=ap->eqcParent) {
    //    if(doprint) printf("    FIND: %i %i\n", ap, ap->eqcParent);
    
    /**
    if(ap->eqcParent->eqcAlgnPos!=NULL) {
      if( (ap->eqcAlgnPos!=NULL) && *ap->eqcParent->eqcAlgnPos < *ap->eqcAlgnPos) {
	//	if(doprint)  printf("    1.1 ALGNPOS: %i %i\n", ap, ap->eqcParent);
        *ap->eqcParent->eqcAlgnPos = *ap->eqcAlgnPos;
      } 
    } else {
      //	ap->eqcParent->eqcAlgnPos = ap->eqcAlgnPos;
    }
    */
    ap->eqcParent = _find_eqc(ap->eqcParent);
  }
  //  if(doprint)   printf("    1.ALGNPOS: %i %i\n", ap, ap->eqcParent);
  //  if(doprint)   printf("    2.ALGNPOS: %i %i\n", ap, ap->eqcParent);
  return ap->eqcParent;
}

/**
 * returnes the representative of the equivalence class
 */
struct algn_pos *find_eqc(struct algn_pos **ap, int seqnum, int pos) {
  //if(1) printf("%i %i %i\n", ap, seqnum,pos);
  struct algn_pos *tap = &ap[seqnum][pos];
  struct algn_pos *eq_ap;
  eq_ap = _find_eqc(tap);
  struct diag_cont *old_dgc;
  
  //  if(eq_ap->eqcAlgnPos != tap->eqcAlgnPos) {
  if(eq_ap != tap) {
    
    /*
    if(tap->state & para->STATE_ORPHANE) {
      printf(" ALARM ORPHANE %i %i\n", eq_ap->state, tap->state);
    }
    if(eq_ap->state & para->STATE_ORPHANE) {
      printf("    ALARM ORPHANE %i %i\n", eq_ap->state, tap->state);
    }
    */
    //if((tap->eqcAlgnPos!=NULL) && (*tap->eqcAlgnPos > *eq_ap->eqcAlgnPos))
    //  *eq_ap->eqcAlgnPos = *tap->eqcAlgnPos;

    //    if(eq_ap->eqcAlgnPos != tap->eqcAlgnPos) 
    if( (!(tap->state & para->STATE_INHERITED))) { //&& !oldparentIsInherited) {
      if(tap->eqcAlgnPos !=NULL) {
	//if(pos==175) printf("free eqcAlgnPos: %i\n", tap->eqcAlgnPos);
	//balance -= sizeof(int);
	free(tap->eqcAlgnPos);
	tap->eqcAlgnPos = NULL;
	//freess += sizeof(int);
      }

      if(tap->proceed !=NULL) {
	//printf("free proceed: %i\n", tap->proceed);
	//balance -= sizeof(char);
	free(tap->proceed);
	//freess += sizeof(char);
	//printf("after free proceed: %i\n", tap->proceed);
	tap->proceed = NULL;
      }
      
      //if(tap->predFPos>=0) 
      if( (eq_ap->predF != tap->predF) && (1 || (tap->predFPos<0)) && (tap->predF!=NULL)){
	//printf("          free predF: %i %i %i %i\n", tap->predF, tap->isInherited, seqnum, pos);
	//balance -= sizeof(int)*sslen;
	  //printf (" 3. free: %i %i %i\n",allocs, frees, allocs-frees);
	//freess += sizeof(int)*sslen;
	free(tap->predF);
	//printf("          after free predF: %i\n", tap->predF);
	tap->predF=NULL;
      }
      //if(tap->succFPos>=0) 
      if((eq_ap->succF != tap->succF) && (1 || (tap->succFPos<0)) && (tap->succF!=NULL)) {
	//printf("free succ: %i\n", tap->succF);
	//balance -= sizeof(int)*sslen;
	  //printf (" 4. free: %i %i %i\n",allocs, frees, allocs-frees);
	//freess += sizeof(int)*sslen;
	free(tap->succF);
	//printf("after free succ: %i\n", tap->succF);
	tap->succF=NULL;
      }
    } /*else {
      if(eq_ap->state & para->STATE_INHERITED) {
	printf(" inherited alarm!\n");
      }
      if(eq_ap->state & para->STATE_ORPHANE) {
	printf(" orphane alarm!\n");
      }
      }*/
    old_dgc = tap->dg_cont;
    *tap = *eq_ap;
    tap->dg_cont = old_dgc;
    tap->row = seqnum;
    tap->col = pos;
    tap->state = (tap->state | para->STATE_INHERITED);
  }
  //if(seqnum==0 &&pos==175)  printf("after !=\n");
  // tap = eq_ap;

  struct algn_pos *ttap;
  if(tap->predFPos>=0 ) {
    //if(seqnum==0 && pos==175) printf("          alarm predF: %i %i %i %i\n", tap->predF, tap->predFPos, seqnum, pos);
    //    printf ("PRE Pos %i %i \n", tap->predFPos, pos);
    if( (tap->predFPos==pos)|| !(tap->state & para->STATE_ORPHANE)) {
      printf("pred ALARM %i %i\n", tap->predFPos, pos);
      exit(99);
    }
    ttap=find_eqc(ap, seqnum, tap->predFPos);
    tap->predF = ttap->predF;
  }
  if(tap->succFPos>=0) {
    //if(seqnum==0 && pos==175) printf("          alarm succF: %i %i %i %i\n", tap->succF, tap->succFPos, seqnum, pos);
    //printf ("2. PRE Pos %i %i \n", tap->predFPos, tap->succFPos);
    if( (tap->succFPos==pos)|| !(tap->state & para->STATE_ORPHANE)) {
      printf("succ ALARM %i %i\n", tap->succFPos, pos);
      exit(99);
    }
    ttap = find_eqc(ap, seqnum, tap->succFPos);
    tap->succF = ttap->succF;
  }
  //if(seqnum==0 && pos==175) printf(" end qgc\n");
  return tap;
}

/**
 * copy alignment
 *
 * doDgc = 0: ignore all backlogdiags and position dg_cont's
 * doDgc = 1: free the target backlogdiags and dg_conts'
 * doDgc = 2: same as 1 but also copy the original backlog diags and dg_conts to the target
struct alignment* copy_alignment( struct alignment *o_algn, struct alignment *algn, char doDgc) {

  struct seq_col *scol = o_algn->scol;
  int slen = scol->length;
  int i,j;
  struct algn_pos *apos;
  struct algn_pos *o_apos;
  struct algn_pos *tpos;
  struct seq *sq;
  struct diag_cont *ptdgc, *tdgc, *o_tdgc;

  if(doDgc>0) {
    tdgc = algn->backlog_diags;
    while(tdgc!=NULL) {
      algn->backlog_diags = algn->backlog_diags->next;
      free(tdgc);
      tdgc = algn->backlog_diags;
    }
    
    if(doDgc>1) {
      o_tdgc = o_algn->backlog_diags;
      ptdgc = NULL;
      while(o_tdgc!=NULL) {
	tdgc = malloc(sizeof(struct diag_cont));
	*tdgc = *o_tdgc;
	if(ptdgc == NULL) {
	  algn->backlog_diags = tdgc;
	} else {
	  ptdgc->next = tdgc;
	}
	ptdgc = tdgc;
      }
    }
  }


  memcpy(algn->seq_is_orphane, o_algn->seq_is_orphane, sizeof(char)*slen);
  algn->total_weight = o_algn->total_weight;
  //  printf(" enter copy\n");
  for(i=0;i<slen;i++) {
    sq = &(scol->seqs[i]);
    for(j=0;j<sq->length;j++) {
      apos = &(algn->algn[i][j]);
      o_apos = &(o_algn->algn[i][j]);

      if(doDgc>0) {
	tdgc = apos->dg_cont;
	while(tdgc!=NULL) {
	  apos->dg_cont = apos->dg_cont->next;
	  free(tdgc);
	  tdgc = apos->dg_cont;
	}
	
	if(doDgc>1) {
	  o_tdgc = o_apos->dg_cont;
	  ptdgc = NULL;
	  while(o_tdgc!=NULL) {
	    tdgc = malloc(sizeof(struct diag_cont));
	    *tdgc = *o_tdgc;
	    if(ptdgc == NULL) {
	      apos->dg_cont = tdgc;
	    } else {
	      ptdgc->next = tdgc;
	    }
	    ptdgc = tdgc;
	  }
	}
      }
      
      if(! (apos->state & para->STATE_ORPHANE)) {
	if(o_apos->state & para->STATE_ORPHANE) {
	  if(!(apos->state & para->STATE_INHERITED)) {
	    //	    printf (" free1\n")
	    //frees += sizeof(int)*slen*2;
	    //printf (" 1. frees: %i %i %i\n",allocs, frees, allocs-frees);

	    free(apos->predF);
	    free(apos->succF);
	    //balance -= 2*sizeof(int)*slen;
	  }
	  apos->predF=NULL;
	  apos->succF=NULL;
	} else {
	  if(apos->state & para->STATE_INHERITED) {
	    if(! (o_apos->state & para->STATE_INHERITED)) {
	      //printf (" 1. malloc: %i %i %i\n",allocs, frees, allocs-frees);
	      apos->predF=malloc(sizeof(int)*slen);
	      apos->succF=malloc(sizeof(int)*slen);
	      //allocss += sizeof(int)*slen*2;
	      //balance += 2*sizeof(int)*slen;
	      if( (apos->predF==NULL) || (apos->succF==NULL)) error("copy_alignment(): (1) Out of memory !");
	      memcpy(apos->predF, o_apos->predF, sizeof(int)*slen);
	      memcpy(apos->succF, o_apos->succF, sizeof(int)*slen);
	    } else {
	      apos->predF=NULL;
	      apos->succF=NULL;
	    }
	  } else {
	    if(! (o_apos->state & para->STATE_INHERITED)) {
	      memcpy(apos->predF, o_apos->predF, sizeof(int)*slen);
	      memcpy(apos->succF, o_apos->succF, sizeof(int)*slen);
	    } else {
	      //	      printf (" free2\n")
	      //frees += sizeof(int)*slen*2;
	      //printf (" 2. frees: %i %i %i\n",allocs, frees, allocs-frees);
	      free(apos->predF);
	      free(apos->succF);
	      //freess += sizeof(int)*slen*2;

	      //balance -= 2*sizeof(int)*slen;
	      apos->predF=NULL;
	      apos->succF=NULL;
	    }
	  }
	}
      } else {
	if( !(o_apos->state & para->STATE_ORPHANE)) {
	  if(o_apos->state & para->STATE_INHERITED) {
	    apos->predF=NULL;
	    apos->succF=NULL;
	  } else {
	    //allocs += sizeof(int)*slen*2;
	    //printf (" 2. malloc: %i %i %i\n",allocs, frees, allocs-frees);
	    apos->predF=malloc(sizeof(int)*slen);
	    apos->succF=malloc(sizeof(int)*slen);
	    //allocss += sizeof(int)*slen*2;
	    //balance += 2*sizeof(int)*slen;
	    if( (apos->predF==NULL) || (apos->succF==NULL)) error("copy_alignment(): (2) Out of memory !");
	    memcpy(apos->predF, o_apos->predF, sizeof(int)*slen);
	    memcpy(apos->succF, o_apos->succF, sizeof(int)*slen);
	  }
	} else {
	  apos->predF=NULL;
	  apos->succF=NULL;
	}
      }
 

      apos->predFPos = o_apos->predFPos;
      apos->succFPos = o_apos->succFPos;


      //printf(" before %i %i\n", i,j);
      apos->eqcRank = o_apos->eqcRank;
      tpos = o_apos->eqcParent;
      apos->eqcParent = &(algn->algn[tpos->row][tpos->col]);
      if( (apos->state & para->STATE_INHERITED)) {
	if( !(o_apos->state & para->STATE_INHERITED)) {
	  //allocs += sizeof(int)+sizeof(char);
	  apos->eqcAlgnPos = malloc(sizeof(int));
	  *(apos->eqcAlgnPos)=j;
	  apos->proceed = malloc(sizeof(char));
	  *(apos->proceed)=0;
	  //allocss += sizeof(int)+sizeof(char);
	}
      } else {
	if((o_apos->state & para->STATE_INHERITED)) {
	  //frees += sizeof(int)+sizeof(char);
	  free(apos->eqcAlgnPos);
	  free(apos->proceed);
	  apos->eqcAlgnPos=NULL;
	  apos->proceed=NULL;
	  //balance -= sizeof(int)+sizeof(char);
	}
      }
      //printf(" after %i %i\n", i,j);
      
      apos->state = o_apos->state;
    }
  }

  //  printf(" leave copy\n");
  return algn;
}
 
 */

/**
 * adds the given diagional to the given alignment and updates the
 * datastructure (i.e. frontiers). The given diag must be consistent
 * to the given alignment !
 */
char align_diag(struct alignment *algn, struct scr_matrix *smatrix, struct diag* dg) {

  char alignedSomething = 0;
  int i,j,k;
  char al;
  if(dg->multi_dg) {
    //return 0;
    for(i=0;i<dg->multi_length;i++) {
      if(dg->multi_cont[i]!=NULL) {
	//printf(" before %f %f\n", dg->total_weight, dg->multi_cont[i]->weight);
	al = align_diag(algn,smatrix, dg->multi_cont[i]);
	if(al) algn->total_weight -= dg->multi_cont[i]->total_weight;
	alignedSomething = alignedSomething || al;
	
	//printf(" after \n");
      }
    }
    return alignedSomething;
  }
  
  //if((dg->length==0) ) return 0;
  if((dg->length==0) || (!dg->meetsThreshold && !dg->anchor)) return 0;
  /*
    printf(" dg %i %i %i\n", dg, dg->multi_dg, dg->length);
  char adapted = adapt_diag(algn, smatrix, dg);
  if(adapted) {
    printf(" dg %i %i %i\n", dg, dg->multi_dg, dg->length);
    error(" inconsistent diag!\n");
  }
  */

  struct seq_col *scol = algn->scol;
  int s1 = dg->seq_p1.num;
  int s2 = dg->seq_p2.num;

  //printf("%i %i\n", s1,s2);

  //char o1 = algn->seq_is_orphane[s1];
  //char o2 = algn->seq_is_orphane[s2];

  int length = dg->length;
  int sp1 = dg->seq_p1.startpos;
  int sp2 = dg->seq_p2.startpos;

  struct algn_pos **ap=algn->algn;
  struct algn_pos *tp;
  
  struct algn_pos *apos1, *apos2, *tpos;

  int p1, p2,plen;
  int p;
  struct seq *sq = scol->seqs;
  int s, slen = scol->length;
  int *oldpredF, *oldsuccF, *otherpredF, *othersuccF;
  int pF,sF;
  
  int *c2n = smatrix->char2num;
  int *sdata = smatrix ->data;
  char *data1 = dg->seq_p1.sq->data;
  char *data2 = dg->seq_p2.sq->data;
  int smatrixlen = smatrix->length;
  int a1,a2;
  int score1;
  
  char seenNonOrphane;
  int opos;
  int ms;
  char skip = 0;
  char found;

  struct diag_cont *dgc, *dgc_next;
  double ttim;
  /*
  int nextpos[slen];
  int firstpos;
  int oldpos;
  char firstRound;
  */
  for(i=0;i<length;i++) {
    p1 = sp1+i; p2 = sp2+i;
    //    printf("pos %i %i %i %i\n",s1,s2,p1,p2);
    //printf(" %i \n", scol->length);
    skip = 0;
    
    //if(sp1==30) printf("%i %i %i %i %i \n", p1,p2,dg->length, dg->seq_p1.sq->length, dg->seq_p2.sq->length);
    /*
    //    if(dg->onlyOverThres==1) {
      a1 = c2n[data1[p1]];
      a2 = c2n[data2[p2]];
      score1 = sdata[smatrixlen*a1+a2];
      
      if( (score1<=3) && (i>length/3)) {
	skip = 1;
      }
      //}
      */
    //    printf(" %i %i %i\n",p1,p2,skip);
	//} else {
    /*
    a1 = c2n[data1[p1]];
    a2 = c2n[data2[p2]];
    score1 = sdata[smatrixlen*a1+a2];
    
    if(score1<=4) {
      skip = 1;
    }
    */
    //if( (i==(length-1)) ){

    //ttim = clock();
    //tim += (clock()-ttim)/CLOCKS_PER_SEC;
    // printf(" tim %f\n", tim);


      //printf(" diag n1=%i s1=%i n2=%i s2=%i l=%i\n", s1,sp1,s2,sp2,length);
      //allocss += 2*sizeof(struct diag_cont );
      //}
    // TODO: tune ration is 1:45 (minor)
    /*
    dgc_next = ap[s1][p1].dg_cont;
    dgc = malloc(sizeof(struct diag_cont));
    dgc->dg = dg;
    dgc->next = dgc_next;
    ap[s1][p1].dg_cont = dgc;
    //if((s1==0) && (p1==133)) printf(" dddgc %i %i\n", dgc, dgc->next);
    
    dgc_next = ap[s2][p2].dg_cont;
    dgc = malloc(sizeof(struct diag_cont));
    dgc->dg = dg;
    dgc->next = dgc_next;
    ap[s2][p2].dg_cont = dgc;
    */


    if(! skip) {
      //printf("apos1 %i %i\n", s1, p1);
      apos1 = find_eqc(ap, s1,p1);//&(ap[s1][p1]);
      //printf("apos2 %i %i %i\n", s2, p2, scol->seqs[s2].length);
      apos2 = find_eqc(ap, s2,p2); //&(ap[s2][p2]);
      //printf("after apos2 %i %i\n", s2, p2);

      oldpredF = apos1->predF;
      oldsuccF = apos1->succF;
      if(oldpredF!=NULL && oldsuccF!=NULL)
	if(oldpredF[s2]==p2 && oldsuccF[s2]==p2) skip = 1;

      if(para->DEBUG>4) {
	if(!skip) {
	  if(oldpredF!=NULL) if(oldpredF[s2]>=p2) {
	    printf(" Incons1 %i %i %i\n", s2, p2,oldpredF[p2]);
	    error(" ERROR");
	  }
	  if(oldsuccF!=NULL) if(oldsuccF[s2]<=p2) {
	    printf(" Incons2 %i %i %i\n",s2, p2,oldsuccF[p2]);
	    error(" ERROR");
	  }
	  oldpredF=apos2->predF;
	  oldsuccF=apos2->succF;
	  if(oldpredF!=NULL) if(oldpredF[s1]>=p1) {
	    printf(" Incons3 %i %i %i\n", s1, p1,oldpredF[p1]);
	    error(" ERROR");
	  }
	  if(oldsuccF!=NULL) if(oldsuccF[s1]<=p1) {
	    printf(" Incons4 %i %i %i\n",s1, p1,oldsuccF[p1]);
	    error(" ERROR");
	  }
	}
      }
    }
      //}

    if(! skip) {
      alignedSomething = 1;
      for(k=0;k<2;k++) {
	tpos = (k==0 ? apos1 : apos2);
	//printf("tpos %i\n", tpos);
	oldpredF = tpos->predF;
	oldsuccF = tpos->succF;
	otherpredF = (k==0 ? apos2->predF : apos1->predF);
	othersuccF = (k==0 ? apos2->succF : apos1->succF);
	//printf("pre isorphane %i\n", tpos);
	if(tpos->state & para->STATE_ORPHANE) {
	  //if(scol->length>21) printf("isorphane %i %i\n", tpos,sizeof(int)*scol->length);
	  //printf("step 1\n");
	  // printf (" 3. malloc: %i\n",sizeof(int)*slen*2);
	  //allocs += sizeof(int)*slen*2;
	  //	  printf (" 3. malloc: %i %i %i\n",allocs, frees, allocs-frees);
	  //balance += sizeof(int)*slen*2;
	  tpos->predF = malloc(sizeof(int)*scol->length);
	  //if(k==0) printf("           apos1->predF = %i\n",tpos->predF),
	  //printf("pre succForphane %i %i\n", tpos->succF, scol->length);
	  //printf(" step 2\n");
	  tpos->succF = malloc(sizeof(int)*scol->length);
	  //allocss += 2*sizeof(int )*scol->length;
	  //printf("  step 3\n");
	  //printf("succForphane %i\n", tpos);
	  if( (tpos->predF==NULL) || (tpos->succF==NULL)) error("align_diag(): (1) Out of memory !");
	  
	}
	//      printf("init loop %i\n", tpos);
	
	for(j=0;j<scol->length;j++) {
	  pF = -1;
	  sF = sq[j].length;
	  
	  
	  // propagate predF and succF
	  if(oldpredF!=NULL && oldpredF[j]>pF) pF = oldpredF[j];
	  //if(j==0) printf("1 pf=%i\n", pF);
	  if(otherpredF!=NULL && otherpredF[j]> pF) pF = otherpredF[j];
	  //if(j==0) printf("2 pf=%i\n", pF);
	  
	  if(oldsuccF!=NULL && oldsuccF[j]<sF) sF = oldsuccF[j];
	  if(othersuccF!=NULL && othersuccF[j]<sF) sF = othersuccF[j];
	  //}
	  if(j==s1) {
	    pF  =  p1;
	    sF  =  p1;
	  } else if(j==s2) {
	    pF  =  p2;
	    sF  =  p2;
	  }
	  /*
	    if(pF > sF) {
	    if(oldpredF!=NULL) printf(" PRE 1. ALARM oldpredF: %i \n",oldpredF[j] );
	    if(oldsuccF!=NULL) printf(" PRE 1. ALARM oldsuccF: %i \n",oldsuccF[j] );
	    if(otherpredF!=NULL) printf(" PRE 1. ALARM otherpredF: %i \n",otherpredF[j] );
	    if(othersuccF!=NULL) printf(" PRE 1. ALARM othersuccF: %i \n",othersuccF[j] );
	    
	    printf("1. ALARM j=%i %i %i %i %i %i %i %i\n",
	    j,
	    (k==0 ? s1 : s2), 
	    (k==0 ? p1 : p2),
	    (k==0 ? p2 : p1), 
	    oldpredF!=NULL ? oldpredF[k==0 ? s2 : s1] : -2, 
	    oldsuccF!=NULL ? oldsuccF[k==0 ? s2 : s1]: 99999, 
	    tpos->predFPos, 
	    tpos->succFPos);
	    exit(1);
	    }
	  */
	  //if(pF==sF && pF==0 && j==0) printf("pf=%i sf=%i j=%i s1=%i s2=%i p1=%i p2=%i iso=%i opf=%i osf=%i otpf=%i otsf=%i\n", pF,sF,j,s1,s2,p1,p2,tpos->isOrphane,oldpredF, oldsuccF, otherpredF, othersuccF);
	  tpos->predF[j]=pF;
	  tpos->succF[j]=sF;
	}
	//if(s1==0 && p1==0) printf(" SET IT 1\n");
	//if(s2==0 && p2==0) printf(" SET IT 2\n");
	//if(tpos->state & para->STATE_ORPHANE) 
	tpos->state = tpos->state & (tpos->state ^ para->STATE_ORPHANE);
	tpos->predFPos = -1;
	tpos->succFPos = -1;
      }
      //printf("end pre/succ %i\n", tpos);
      apos1->predF[s1]=p1;
      apos1->succF[s1]=p1;
      apos2->predF[s2]=p2;
      apos2->succF[s2]=p2;
      
      apos1->predF[s2]=p2;
      apos1->succF[s2]=p2;
      apos2->predF[s1]=p1;
      apos2->succF[s1]=p1;
      
      
      
      if(apos2->eqcRank< apos1->eqcRank) {
	apos2->eqcParent = apos1;
      } else {
	apos1->eqcParent = apos2;
	if(apos2->eqcRank==apos1->eqcRank) 
	  apos2->eqcRank++;
      }

      //printf("end ranking %i, %i    %i\n", s1,p1,apos1->predF);
      apos1 = find_eqc(ap, s1,p1);//&(ap[s1][p1]);
      //printf("end first egc %i\n", tpos);
      apos2 = find_eqc(ap, s2,p2); //&(ap[s2][p2]);
      //printf("end second egc %i\n", tpos);
      
      // update the other affected sites
      for(ms=0;ms<slen;ms++) {
	// spaeter ueberlegen: if((ms!=s1 && ms!=s2) || ( (ms==s1 || ms==s2) && (i==0 || i== (length-1) ))) {
	//	if(apos1->predF[ms]==apos1->succF[ms]) {
      
	p = apos1->predF[ms]; // -( (apos1->predF[ms]==apos1->succF[ms]) ? 1 : 0); // GOGOGOGO 
	opos=apos1->predF[ms];
	//    printf("SUPERPRE WHILE %i %i %i\n",k,s1,s2);
	tp = ap[ms];
	//printf("PRE WHILE %i\n",k);
	seenNonOrphane=0;
	found = 1;
	//firstRound = 1;
	//firstpos = 0;
	while( (p>=0)&& found) { // && tp[p].isOrphane) {
	  //printf(" WHILE %i\n",p); 
	  if(tp[p].state & para->STATE_ORPHANE) {
	    if( (! seenNonOrphane)) {
	      //	      if(ms==0 && p==0) printf("setting s1=%i p1=%i iso=%i ms=%i p=%i opos=%i\n",s1,p1, tp[p].isOrphane, ms,p,opos);
	      tp[p].succFPos = opos;
	    } else p = tp[p].predFPos+1;
	  } else {
	    if( (p!=opos) || (apos1->succF[ms]!=apos1->predF[ms])) 
	      seenNonOrphane = 1;
	    //if(p==442) printf("  pre find %i %i %i\n",s1, ms,p);
	    tpos = find_eqc(ap, ms,p);//&(ap[s1][p1]);
	    //if(p==442)printf("  post find %i %i\n",ms,p);
	    if(! (tpos->state & para->STATE_ORPHANE) && ! (tpos->state & para->STATE_INHERITED)) {
	      //printf(" %i %i %i %i\n",ms,p,s1,p1);
	      if(seenNonOrphane) found = 0;
	      for(s=0;s<slen;s++) {
		if(tpos->succF[s]>apos1->succF[s]) {
		  tpos->succF[s]=apos1->succF[s];
		  found = 1;
		}
	      }
	      /*
	      oldpos = firstpos;
	      s = firstpos;
	      while(s<slen) {
		if(firstRound) { 
		  nextpos[s]=s+1;
		} 
		if(tpos->succF[s]>apos1->succF[s]) {
		  tpos->succF[s]=apos1->succF[s];
		  found = 1;
		  oldpos = s;
		} else {
		  if(! firstRound) {
		    nextpos[oldpos] = nextpos[s];
		    if(oldpos==firstpos) {
		      firstpos=nextpos[oldpos];
		      oldpos = firstpos;
		    }
		  }
		}
		if(firstRound) {
		  oldpos = s;
		  s++;
		} else {
		  s = nextpos[s];
		}
	      }
	      */
	    }
	  }
	  //firstRound = 0;
	  p--;
	}
	
	//printf("END WHILE\n");
	
	//printf(" PRE 2 %i %i %i\n", apos1->predF, apos1->succF,ms );
	p = apos1->succF[ms]; // +( (apos1->predF[ms]==apos1->succF[ms]) ? 1 : 0); // GOGOGOGO ;
	opos= apos1->succF[ms];
	plen = scol->seqs[ms].length;
	seenNonOrphane=0;
	//printf("2. PRE WHILE %i\n",k);
	// if(opos>=plen) opos = -1;
	found = 1;
	while( (p<plen ) && found) {
	  //      step = (int)(p-tep);
	  if(tp[p].state & para->STATE_ORPHANE) {
	    if( (! seenNonOrphane)) {
	      /*
	      if(p==opos) {
		printf("ALARM set predFPos s1=%i s2=%i p1=%i p2=%i ms=%i sF=%i pF=%i p=%i opos=%i iso=%i %i %i %i\n",
                   s1,s2,p1,p2,ms, apos1->succF[ms],
                   apos1->predF[ms],p,opos, tp[p].state,(int) &tp[p],(int) apos1, (int) apos2);
		exit(98);
	      }
	    */
	      tp[p].predFPos = opos;
	    } else {
	      if(tp[p].succFPos < 0) 
		p = plen+1;
	      else
		p = tp[p].succFPos-1;
	    }

	  } else {
	    if( (p!=opos)|| (apos1->succF[ms]!=apos1->predF[ms])) 
	      seenNonOrphane = 1;
	    //printf("  pre find %i %i\n",ms,p);
	    tpos = find_eqc(ap, ms,p);//&(ap[s1][p1]);
	    //printf("  end find\n");
	    if(! (tpos->state & para->STATE_ORPHANE)  && !(tpos->state & para->STATE_INHERITED)) {
	      if(seenNonOrphane) found = 0;
	      for(s=0;s<slen;s++) {
		//		printf("   s=%i slen=%i tpos->predF=%i\n",s,slen,tpos->predF);
		
		if( tpos->predF[s]<apos1->predF[s]) {
		  //printf("   inner before s=%i slen=%i\n",s,slen);
		  tpos->predF[s]=apos1->predF[s];
		  //printf("   inner after s=%i slen=%i\n",s,slen);
		  found = 1;
		}
	      }
	    }
	  }
	  p++;
	}
	//printf("2. END WHILE %i\n",k);
      }

    }
  }
  
  
  // printf("s1: %i s2: %i\n", algn->seq_is_orphane[s2],s2);
  algn->seq_is_orphane[s1]=0;
  algn->seq_is_orphane[s2]=0;

  int maxposs1 = dg->seq_p1.startpos + dg->length-1;
  int maxposs2 = dg->seq_p2.startpos + dg->length-1;
  if(dg->seq_p1.sq->max_seen < maxposs1)
    dg->seq_p1.sq->max_seen = maxposs1;

  if(dg->seq_p2.sq->max_seen < maxposs2)
    dg->seq_p2.sq->max_seen = maxposs2;

  //if(alignedSomething) {
  /*
    if(algn->aligned_diags_amount >= algn->max_aligned_diags_amount) {
      algn->max_aligned_diags_amount = algn->aligned_diags_amount + 16;
      algn->aligned_diags = realloc(algn->aligned_diags, 
				    sizeof(struct diag*)*algn->max_aligned_diags_amount);
      if(algn->aligned_diags==NULL) error(" align_diag(): Out of Memory!");
    }
    algn->aligned_diags[algn->aligned_diags_amount++] = dg;
  */
    /*
  } else {
    
    dgc = malloc(sizeof(struct diag_cont));
    dgc->dg = dg;
    dgc->next = NULL;
    algn->backlog_diags = enter_sorted(algn->backlog_diags, dgc);
    
  }
    */
  //printf(" diaglen: %i\n", algn->aligned_diags_amount);
    if(! dg->anchor) algn->total_weight += dg->weight;//*dg->weight_fac;
    return(alignedSomething);
}





/**
 *
 * prepares the alignment: calculates the position number of each residue 
 *  
 */
void prepare_alignment(struct alignment *algn) {
  struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  unsigned int i,j,s,ts, hasmore, max;
  int *predF, *succF;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;
  int  tproc[slen];
  //  char proceed[slen];

  for(i=0;i<slen;i++) {
    tproc[i]=0;
  }

  //
  // prepare block
  //
  struct algn_pos *ap1;
  hasmore = 1;
  char alarmHasProceed;

  for(j=0;hasmore;j++) {
    hasmore = 0;
    //memset(proceed, 0, slen*sizeof(char));
    //    printf("Position: %i\n", j);
    //    for(k=0;k<2;k++) {
    for(s=0;s<slen;s++) {
      sq = &scol->seqs[s];
      if(tproc[s]<sq->length) {
        ap1 = find_eqc(ap,s,tproc[s]);
	*ap1->proceed = 1;
	
      }
    }
    for(s=0;s<slen;s++) {
      sq = &scol->seqs[s];
      if(tproc[s]<sq->length) {
	//printf(" DO IT %i %i %i %i\n",j,s,tproc[s], *ap1->eqcAlgnPos);
        ap1 = find_eqc(ap,s,tproc[s]);
//		printf("alig.c ap1 = %d\tap = %d\n",ap1,ap);
	
	if(j>=*ap1->eqcAlgnPos) {
	  predF = ap1->predF;
	  succF = ap1->succF;
	  
	  *ap1->eqcAlgnPos=j;// *tap1->eqcAlgnPos;
	  if(predF!=NULL) {
	    for(ts=0;ts<slen;ts++) {
	      //printf(" MIST %i %i %i %i %i %i %i\n",j,s,ts,tproc[s], tproc[ts], predF[ts], succF!=NULL ? succF[ts]: 99999);
	      if( (tproc[ts]<predF[ts]) && !(tproc[ts]==predF[ts] && succF!=NULL && succF[ts]==predF[ts])) {
		*ap1->eqcAlgnPos=j+1;// *tap1->eqcAlgnPos;
		//printf(" 2. MIST %i %i %i %i %i %i %i\n",j,s,ts,tproc[s], tproc[ts], predF[ts], succF!=NULL ? succF[ts]: 99999);
		*ap1->proceed = 0;
		break;
	      } else {
		/*
		 *ap1->eqcAlgnPos=j;// *tap1->eqcAlgnPos;
		 *ap1->proceed = 1;
		 */
	      }
	    }
	  }
	} else {
	  *ap1->proceed = 0;
	}
	/*
	  if(j>=143) {
	  printf("ALARM j=%i s=%i tproc[s]=%i algnPos=%i\n",j,s,tproc[s],*ap[s][tproc[s]].eqcAlgnPos);
	  }
	*/
      }
    }
    alarmHasProceed=0;
    for(s=0;s<slen;s++) {
      sq = &scol->seqs[s];
      if(tproc[s]<sq->length) {
	ap1 = find_eqc(ap,s,tproc[s]);
	if(*ap1->proceed) {
	  alarmHasProceed = 1;
	  tproc[s]++;
	}
      }
      if(tproc[s]<sq->length) hasmore = 1;
      //printf("%i %i\n", sq->length,tproc[s]);
    }
    if(! alarmHasProceed && hasmore) {
      printf("IO ALARM! %i\n",j);
      exit(1);
      hasmore=0;
    }
    if(!hasmore) max = j+1;
  }
  algn->max_pos= max;
}



