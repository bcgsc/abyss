#include <time.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <string.h> 
#include <ctype.h>
 
#include "parameters.h" 
#include "struct.h" 
#include "io.h"
 
extern void error(char *message); 
extern void merror(char *msg1, char *msg2); 
 
extern struct algn_pos *find_eqc(struct algn_pos **ap, int seqnum, int pos); 
extern void prepare_alignment(struct alignment *algn);
 
/** 
 * 
 * diag.c: Creation of diagonals & calculation of scores and weights  
 * 
 */ 
 
 
/** 
 * factory method that creates a diagonal from the given sequence parts 
 * 
 * The pointer returned has to be deallocted explicitely from memory. 
 */ 
struct diag* create_diag(int n1, struct seq* sq1, unsigned int sp1,  
			 int n2, struct seq* sq2, unsigned int sp2, 
			 int dlength) { 
  struct diag* dg = malloc(sizeof(struct diag)); 
  if(dg==NULL) error("create_diag(): Out of Memory !");
 
  if(sq1->length < sp1+dlength) { 
    printf(" startpos=%i diaglength=%i seqlength=%i\n",sp1,dlength,sq1->length);
    merror("create_diag(): startpos+diaglength exceeds sequence length in diag ", sq1->name); 
  } 
  if(sq2->length < sp2+dlength) { 
    printf(" startpos=%i diaglength=%i seqlength=%i\n",sp2,dlength,sq2->length);
    merror("create_diag(): startpos+diaglength exceeds sequence length in diag ", sq2->name); 
  } 
 
  dg->seq_p1.num = n1; 
  dg->seq_p1.sq = sq1; 
  dg->seq_p1.startpos = sp1; 
  //dg->seq_p1.leftmargin = sp1;
  //dg->seq_p1.rightmargin = sq1->length - dlength - sp1;

  dg->seq_p2.num = n2; 
  dg->seq_p2.sq = sq2; 
  dg->seq_p2.startpos = sp2; 
  //dg->seq_p1.leftmargin = sp2;
  //dg->seq_p1.rightmargin = sq2->length - dlength - sp2;
 
  dg->pred_diag = NULL; 
  dg->col_pred_diag = NULL; 
  dg->length = dlength; 
  //dg->onlyOverThres= 0; 
  dg->score= -1; 
  dg->orig_score= -1; 
  dg->weight = 0.0; 
  dg->weight_sum = 0.0; 
  dg->ov_weight = 0.0; 
  dg->weight_fac = 1.0; 
  dg->anchor = 0;
  dg->marked = 0;
  dg->multi_dg = 0;
  dg->multi_cont = NULL;
  dg->multi_length = 0;
  dg->meetsThreshold = 0;

  dg->neighbours = NULL;
  dg->degree = 0;
  dg->max_degree= 0;
  return dg; 
} 
 
/** 
 * frees the memory of the given diagonal and the included seq_parts 
 */ 
void free_diag(struct diag* dg) { 
  if(dg->multi_dg) {
    int i;
    for(i=0;i<dg->multi_length;i++) {
      if(dg->multi_cont[i]!=NULL) free_diag(dg->multi_cont[i]);
    }
    free(dg->multi_cont);
  }
  free(dg); 
} 

/**
 * frees the (sub) tree of the given gt_node
 */
void free_gt_node(struct gt_node *gtn) {
  if(! gtn->isLeaf) {
    free_gt_node(gtn->succ1);
    free_gt_node(gtn->succ2);
  }
  free(gtn->seq_num);
  free(gtn);
}
 
/** 
 * calculuates "m over n" 
 */ 
unsigned long binomial(long m, long n) { 
  double result=1.0; 
  long i; 
  for(i=0;i<n;i++) { 
    result *= ((double)(m-i))/(double)(n-i); 
  } 
  return (unsigned long)result; 
} 
 
/** 
 * creates temporary probability distribution 
 */ 
long double **create_tmp_pdist(struct prob_dist *pdist) { 
  int length = pdist->max_dlen; 
  struct scr_matrix *smatrix = pdist->smatrix; 
 
  long double **dist = calloc(length+1, sizeof(long double *)); 
  if(dist==NULL) error("create_tmp_pdist(): (1) Out of memory when allocating data !"); 
 
  int i; 
  long mxscr, sm_max_scr=smatrix->max_score; 
  for(i=0;i<=length;i++) { 
    mxscr = i *sm_max_scr; 
    dist[i] = calloc(mxscr+1, sizeof(long double )); 
    if(dist[i]==NULL) error("create_tmp_pdist(): (3) Out of memory at iteration" ); 
  } 
  return dist; 
} 
 
/** 
 * frees temporary probability distribution 
 */ 
void free_tmp_pdist(long double **dist, int length) { 
  int i; 
  for(i=0;i<=length;i++) { 
    free(dist[i]); 
  } 
  free(dist); 
} 
 
void fill_tmp_pdist(struct prob_dist *pdist, long double **tmp_dist, int slen1, int slen2) { 
  unsigned int length = pdist->max_dlen; 
  struct scr_matrix * smatrix = pdist->smatrix; 
 
 
  unsigned int i; 
  long mxscr, sm_max_scr=smatrix->max_score,scr; 
 
  long double factor, np, np2,prob; 
  long double seq_factor= (((long double)slen1)*(slen2)); 
 
  for(i=1;i<=length;i++) { 
    mxscr = i *sm_max_scr; 
     
    factor = (long double)(seq_factor)/(long double)(4.0*i*i); // original ! 
     
    for(scr=0;scr<=mxscr;scr++) { 
      prob = pdist->data[i][scr]; 
 
      np2 =prob * (factor); 
      if(np2>=para->DIAG_CALC_WEIGHT_THRESHOLD) { // curent 
	np = (long double)1.0- pow(1.0-prob,factor); // current 
      } else { 
	np = np2; 
      } 
      tmp_dist[i][scr] = -log(np); 
    } 
  } 
} 
 
/** 
 * calculates the score of the given diag by using the given score matrix. the  
 * resulting score is stored within the diag 
 * omitScore = -1: score calculation but weight interpolation with seqlen = 100 
 * omitScore = 0:  normal 
 * omitScore = 1:  no score calculation 
 */ 
static void real_calc_weight(struct diag* dg, struct scr_matrix* smatrix,  
		 struct prob_dist *pdist, char omitScore, long double **tmp_dist, struct alignment *algn ) { 
   
  if(dg->multi_dg) {
    int i;
    dg->weight = 0.0;
    dg->meetsThreshold = 0;
    dg->length = 0;
    for(i=0;i<dg->multi_length;i++) {
      if(dg->multi_cont[i]!=NULL) {
	//printf("   before real calc %i %i %f\n", dg->length, dg->score, dg->weight);
	real_calc_weight(dg->multi_cont[i],smatrix, pdist, omitScore, tmp_dist, algn);
	if(dg->multi_cont[i]->length > dg->length) 
	  dg->length = dg->multi_cont[i]->length;

	//printf("   real calc %i %i %f\n", dg->length, dg->score,dg->weight);
	dg->weight += dg->multi_cont[i]->weight;
	dg->meetsThreshold = dg->meetsThreshold | dg->multi_cont[i]->meetsThreshold;
	//printf(" readjusted %i %i %i %i %i %i %f %f %i\n", dg, dg->multi_cont[i],dg->multi_cont[i], dg->multi_dg, dg->multi_length, dg->weight, dg->multi_cont[i]->weight, dg->length);
	
      }
    }
    if(dg->length<=0) dg->meetsThreshold = 0;
    dg->total_weight = dg->weight;
    return;
  }

  if(dg->length==0) { 
    dg->score = 0; 
    dg->weight = 0.0; 
    dg->meetsThreshold = 0; 
    return; 
  } 
   
  unsigned int len = dg->length; 
 
  int pos; 
  long double np=0.0,np2; 
   
  if(omitScore<=0) { 
    unsigned int sp1=dg->seq_p1.startpos; 
    unsigned int sp2=dg->seq_p2.startpos; 
    char *data1 = dg->seq_p1.sq->data; 
    char *data2 = dg->seq_p2.sq->data; 
    int a1, a2; 
    int *c2n = smatrix->char2num; 
    int *sdata = smatrix ->data; 
    int slen = smatrix->length; 
 
    dg->score = 0; 
    for(pos=0;pos<len;pos++) { 
      a1 = c2n[(int) data1[sp1+pos]]; 
      a2 = c2n[(int) data2[sp2+pos]]; 
      dg->score+=(long)sdata[slen*a1+a2]; 
    } 
  } 
   
  long double prob; 
 
  long double factor; 
 
  dg->meetsThreshold = 0; 
   
   
  if(dg->score <= pdist->smatrix->avg_sim_score*dg->length) { 
    dg->total_weight = 0.0; 
    dg->ov_weight = 0.0; 
    dg->weight_fac = 1.0; 
    dg->weight = 0.0; 
    return; 
  } 
   
  // interpolate only for splitted diags that fall beyond the threshold
  /*
  if( (dg->orig_score>0)) {
    //if( (dg->weight<=0.0) && (dg->orig_score>0)) {
    dg->weight = ((double)dg->score)/((double)dg->orig_score)*dg->weight_sum;
    //    printf(" ws %.20f %.20f \n", dg->weight, dg->weight_sum);
  //if(dg->weight_sum>0.0) dg->weight = dg->weight_sum;
  }else */ 
  if(len<=pdist->max_dlen) { 
 
    if(tmp_dist==NULL) { 
      prob  =pdist->data[len][dg->score]; 
 
      /*
      if(omitScore>=0) { 
	factor = (long double)((len + dg->seq_p1.leftmargin + dg->seq_p1.rightmargin )* (len + dg->seq_p2.leftmargin + dg->seq_p2.rightmargin ))/(long double)(4.0*len*len); // original 
      } else { 
	factor = (long double)(10000.0)/(long double)(4.0*len*len); 
      } 
      */
      //printf(" %i %i\n", dg->seq_p1.sq->length, len+dg->seq_p1.leftmargin+dg->seq_p1.rightmargin);
      
      if(omitScore>=0) { 
	factor = (long double)((dg->seq_p1.sq->length )* (dg->seq_p2.sq->length ))/(long double)(4.0*len*len); // original 
      } else { 
	factor = (long double)(10000.0)/(long double)(4.0*len*len); 
      } 
      
      np2 =prob * (factor); 
      if(np2>=para->DIAG_CALC_WEIGHT_THRESHOLD) {  
	np = (long double)1.0- pow(1.0-prob,factor); // current 
      } else { 
	np = np2; 
      } 
      dg->weight = -log(np);  
    } else { 
      dg->weight = tmp_dist[len][dg->score]; 
    } 
  }

  dg->total_weight = (dg->weight);//+dg->ov_weight);//* dg->weight_fac; 
  if( (dg->length >= para->DIAG_MIN_LENGTH) && (dg->weight > para->DIAG_THRESHOLD_WEIGHT)) { 
    dg->meetsThreshold = 1; 
  } else { 
  } 
} 
 
void calc_weight(struct diag* dg, struct scr_matrix* smatrix,  
		 struct prob_dist *pdist) { 
  real_calc_weight(dg, smatrix, pdist, 0,NULL,NULL); 
} 
 
 
 
/** 
 * calculates the overlap weight for the given diag 
 */ 
void calc_ov_weight(struct diag* dg, struct diag_col *dcol, struct scr_matrix* smatrix,  
		    struct prob_dist *pdist) { 
  int sn1 = dg->seq_p1.num; 
  int sn2 = dg->seq_p2.num; 
  int snt, sn; 
  struct seq *seq1 = dg->seq_p1.sq; 
  struct seq *seq2 = dg->seq_p1.sq; 
  struct diag* tdg = create_diag(sn1, seq1,0, 
				   1, seq2,0,0);   
  struct diag* dg2; 
  struct simple_diag_col *sdcol; 
 
  int i,j, slen=dcol->seq_amount, dlen; 
  int sp1 = dg->seq_p2.startpos,tsp1; 
  int tep1,t; 
  double w; 
  struct seq_part *seq_p, *d_seq_p1, *d_seq_p2; 
  dg->ov_weight = 0.0; 
  int tstartpos; 
  if(dg->length >0) { 
    for(sn=0;sn<2;sn++) { 
      tstartpos = (sn==0 ? dg->seq_p1.startpos : dg->seq_p2.startpos); 
      tdg->seq_p1.sq = (sn==0 ? seq1 : seq2); 
      tdg->seq_p1.num = (sn==0 ? sn1 : sn2);; 
      // other way 
      seq_p = (sn==0 ? &dg->seq_p2 : &dg->seq_p1); 
      tsp1 = seq_p->startpos; 
      tep1 = seq_p->startpos+dg->length-1; 
 
      snt = (sn==0 ? sn2 : sn1); 
      //printf(" slen %i\n", slen);
      for(i=0; i<slen;i++) { 
	if((i!=sn1) && (i!=sn2)) { 
	  //printf("OV %i %i!\n",i, slen); 
	   
	  sdcol = dcol->diag_matrix[(snt < i) ? (snt*slen + i) : (i*slen+snt)]; 
	  
	  tdg->seq_p2.num=i; 
	  dlen = sdcol->length; 
	  for(j=0;j<dlen;j++) { 
	    dg2 = sdcol->data[j]; 
	    if(snt<i) { 
	      d_seq_p1 = &dg2->seq_p1; 
	      d_seq_p2 = &dg2->seq_p2; 
	    } else { 
	      d_seq_p1 = &dg2->seq_p2; 
	      d_seq_p2 = &dg2->seq_p1; 
	    } 
	    if(j==0) { 
	      tdg->seq_p2.sq = d_seq_p2->sq; 
	    } 
	    if(dg2->length >0) { 
	      if(d_seq_p1->startpos>tsp1) tsp1 = d_seq_p1->startpos; 
	      t=d_seq_p1->startpos+dg2->length-1; 
	      if(t<tep1) tep1 = t; 
	      if(tsp1<=tep1) { 
		//tdg->seq_p2.sq=dg2->seq_p2.sq; 
		tdg->seq_p1.startpos =  tstartpos + tsp1- sp1; 
		tdg->seq_p2.startpos = d_seq_p2->startpos + tsp1- d_seq_p1->startpos; 
		 
		tdg->length = tep1-tsp1+1; 
		//real_calc_weight(tdg, smatrix, pdist,-1,NULL,NULL); 
		real_calc_weight(tdg, smatrix, pdist,0,NULL,NULL); 
		if(tdg->meetsThreshold ) { 
		  w = tdg->weight; 
		  //printf("add %.20f\n",w); 
		  dg->ov_weight += w/2.0; 
		  //dg2->ov_weight += w; 
		} 
	      } 
	    } 
	  } 
	} 
      } 
    } 
  } 
  dg->total_weight = (dg->weight);//+dg->ov_weight)*dg->weight_fac;// + dg->ov_weight; 
  free_diag(tdg); 
} 
 
 
/** 
 * creates the collection of all diags  
 * 
 * The pointer returned (and the ones included in the struct)  
 * has to be deallocted explicitely from memory. 
 */ 
struct diag_col* create_diag_col(int seq_amount) { 
  struct diag_col *dc = calloc(1, sizeof(struct diag_col)); 
  if(dc==NULL) error("create_diag_col(): (1) Out of memory !"); 
 
  //printf("go for k\n"); 
  dc->diag_matrix = malloc(seq_amount*seq_amount* 
			   sizeof(struct simple_diag_col *)); 
  //printf("2go for k\n"); 
  dc->gt_root = NULL;
  if(dc->diag_matrix==NULL) error("create_diag_col(): (2) Out of memory !"); 
  return dc; 
} 
 
/** 
 * frees a diagcol and all data included in it 
 */ 
void free_diag_col(struct diag_col* dcol) { 
  int s1,s2,sl=dcol->seq_amount; 
  // printf("damount: %i\n", dcol->diag_amount); 
  //printf("--------------------------------\n"); 
  for(s1 =0;s1<dcol->diag_amount;s1++) { 
    //    print_diag(dcol->diags[s1]); 
    //printf(" NO FREEEEEEEEEE %i %i\n", s1,dcol->diags[s1]); 
    free_diag(dcol->diags[s1]); 
  } 
  free(dcol->diags); 
  for(s1=0;s1<sl;s1++) { 
    for(s2=s1+1;s2<sl;s2++) { 
      free(dcol->diag_matrix[s1+sl*s2]->data); 
      free(dcol->diag_matrix[s1+sl*s2]); 
    } 
  } 
  free(dcol->diag_matrix); 
  if(dcol->gt_root!=NULL) free_gt_node(dcol->gt_root);
  free(dcol); 
} 
 

/** 
 * finds all relevant multi diags with respect to the given alignment
 * and two groups of sequences that come from the guided method
 * 
 * The pointer returned (and the ones included in the struct)  
 * has to be deallocted explicitely from memory. 
 */ 
inline struct simple_diag_col* find_diags_guided(struct scr_matrix *smatrix,  
				struct prob_dist *pdist, struct gt_node* n1,  
				struct gt_node* n2, struct alignment *algn, 
				double thres_weight,
				int*** diag_info) { 
  
  struct simple_diag_col* dcol = calloc(1, sizeof(struct simple_diag_col)); 
  if(dcol==NULL) error("find_diags_guided(): (1) Out of memory !"); 

  prepare_alignment(algn);

  

  struct seq_col *scol = algn->scol;
  struct seq *seq1 = &(scol->seqs[n1->seq_num[0]]);
  struct seq *seq2 = &(scol->seqs[n2->seq_num[0]]);
  int slen = scol->length;

  int ni1, ni2;
  struct seq *seq_n1, *seq_n2;

  unsigned int size = 16; 
  int length = 0; 
  struct diag **data = calloc(size, sizeof(struct diag* )); 
  //   printf("go for k\n"); 
  if(data==NULL) error("find_diags_guided(): (2) Out of memory !"); 
   
  struct algn_pos *tap;// = find_eqc(algn->algn, seq1->num, seq1->length-1); 
  long slen1 = algn->max_pos;
  long slen2 = slen1; 


  unsigned int max_dlen = pdist->max_dlen; 

  struct diag *sdg, *tsdg;
  struct diag *dg = create_diag(seq1->num, seq1,0, 
				seq2->num, seq2,0,0);   

  struct diag* tdg,*tdg2; 
 
  int  i,j,k,l,kposi,kposj,jmin=0,jmax=slen2; 

 
  int sn1 = seq1->num; 
  int sn2 = seq2->num; 

  int *c2n = smatrix->char2num; 
  int *sdata = smatrix ->data; 
  char *data1;// = seq1->data; 
  char *data2;// = seq2->data; 
  int smatrixlen = smatrix->length; 
  int a1,a2; 
 
  int maxslen = slen1; 
  if(slen2>maxslen) maxslen = slen2; 

  //long double **tmp_dist = create_tmp_pdist(pdist);
  //fill_tmp_pdist(pdist, tmp_dist, slen1, slen1);
  
 
  int max_pool = (slen1+slen2-1); 
  struct diag **diag_col = malloc(sizeof(struct diag*)*(slen1+1)); 
  struct diag **diag_row = malloc(sizeof(struct diag*)*(slen2+1)); 
  struct diag **pool_diags=malloc(sizeof(struct diag *)*max_pool); 
  int pooled = 0; 
  double thres_sim_score =para->PROT_SIM_SCORE_THRESHOLD; 
  char hasAli = (algn!=NULL); 
  double diag_thres = para->DIAG_THRESHOLD_WEIGHT; 
  double avg_sim_score = pdist->smatrix->avg_sim_score; 
   
  int maxd,maxd2,cons_maxd; 
  int score1=0,score2 = 0; 
  char prevail;
  int extended,passed;//,actuals; 
  int passed2;
  int max_under = 0;//para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS;
  int max_under_inner =0;//2*para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS;
  
  double dg_thres_weight = thres_weight;// * n1->seq_num_length*n2->seq_num_length;

  //int min_motives = n1->seq_num_length;
  //if(n2->seq_num_length > min_motives) min_motives = n2->seq_num_length;
  //min_motives++;
  int min_motives = sqrt(n1->seq_num_length * n2->seq_num_length);
  if(min_motives<2) min_motives = 2;
  //  if(min_motives>(n1->seq_num_length*n2->seq_num_length)) min_motives = (n1->seq_num_length * n2->seq_num_length);

  char stst;
  char mstatus;
  int slensq = slen*slen;
  int sqind;
  char startstop[slensq];
  char startstop_pre[slen];
  int c_klen[slensq];
  int c_kscore[slensq];

  int rev_pos[slen][slen1];
  struct gt_node *tgtn;


  // build reverse algn pos 
  
  int iterlen = 0;
  int iter[slen*slen];
  for(ni1=0;ni1<n1->seq_num_length;ni1++) {
    sn1 = n1->seq_num[ni1];
    for(ni2=0;ni2<n2->seq_num_length;ni2++) {
      //printf("   begin MARKER MARKER MARKER %i %i %i\n",i,j,k);
      sn2 = n2->seq_num[ni2];
      sqind = sn1*slen+sn2;
      iter[iterlen++] = sqind;
    }
  }
  

  dg->multi_dg = 1;
  dg->marked = 1;
  dg->multi_length = slensq;
  dg->multi_cont = malloc(sizeof(struct diag *)*slensq);
  memset(dg->multi_cont, 0, sizeof(struct diag *)*slensq);
  for(i=0;i<slensq;i++) {
    dg->multi_cont[i]=NULL;
  }

  maxd = n1->seq_num_length + n2->seq_num_length;
  for(i=0;i<maxd;i++) {
    k = i;
    tgtn = n1;
    if(i>=n1->seq_num_length) {
      k-=n1->seq_num_length;
      tgtn = n2;
    }
    ni1  = tgtn->seq_num[k];
    seq_n1 = &(scol->seqs[ni1]);
    for(j=0;j<slen1;j++) {
      rev_pos[ni1][j]=-1;
    }
    for(j=0;j<seq_n1->length;j++) {
      tap = find_eqc(algn->algn, seq_n1->num, j);     
      rev_pos[ni1][*(tap->eqcAlgnPos)] = j;
    }
  }

  // DIALIGN  
  for(k=0;k<=slen1;k++) { 
    //printf(" begin %i %i %i %i\n",k,slen, seq1->num, seq2->num);
    diag_col[k]=create_diag(seq1->num, seq1,0, 
			    seq2->num, seq2,0,0);  
    //printf(" middle\n");
    diag_col[k]->multi_dg = 1;
    diag_col[k]->marked = 1;
    diag_col[k]->multi_length = 0;
     
    //diag_col[k]->length = 1; 
    diag_col[k]->weight_sum = 0.0; 
    diag_col[k]->weight = 0.0; 
    diag_col[k]->meetsThreshold = 0; 
    pool_diags[pooled] = diag_col[k]; 
    pool_diags[pooled]->pool_pos = pooled; 
    pooled++; 
    //printf(" end %i %i\n",k,slen);
  } 
  for(k=0;k<=slen2;k++) { 
    diag_row[k] = diag_col[0]; 
  } 
 
  double old_weight; 
 
  if(max_dlen> slen1/2.0) max_dlen = slen1/2.0; 
  if(max_dlen> slen2/2.0) max_dlen = slen2/2.0; 
 
  for(i=0;i<=slen1;i++) { 
 
    // merge row/col 
    //printf("before \n");
    if(i>0) { 
      tdg = diag_col[i]; 
      while(tdg!=NULL) { 
	kposi = tdg->seq_p2.startpos+tdg->length; 
	//printf("  kposi %i %i %i\n", diag_row[kposi], kposi, slen2);
	if(tdg->weight_sum > diag_row[kposi]->weight_sum) { 
	  diag_row[kposi] = tdg; 
	  prevail = 1; 
	} else { 
	  prevail = 0; 
	} 
	tdg2 = tdg; 
	tdg = tdg->col_pred_diag; 
	if(! prevail) { 
	  pool_diags[tdg2->pool_pos]=NULL; 
	  //printf(" BEFORE\n");
	  free_diag(tdg2); 
	  //printf(" AFTER\n");
	} 
      } 
    } 
    //printf("after \n");

    for(j=1;j<=slen2;j++) { 
      if(diag_row[j-1]->weight_sum > diag_row[j]->weight_sum) { 
	diag_row[j] = diag_row[j-1]; 
      } 
    } 

    if(i==slen1) break; 

    for(j=jmin;j<jmax;j++) { 
      old_weight = 0.0;
      memset(startstop, 0, sizeof(char)*slen*slen);
      memset(startstop_pre, 0, sizeof(char)*slen);
      maxd = max_dlen;
      if( (i+maxd)>slen1) maxd = slen1-i;
      if( (j+maxd)>slen1) maxd = slen1-j;
      /*
      for(l=0;l<iterlen;l++) {
	sdg = dg->multi_cont[iter[l]];
	if(sdg!=NULL) {
	  free_diag(sdg);
	  dg->multi_cont[iter[l]]=NULL;
	}
      }
      */
      dg->seq_p1.startpos = i;
      dg->seq_p2.startpos = j;

      passed = 0;
      passed2 = 0;
      //      actuals = 0;
      for(k=1;k<=maxd;k++) { 
	//printf(" %i %i %i %i\n",i,j,k,slen1);
	kposi = i+k-1; 
	kposj = j+k-1; 
	maxd2 = maxd;
	maxd = 0;
	//printf(" k= %i\n",k);
	extended = 0;
	//mstatus = 0;

	for(ni1=0;ni1<n1->seq_num_length;ni1++) {
	  sn1 = n1->seq_num[ni1];
	  seq_n1 = &(scol->seqs[sn1]);

	  if( (rev_pos[sn1][kposi]>=0)&& (startstop_pre[sn1]!=2) )  {
	    data1 = seq_n1->data;
	    startstop_pre[sn1]=2;
	    for(ni2=0;ni2<n2->seq_num_length;ni2++) {
	      //printf("   begin MARKER MARKER MARKER %i %i %i\n",i,j,k);
	      sn2 = n2->seq_num[ni2];
	      seq_n2 = &(scol->seqs[sn2]);
	      sqind = sn1*slen+sn2;
	      stst = startstop[sqind];
	      //if(stst==1) mstatus += stst;

	      //printf("   end MARKER MARKER MARKER %i %i\n",sn2, kpos);
	      if( (rev_pos[sn2][kposj]>=0) && (stst!=2) 
		  && (diag_info[ni1][ni2][rev_pos[sn1][kposi]]==rev_pos[sn2][kposj])
		  ) {

		dg->length = k; 

		data2 = seq_n2->data;

		a1 = c2n[(int) data1[rev_pos[sn1][kposi]]]; 
		a2 = c2n[(int) data2[rev_pos[sn2][kposj]]]; 
		score2 = sdata[smatrixlen*a1+a2]; 
		sdg = NULL;
		if( (stst==0) && (score2 >= thres_sim_score)) {
		  //if(rev_pos[sn1][kpos]>1000) printf(" %i %i %i\n", sn1,kpos, rev_pos[sn1][kpos]);
		  sdg = dg->multi_cont[sqind];
		  if(sdg==NULL) {
		    sdg = create_diag(sn1, seq_n1, rev_pos[sn1][kposi],
				      sn2, seq_n2, rev_pos[sn2][kposj],0);
		  } else {
		    sdg->seq_p1.num = sn1;
		    sdg->seq_p1.sq = seq_n1;
		    sdg->seq_p1.startpos = rev_pos[sn1][kposi];
		    sdg->seq_p2.num = sn2;
		    sdg->seq_p2.sq = seq_n2;
		    sdg->seq_p2.startpos = rev_pos[sn2][kposj];
		    sdg->length = 0;
		  }
		  sdg->score = 0;
		  dg->multi_cont[sqind] = sdg;
		  startstop[sqind] = 1;
		  passed++;
		  c_klen[sqind] = 0;
		  c_kscore[sqind] = 0;
		} else if(stst==1) {
		  sdg = dg->multi_cont[sqind];
		} else {
		  //startstop[sqind]=2;
		  if(k>max_under) startstop[sqind] = 2;
		  sdg == NULL;
		}
		
		if(sdg!=NULL) {
		  maxd = maxd2; // reallow next k
		  //printf("  length compare %i %i %i %i\n",sqind, sdg->length+1,k,sdg);
		  c_klen[sqind]++;
		  c_kscore[sqind]+=score2;

		  // if( (sdg->seq_p1.startpos==20) && (sdg->seq_p2.startpos==2) ) printf(" %i %i %i %i %i\n", k, c_klen[sqind], sdg->seq_p1.startpos, sdg->seq_p2.startpos);

		  
		  if(((sdg->score+c_kscore[sqind]) < (avg_sim_score*(sdg->length+c_klen[sqind])))) { 
		    if( (sdg->length+c_klen[sqind])>max_under_inner) startstop[sqind]=2;
		  } 
		  /*
		  if( (c_klen[sqind]>=1) &&  
		      (c_kscore[sqind]< (para->PROT_DIAG_AVG_SCORE_THRESHOLD*c_klen[sqind]))) { 
		    passed--;
		    startstop[sqind]=2;
		  }    
		  */
		  
		  if( (c_klen[sqind]>=para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS) &&  
		      (c_kscore[sqind]< (para->PROT_DIAG_AVG_SCORE_THRESHOLD*c_klen[sqind]))) { 
		    if ( ((sdg->length+c_klen[sqind])>para->PROT_DIAG_MIN_LENGTH_THRESHOLD)) { 
		      passed--;
		      startstop[sqind]=2;
		    }
		  }    
		  
		  //printf(" diag length %i %i %i %i %i\n", sn1, sn2, i,j, sdg->length);
		  if( (score2 >= thres_sim_score) && (startstop[sqind]==1)) {
		    //printf("  very before score %i\n", sdg->score);
		    sdg->score += c_kscore[sqind];
		    sdg->length += c_klen[sqind];
		    //printf(" begin CONSIDER SDG %i %i %i %f %f\n",sn1,sn2,sdg->length,sdg->weight, diag_thres);
		    //sdg->weight = tmp_dist[sdg->length][sdg->score]; 
		    calc_weight(sdg, smatrix, pdist);
		    //printf(" before score %i %i  %.20f %i %i\n", sdg->score, sdg->length, sdg->weight, c_kscore[sqind], c_klen[sqind]);
		    /*
		    score1 = sdg->score;
		    if(score1!=sdg->score) {
		      print_diag(sdg);
		      printf("  score %i %i %i %.20f sc %i %i len %i %i %i \n", score1, sdg->score, sdg->length, sdg->weight, c_kscore[sqind],score2, c_klen[sqind],a1,a2);
		      error(" score changed!\n");
		    }
		    */
		    
		    //printf("  after score %i %i %.20f\n", sdg->score, sdg->length, sdg->weight);
		    //printf(" END\n");
		    //if(sdg->weight>0.0) printf(" end CONSIDER SDG %i %i %i %i %f %f\n",sn1, sn2, sdg->length, sdg->score, sdg->weight, diag_thres);
		    c_klen[sqind]=0;
		    c_kscore[sqind]=0;
		    sdg->meetsThreshold = (sdg->weight>diag_thres ? 1 : 0); 
		    extended++;
		  }
		}
	      } else {
		if(startstop[sqind]==1) {
		  startstop[sqind]=2;
		  passed--;
		}
	      }
	      if(startstop[sqind]!=2) startstop_pre[sn1]=0;
	    }
	  } else {
	    for(ni2=0;ni2<n2->seq_num_length;ni2++) {
	      sn2 = n2->seq_num[ni2];
	      sqind = sn1*slen+sn2;
	      if(startstop[sqind]==1) {
		startstop[sqind]=2;
		passed--;
	      }
	    }
	  }
	}

	// if no inner starting diag break
	if( (k==1) && !extended) maxd = 0;
	// if only one inner diag and proceeded max_under further then break 
	if( (passed2>=0) && (passed<min_motives) && ( (k-passed2)>=max_under)) {
	  //printf(" break %i\n");
	  maxd = 0;
	}
	// find the last k  where >=2 inner diags where active
	if(passed >= min_motives) passed2 = k;

	if( ((extended) && (passed>=min_motives))) { // && (passed>1)){// && (mstatus>1)) {
	  dg->weight = 0.0;
	  extended = 0;
	  for(l=0;l<iterlen;l++) {
	    sdg = dg->multi_cont[iter[l]];
	    if( (sdg!=NULL) && (startstop[iter[l]]>0) && (sdg->weight >= thres_weight)) {
	      extended++;
	      dg->weight += sdg->weight;
	    }
	  }
	  dg->meetsThreshold = (dg->weight>diag_thres ? 1 : 0); 
	    
	  if(dg->meetsThreshold && (extended >=min_motives) && (dg->weight>=old_weight) && (dg->weight>=dg_thres_weight)) {
	    old_weight = dg->weight;
	    
	    if(max_pool<=pooled) { 
	      max_pool += maxslen; 
	      pool_diags = realloc(pool_diags, sizeof(struct diag*)*max_pool); 
	      //printf(" pool size %i %.20f %.20f\n", max_pool, dg->weight, old_weight);
	    } 
	    tdg = malloc(sizeof(struct diag)); 
	    pool_diags[pooled] = tdg; 
	    dg->pool_pos = pooled; 
	    pooled++; 
	    *tdg = *dg; 
	    tdg->pred_diag = diag_row[j]; 
	    tdg->weight_sum =  diag_row[j]->weight_sum+tdg->weight; 
	    tdg->col_pred_diag = diag_col[i+k]; 
	    
	    diag_col[i+k] = tdg; 

	    dg->multi_cont = malloc(sizeof(struct diag *)*slensq);
	    memset(dg->multi_cont, 0, sizeof(struct diag *)*slensq);
	    for(l=0;l<iterlen;l++) {
	      
	      tsdg = tdg->multi_cont[iter[l]];
	      if((tsdg!=NULL)) { // && (startstop[iter[l]]==1)){
		sdg = malloc(sizeof(struct diag)); 

		if(startstop[iter[l]]==0) {
		  tsdg->length = 0;
		  tsdg->score = 0;
		  tsdg->weight = 0.0;
		  tsdg->total_weight = 0.0;
		}

		*sdg = *tsdg;
		dg->multi_cont[iter[l]]=sdg;
	      } else {
		dg->multi_cont[iter[l]]=NULL;
	      }
	      
	      //dg->multi_cont[iter[l]]=NULL;
	    }
	  }
	}
	//if(k > 20) printf("  maxk= %i\n",k);

      }
    }
  }
  //printf(" after main %i\n",k);
  tdg = diag_row[slen2]; 
  dcol->total_weight = 0; 
  double lencomp = (log(slen1)+log(slen2)); 
  length = 0;
  
  while((tdg!=NULL)) { 
    //if (tdg->weight <=0.0) break; 
    if(tdg->meetsThreshold) { 
      //      printf(" add tdg %i %i %i\n", tdg, tdg->length, tdg->meetsThreshold);
      dcol->total_weight += tdg->weight+lencomp; 
 
      data[length] = tdg; 
      l = 0;
      tdg->weight = 0.0;
      for(i=0;i<slensq;i++) {
	sdg = tdg->multi_cont[i];
	if(sdg!=NULL) {
	  calc_weight(sdg, smatrix, pdist);
	  
	  if ((sdg->length>0) && (sdg->meetsThreshold) && (sdg->weight > thres_weight)){
	    tdg->multi_cont[l++] = sdg;
	    tdg->weight += sdg->weight;
	    //printf ("  tdg %i %f %f\n", tdg, tdg->weight, sdg->weight);
	  }
	}
      }
      tdg->multi_length = l;
      if( (tdg->multi_length>=min_motives) && (tdg->weight >= dg_thres_weight)) {
	tdg->marked = 1;
	tdg->total_weight = tdg->weight;
	tdg->weight_sum = -1.0; 
	length++; 
      } else {
	// TODO free ??
      }

      if(length >= size) { 
	size += 64; 
	data = realloc(data, sizeof(struct diag *)*size); 
	if(data==NULL) error("find_diags(): (3) Out of memory !"); 
      } 
    } 
    tdg = tdg->pred_diag; 
  } 
 
  for(k=0;k<pooled;k++) { 
    if(pool_diags[k]!=NULL) 
      if(pool_diags[k]->weight_sum>-1.0) { 
	free_diag(pool_diags[k]); 
      } 
  } 
   
  free(pool_diags); 
  free(diag_col); 
  free(diag_row); 
  free_diag(dg); 
  //free_tmp_pdist(tmp_dist, pdist->max_dlen);

  dcol->length = length; 
 
  data = realloc(data, sizeof(struct diag *)*length); 
  dcol->data = data; 

  return dcol;
}

/** 
 * finds all relevant diags by the DIALIGN METHOD  
 * 
 * The pointer returned (and the ones included in the struct)  
 * has to be deallocted explicitely from memory. 
 */ 
struct simple_diag_col* find_diags_dialign(struct scr_matrix *smatrix,  
				struct prob_dist *pdist, struct seq* seq1,  
				struct seq* seq2, struct alignment *algn, 
				 long double **tmp_dist, int round) { 
  struct simple_diag_col* dcol = calloc(1, sizeof(struct simple_diag_col)); 
  if(dcol==NULL) error("find_diags_dialign(): (1) Out of memory !"); 
   
  unsigned int size = 16; 
  int length = 0; 
  struct diag **data = calloc(size, sizeof(struct diag* )); 
  //   printf("go for k\n"); 
  if(data==NULL) error("find_diags_dialign(): (2) Out of memory !"); 
   
  long slen1 = seq1->length; 
  long slen2 = seq2->length; 
  unsigned int max_dlen = pdist->max_dlen; 
 
  struct diag *dg = create_diag(seq1->num, seq1,0, 
				seq2->num, seq2,0,0);   
  struct diag* tdg,*tdg2; 
 
  int  i,j,k,kpos,jmin=0,jmax=slen2; 
 
  int sn1 = seq1->num; 
  int sn2 = seq2->num; 

  int *c2n = smatrix->char2num; 
  int *sdata = smatrix ->data; 
  char *data1 = seq1->data; 
  char *data2 = seq2->data; 
  int smatrixlen = smatrix->length; 
  int a1,a2; 
  int c_klen,c_kscore; 
 
  int maxslen = slen1; 
  if(slen2>maxslen) maxslen = slen2; 
 
  int max_pool = (slen1+slen2-1); 
  struct diag **diag_col = malloc(sizeof(struct diag*)*(slen1+1)); 
  struct diag **diag_row = malloc(sizeof(struct diag*)*(slen2+1)); 
  struct diag **pool_diags=malloc(sizeof(struct diag *)*max_pool); 
  int pooled = 0; 
  double thres_sim_score =para->PROT_SIM_SCORE_THRESHOLD; 
  char hasAli = (round>1); //(algn!=NULL); 
  struct algn_pos **ap,*tap; 
  double diag_thres = para->DIAG_THRESHOLD_WEIGHT; 
  double avg_sim_score = pdist->smatrix->avg_sim_score; 
   
  int maxd,maxd2,cons_maxd; 
  double score1=0,score2 = 0; 
  char prevail; 
 
  if(hasAli || para->DO_ANCHOR) { 
    ap = algn->algn; 
  } 

  // DIALIGN  
  for(k=0;k<=slen1;k++) { 
    diag_col[k]=create_diag(seq1->num, seq1,0, 
			    seq2->num, seq2,0,1);  
     
    //diag_col[k]->length = 1; 
    diag_col[k]->weight_sum = 0.0; 
    diag_col[k]->weight = 0.0; 
    pool_diags[pooled] = diag_col[k]; 
    pool_diags[pooled]->pool_pos = pooled; 
    pooled++; 
  } 
  for(k=0;k<=slen2;k++) { 
    diag_row[k] = diag_col[0]; 
  } 
 
  double old_weight; 
 
  if(max_dlen> slen1/2.0) max_dlen = slen1/2.0; 
  if(max_dlen> slen2/2.0) max_dlen = slen2/2.0; 
 
  for(i=0;i<=slen1;i++) { 
 
    // merge row/col 
    if(i>0) { 
      tdg = diag_col[i]; 
      while(tdg!=NULL) { 
	kpos = tdg->seq_p2.startpos+tdg->length; 
	if(tdg->weight_sum > diag_row[kpos]->weight_sum) { 
	  diag_row[kpos] = tdg; 
	  prevail = 1; 
	} else { 
	  prevail = 0; 
	} 
	tdg2 = tdg; 
	tdg = tdg->col_pred_diag; 
	if(! prevail) { 
	  pool_diags[tdg2->pool_pos]=NULL; 
	  free_diag(tdg2); 
	} 
      } 
    } 
    for(j=1;j<=slen2;j++) { 
      if(diag_row[j-1]->weight_sum > diag_row[j]->weight_sum) { 
	diag_row[j] = diag_row[j-1]; 
      } 
    } 
    if(i==slen1) break; 
    if(hasAli || para->DO_ANCHOR) { 
      tap = find_eqc(ap, sn1, i); 
       
      if(tap->predF!=NULL) { 
	jmin = tap->predF[sn2]+1; 
      } else { 
	jmin = 0; 
      } 
       
      if(tap->succF!=NULL) { 
	jmax = tap->succF[sn2]; 
      }else { 
	jmax = slen2; 
      } 
       
      if(jmin<0) jmin = 0; 
    }  

    if(para->DO_ANCHOR) {
      if( (jmin-1)==jmax) {
	jmin--;
	jmax++;
      }
    }

    for(j=jmin;j<jmax;j++) { 
       
      if(i<slen1 && j<slen2) { 
	a1 = c2n[(int) data1[i]]; 
	a2 = c2n[(int) data2[j]]; 
	score1 = sdata[smatrixlen*a1+a2]; 
      } else { 
	score1 = 0; 
      } 
       
      if(score1>=thres_sim_score) { 
	maxd = slen1 - i; 
	maxd2 = slen2 - j; 
	if(maxd >maxd2) maxd = maxd2; 
	if(maxd > max_dlen) maxd = max_dlen; 
	 
	dg->seq_p1.startpos = i; 
	dg->seq_p2.startpos = j; 
	dg->score = score1; 
 
	cons_maxd = maxd+1; 
	old_weight = 0.0; 
	 
	c_klen = 0; 
	c_kscore = 0; 
	 
	for(k=1;k<=maxd;k++) { 
	  dg->length = k; 
	  kpos = i+k; 
	  if(hasAli || para->DO_ANCHOR) { 
	    a1 = i+k-1; 
	    a2 = j+k-1; 
	    tap = find_eqc(ap, sn1, a1); 

	    if (! ( (para->DO_ANCHOR) && (tap->predF!=NULL) && (tap->succF!=NULL) && 
		    (tap->predF[sn2]==tap->succF[sn2]) &&
		    (tap->predF[sn2]==a2)
		  )
		) {
	      
	      if(tap->predF!=NULL) { 
		if( (tap->predF[sn2] - a2)>0) break; 
		/*
		  if(tap->predF[sn2]==a2) {
		  if( (tap->succF==NULL) || (tap->predF[sn2]!=tap->succF[sn2])) break;
		  }
		*/
	      }  
	      if(tap->succF!=NULL) { 
		if((a2 - tap->succF[sn2])>0) break; 
		/*
		  if(tap->succF[sn2]==a2) {
		  if( (tap->predF==NULL) || (tap->predF[sn2]!=tap->succF[sn2])) break;
		  }
		*/
	      } 
	    }
	  }  
	   
	  if(k>1) { 
	    a1 = c2n[(int) data1[kpos-1]]; 
	    a2 = c2n[(int) data2[j+k-1]]; 
	    score2 = sdata[smatrixlen*a1+a2]; 
	    dg->score += score2; 
	  } else { 
	    score2 = score1; 
	  } 
 
	   
	  if( dg->score < (avg_sim_score*(double)k)) {
 	    break; 
	  } 

	  c_klen++; 
	  c_kscore+=score2; 
	  
	  if( (c_klen>=para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS) &&  
	      (c_kscore< (para->PROT_DIAG_AVG_SCORE_THRESHOLD*c_klen))) { 
	    if ( (k>para->PROT_DIAG_MIN_LENGTH_THRESHOLD)) { 
	      break; 
	    } else { 
	      if(maxd>para->PROT_DIAG_MIN_LENGTH_THRESHOLD) maxd = para->PROT_DIAG_MIN_LENGTH_THRESHOLD; 
	    } 
	  }   
 
	  if(score2 >= thres_sim_score) { 
	    c_klen=0; 
	    c_kscore=0; 
	    if(1) { 
	      if(!hasAli) { 
		dg->weight = tmp_dist[k][dg->score]; 
		dg->meetsThreshold = (dg->weight>diag_thres ? 1 : 0); 
	      } else { 
		real_calc_weight(dg, smatrix, pdist,1,tmp_dist,algn);	     
	      }	       
	      if(dg->meetsThreshold && (dg->weight>=old_weight)) { 
		old_weight = dg->weight; 
		if(max_pool<=pooled) { 
		  max_pool += maxslen; 
		  pool_diags = realloc(pool_diags, sizeof(struct diag*)*max_pool); 
		  //printf(" old pool size %i\n", max_pool);
		} 
		tdg = malloc(sizeof(struct diag)); 
		pool_diags[pooled] = tdg; 
		dg->pool_pos = pooled; 
		pooled++; 
		*tdg = *dg; 
		tdg->pred_diag = diag_row[j]; 
		tdg->weight_sum =  diag_row[j]->weight_sum+tdg->weight; 
		tdg->col_pred_diag = diag_col[kpos]; 
 
		diag_col[kpos] = tdg; 
	      } 
	    } 
	  } 
	} 
      } 
       
    } 
  } 
   
  tdg = diag_row[slen2]; 
  dcol->total_weight = 0; 
  double lencomp = (log(slen1)+log(slen2)); 
 
  while((tdg!=NULL)) { 
    if (tdg->weight <=0.0) break; 
    if(1) { 
      dcol->total_weight += tdg->weight+lencomp; 
 
      data[length] = tdg; 
      tdg->weight_sum = -1.0; 
      length++; 
      if(length >= size) { 
	size += 64; 
	data = realloc(data, sizeof(struct diag *)*size); 
	if(data==NULL) error("find_diags(): (3) Out of memory !"); 
      } 
    } 
    tdg = tdg->pred_diag; 
  } 
 
  for(k=0;k<pooled;k++) { 
    if(pool_diags[k]!=NULL) 
      if(pool_diags[k]->weight_sum>-1.0) { 
	free(pool_diags[k]); 
      } 
  } 
   
  free(pool_diags); 
  free(diag_col); 
  free(diag_row); 
  free_diag(dg); 
  dcol->length = length; 
 
  data = realloc(data, sizeof(struct diag *)*length); 
  dcol->data = data; 
 
  if(para->DEBUG>5) { 
    for(i=0;i<length;i++) { 
      print_diag(data[i]); 
      printf("\n"); 
    } 
  } 
 
  return dcol; 
} 

 
/** 
 * finds all relevant diags by dynamic programming on diagonal stripes 
 * 
 * The pointer returned (and the ones included in the struct)  
 * has to be deallocted explicitely from memory. 
 */ 
inline struct simple_diag_col* find_diags_dyn(struct scr_matrix *smatrix,  
				struct prob_dist *pdist, struct seq* seq1,  
				struct seq* seq2, struct alignment *algn, 
				 long double **tmp_dist) { 
 
  struct simple_diag_col* dcol = calloc(1, sizeof(struct simple_diag_col)); 
  if(dcol==NULL) error("find_diags_dyn(): (1) Out of memory !"); 
   
  unsigned int size = 64; 
  int  l, k,lastk, maxl; 
  int length = 0; 
  struct diag **data = calloc(size, sizeof(struct diag *)); 
  if(data==NULL) error("find_diags_dyn(): (2) Out of memory !"); 
   
  int slen1 = seq1->length; 
  int slen2 = seq2->length; 
  unsigned int max_dlen = pdist->max_dlen; 
 
  struct diag* dg = create_diag(seq1->num, seq1,0, 
				seq2->num, seq2,0,0);   
  struct diag* tdg; 
 
  int i,j,d; 
 
  int sn1 = seq1->num; 
  int sn2 = seq2->num; 
  //    printf("%i\n",slen2);  
  int *c2n = smatrix->char2num; 
  int *sdata = smatrix ->data; 
  char *data1 = seq1->data; 
  char *data2 = seq2->data; 
  int slen = smatrix->length; 
 
  int a1,a2; 
  int score1=0,score2 = 0; 
 
  int maxslen = slen1; 
  if(slen2>maxslen) maxslen = slen2; 
  double avslen = ((double)slen1+slen2)/2.0; 
 
  int delta; 
  int sim_thr_pred_pos[maxslen]; 
  //int sim_thr_succ_pos[maxslen]; 
  int scores[maxslen]; 
  long score_sum[maxslen]; 
  long s_sum; 
  int old_thr_pos; 
 
  //double *dyn_weight = calloc(maxslen, sizeof(double)); 
  double weight; 
  struct diag **dyn_diags=malloc(sizeof(struct diag *)*maxslen); 
  int max_pool = maxslen; 
  struct diag **pool_diags=malloc(sizeof(struct diag *)*max_pool); 
  int pooled = 0; 
 
  int thres_sim_score = para->PROT_SIM_SCORE_THRESHOLD; 
 
  int kscore, klen; 
 
  char hasAli = (algn!=NULL); 
 
  double diag_thres = para->DIAG_THRESHOLD_WEIGHT; 
  double avg_sim_score = pdist->smatrix->avg_sim_score; 
 
  struct algn_pos **ap,*tap; 
  if(hasAli) { 
    ap = algn->algn; 
    thres_sim_score =thres_sim_score/2; 
    if(thres_sim_score<4) thres_sim_score = 4; 
    diag_thres = 0.0; 
  } 
 
 
  for(d=-slen1+1;d<slen2;d++) { 
    //printf(" d=%i\n",d);  
     
    if(d<=0) { 
      i = - d; 
      j=0; 
      maxl = slen1-i; 
      if(slen2<maxl) maxl = slen2; 
    } else { 
      i=0; 
      j=d; 
      maxl = slen2 -j; 
      if(slen1<maxl) maxl = slen1; 
    } 
 
    // prepare values 
    old_thr_pos=-1; 
    s_sum=0; 
    for(k=0;k<maxl;k++) { 
      // hier: hasAlipredF/SuccF abfragen !!!! 
      if(hasAli) { 
	a1 = i+k; 
	a2 = j+k; 
	//printf("pre %i\n", a1); 
	tap = find_eqc(ap, sn1, a1); 
	//printf("after %i\n", a1); 
	if(tap->predF!=NULL) { 
	  if( (tap->predF[sn2] - a2)>0) break; 
	}  
	if(tap->succF!=NULL) { 
	  if((a2 - tap->succF[sn2])>0) break; 
	} 
      }  
 
      a1 = c2n[(int) data1[i+k]]; 
      a2 = c2n[(int) data2[j+k]]; 
      score1 = sdata[slen*a1+a2]; 
      scores[k] = score1; 
      s_sum+= score1; 
      score_sum[k] = s_sum; 
      if(score1 < thres_sim_score) { 
	sim_thr_pred_pos[k] = old_thr_pos; 
      } else { 
	//if(old_thr_pos>=0) sim_thr_succ_pos[old_thr_pos] = k; 
	old_thr_pos = k; 
      } 
    } 
    maxl = k; 
    //if(old_thr_pos>=0) sim_thr_succ_pos[old_thr_pos] = maxl; 
 
    dyn_diags[0] = create_diag(seq1->num, seq1,0, 
				seq2->num, seq2,0,0);  
    dyn_diags[0]->weight_sum = 0.0; 
    dyn_diags[0]->weight = 0.0; 
    pool_diags[0] = dyn_diags[0]; 
    pooled = 1; 
    lastk=0; 
 
    for(k=0;k<maxl;k++) { 
      //printf("process %i %i\n", k,maxl); 
 
      if(k>0) { 
	dyn_diags[k] = dyn_diags[lastk]; 
	//dyn_weight[k] = dyn_weight[lastk]; 
      } 
      if(hasAli) { 
	a1 = i+k; 
	a2 = j+k; 
	//printf("pre %i\n", a1); 
	tap = find_eqc(ap, sn1, a1); 
	//printf("after %i\n", a1); 
	if(tap->predF!=NULL) { 
	  if( (tap->predF[sn2] - a2)>0) break; 
	}  
	if(tap->succF!=NULL) { 
	  if((a2 - tap->succF[sn2])>0) break; 
	} 
      }  
 
      score1 = scores[k]; 
      if(score1>=thres_sim_score) { 
 
	for(l=para->DIAG_MIN_LENGTH;l<=max_dlen; l++) { 
	  delta = k-l+1; 
	  dg->seq_p1.startpos = i+delta; 
	  dg->seq_p2.startpos = j+delta; 
 
	  kscore = 0; 
	  klen = 0; 
 
	  if((dg->seq_p1.startpos<0) || (dg->seq_p2.startpos<0)) { 
	    break; 
	  } else { 
 
 
	    dg->length = l; 
	    //printf("%i %i \n", i-l+1, j-l+1); 
 
	    score2 = scores[delta]; 
	    klen++; 
	    kscore += score2; 
 
	    if( (kscore < avg_sim_score*klen)) { // && (dg->score>1.2*k*avg_sim_score)) { 
	       
	      /* experiment */ 
	      if( ( (k>=0.2*avslen) || (k>=20) || ( (i-j)>0.3*avslen)) && klen>=4) break; 
	    } 
	   
	    if(kscore >= avg_sim_score*klen ) { 
	      //if(score2 >= thres_sim_score) { 
	      kscore = 0; 
	      klen = 0; 
	      dg->score = score_sum[k] - (delta > 0 ? score_sum[delta-1] : 0); 
	      if(dg->score <= avg_sim_score*dg->length) break; 
	      if(!hasAli) { 
		dg->weight = tmp_dist[dg->length][dg->score]; 
		dg->meetsThreshold = (dg->weight>diag_thres ? 1 : 0); 
	      } else { 
		real_calc_weight(dg, smatrix, pdist,1,tmp_dist,algn);	     
	      } 
	       
	      if(dg->meetsThreshold) { 
		if(max_pool<=pooled) { 
		  max_pool += maxslen; 
		  pool_diags = realloc(pool_diags, sizeof(struct diag*)*max_pool); 
		} 
		if(delta==0) { 
		  if(dg->weight > dyn_diags[k]->weight_sum) { 
		    dg->weight_sum = dg->weight; 
		    dyn_diags[k] = malloc(sizeof(struct diag));  
		    pool_diags[pooled] = dyn_diags[k]; 
		    pooled++; 
		    *dyn_diags[k] = *dg; 
		    dyn_diags[k]->pred_diag = NULL; 
		  } 
		} else { 
		  weight = dg->weight + dyn_diags[delta-1]->weight_sum; 
		  if( (weight) >= dyn_diags[k]->weight_sum) { 
		    dg->weight_sum = weight; 
		    dyn_diags[k] = malloc(sizeof(struct diag));  
		    pool_diags[pooled] = dyn_diags[k]; 
		    pooled++; 
		    *dyn_diags[k] = *dg; 
		    if(dyn_diags[delta-1]->weight >0) { 
		      dyn_diags[k]->pred_diag = dyn_diags[delta-1]; 
		    } else { 
		      dyn_diags[k]->pred_diag = NULL; 
		    } 
		  } 
		} 
 
		lastk = k; 
	      } 
	       
	    } else { 
	      l += (delta - sim_thr_pred_pos[delta])-1; 
	    } 
	  } 
	} 
      } 
    } 
    tdg = dyn_diags[lastk]; 
    while((tdg!=NULL)) { 
      if (tdg->weight <=0.0) break; 
       
      data[length] = tdg; 
      tdg->weight_sum = -1.0; 
      length++; 
      if(length >= size) { 
	size += 64; 
	data = realloc(data, sizeof(struct diag *)*size); 
	if(data==NULL) error("find_diags(): (3) Out of memory !"); 
      } 
      tdg = tdg->pred_diag; 
    } 


    for(k=0;k<pooled;k++) { 
      if(pool_diags[k]->weight_sum>-1.0) free(pool_diags[k]); 
    } 
  } 
 
  data = realloc(data, sizeof(struct diag *)*length); 
  free(pool_diags); 
  free(dyn_diags); 
  free_diag(dg); 
  dcol->length = length; 
  dcol->data = data; 
 
  if(para->DEBUG>5) { 
    for(i=0;i<length;i++) { 
      print_diag(data[i]); 
      printf("\n"); 
    } 
  } 
 
  return dcol; 
} 
 
/**
 * changes all_diags such that each diag is splitted, whenever there is 
 * another diag that has parts of it in common
void split_diags(struct seq_col *in_seq_col, struct diag_col *all_diags) { 

  // note that we re-use *pred_diag for the successor diag here to save memory !!!


  struct diag *dg;
  struct diag *tdg; 
  struct diag *pdg; 
  struct diag *sdg;
  struct diag **data;
  struct diag **n_data;
  struct simple_diag_col *sdcol;

  unsigned int s1, s2, sl = in_seq_col->length;
  unsigned int dl, d, si,k;
  unsigned int odgpos1[sl];
  unsigned int odgpos2[sl];
  int min,max,xpos,txpos1, txpos2;
  unsigned int odglen;
  char next;
  char hasSplit=1;
  while(hasSplit) {
    //printf("hasSplit %i\n",hasSplit);
    hasSplit = 0;
    for(s1=0;s1<sl;s1++) { 
      for(s2=s1+1;s2<sl;s2++) { 
	sdcol = all_diags->diag_matrix[s1+sl*s2];
	dl = sdcol->length;
	data = sdcol->data;
	memset(odgpos1,0,sl*sizeof(int));
	memset(odgpos2,0,sl*sizeof(int));
	pdg=NULL;
	for(d=0;d<dl;d++) {
	  dg = data[d];
	  dg->weight_sum = dg->weight;
	  dg->orig_score = dg->score;
	  if(pdg!=NULL) pdg->pred_diag = dg;
	  pdg = dg;
	  dg->pred_diag=NULL;
	  
	  next=0;
	  while(! next) {
	    next = 1;
	    xpos = -1;
	    
	    // round 1
	    //print_diag(dg);
	    //printf(" %i %i min=%i max=%i\n", s1,s2, min,max);
	    min = dg->seq_p1.startpos;
	    max = dg->seq_p1.startpos+dg->length-1;
	    for(si=0;si<sl;si++) {
	      if( (si==s1) || (si==s2)) continue;
	      
	      // if( (s1==0) && (s2==3) && (si==2)) printf("%i %i  min=%i max=%i\n", s1,s2,si,min,max);
	      if(s1<si){
		odglen =all_diags->diag_matrix[s1+sl*si]->length;
		n_data = all_diags->diag_matrix[s1+sl*si]->data;
	      } else {
		odglen = all_diags->diag_matrix[si+sl*s1]->length;
		n_data = all_diags->diag_matrix[si+sl*s1]->data;
	      }
	      //printf(" %i %i %i\n", si, odglen, odgpos1[si]);
	      while( (odgpos1[si]<odglen)) {
		sdg = n_data[odgpos1[si]];
		txpos1 = (s1==sdg->seq_p2.num) ? sdg->seq_p2.startpos : sdg->seq_p1.startpos;
		txpos2 = txpos1 + sdg->length-1;
		//if( (s1==0) && (s2==2) && (si==3)) printf(" pppp %i %i %i %i\n", min,max,txpos1, txpos2);
		if(txpos2<min) {
		  odgpos1[si]++;
		} else {
		  txpos2++;
		  if( ( (txpos1>min) && (txpos1<=max) && ((txpos1<xpos) || (xpos<0)))) xpos=txpos1;
		  if( (txpos2<=max) && ((txpos2<xpos) || (xpos<0))) xpos=txpos2;
		  //printf("sss %i %i %i %i %i xp=%i\n", odgpos1[si], min,max,txpos1,txpos2,xpos);
		  break;
		}
	      }
	    }
	    
	    // round 2
	    min = dg->seq_p2.startpos;
	    max = dg->seq_p2.startpos+dg->length-1;
	    if(xpos>=0) xpos = dg->seq_p2.startpos + (xpos-dg->seq_p1.startpos);
	    for(si=0;si<sl;si++) {
	      if( (si==s1) || (si==s2)) continue;
	      //if(si==s2) continue;
	      
	      if(s2<si){
		odglen = all_diags->diag_matrix[s2+sl*si]->length;
		n_data = all_diags->diag_matrix[s2+sl*si]->data;
	      } else {
		odglen = all_diags->diag_matrix[si+sl*s2]->length;
		n_data = all_diags->diag_matrix[si+sl*s2]->data;
	      }
	      while( (odgpos2[si]<odglen)) {
		//printf(" before split %i %i %i %i %i\n", si, s1, s2, odgpos2[si],odglen);
		sdg = n_data[odgpos2[si]];
		//printf(" after split %i\n", sdg);
		txpos1 = (s2==sdg->seq_p2.num) ? sdg->seq_p2.startpos : sdg->seq_p1.startpos;
		txpos2 = txpos1 + sdg->length-1;
		if(txpos2<min) {
		  odgpos2[si]++;
		} else {
		  txpos2++;
		  if( (txpos1>min) && (txpos1<=max) && ((txpos1<xpos) || (xpos<0))) xpos=txpos1;
		  if( (txpos2<=max) && ((txpos2<xpos) || (xpos<0))) xpos=txpos2;
		  break;
		}
	      }
	    }
	    
	    // split diag ?
	    if( (xpos>=0) && (xpos<=max) && (xpos>min)) {
	      hasSplit = 1;
	      next = 0;
	      //print_diag(dg);
	      tdg = create_diag(dg->seq_p1.num, dg->seq_p1.sq, xpos-min+dg->seq_p1.startpos,
				dg->seq_p2.num, dg->seq_p2.sq, xpos,
				max-xpos+1);
	      tdg->weight_sum = dg->weight_sum;
	      tdg->orig_score = dg->score;
	      dg->length = xpos - min;
	      //print_diag(dg);
	      //print_diag(tdg);
	      //if( (s1==0) && (s2==2)) printf(" SPLIT %i %i %i %i %i !\n", dg, dg->seq_p1.startpos, dg->length, tdg->seq_p1.startpos, tdg->length);
	      //printf(" 2SPLIT %i !\n", dg->length);
	      sdcol->length++;
	      dg->pred_diag = tdg;
	      tdg->pred_diag=NULL;
	      dg = tdg;
	      pdg = dg;
	    } else {
	      //if(s1==0) printf(" ELSE %i %i %i  !\n", dg, dg->seq_p1.startpos, dg->length);
	      
	    }
	  }
	}
	dl = sdcol->length;
	if(dl>0) dg=sdcol->data[0];
	//printf(" before realloc %i\n", sdcol->length);
	sdcol->data=realloc(sdcol->data, sdcol->length*sizeof(struct diag*));
	//printf(" after realloc %i\n", sdcol->length);
	
	if( (dl>0)&& (sdcol->data==NULL)) {
	  error(" split_diags: Out of Memory - reallocating!\n");
	}
	for(k=0;k<dl;k++) {
	  //printf (" %i %i %i %i\n",s1,s2,k,dg);
	  if(dg==NULL) {
	    error("split_diags: ALARM dg=NULL !\n");
	  }
	  //printf(" diag %i %i %i %i\n",dg->seq_p1.num, dg->seq_p2.num, dg->seq_p1.startpos, dg->length);
	  sdcol->data[k] = dg;
	  dg = dg->pred_diag;
	}
      }
    }
  }
}

 */


/**
 * builds the upgma guide tree in the given diag_col
 */
void build_guide_tree(struct diag_col *dcol) {
  int slen = dcol->seq_amount;
  double weights[slen][slen];
  struct gt_node *nodes[slen];
  struct gt_node *gtn, *gtn1, *gtn2;
  int max1=0, max2=1;

  int i,j,k;

  for(i=0;i<slen;i++) {
    gtn = malloc(sizeof(struct gt_node ));
    gtn->isLeaf = 1;
    gtn->seq_num = malloc(sizeof(int)*1);
    gtn->seq_num[0] = i;
    gtn->seq_num_length = 1;
    gtn->succ1 = NULL;
    gtn->succ2 = NULL;
    nodes[i] = gtn;
    for(j=i+1;j<slen;j++) {
      weights[i][j] = dcol->diag_matrix[slen*i+j]->total_weight;
      weights[j][i] = weights[i][j];
      if(weights[i][j] > weights[max1][max2]) {
	max1 = i;
	max2 = j;
      }
    }
  }
  
  for(k=0;k<(slen-1);k++) {
    gtn1 = nodes[max1];
    gtn2 = nodes[max2];
    gtn = malloc(sizeof(struct gt_node ));
    gtn->isLeaf = 0;
    gtn->seq_num_length = gtn1->seq_num_length + gtn2->seq_num_length;
    gtn->seq_num = malloc(sizeof(int)*gtn->seq_num_length);

    for(i=0;i<gtn1->seq_num_length;i++) {
      gtn->seq_num[i] = gtn1->seq_num[i];
    }
    for(i=0;i<gtn2->seq_num_length;i++) {
      gtn->seq_num[gtn1->seq_num_length+i] = gtn2->seq_num[i];
    }

    gtn->succ1 = gtn1;
    gtn->succ2 = gtn2;
    nodes[max1] = gtn;
    nodes[max2] = NULL;
    for(i=0;i<slen;i++) {
      if( (i!=max1) && (i!=max2)) {
	weights[i][max1] = 0.1*0.5*(weights[i][max1]+weights[i][max2]) + 0.9* 
	  ( (weights[i][max1] > weights[i][max2]) ? weights[i][max1] : weights[i][max2]);
	//weights[i][max1] = 1.0*(weights[i][max1]+weights[i][max2]) + 0.0*                          ( (weights[i][max1] > weights[i][max2]) ? weights[i][max1] : weights[i][max2]);
	weights[max1][i] = weights[i][max1];
      }
    }
    max1 = -1;
    max2 = -1;
    for(i=0;i<slen;i++) {
      for(j=i+1;j<slen;j++) {
	if( (nodes[i]!=NULL) && (nodes[j]!=NULL)) {
	  if( (max1<0) || (weights[max1][max2]< weights[i][j])) {
	    max1=i; 
	    max2=j;
	  }
	}
      }
    }
  }
  
  dcol->gt_root = nodes[0];
  /*
  printf(" root %i\n", nodes[0]);
  printf(" left1 %i\n", nodes[0]->succ1->succ1->seq_num);
  printf(" left2 %i\n", nodes[0]->succ1->succ2->seq_num);
  printf(" right %i\n", nodes[0]->succ2->seq_num);
  */
}


/** 
 * Finds all diags of each pair of sequences in in_seq_col by using 
 * the function above 
 * 
 * The pointer returned (and the ones included in the struct)  
 * has to be deallocted explicitely from memory. 
 */ 
struct diag_col *find_all_diags(struct scr_matrix *smatrix,  
				struct prob_dist *pdist, 
				struct seq_col *in_seq_col, struct alignment *algn, int round) { 
  unsigned int s1, s2, rs2, sl = in_seq_col->length, sp, ap; 
  struct diag_col *all_diags = create_diag_col(sl); 
  struct simple_diag_col *sdcol; 
 
  unsigned int diag_amount = 0; 
  struct diag *dg; 
 
  char hasAli = (round >1);//(algn!=NULL); 
 
  long double **tmp_dist = NULL; 
  if(!hasAli) tmp_dist = create_tmp_pdist(pdist); 
 
  int s2max = sl; 
  int s2width =(int) sqrt(sl); 
 
  double total=0.0; 
  //double imp[sl]; 
  //  for(s1=0;s1<sl;s1++) { 
  //  imp[s1] = 0.0; 
  //} 
  double totala[sl]; 
  memset(totala,0,sizeof(double)*sl); 
  for(s1=0;s1<sl;s1++) { 
    if(para->FAST_PAIRWISE_ALIGNMENT && s2width+1<sl) { 
      s2max = s1+s2width+1; 
      //printf("%i %i\n",s1, s2max); 
    } 
    //printf("before enter %i\n", s1); 
    for(s2=s1+1;s2<s2max;s2++) { 
      rs2 = s2 % sl; 
      if(!hasAli) { 
	fill_tmp_pdist(pdist,tmp_dist,in_seq_col->seqs[s1].length,in_seq_col->seqs[rs2].length ); 
      } 
      if(para->DEBUG>5) printf("%i %i\n", s1,s2); 
      //time1 = clock(); 
      //sdcol=find_diags_dyn(smatrix, pdist, &in_seq_col->seqs[s1], 
      //		 &in_seq_col->seqs[s2],algn,tmp_dist); 
       
      /* 
      doAlign = 1; 
      if(hasAli) { 
	if(algn->redo_seqs[s1*sl+s2]==0)  
	  doAlign = 0; 
      } 
       
      if(doAlign) { 
      */ 
      if(in_seq_col->seqs[s1].length > 0 && in_seq_col->seqs[s2].length > 0) { 
	//if(para->DEBUG>1) printf(" %i %i %i\n", s1, rs2,sl-1); 
	//	printf("find diags %i %i\n",s1,s2); 
        sdcol=find_diags_dialign(smatrix, pdist, &in_seq_col->seqs[s1], 
				 &in_seq_col->seqs[rs2],algn,tmp_dist, round); 
	//	imp[s1] += sdcol->total_weight; 
	//imp[s2] += sdcol->total_weight; 
	total += sdcol->total_weight; 
	//totala[s1] +=sdcol->total_weight;  
	//totala[s2] +=sdcol->total_weight;  
	//printf(" num of diags:%i\n ", sdcol->length); 
	/* 
	  } else { 
	  sdcol = calloc(1, sizeof(struct simple_diag_col)); 
	  sdcol->length = 0; 
	  } 
	*/ 
	//printf("%i %i %f\n", s1,s2, (clock()-time1)/CLOCKS_PER_SEC); 
	 
	all_diags->diag_matrix[s1+sl*rs2] = sdcol; 
	all_diags->diag_matrix[rs2+sl*s1] = sdcol; 
	diag_amount += sdcol->length; 
      } 
    } 
  } 
  if(!hasAli) free_tmp_pdist(tmp_dist, pdist->max_dlen); 
  all_diags->diags= calloc(diag_amount, sizeof(struct diag*)); 
  if(all_diags->diags==NULL) error("find_all_diags(): (1) Out of memory !"); 
 
  ap=0; 
  for(s1=0;s1<sl;s1++) { 
    for(s2=s1+1;s2<sl;s2++) { 
      sdcol=all_diags->diag_matrix[s1+sl*s2]; 
      if(sdcol!=NULL) { 
	for(sp=0;sp<sdcol->length;sp++) { 
	  //	  if(hasAli || (sdcol->data[sp]->weight >0.01*(sdcol->total_weight/sdcol->length))) { 
	  sdcol->data[sp]->pred_diag = NULL;
	  all_diags->diags[ap]=sdcol->data[sp]; 
	  ap++; 
	  //} else { 
	  //  free(sdcol->data[sp]); 
	  //diag_amount--; 
	  // } 
	} 
      } 
    } 
  } 
   
  all_diags->seq_amount = sl; 
   
  for(s1=0;s1<sl;s1++) { 
    for(s2=s1+1;s2<sl;s2++) { 
      sdcol=all_diags->diag_matrix[sl*s1+s2]; 
      if(sdcol!=NULL) { 
	sdcol->weight_fac =pow(sdcol->total_weight/total,2.0); 
	for(sp=0;sp<sdcol->length;sp++) { 
	  dg = sdcol->data[sp]; 
	  //	  if(hasAli) print_diag(dg); 
	  dg->weight_fac = sdcol->weight_fac; 
	  //	  dg->ov_weight = sdcol->total_weight; 
	  /* 
	  if(1 || !hasAli) { 
	    dg->weight_fac = sdcol->total_weight*totala[s1]/(sl-1)*sdcol->total_weight*totala[s2]/(sl-1); 
	  } 
	*/ 
	  //dg->weight_fac =pow(sdcol->total_weight*(sl-1),2.0)/(totala[s1]*totala[s2]); 
 
	  if(!hasAli) {
	    if(para->DO_OVERLAP) {
	      //printf(" do overlap\n");
	      //dg->weight_fac = 1.0;
	      calc_ov_weight(dg,all_diags, smatrix,pdist); 
	    }
	    dg->total_weight = (dg->weight);//+dg->ov_weight);// *dg->weight_fac; 
	  } else {
	    // changed in TX 1.0.0
	    if(para->FAST_MODE) dg->weight_fac = 1.0;
	    dg->total_weight = (dg->weight);//+dg->ov_weight);// *dg->weight_fac; 
	  }
	  //if(sp==sdcol->length-1) sdcol->total_weight *= dg->weight_fac;
	} 
      } 
    } 
  } 
   
  all_diags->diag_amount = diag_amount; 

  if(! hasAli) build_guide_tree(all_diags);
  return all_diags; 
} 


 
/** 
 * Finds all diags of each pair of sequences in in_seq_col by using 
 * the function above 
 * 
 * The pointer returned (and the ones included in the struct)  
 * has to be deallocted explicitely from memory. 
 */ 
struct diag_col *old_find_all_diags(struct scr_matrix *smatrix,  
				    struct prob_dist *pdist, 
				    struct seq_col *in_seq_col, 
				    struct alignment *algn,
				    int round ) { 
  unsigned int s1, s2, rs2, sl = in_seq_col->length, sp, ap; 
  struct diag_col *all_diags = create_diag_col(sl); 
  struct simple_diag_col *sdcol; 
 
  unsigned int diag_amount = 0; 
  struct diag *dg; 
 
  char hasAli = (round >1);//(algn!=NULL); 
 
  long double **tmp_dist = NULL; 
  if(!hasAli) tmp_dist = create_tmp_pdist(pdist); 
 
  int s2max = sl; 
  int s2width =(int) sqrt(sl); 
 
  double total=0.0; 
  //double imp[sl]; 
  //  for(s1=0;s1<sl;s1++) { 
  //  imp[s1] = 0.0; 
  //} 
  //double totala[sl]; 
  //memset(totala,0,sizeof(double)*sl); 
  for(s1=0;s1<sl;s1++) { 
    if(para->FAST_PAIRWISE_ALIGNMENT && s2width+1<sl) { 
      s2max = s1+s2width+1; 
      //printf("%i %i\n",s1, s2max); 
    } 
    //printf("before enter %i\n", s1); 
    for(s2=s1+1;s2<s2max;s2++) { 
      rs2 = s2 % sl; 
      if(!hasAli) { 
	fill_tmp_pdist(pdist,tmp_dist,in_seq_col->seqs[s1].length,in_seq_col->seqs[rs2].length ); 
      } 
      if(para->DEBUG>5) printf("%i %i\n", s1,s2); 
      //time1 = clock(); 
      //sdcol=find_diags_dyn(smatrix, pdist, &in_seq_col->seqs[s1], 
      //		 &in_seq_col->seqs[s2],algn,tmp_dist); 
       
      /* 
      doAlign = 1; 
      if(hasAli) { 
	if(algn->redo_seqs[s1*sl+s2]==0)  
	  doAlign = 0; 
      } 
       
      if(doAlign) { 
      */ 
      if(in_seq_col->seqs[s1].length > 0 && in_seq_col->seqs[s2].length > 0) { 
	//if(para->DEBUG>1) printf(" %i %i %i\n", s1, rs2,sl-1); 
	//	printf("find diags %i %i\n",s1,s2); 
        sdcol=find_diags_dialign(smatrix, pdist, &in_seq_col->seqs[s1], 
				 &in_seq_col->seqs[rs2],algn,tmp_dist, round); 
	//	imp[s1] += sdcol->total_weight; 
	//imp[s2] += sdcol->total_weight; 
	total += sdcol->total_weight; 
	//totala[s1] +=sdcol->total_weight;  
	//totala[s2] +=sdcol->total_weight;  
	//printf(" num of diags:%i\n ", sdcol->length); 
	/* 
	  } else { 
	  sdcol = calloc(1, sizeof(struct simple_diag_col)); 
	  sdcol->length = 0; 
	  } 
	*/ 
	//printf("%i %i %f\n", s1,s2, (clock()-time1)/CLOCKS_PER_SEC); 
	 
	all_diags->diag_matrix[s1+sl*rs2] = sdcol; 
	all_diags->diag_matrix[rs2+sl*s1] = sdcol; 
	diag_amount += sdcol->length; 
      } 
    } 
  } 

  
  // new: call of split function
  //printf("before split\n");
  /*
  if(! hasAli) {
    diag_amount = 0;
    split_diags(in_seq_col, all_diags);
    //total=0.0;
    for(s1=0;s1<sl;s1++) { 
      for(s2=s1+1;s2<sl;s2++) { 
	sdcol=all_diags->diag_matrix[s1+sl*s2]; 
	if(sdcol!=NULL) {
	  //sdcol->total_weight = 0.0;
	  for(sp=0;sp<sdcol->length;sp++) { 
	    dg = sdcol->data[sp]; 
	    real_calc_weight(dg, smatrix, pdist, 0, NULL,algn);
	    //sdcol->total_weight += dg->weight;
	    //total += dg->weight;
	  }
	  diag_amount +=sdcol->length;
	}
      }
    }
  }
  */
  //all_diags->diag_amount = diag_amount;
  //  printf("after split\n");
  //double max = 0.0;


  all_diags->diags= calloc(diag_amount, sizeof(struct diag*)); 
  if(all_diags->diags==NULL) error("find_all_diags(): (1) Out of memory !"); 

  double sdtotal;
 
  ap=0; 
  for(s1=0;s1<sl;s1++) { 
    for(s2=s1+1;s2<sl;s2++) { 
      sdcol=all_diags->diag_matrix[s1+sl*s2]; 
      if(sdcol!=NULL) { 
	for(sp=0;sp<sdcol->length;sp++) { 
	  //	  if(hasAli || (sdcol->data[sp]->weight >0.01*(sdcol->total_weight/sdcol->length))) { 
	    all_diags->diags[ap]=sdcol->data[sp]; 
	    ap++; 
	    //} else { 
	    //  free(sdcol->data[sp]); 
	    //diag_amount--; 
	    // } 
	} 
      } 
    } 
  } 
   
  all_diags->seq_amount = sl; 
   
  all_diags->total_weight = 0.0;
  for(s1=0;s1<sl;s1++) { 
    for(s2=s1+1;s2<sl;s2++) { 
      sdcol=all_diags->diag_matrix[sl*s1+s2]; 
      sdcol->total_weight = 0.0;
      sdtotal = 0.0;
      if(sdcol!=NULL) { 
	for(sp=0;sp<sdcol->length;sp++) { 
	  dg = sdcol->data[sp]; 
	  //	  if(hasAli) print_diag(dg); 
	  dg->weight_fac = pow(sdcol->total_weight/total,2.0); 
	  //dg->weight_fac = 1.0;

	  //	  dg->ov_weight = sdcol->total_weight; 
	  /* 
	  if(1 || !hasAli) { 
	    dg->weight_fac = sdcol->total_weight*totala[s1]/(sl-1)*sdcol->total_weight*totala[s2]/(sl-1); 
	  } 
	*/ 
	  //	  dg->weight_fac =pow(sdcol->total_weight*(sl-1),2.0)/(totala[s1]*totala[s2]); 
 
	  if(!hasAli) {
	    if(para->DO_OVERLAP) {
	      //dg->weight_fac = 1.0;
	      //error("calc ov\n");
	      calc_ov_weight(dg,all_diags, smatrix,pdist); 
	    }
	    //dg->total_weight = (dg->weight+dg->ov_weight) *dg->weight_fac; 
	  } else {
	    dg->weight_fac = 1.0;
	  }
	  dg->total_weight = (dg->weight+dg->ov_weight);//*dg->weight_fac; 
	  //if(dg->total_weight > max) max = dg->total_weight;
	  //dg->weight = dg->score;
	  //	  print_diag(dg);
	  sdtotal += dg->weight;
	  all_diags->total_weight += dg->total_weight;
	}
	sdcol->total_weight = sdtotal;
      }
    } 
  } 

  all_diags->average_weight = all_diags->total_weight / all_diags->diag_amount;
  //printf(" max %f\n", max);
  if(!hasAli) free_tmp_pdist(tmp_dist, pdist->max_dlen); 
  all_diags->diag_amount = diag_amount; 

  if(! hasAli) build_guide_tree(all_diags);
  return all_diags; 
} 
