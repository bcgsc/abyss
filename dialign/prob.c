#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "struct.h"
#include "parameters.h"

extern void error(char *message);
extern void merror(char *msg1, char *msg2);

// io.c#
extern void print_diag(struct diag* aDiag);

// diag.c
//extern struct seq_part* create_seq_part(int num, struct seq* aSeq, 
//					long startpos);
extern struct diag* create_diag(struct seq_part* part1, struct seq_part* part2,
				int dlength);

/**
 *
 * prob.c: Probability measurement (Creation  of random sequences, diags,...)
 *
 * 2003-10-14  A.R.Subramanian
 *             (Initial)
 */


/**
 * auxiliary method that fills random data into the sequence (see next func.)
 */
void fill_random_seq(struct seq* sq, struct scr_matrix *smatrix) {
  long length = sq->length;
  char *data = sq->data;

  int *n2c = smatrix->num2char;
  float slen = (float)smatrix->length;
  int pos, rnd;
  for(pos=0; pos<length;pos++) {
    do {
      rnd = (int)( (slen*random())/(RAND_MAX+1.0));
      data[pos]=(char)(n2c[rnd]);
    } while(n2c[rnd]=='X' || (n2c[rnd]<'A') || (n2c[rnd]>'Z'));
  }

}

/**
 * create random squence of the given length using the characters
 * of the given score matrix
 *
 * The pointer returned (and the ones included in the struct) 
 * has to be deallocted explicitely from memory.
 */
struct seq* create_random_seq(struct scr_matrix *smatrix, int length) {
  struct seq* sq = calloc(1, sizeof(struct seq));
  sq->length=length;
  sq->data = calloc(length, sizeof(char));
  sq->name = calloc(strlen("random")+1, sizeof(char));
  strcpy(sq->name, "random");
  fill_random_seq(sq, smatrix);
  return sq;
}



/**
 * calculates the score distribution up to the given maximum diaglen
 *
 * The pointer returned and the included dist pointer has to be deallocted 
 * explicitely from memory.
 */
struct prob_dist* calc_score_dist(struct scr_matrix *smatrix, int mxdlen) { 
  long sm_max_scr = smatrix->max_score;
  long maxdlen = mxdlen;
  struct prob_dist *sdist = calloc(1, sizeof(struct prob_dist));

  long **edist = (calloc(maxdlen+1, sizeof(long  *)));
  long **tdist = (calloc(maxdlen+1, sizeof(long  *)));
  long double **dist = (sdist->data=calloc(maxdlen+1, sizeof(long double *)));

  if(sdist==NULL || dist==NULL || edist==NULL|| tdist==NULL) error("calc_score_dist(): Out of memory !");
  sdist->max_dlen = maxdlen;
  int *sdata = smatrix->data;
  int *smdist = smatrix->dist;

  sdist->smatrix = smatrix;
  sdist->max_dlen = maxdlen;

  int reduce = 3; // special characters '?', '#', '%', we don't account for
  if(smatrix->char2num['X']>0) reduce++; // don't account for 'X'

  long double square = (smatrix->length-reduce) * (smatrix->length-reduce);
  long scr;
  unsigned long i,j, scr2, mxscr, omxscr;

  for(i=1;i<=maxdlen;i++) {
    mxscr = i*sm_max_scr;
    dist[i] = calloc(mxscr+1, sizeof(long double ));
    edist[i] = calloc(mxscr+1, sizeof(long  ));
    tdist[i] = calloc(mxscr+1, sizeof(long  ));
    if(dist[i]==NULL) error("calc_score_dist(): Out of memory at iteration" );
    if(i==1) {
      for(j=0;j<=mxscr;j++) {
	dist[i][j] = smdist[j]; //static, wie oft kommt ein score j in der scr-Matrix vor, nur für Fragmentlänge 1, ansonsten initialisiere mit Null
	edist[i][j]=0;
	tdist[i][j]=0;
      }
    } else {
      for(scr=0;scr<=mxscr;scr++) {
		dist[i][scr]=0.0;
		edist[i][scr]=0;
		tdist[i][scr]=0;
		omxscr = (i-1)*sm_max_scr;
		for(j=0;j<=sm_max_scr;j++) {
		  if((scr-j)<=omxscr) dist[i][scr] += dist[1][j]*dist[i-1][scr-j]; // statisch wenn größer eins
		}
      }
    }
  }
  
  struct seq *sq1, *sq2;
  long max_experiment=550000, ex, found;
  struct diag dg;
  char *data1;
  char *data2;
  int *c2n = smatrix->char2num;
  int a1, a2,pos,dpos1,dpos2;
  long double experiment,oldprob,opf,opf2;
  int dlen = 100;
  //  int maxdlen=100;
  long double factor = (long double)dlen*dlen;

  //long double expect;
  for(i=1;i<=maxdlen;i++) {
//    printf(" Diag of length %i\n", i);
    mxscr = i*sm_max_scr;
    //expect =0.0;
    dlen = i;
    factor = (long double)dlen*dlen;
    for(scr=0;scr<=mxscr;scr++) {
      //expect += ((long double)scr) * dist[i][scr]*pow(1.0/square,i);;
      //printf("  max %i\n", mxscr);
      for(scr2=scr+1;scr2<=mxscr;scr2++) {
	dist[i][scr] += dist[i][scr2];
      }

      dist[i][scr] = dist[i][scr]*pow(1.0/square,i);
    }
  }
  // EXPERIMENTS:

  srandom((int)time(NULL));
  sq1 = create_random_seq(smatrix, 2*maxdlen);
  sq2 = create_random_seq(smatrix, 2*maxdlen);
  dg.seq_p1.sq = sq1;
  dg.seq_p2.sq = sq2;
  dg.length =2*maxdlen;
  //printf(" pre random\n");
  //	for(ex=0;(ex<max_experiment || found==0);ex++) {
  char seenmax=0;
  for(ex=0;(ex<max_experiment);ex++) {
    if( (ex%1000)==0) fprintf(stderr," Experiment %li\n", ex);
    fill_random_seq(sq1, smatrix);
    fill_random_seq(sq2, smatrix);
    data1 = sq1->data;
    data2 = sq2->data;
    for(dpos1=0;dpos1<=maxdlen;dpos1++) {
      for(dpos2=0;dpos2<=maxdlen;dpos2++) {
		dg.seq_p1.startpos = dpos1;//(int)( (100*random())/(RAND_MAX+1.0));
		dg.seq_p2.startpos = dpos2; //(int)( (100*random())/(RAND_MAX+1.0));
		dg.score = 0;
		for(pos=0;pos<maxdlen;pos++) {
		  a1 = c2n[(int) data1[dg.seq_p1.startpos+pos]];
		  a2 = c2n[(int) data2[dg.seq_p2.startpos+pos]];
		  dg.score+=(long)sdata[smatrix->length*a1+a2];
		  if(dpos1<=pos+1 && dpos2<=pos+1) // t6
		    tdist[pos+1][dg.score]=1; // kann mehrmals dasselbe sein?? kann immer wieder überschrieben werden, warum?
		}
      }
    }
    for(i=1;i<=maxdlen;i++) {
      //printf(" Diag of length %i\n", i);
      mxscr = i*sm_max_scr;
      seenmax =0;
      for(scr=mxscr;scr>=0;scr--) {
		if(!seenmax) {
	  		if(tdist[i][scr]==1) {
	    		edist[i][scr]=edist[i][scr]+1; //für jedes Exp. maximalen score der jeweiligen Fragmentlänge hochzählen
	    		seenmax = 1;
	  		}
		} 
		tdist[i][scr]=0; // alles wieder auf Null setzen !
     }
    }
  }
  free(sq1);
  free(sq2);


  for(i=1;i<=maxdlen;i++) {
    //    printf(" Diag of length %i\n", i);
    mxscr = i*sm_max_scr;
    //expect =0.0;
    dlen = i;
    factor = (long double)(dlen+1)*(dlen+1); // t6
    //factor = (long double)10201.0;
    for(scr=0;scr<=mxscr;scr++) {
      for(scr2=scr+1;scr2<=mxscr;scr2++) {
	edist[i][scr] += edist[i][scr2]; // da in scr2 ja auch schon der scr1 beinhaltet ist; durch aufsummieren bekommen wir die Anzahl, wie viele Fragmente der Länge i min. den score scr haben (oder mehr)
      }

      found = edist[i][scr]; // random exp.
      oldprob = dist[i][scr]; // statisch
      opf2 = (long double)oldprob*factor;
      opf = (long double)1.0-pow(1.0-oldprob,factor);
      experiment =((long double)found)/ ( ((long double)(max_experiment)));
      //if(switcher==0 && oldprob<0.000001) { // old 0.00001 t3

      //if(found==0) 
      //printf ("%Le %Le %Le %Le %Le %Le %i\n",experiment,opf,opf2,oldprob,factor,square,found);
      
      if(scr>0 && (found==0 || (experiment>opf || experiment>opf2 || experiment<10.0/max_experiment))) {
		if(  (opf<=dist[i][scr-1]) && ( (opf > opf2) || (opf2 > dist[i][scr-1]))) {
		  experiment = opf;
		}
		if( ( (opf2 <= dist[i][scr-1]))) {
	  		experiment = opf2;
		}
      }
      if(experiment == 0.0) {
		experiment = dist[i][scr-1];
      }

      dist[i][scr] = experiment;

      //if(isnan(dist[i][scr]) || isinf(dist[i][scr])) dist[i][scr]=1.0;
      if(para->DEBUG>1)printf("%li %li %Le\n", i, scr,dist[i][scr] );
    }
    //    printf("%i %Le\n", i, expect);
  }
  
  return sdist;
}

 
