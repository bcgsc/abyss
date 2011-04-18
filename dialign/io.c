/**
 *
 * io.c: IO operations, e.g. read input files, console output etc.
 *
 * A.R.Subramanian
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
#include "io.h"


// alig.c
extern struct algn_pos *find_eqc(struct algn_pos **ap, int seqnum, int pos);
extern void prepare_alignment(struct alignment *algn);


extern struct diag* create_diag(int n1, struct seq* sq1, unsigned int sp1,  
				int n2, struct seq* sq2, unsigned int sp2, 
				int dlength);

extern int errno;


/**
 * print version
 */
void version() {
 
  printf("  This is DIALIGN-TX Version %s - A Multiple Sequence alignment program.\n",para->VERSION );
  printf("             Author: Amarendran R. Subramanian, 2004-2008 \n");
  printf("                        subraman@informatik.uni-tuebingen.de\n\n");
  printf("                    \n"); 
  printf("   Research work using DIALIGN-TX should cite:\n\n");
  printf("   DIALIGN-TX: improvement of the segment-based approach for multiple\n");
  printf("   sequence alignment by combining greedy and progressive alignment strategies\n");
  printf("   Amarendran R. Subramanian, Michael Kaufmann, Burkhard Morgenstern,\n");
  printf("   Algorithms for Molecular Biology 3:6, 2008\n\n");
  printf("   DIALIGN-T: An improved algorithm for segment-based multiple sequence alignment\n");
  printf("   Amarendran R. Subramanian, Jan Weyer-Menkhoff, Michael Kaufmann,\n");
  printf("   Burkhard Morgenstern, BMC Bioinformatics 6:66, 2005\n");
  //  printf("                    Special Thanks to:\n");
  //printf("                       Burkhard Morgenstern\n");
  //printf("                       Michael Kaufmann\n"); 
  //printf("                       David Mathog\n"); 
  //printf("                       Numereous DIALIGN-T/X Users\n"); 
}

/**
 * print error message and exit
 */
void error(char *message)
{
  printf("ERROR: %s\n", message);
//  if(errno) perror("perror()");
  exit(1);
}

/**
 * print error message and exit
 */
void merror(char *msg1, char *msg2)
{
  printf("ERROR: %s %s\n", msg1, msg2);
//  if(errno) perror("perror()");
  exit(1);
}

/**
 * strips off leading whitespace characters
 */
void strip_leading_ws( char *str ) {
  int s,d;
  for(s=d=0; str[s]==' ' || str[s]=='\t'; s++){}
  for(;;s++,d++){ 
    str[d]=str[s];   
    if( str[s] == '\0')break;
  }
}

/**
 * builds a pathname from a dir-name and a filename.
 *
 * The pointer returned (and the ones included in the struct) 
 * has to be deallocted explicitely from memory.
 */
char* build_pathname(char *dir, char *file) {
  int dirslen = strlen(dir);

  char *pathn = calloc(dirslen+strlen(file)+2, sizeof(char));
  if(pathn==NULL) error("build_pathname(): out of memory !");
  
  strcpy(pathn, dir);
  if(dirslen>0 && dir[dirslen-1]!='/') strcat(pathn, "/");
  strcat(pathn, file);
  //  if(para->DEBUG)
  //    printf("DEBUG build_pathname(): Created filename: %s from dir=%s and file=%s\n", pathn, dir, file);

  return pathn;
}


/**
 * prints a sequence
 */
void print_seq(struct seq * aSeq) {
  int row;
  int slen = aSeq->length;
  int maxrow = slen/para->PRINT_SEQ_LINE_LENGTH;
  int row_residue = slen % para->PRINT_SEQ_LINE_LENGTH;

  printf("Sequence: %s\nLength: %i\n", aSeq->name, slen);


  char line[para->PRINT_SEQ_LINE_LENGTH+1]; 
  
  for(row=0;row <=maxrow; row++) {
    if(row<maxrow) {
      strncpy(line, &(aSeq->data[row*para->PRINT_SEQ_LINE_LENGTH]),
	      para->PRINT_SEQ_LINE_LENGTH);
      line[para->PRINT_SEQ_LINE_LENGTH]='\0';
    } else{
      if(row_residue==0) break;
      strncpy(line, &(aSeq->data[row*para->PRINT_SEQ_LINE_LENGTH]),
	      row_residue);
      line[row_residue]='\0';
    }
    printf("%s\n", line);
  }

}

/**
 * prints a diagional
 */
void print_diag(struct diag *aDiag) {
  int row;
  int slen = aDiag->length;
  int maxrow = slen/para->PRINT_SEQ_LINE_LENGTH;
  int row_residue = slen % para->PRINT_SEQ_LINE_LENGTH;

  printf("Diag: %s\n      %s\nLength: %i startpos1: %i startpos2: %i\n", aDiag->seq_p1.sq->name, 
	 aDiag->seq_p2.sq->name, slen, aDiag->seq_p1.startpos,aDiag->seq_p2.startpos );

  char *data1 =aDiag->seq_p1.sq->data;
  char *data2 =aDiag->seq_p2.sq->data;
  char *data;
  
  int startpos1 = aDiag->seq_p1.startpos;
  int startpos2 = aDiag->seq_p2.startpos;
  int snum, startpos;

  char line[para->PRINT_SEQ_LINE_LENGTH+1]; 

  for(row=0;row <=maxrow; row++) {
    for (snum=0;snum<2;snum++) {
      startpos = (snum==0 ? startpos1 : startpos2);
      data = (snum==0 ? data1 : data2);
      if(row<maxrow) {
	strncpy(line, &(data[startpos+row*para->PRINT_SEQ_LINE_LENGTH]),
		para->PRINT_SEQ_LINE_LENGTH);
	line[para->PRINT_SEQ_LINE_LENGTH]='\0';
      } else{
	if(row_residue==0) break;
	strncpy(line, &(data[startpos+row*para->PRINT_SEQ_LINE_LENGTH]),
		row_residue);
	line[row_residue]='\0';
      }
      printf("%s\n", line);
    }
    if(row<maxrow) printf("\n");
  }
  printf("Score: %li Weight: %e \n", aDiag->score, aDiag->weight);
}

/**
 * prints a score matrix
 */
void print_scr_matrix(struct scr_matrix *aSmatrix) {
  int len =aSmatrix->length;
  printf("Length: %i, maximal score: %i\n", len, aSmatrix->max_score);

  int r, c;
  for (r=-1; r<len;r++) {
    for (c=0; c<=len;c++) {
      if(r==-1) {
	if(c<len) printf("  %c ", aSmatrix->num2char[c]);
      } else if(c<len) {
	printf("%3i ", aSmatrix->data[len*r+c]);
      } else {
	printf(" %c ", aSmatrix->num2char[r]);
      }
    }
    printf("\n");
  }
}

/**
 * reads score matrix from the file
 * indicated by filename parameter.
 *
 * The pointer returned (and the ones included in the struct) 
 * has to be deallocted explicitely from memory.
 */
struct scr_matrix* read_scr_matrix(char *filename) {

  if(para->DEBUG)
    printf("DEBUG read_scr_matrix(): Processing input file: %s\n", filename);

  int length = 0;

  struct scr_matrix *smatrix = calloc(1, sizeof(struct scr_matrix));
  int *char2num = (smatrix->char2num =calloc(256, sizeof(int)));
  int *num2char = (smatrix->num2char =calloc(256, sizeof(int)));

  if(smatrix==NULL || char2num==NULL) error("read_scr_matrix(): Out of memory !");
  
  FILE *fp;

  char rline[100];

  if( (fp = fopen( filename , "r")) == NULL) { 
    merror("read_scr_matrix(): Cannot open input score matrix file", filename );
  }
  
  int sret;
  char amino;
  while( (sret=fscanf(fp, " %c ", &amino))!=EOF && sret>0) {
    if(amino=='#') {
      break;
    } else {
      num2char[length]=amino;
      char2num[(int)amino]=length++;
    }
  }

  num2char[length]='?';
  char2num[(int)'?']=length++;
  num2char[length]='#';
  char2num[(int)'#']=length++;
  num2char[length]='$';
  char2num[(int)'$']=length++;
  int additional=3;

  if(sret==EOF) merror("read_scr_matrix(): Unexpected end of file ",filename);
  if(length==0) merror("read_scr_matrix(): Invalid format of file ",filename);

  smatrix->length = length; 
  int *data = (smatrix->data = calloc(length*length, sizeof(int)));
  if(data==NULL) error("read_scr_matrix(): Out of memory when allocating data !");

  // read the matrix entries
  int r,c;
  int is;//,tis;
  long double frac = (long double)2.0/(long double)(length*length);
  long double avg_score = 0.0;
// generate Array from Matrix 
  for( r=0; r<length; r++) {
    for( c=r; c<length; c++) {
      // check whether it is a regular acid or a special character like '$',...
      if( (r<length-additional) && (c<length-additional)) {
		fscanf( fp, "%i", &is);
      } else {
		is = 0;
      }
      //tis = is;
      is +=para->SCR_MATRIX_ADD;
      if(smatrix->max_score<is) smatrix->max_score =is;
      avg_score += (long double)is*frac;
      data[length*r+c] = is;
      // ensure symmetry of the weight matrix
      data[length*c+r] = is;
    }
    fscanf(fp, "%s\n", rline);
  }
  fclose(fp);

// Array is a symmetric Matrix
/************************************ vvv Vorsicht !!! ***********/
  smatrix->avg_sim_score = para->PROT_SIM_SCORE_THRESHOLD;
  int ms = smatrix->max_score;
  
  int *dist = (smatrix->dist=calloc(ms+1, sizeof(int)));
  for(r=0;r<ms;r++) {
    dist[r]=0;
  }
// count how often score x appears and save in dist[x]
  for(r=0;r<length;r++) {
    for(c=0;c<length;c++) {
      if(num2char[r]!='X' && (num2char[r] >='A') && (num2char[r]<='Z'))
	if(num2char[c]!='X' && (num2char[c] >='A') && (num2char[c]<='Z'))
	  dist[data[length*r+c]]++;
    }
  }
  return smatrix;
}

/**
 * reads the probability distribution for diagonal lengths from the file
 * indicated by filename parameter.
 *
 * The pointer returned (and the ones included in the struct) 
 * has to be deallocted explicitely from memory.
 */
struct prob_dist* read_diag_prob_dist(struct scr_matrix* smatrix, char *filename) {

  if(para->DEBUG>1)
    printf("DEBUG read_diag_prob_dist(): Processing input file: %s\n", filename);

  int length = 0;

  struct prob_dist *pdist = calloc(1, sizeof(struct prob_dist));
  pdist->smatrix = smatrix;

  if(pdist==NULL) error("read_diag_prob_dist(): Out of memory !");
  
  FILE *fp;


  if( (fp = fopen( filename , "r")) == NULL) { 
    merror("read_diag_prob_dist(): Cannot open input file", filename );
  }
  
  int sret;
  sret=fscanf(fp, "%i\n", &length);
  if(sret<=0)     
    merror("read_diag_prob_dist(): Invalid format in file", filename );

  if(length==0) merror("read_scr_matrix(): Invalid format of file ",filename);

  //  fscanf(fp, "%s\n", rline);
  //  printf("rline:%s %i\n",rline, length);
  
  //  length=40;
  
  pdist->max_dlen = length; 
  pdist->data = calloc(length+1, sizeof(long double *));
  long double **dist =pdist->data;
  if(dist==NULL) error("read_diag_prob_dist(): (1) Out of memory when allocating data !");

  pdist->log_data = calloc(length+1, sizeof(double *));
  double **log_dist =pdist->log_data;
  if(log_dist==NULL) error("read_diag_prob_dist(): (1.1) Out of memory when allocating data !");


  // read the entries
  unsigned long i, scr, mxscr, sm_max_scr=smatrix->max_score;
  unsigned long ti, tscr; 
  long double weight;
  long size=0;

  for( i=1; i<=length; i++) {
    mxscr = i*sm_max_scr;
    size += mxscr+1;
    dist[i] = calloc(mxscr+1, sizeof(long double ));
    log_dist[i] = calloc(mxscr+1, sizeof(long double ));
    if(dist[i]==NULL) error("read_diag_prob_dist(): (3) Out of memory at iteration" );
    for(scr=0;scr<=mxscr;scr++) {
      dist[i][scr]=1.0;
    }
    for(scr=0;scr<=mxscr;scr++) {
      dist[i][scr]=1.0;
      fscanf( fp, "%li %li %Le\n", &ti,&tscr,&weight );
      //if(i!=ti || tscr!=scr) merror("read_scr_matrix(): (4) Invalid format of file ",filename);
      scr = tscr;
      if(weight==0.0) weight = 1.0;
      dist[i][scr]=weight;
      log_dist[i][scr]=-log(weight);
      if(para->DEBUG>5)printf("%li %li %Le\n", i, scr,dist[i][scr] );
    }
  }
  if(para->DEBUG >1) printf("DEBUG: PROB DIST SIZE: %li bytes\n", size*sizeof(long double)+length*sizeof(long double *));

  //pdist->avg_sim_score = 4; 

  return pdist;
}


/**
 * reads sequence collection (seq_col) from the file
 * indicated by filename parameter.
 *
 * The pointer returned (and the ones included in the struct) 
 * has to be deallocted explicitely from memory.
 */
struct seq_col* read_fasta(char *filename) {

  if(para->DEBUG)
    printf("DEBUG read_seq_col(): Processing input file: %s\n", filename);

  struct seq_col *scol = calloc(1, sizeof(struct seq_col));
  struct seq* seqs = (scol->seqs = calloc(para->MAX_SEQ_AMOUNT, sizeof(struct seq)));
  struct seq* seq;

  if(scol==NULL || seqs==NULL) error("read_fasta(): Out of memory !");

  FILE *fp;

  char rline[para->MAX_FASTA_LINE_LENGTH];

  if( (fp = fopen( filename , "r")) == NULL) { 
    merror("read_fasta(): Cannot open input FASTA file", filename );
  }
  
  scol->length = 0;
  char *data; char ch;
  unsigned int *slen;
  int data_maxlen;
  int c, rlen;
  int valid =0;
  char long_name=0;
  int name_length;

  while( fgets( rline , para->MAX_FASTA_LINE_LENGTH+1, fp ) != NULL ) {
    //    if(para->DEBUG)
    //      printf(rline);

    strip_leading_ws(rline);

    rlen = strlen(rline);

    if(rline[0] == '>') {
      valid = 1;
      long_name=1;
		
      seq = &(scol->seqs[scol->length]);
      seq->max_seen = 0;
      seq->name = calloc(rlen, sizeof(char));
      seq->num = scol->length;
	  seq->orf_frame=0;
	  seq->crick_strand=0;
      strncpy(seq->name, &(rline[1]), rlen-2);
      slen = (unsigned int *) &(seq->length);

      *slen = 0;

      data_maxlen = 1024;
      seq->data = (data = calloc(data_maxlen, sizeof(char)));
      if(data==NULL) error("read_fasta(): Out of memory: seq->data alloc");

      scol->length++;
	  if(rlen < para->MAX_FASTA_LINE_LENGTH) 
		{
			long_name = 0; // if relen =  then line is longer	
		}

    } 
	else if(long_name)
	{
		long_name=1;
		if(rlen < para->MAX_FASTA_LINE_LENGTH) 
		{
			long_name = 0; // if relen =  then line is longer
			rline[rlen-1]='\0';
		}

		
		name_length = strlen(seq->name);
		if(NULL == (seq->name = realloc(seq->name,rlen + name_length) ) )
		{
			error("read_fasta(): Out of memory: seq->data realloc");
		}
		
		strcat(seq->name, rline);
		
	}
	else {
      for( c=0;c<rlen;c++) {
	ch = rline[c];
	if(ch=='?' || ch=='$' || ch=='#') ch='X';
	if( (ch >= 65 && ch < 91) || (ch >= 97 && ch < 123)) {
	  if(! valid) merror("read_fasta(): File is not in FASTA format:", filename);
	  
	  // realloc memory when necessary
	  if(*slen >= data_maxlen) 
	  {
	    data_maxlen += 1024;
	    seq->data = (data = realloc(data, data_maxlen*sizeof(char)));
	    if(data==NULL) error("read_fasta(): Out of memory: seq->data alloc");
	  }
	  data[(*slen)++] = ((ch >= 97) ? toupper(ch) : ch);
	}
      }
    }
  }
	
  int avg=0;
  for(c=0;c<scol->length;c++) {
    avg += scol->seqs[c].length;
  }
  scol->avg_length = avg/scol->length;

  if(para->DEBUG)
    printf("\n");
  
  fclose(fp);
  return scol;
}

/**
 * reads the given anchor file and returns the accordingly created
 * simple_diag_col structure pointer
 */
struct simple_diag_col* read_anchors(char *filename, struct seq_col* scol) {

  if(para->DEBUG)
    printf("DEBUG read_anchors(): Processing anchor file: %s\n", filename);

  struct simple_diag_col *sdcol = malloc(sizeof(struct simple_diag_col));

  FILE *fp;

  if( (fp = fopen( filename , "r")) == NULL) { 
    merror("read_anchros(): Cannot open input anchor file", filename );
  }

  //fscanf( fp, "%li %li %Le\n", &ti,&tscr,&weight );

  long int s1,s2,sp1,sp2,len;
  double score;
  
  int alloc_size = 64;
  sdcol->data = malloc(sizeof (struct diag*)*alloc_size);
  sdcol->length=0;

  while( fscanf(fp,"%li %li %li %li %li %le\n",&s1,&s2,&sp1,&sp2,&len,&score ) == 6) {
    if(sdcol->length >= alloc_size) {
      alloc_size+=16;
      sdcol->data = realloc(sdcol->data,sizeof (struct diag*)*alloc_size);
    }
    sdcol->data[sdcol->length]= create_diag(s1-1,&(scol->seqs[s1-1]),sp1-1,
					    s2-1,&(scol->seqs[s2-1]),sp2-1,len);
    sdcol->data[sdcol->length]->anchor = 1;
    sdcol->data[sdcol->length]->meetsThreshold = 1;
    //printf(" total weight %e\n", score);
    sdcol->data[sdcol->length++]->total_weight = score;
  }

  if(para->DEBUG) printf("DEBUG  read_anchors(): Read %i anchors from file\n", sdcol->length);
  fclose(fp);
  return sdcol;
}


/**
 * prints the given alignment  in a simple way. 
 */
void simple_print_alignment_default(struct alignment *algn) {
  struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  unsigned int i,j,s,pos,max,tmax;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;
  int proc[slen];
  //  char proceed[slen];

  for(i=0;i<slen;i++) {
    proc[i]=0;
  }

  
  prepare_alignment(algn);
  max = algn->max_pos;
  struct algn_pos *ap1;

  //
  // print block
  // WARNING!  print_info reallocates memory, this loses the algn pointer
  //   since the new one is not stored when it is returned here, so it cannot
  //   be freed later. In the present version the leak is relatively minor and
  //   should not crash anything.
  //
	printf("%s",print_info(algn));

  printf("\n   ALIGNMENT OUTPUT:\n");
  printf("   -----------------\n\n");
  //  printf("%i\n", max);
  for(pos=0;pos<max;pos+=50) {
    tmax = pos+50;
    //printf("tmax: %i\n", tmax);
    
    if(tmax>max) tmax = max;
    for(s=0;s<slen;s++) {
      sq = &scol->seqs[s];
      for(j=pos;j<tmax;j++) {
	
	if( (j%50)==0) {
	  printf("%12.12s%10i   ", sq->name, proc[s]+1);
	} else if( (j%10)==0) {
	  printf(" ");
	}
	if(proc[s]<sq->length) {
	  ap1 = find_eqc(ap,s,proc[s]);
	  if( (*ap1->eqcAlgnPos) < j) {
	    printf ("\nALARM %i %i %i %i\n", s,j, proc[s], *ap1->eqcAlgnPos);
	  }
	  //	  if(proc[0]==244 && s==0) printf("\nMOVE FORWARD: %i %i %i\n",j, *ap[0][244].eqcAlgnPos, *ap[1][241].eqcAlgnPos);
	  
	  if( (*ap1->eqcAlgnPos) == j) {
	    if(ap1->state & para->STATE_ORPHANE) {
	      printf("%c", tolower(sq->data[proc[s]]));
	    } else {
	      printf("%c", toupper(sq->data[proc[s]]));
	    }
	    proc[s]++;
	  } else {
	    printf("-");
	    //printf("%i",*ap[s][proc[s]].maxpos);
	  }
	} else {
	  printf("-");
	  //printf("%i",*ap[s][proc[s]].maxpos);
	  //	  printf("\n%i %i %i\n", s, j,proc[s]);
	}
      }
      printf("\n");
    }
    printf("\n");
    printf("\n");
  }
}

void simple_print_alignment_dna_retranslate(struct alignment *algn) 
{
	char *tmp = "000";
  struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  unsigned int i,j,s,pos,max,tmax;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;
  int proc[slen];
  //  char proceed[slen];

  for(i=0;i<slen;i++) {
    proc[i]=0;
  }

  prepare_alignment(algn);
  max = algn->max_pos;
  struct algn_pos *ap1;

  //
  // print block
  //
	printf("%s",print_info(algn));
  printf("\n   ALIGNMENT OUTPUT:\n");
  printf("   -----------------\n\n");
  //  printf("%i\n", max);
  for(pos=0;pos<max;pos+=16) {
    tmax = pos+16;
    //printf("tmax: %i\n", tmax);
    
    if(tmax>max) tmax = max;
    for(s=0;s<slen;s++) {
      sq = &scol->seqs[s];
      for(j=pos;j<tmax;j++) {
	
	if( (j%16)==0) {
	  printf("%12.12s%10i   ", sq->name, ((proc[s])*3)+1);
	} else if( (j%4)==0) {
	  printf(" ");
	}
	if(proc[s]<sq->length) {
	  ap1 = find_eqc(ap,s,proc[s]);
	  if( (*ap1->eqcAlgnPos) < j) {
	    printf ("\nALARM %i %i %i %i\n", s,j, proc[s], *ap1->eqcAlgnPos);
	  }
	  //	  if(proc[0]==244 && s==0) printf("\nMOVE FORWARD: %i %i %i\n",j, *ap[0][244].eqcAlgnPos, *ap[1][241].eqcAlgnPos);
	  
	  if( (*ap1->eqcAlgnPos) == j) {

		tmp = retranslate(sq->dna_num[proc[s]]);
	     if(ap1->state & para->STATE_ORPHANE) {
			printf("%c%c%c", tolower(tmp[0]),tolower(tmp[1]), tolower(tmp[2]) );
	    } else {
			printf("%c%c%c", toupper(tmp[0]), toupper(tmp[1]), toupper(tmp[2]) );
	    }
	    proc[s]++;
	  } else {
	    printf("---");
	    //printf("%i",*ap[s][proc[s]].maxpos);
	  }
	} else {
	  printf("---");
	  //printf("%i",*ap[s][proc[s]].maxpos);
	  //	  printf("\n%i %i %i\n", s, j,proc[s]);
	}
      }
      printf("\n");
    }
    printf("\n");
    printf("\n");
  }

}


/**
 * prints the given alignment in fasta format 
 * to the given file. 
 */
void fasta_print_alignment_default(struct alignment *algn, char *filename) {
  struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  unsigned int j,s,proc, max;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;

  prepare_alignment(algn);
  max = algn->max_pos;
  struct algn_pos *ap1;

  FILE *fp;
  if( (fp = fopen( filename , "w")) == NULL) { 
    merror("fasta_print_alignment(): Cannot open  file", filename );
  }
//  fprintf(fp,"%s",print_info(algn));
  max = algn->max_pos;
  for(s=0;s<slen;s++) {
    sq = &(scol->seqs[s]);
    fprintf(fp, ">%s",sq->name);
    proc = 0;
    for(j=0;j<max;j++) {
      if(proc <sq->length) {
	if( (j%60)==0) fprintf(fp,"\n");
	ap1 = find_eqc(ap,s,proc);
	if(*ap1->eqcAlgnPos==j) {
	  if(ap1->state & para->STATE_ORPHANE) {
	    fprintf(fp, "%c", tolower(sq->data[proc]));
	  } else {
	    fprintf(fp, "%c", toupper(sq->data[proc]));
	  }
	  proc++;
	} else {
	  fprintf(fp,"-");
	}
      } else {
		if( (j%60)==0) fprintf(fp,"\n");
	fprintf(fp,"-");
      }
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void fasta_print_alignment_dna_retranslate(struct alignment *algn, char *filename) 
{
	char *tmp = "000";
	struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  unsigned int j,s,proc, max;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;

  prepare_alignment(algn);
  max = algn->max_pos;
  struct algn_pos *ap1;

  FILE *fp;
  if( (fp = fopen( filename , "w")) == NULL) { 
    merror("fasta_print_alignment(): Cannot open  file", filename );
  }
//  fprintf(fp,"%s",print_info(algn));
  max = algn->max_pos;
  for(s=0;s<slen;s++) 
	{
    sq = &(scol->seqs[s]);
    fprintf(fp, ">%s",sq->name);
    proc = 0;
    	for(j=0;j<max;j++) 
		{
    		if(proc <sq->length) 
			{
				if( (j%20)==0) fprintf(fp,"\n");
				ap1 = find_eqc(ap,s,proc);
				if(*ap1->eqcAlgnPos==j) 
				{
					tmp = retranslate(sq->dna_num[proc]);
//					printf(fp,"%c,\n",sq->dna_num[proc]);
	 				if(ap1->state & para->STATE_ORPHANE) {
 						fprintf(fp, "%c%c%c", tolower(tmp[0]),tolower(tmp[1]), tolower(tmp[2]));
	 				}
					else 
					{
	   					fprintf(fp, "%c%c%c", toupper(tmp[0]), toupper(tmp[1]), toupper(tmp[2]));
	 				}
					proc++;
				} 
				else 
				{
	  				fprintf(fp,"---");
				}
			}
	 		else 
			{
				if( (j%20)==0) fprintf(fp,"\n");
				fprintf(fp,"---");
			}
		}
   		 	fprintf(fp,"\n");
	}
  	fclose(fp);
}

char* print_info(struct alignment *algn)
{
	
	int i;
	char *output;
	char *line, *line2;

	if(NULL == ( output = (calloc(63,sizeof(char)))))
	{
		error("print_info(): Out of memory !");
	}
	strcat(output, "\n");
	for(i=0;i!=60; ++i)
	{
		strcat(output, "*");
	}
	strcat(output, "\n");
	if (para->DNA_TRANSLATION) {
		if(para->FIND_ORF){
			if(!para->ORF_FRAME){
//				-L :
				if(NULL == ( line = (calloc(62, sizeof(char)))))
				{
					error("print_info(): Out of memory !");
				}
				if(NULL == ( output = (realloc(output,strlen(output)+15*61))))
				{
					error("print_info(): Out of memory !");
				}
				line = "Multiple Sequence Alignment (with translation)";
				line = output_line(line);
				strcat(output, line);
				line = "Input sequences in DNA";
				line = output_line(line);
				strcat(output, line);
				if(!para->OUTPUT) line = "Alignment output in aminoacids";
				else line = "Alignment output in DNA";
				line = output_line(line);
				strcat(output, line);
				line = "Sequences translated into aminoacids";
				line = output_line(line);
				strcat(output, line);
				line = "Only longest open reading frames aligned";
				line = output_line(line);
				strcat(output, line);
				line = "Sequence lengths cut mod 3 = 0";
				line = output_line(line);
				strcat(output, line);
				strcat(output,blank_line());
				line = "reading frame 1 :      123 123 123 123 ...";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 2 :    X 123 123 123 123 ...";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 3 :   XX 123 123 123 123 ...";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 4 :  ... 321 321 321 321 XX";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 5 :  ... 321 321 321 321 X";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 6 :  ... 321 321 321 321";
				line = output_line_left(line);
				strcat(output, line);
				strcat(output,blank_line());
				strcat(output,blank_line());
				for (i = 0; i != algn->scol->length; ++i)
				{
					int k;
					char *tmp;
					int tmp2;
					tmp = " :    reading frame  =  ";
					line[0] = algn->scol->seqs[i].name[0];
					
					for(k=1; k!=12 ; ++k)
					{
						line[k]=algn->scol->seqs[i].name[k];
						tmp2=algn->scol->seqs[i].orf_frame+48;
					}

					line[12]='\0';
					strcat(line, tmp);
					line[strlen(line)-1]=tmp2;
					line[strlen(line)]='\0';
					if(NULL == ( output = (realloc(output,strlen(output)+62))))
					{
						error("print_info(): Out of memory !");
					}
					line = output_line(line);
					strcat(output, line);

				} 
				free(line);	
			}

			else{
//				-O :
				if(NULL == ( line = (calloc(62, sizeof(char)))))
				{
					error("print_info(): Out of memory !");
				}
				if(NULL == ( output = (realloc(output,strlen(output)+15*61))))
				{
					error("print_info(): Out of memory !");
				}
				line = "Multiple Sequence Alignment (with translation)";
				line = output_line(line);
				strcat(output, line);
				line = "Input sequences in DNA";
				line = output_line(line);
				strcat(output, line);
				if(!para->OUTPUT) line = "Alignment output in aminoacids";
				else line = "Alignment output in DNA";
				line = output_line(line);
				strcat(output, line);
				line = "Sequences translated into aminoacids";
				line = output_line(line);
				strcat(output, line);
				line = "reading frames found due to longest ORF";
				line = output_line(line);
				strcat(output, line);
				line = "Sequence lengths cut mod 3 = 0";
				line = output_line(line);
				strcat(output, line);
				strcat(output,blank_line());
				line = "reading frame 1 :      123 123 123 123 ...";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 2 :    X 123 123 123 123 ...";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 3 :   XX 123 123 123 123 ...";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 4 :  ... 321 321 321 321 XX";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 5 :  ... 321 321 321 321 X";
				line = output_line_left(line);
				strcat(output, line);
				line = "reading frame 6 :  ... 321 321 321 321";
				line = output_line_left(line);
				strcat(output, line);
				strcat(output,blank_line());
				strcat(output,blank_line());
				for (i = 0; i != algn->scol->length; ++i)
				{
					int k;
					char *tmp;
					int tmp2;
					tmp = " :    reading frame  =  ";
					line[0] = algn->scol->seqs[i].name[0];
					
					for(k=1; k!=12 ; ++k)
					{
						line[k]=algn->scol->seqs[i].name[k];
						tmp2=algn->scol->seqs[i].orf_frame+48;
					}

					line[12]='\0';
					strcat(line, tmp);
					line[strlen(line)-1]=tmp2;
					line[strlen(line)]='\0';
					if(NULL == ( output = (realloc(output,strlen(output)+62))))
					{
						error("print_info(): Out of memory !");
					}
					line = output_line(line);
					strcat(output, line);

				} 
				free(line);
			}
		}
		else{
//			-T :
			if(NULL == ( line = (calloc(62, sizeof(char)))))
			{
				error("print_info(): Out of memory !");
			}
			if(NULL == ( output = (realloc(output,strlen(output)+5*61))))
			{
				error("print_info(): Out of memory !");
			}
			line = "Multiple Sequence Alignment (with translation)";
			line = output_line(line);
			strcat(output, line);
			line = "Input sequences in DNA";
			line = output_line(line);
			strcat(output, line);
			if(!para->OUTPUT) line = "Alignment output in aminoacids";
			else line = "Alignment output in DNA";
			line = output_line(line);
			strcat(output, line);
			line = "Sequences translated into aminoacids";
			line = output_line(line);
			strcat(output, line);
			line = "Sequence lengths cut mod 3 = 0";
			line = output_line(line);
			strcat(output, line);
			free(line);
		}
	}		
	else{
// 		-D :
		if(NULL == ( line = (calloc(62, sizeof(char)))))
		{
			error("print_info(): Out of memory !");
		}
		if(NULL == ( line2 = (calloc(62, sizeof(char)))))
		{
			error("print_info(): Out of memory !");
		}
		if(NULL == ( output = (realloc(output,strlen(output)+3*61+1))))
		{
			error("print_info(): Out of memory !");
		}
                (void) strcpy(line,"Multiple Sequence Alignment");
		line = output_line(line);
		strcat(output, line);
		if(para->DNA_PARAMETERS) (void) strcpy(line2,"in DNA");
		else (void) strcpy(line2,"in aminoacids");
		sprintf(line,"Input sequences %s", line2);
		line = output_line(line);
		strcat(output, line);
		if(para->DNA_PARAMETERS) (void) strcpy(line2,"in DNA");
		else (void) strcpy(line2,"in aminoacids");
		sprintf(line,"Alignment output %s", line2);
		line = output_line(line);
		strcat(output, line);
		free(line);
		free(line2);
		
	}
	if(NULL == ( output = (realloc(output,strlen(output)+63))))
	{
		error("print_info(): Out of memory !");
	}
	
	for(i=0;i!=60; ++i)
	{
		strcat(output, "*");
	}
	strcat(output, "\n\n");
	return output;
}


/*SIDE EFFECT, releases original memory held by string */
char* output_line(char *string)
{
	char *tmp;
	if(NULL == ( tmp = (calloc(62, sizeof(char)))))
	{
			error("print_info(): Out of memory !");
	}
	int x, i, y;
	y = 0;
	if (strlen(string) %2 == 0) y = 1;
	x = (60 - strlen(string)) /2 ;
	strcat(tmp, "*");
	for (i = 0; i != x-y; ++i)
	{
		strcat(tmp, " ");
	}
	strcat(tmp,string);
	for (i = 0; i != x-1; ++i)
	{
		strcat(tmp, " ");
	}
	strcat(tmp, "*\n");
        free(string);
	string=&(tmp[0]);
	return string;
}

/*SIDE EFFECT, releases original memory held by string */
char* output_line_left(char *string)
{
	char *tmp;
	if(NULL == ( tmp = (calloc(62, sizeof(char)))))
	{
			error("print_info(): Out of memory !");
	}
	int x, i, y;
	y = 0;
	x = (60 - strlen(string)-3) ;
	strcat(tmp, "* ");
	strcat(tmp,string);
	for (i = 0; i != x; ++i)
	{
		strcat(tmp, " ");
	}
	strcat(tmp, "*\n");
        free(string);
	string=&(tmp[0]);
	return string;
}
char* blank_line()
{
	int i;
	char *string;

	if(NULL == ( string = (calloc(62,sizeof(char)))))
	{
		error("print_info(): Out of memory !");
	}
	strcat(string, "*");
	for (i = 0; i != 58; ++i)
	{
		strcat(string, " ");
	}

	strcat(string, "*\n");
	return string;

}

void print_pdist_matrix(struct prob_dist *sdist, char *filename)
{
	int lh, scr, mxscr;
	FILE *fp;
  	if( (fp = fopen( filename , "w")) == NULL) 
	{ 
    	merror("fasta_print_alignment(): Cannot open  file", filename );
  	}
	fprintf(fp,"%d\n",sdist->max_dlen);
 	for(lh=1; lh!=sdist->max_dlen+1; ++lh)
	{
		mxscr = lh * sdist->smatrix->max_score;
		for(scr=0; scr!=mxscr+1; ++scr) 
		{
			fprintf(fp,"%d %d %Le\n",lh,scr,sdist->data[lh][scr]);
		}
	}

}
