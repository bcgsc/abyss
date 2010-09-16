/**
 *
 * io.h:  
 *
 * 2004-08-30  Dorothea Emig Volker Menrad
 *             
 */

    /************************************************/
    /*                                              */
    /*                  structs                     */
    /*                                              */
    /************************************************/


    /************************************************/
    /*                                              */
    /*              global variable                 */
    /*                                              */
    /************************************************/

    /*********************************************/
    /*                                           */
    /*        functions from io.c                */
    /*                                           */
    /*********************************************/
void version(); // prints version on stdout
void print_scr_matrix(struct scr_matrix* aSmatrix);
struct scr_matrix* read_scr_matrix(char *filename);
void print_seq(struct seq* aSeq);
struct seq_col* read_fasta(char *filename);
void print_diag(struct diag* aDiag);
struct prob_dist* read_diag_prob_dist(struct scr_matrix* smatrix, char *filename);
void simple_print_alignment_default(struct alignment *algn);
void simple_print_alignment_dna_retranslate(struct alignment *algn);
//void simple_print_alignment_dna(struct alignment *algn);
void fasta_print_alignment_default(struct alignment *algn, char *filename);
void fasta_print_alignment_dna_retranslate(struct alignment *algn, char *filename);
//void fasta_print_alignment_dna(struct alignment *algn, char *filename);
char *print_info(struct alignment *align);
char *output_line(char *string);
char *output_line_left(char *string);
char *blank_line();
void print_pdist_matrix(struct prob_dist *sdist,char *filename);
void error(char *message);
void merror(char *msg1, char *msg2);
char* build_pathname(char *dir, char *file);
