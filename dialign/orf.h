/**
 *
 * orf.h:  
 *
 * 2004-08-24  Dorothea Emig Volker Menrad
 *             
 */

    /************************************************/
    /*                                              */
    /*                  structs                     */
    /*                                              */
    /************************************************/


	struct orf
	{
		int length;		// length of found orf
		char *sequence;	// part of the sequence in the orf
		char finish;	
		char *dna_num;	// retranslation
	};


    /************************************************/
    /*                                              */
    /*              global variable                 */
    /*                                              */
    /************************************************/

    /*********************************************/
    /*                                           */
    /*        functions from parameters.c        */
    /*                                           */
    /*********************************************/
struct seq_col* set_longest_orf(struct seq_col *in_seq_col);
struct seq* orf_finder(struct seq *in_seq_col);
//char inverse(char base);
//char translate(char first_base, char second_base, char third_base,char* dna_number_for_retranslate);
