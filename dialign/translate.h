/**
 *
 * translate.h:  
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
    /*        functions from translate.c        */
    /*                                           */
    /*********************************************/

void translate_sequence_collection(struct seq_col *in_sequence_collection);
void translate_sequence(struct seq *in_sequence);
char translate(char first_base, char second_base, char third_base, char *dna_number);
char inverse(char base);
void retranslate_sequence(struct seq *in_sequence);
char* retranslate(char amino);
void translate_sequence_collection_orf_frame(struct seq_col *in_seq_col);
void translate_sequence_collection_default(struct seq_col *in_seq_col);
