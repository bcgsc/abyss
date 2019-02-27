/**
 *
 * parameters.h:  
 *
 * 2004-08-13  Dorothea Emig Volker Menrad
 *             
 */

/************************************************/
/*                                              */
/*                  structs                     */
/*                                              */
/************************************************/

struct parameters
{
	char *VERSION; 
	// 0 no debug statements
	// 1 debugs the current phase of the processing
	// 2 very loquacious debugging
	// 5 hardcore debugging
	int DEBUG; // = 0;

	// maximum amount of input sequences
	int MAX_SEQ_AMOUNT; // = 5000;

	// maximum number of characters per line in a FASTA file
	int MAX_FASTA_LINE_LENGTH; // = 100;

	// maximum amount of characters per line when printing a sequence
	int PRINT_SEQ_LINE_LENGTH; // = 80;

	/*******************************
	 * 
	 * PROTEIN SCORING/WEIGHTING SECTION
	 *
	 *****************************/


	// score matrix (e.g. BLOSUM) file name (in the configuration directory)
	//#define SCR_MATRIX_FILE_NAME "BLOSUM75.scr"
	char *SCR_MATRIX_FILE_NAME; // = "BLOSUM.scr";
	//#define SCR_MATRIX_FILE_NAME "BLOSUM90.scr"

	// defines the minimum weight when the weight is changed to
	//  1-pow(1-prob, factor)
	double DIAG_CALC_WEIGHT_THRESHOLD; // = 0.000000065;
	//#define DIAG_CALC_WEIGHT_THRESHOLD 0.0000002


	// diag prob. distribution  file name (in the configuration directory)
	char *DIAG_PROB_FILE_NAME; // = "BLOSUM.diag_prob_t10";   
	//#define DIAG_PROB_FILE_NAME "BLOSUM.diag_prob_t7"   // current
	//#define DIAG_PROB_FILE_NAME "BLOSUM75.diag_prob_t2"   // 
	//#define DIAG_PROB_FILE_NAME "BLOSUM75.diag_prob_t1"   
	//#define DIAG_PROB_FILE_NAME "BLOSUM90.diag_prob"   
	//#define DIAG_PROB_FILE_NAME "BLOSUM75.diag_prob"   
	//#define DIAG_PROB_FILE_NAME "tp400_prot"

	// add to each score (to prevent negative values)
	int SCR_MATRIX_ADD; // = 0; // BLOSUM(62)
	//#define SCR_MATRIX_ADD 5 // BLOSUM75
	//#define SCR_MATRIX_ADD 6 // BLOSUM90

	/*******************************
	 * 
	 * PROTEIN QUALITY SECTION
	 *
	 *****************************/

	// "even" sim score threshold for protein sequences alignment
	//#define PROT_SIM_SCORE_THRESHOLD 6 //BLOSUM90
	//#define PROT_SIM_SCORE_THRESHOLD 5 //BLOSUM75
	double PROT_SIM_SCORE_THRESHOLD; // = 4; // BLOSUM62

	// maximum number of consecutive positions for frag-window
	int PROT_DIAG_MAX_UNDER_THRESHOLD_POS; // = 4; // old was 4

	// minimum diagonal length for breaking
	double PROT_DIAG_MIN_LENGTH_THRESHOLD; // = 40.0; // old was 40

	// minimal allowed average csore in frag window
	double PROT_DIAG_AVG_SCORE_THRESHOLD; // = 4.0; // BLOSUM62
	//#define PROT_DIAG_AVG_SCORE_THRESHOLD 5.0 // BLOSUM75


	/*******************************
	 * 
	 * GLOBAL QUALITY/SPEED SECTION
	 *
	 *****************************/

	// whether overlap weights are calculated or not
	int DO_OVERLAP; // = 0;

	// minimum diag length
	int DIAG_MIN_LENGTH;// = 1; 

	// diag threshold weight 
	double DIAG_THRESHOLD_WEIGHT;// = -log(0.5);

	// if there are anchors
	char DO_ANCHOR; // = 0;   

	// name of optional anchor file
	char *ANCHOR_FILE_NAME; // = NULL;   

	// sensitivity mode, does not set 
        // DIAG_THRESHOLD_WEIGHT to 0 four rounds >1
	char SENS_MODE;// = 0;

	// fast mode - behaves like dialign t 0.2.2
	char FAST_MODE;// = 0;

	// 1: only use a sqrt(amount_of_seqs) stripe of neighbour sequences to 
	// calculate pairwise alignments
	// 0: all pairwise alignments will be calculated
	int FAST_PAIRWISE_ALIGNMENT;// = 0;

	 char *conf_dir ; 
 	 char *in_file ;
 	 char *out_file ;
     char COMPUTE_PROB;
	int DNA_PARAMETERS; // default Einstellungen für DNA holen
	int DNA_TRANSLATION; // Vergleich auf Proteinebene, input DNA
	int FIND_ORF; // Vergleich auf Proteinebene, input DNA, mit ORF Finder
	int ORF_FRAME; // Vergleich auf Proteinebene, input DNA, nur longest ORF wird aligned
	int OUTPUT; // für DNA = 1, für Aminosäuren = 0


	int STATE_ORPHANE;
	int STATE_INHERITED;
	
};

    /************************************************/
    /*                                              */
    /*              global variable                 */
    /*                                              */
    /************************************************/
extern struct parameters* para;



    /*********************************************/
    /*                                           */
    /*        functions from parameters.c        */
    /*                                           */
    /*********************************************/
//struct parameters*

//initialises the parameters with default values
void init_parameters();

// check, whether there are enough arguments and if options are used
void check_input(int length, char** arguments);

//error message
void wrong_input();

void parameters(int argc, char** argv);
void set_parameters_dna();

