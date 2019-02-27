/**
 *
 * parameters.c:  Read from stdin 
 *
 *             
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <unistd.h>
#include "parameters.h"
#include "struct.h"
#include "io.h"
#ifdef __unix__
	#include <sys/stat.h>
	#include <sys/types.h>
#else
	#include <sys/stat.h>
#endif


extern char *optarg;
extern int optind, opterr, optopt;
struct parameters* para;
/****************************
* PROTEIN DEFAULT VALUES!   *
****************************/
void init_parameters()
{
	if((para =(struct parameters *) malloc(sizeof(struct parameters)) ) == NULL) 
	{
		error("init_parameters(): Out of memory when allocating data !");
	}
	para->VERSION = "1.0.2";
	para->DEBUG = 0;
	para->MAX_SEQ_AMOUNT = 5000;
	para->MAX_FASTA_LINE_LENGTH = 100;
	para->PRINT_SEQ_LINE_LENGTH = 80;
	para->SCR_MATRIX_FILE_NAME="BLOSUM.scr";
	para->DIAG_CALC_WEIGHT_THRESHOLD = 0.000000065;
	para->DIAG_PROB_FILE_NAME="BLOSUM.diag_prob_t10";
	para->SCR_MATRIX_ADD = 0; 
	para->PROT_SIM_SCORE_THRESHOLD = 4.0;
	//para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS = 4;
	//para->PROT_DIAG_MIN_LENGTH_THRESHOLD = 10.0;
	para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS = 4;
	para->PROT_DIAG_MIN_LENGTH_THRESHOLD = 40.0;
	para->PROT_DIAG_AVG_SCORE_THRESHOLD = 4.0;
	para->DO_ANCHOR = 0;
	para->ANCHOR_FILE_NAME = NULL;
	para->DO_OVERLAP = 0;
	para->DIAG_MIN_LENGTH = 1;
	para->FAST_MODE = 0;
	para->SENS_MODE = 0;
	para->DIAG_THRESHOLD_WEIGHT = -log(0.5);
	para->FAST_PAIRWISE_ALIGNMENT = 0;
	para->conf_dir = NULL;
	para->in_file = NULL;
	para->out_file = NULL;
 	para->COMPUTE_PROB=0;
	para->STATE_ORPHANE=1;
	para->STATE_INHERITED=2;
	para->DNA_PARAMETERS = 0;
	para->DNA_TRANSLATION = 0;
	para->FIND_ORF = 0;
	para->ORF_FRAME = 0;
	para->OUTPUT = 0;
}

/****************************
*     DNA DEFAULT VALUES!   *
****************************/
void set_parameters_dna()
{
	para->VERSION = "1.0.2";
	para->DEBUG = 0;
	para->MAX_SEQ_AMOUNT = 5000;
	para->MAX_FASTA_LINE_LENGTH = 100;
	para->PRINT_SEQ_LINE_LENGTH = 80;
	para->SCR_MATRIX_FILE_NAME="dna_matrix.scr";
	para->DIAG_CALC_WEIGHT_THRESHOLD = 0.000000065;
	para->DIAG_PROB_FILE_NAME="dna_diag_prob_100_exp_550000";
	para->SCR_MATRIX_ADD = 0; 
	para->PROT_SIM_SCORE_THRESHOLD = 0.25; 
	para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS = 4; // 1
	para->PROT_DIAG_MIN_LENGTH_THRESHOLD = 40.0;//40.0; // 1 
	para->PROT_DIAG_AVG_SCORE_THRESHOLD = 0.25; //1.0 
	para->DO_ANCHOR = 0;
	para->ANCHOR_FILE_NAME = NULL;
	para->DO_OVERLAP = 0;
	para->DIAG_MIN_LENGTH = 1;
	para->FAST_MODE = 0;
	para->SENS_MODE = 0;
	para->DIAG_THRESHOLD_WEIGHT = -log(0.5);//-log(0.875);
	para->FAST_PAIRWISE_ALIGNMENT = 0; 
	para->conf_dir = NULL;
	para->in_file = NULL;
	para->out_file = NULL;
 	para->COMPUTE_PROB=0;
	para->STATE_ORPHANE=1;
	para->STATE_INHERITED=2;
	para->DNA_PARAMETERS = 1;
	para->DNA_TRANSLATION = 0;
	para->FIND_ORF = 0;
	para->ORF_FRAME = 0;
	para->OUTPUT = 1;
}

void check_input(int argc, char** argv)
{
	int flag_protein_output = 0;
	opterr = 1;
	optind = 0;
	int opt, only_config_dir = 0;

	while((opt = getopt(argc, argv,"PDTLOCFHhA:d:s:a:c:l:m:w:p:v:t:n:g:o:r:u"))!= -1){
	  switch(opt){
		case 'd':
			para->DEBUG = atoi(optarg);
			break;
		case 's':
			para->MAX_SEQ_AMOUNT = atoi(optarg);
			break;
		case 'a':
			para->MAX_FASTA_LINE_LENGTH = atoi(optarg);	
			break;
		case 'c':
			para->PRINT_SEQ_LINE_LENGTH = atoi(optarg);
			break;
		case 'l':
		        para->SENS_MODE = atoi(optarg);
			break;
		case 'A':
			para->ANCHOR_FILE_NAME = optarg;
			para->DO_ANCHOR = 1;
			break;
		case 'm':
			para->SCR_MATRIX_FILE_NAME = optarg;
			break;
		case 'w':
			para->DIAG_CALC_WEIGHT_THRESHOLD = atof(optarg);
			break;
		case 'p':
			para->DIAG_PROB_FILE_NAME = optarg;
			break;
		case 'v':
			para->SCR_MATRIX_ADD = atoi(optarg);
			break;
		case 't':
			para->PROT_SIM_SCORE_THRESHOLD = atoi(optarg);
			break;
		case 'n':
			para->PROT_DIAG_MAX_UNDER_THRESHOLD_POS = atoi(optarg);
			break;
		case 'g':
			para->PROT_DIAG_MIN_LENGTH_THRESHOLD = atof(optarg);
			break;
		case 'r':
			para->PROT_DIAG_AVG_SCORE_THRESHOLD = atof(optarg);
			break;
		case 'o':
			para->DO_OVERLAP = atoi(optarg);
			break;
		case 'u':
			para->DIAG_MIN_LENGTH = atoi(optarg);
			break;
		case 'h':
			wrong_input();
		case 'H':
			wrong_input();
		case 'C':
			only_config_dir = 1;
			para->COMPUTE_PROB = 1;
			break;
		case 'F':
  		        para->DIAG_THRESHOLD_WEIGHT = 0.0;
		        para->FAST_MODE = 1;
			break;
		case 'O':
		        para->DNA_PARAMETERS = 0;
			para->DNA_TRANSLATION = 1;
			para->FIND_ORF = 1;
			para->ORF_FRAME = 1;
			para->OUTPUT = 1;
			break;
		case 'L':
			para->DNA_PARAMETERS = 0;
			para->DNA_TRANSLATION = 1;
			para->FIND_ORF = 1;
			para->ORF_FRAME = 0;
			para->OUTPUT = 1;
			break;
		case 'T':
			para->DNA_PARAMETERS = 0;
			para->DNA_TRANSLATION = 1;
			para->FIND_ORF = 0;
			para->ORF_FRAME = 0;
			para->OUTPUT = 1;
			break;
		case 'D':
			break;
		case 'P':
			flag_protein_output = 1;
			break;
		case '?':
            break;
		}
	}
	
	if(flag_protein_output){
		para->OUTPUT = 0;
	}

	if(argc-optind == 0 ){
		wrong_input();
		error("conf-directory and infile needed !");
	}
	else if(argc-optind > 3){
		wrong_input();
		error("too many arguments -> conf-directory, infile, [outfile] !");
	}
    else{
		struct stat attribut_dir;
    	if(only_config_dir){
			if(stat(argv[optind],&attribut_dir)==-1){
				wrong_input();
				error("conf-directory doesn't exist");
			}
			else {
				if(attribut_dir.st_mode & S_IFDIR){
					para->conf_dir = argv[optind++];
				}
				else {
					wrong_input();
					error("conf-directory is no directory");
				}
			}
		}
		else{
			struct stat attribut_file;
			if(stat(argv[optind],&attribut_dir)==-1){
				wrong_input();
				error("conf-directory doesn't exist");
			}
			else{
				if(attribut_dir.st_mode & S_IFDIR){
					para->conf_dir = argv[optind++];
					if(stat(argv[optind],&attribut_file)==-1){
						wrong_input();
						error("infile doesn't exist");
					}
					else{
						if(attribut_file.st_mode & S_IFREG){
							para->in_file = argv[optind++];
						}
						else{
							wrong_input();
							error("infile isn't a regular file!");
						}
					}
				}
				else{
					wrong_input();
					error("conf-directory is no directory");
				}
			}
			if(argv-optind >0)
				para->out_file = argv[optind];
		}
	}
	return;
}

void wrong_input()
{
	printf("Usage: dialign-t [OPTIONS] <conf-directory> <fasta-file> [<fasta-out-file>]\n");
	printf("\n  -d\tDebug-Mode  [DEFAULT 0]\n");
	printf("   \t\t 0 no debug statements\n");
	printf("   \t\t 1 debugs the current phase of the processing\n");
	printf("   \t\t 2 very loquacious debugging\n");
	printf("   \t\t 5 hardcore debugging\n");
	printf("  -s\tmaximum amount of input sequences [DEFAULT 5000]\n");
	printf("  -a\tmaximum number of characters per line in a FASTA file [DEFAULT 100]\n");
	printf("  -c\tmaximum amount of characters per line when printing a sequence\n     \t[DEFAULT 80]\n");
	printf("  -l\tsensitivity mode, the higher the level the less likely\n");
	printf("    \tspurious random fragments are aligned in local alignments \n     \t[DEFAULT 0]\n");
	printf("   \t\t 0 switched off \n");
	printf("   \t\t 1 level-1, reduced sensitivity\n");
	printf("   \t\t 2 level-2, strongly reduced sensitivity\n");
	printf("  -m\tscore matrix file name (in the configuration directory)\n     \t\t[DEFAULT PROTEIN: BLOSUM.scr]\n \t\t[DEFAULT DNA: dna_matrix.scr]\n");
	printf("  -w\tdefines the minimum weight when the weight formula is changed\n     \tto 1-pow(1-prob, factor) [DEFAULT 0.000000065]\n");
	printf("  -p\tprobability distribution file name (in the configuration\n     \tdirectory) \n \t\t[DEFAULT PROTEIN: BLOSUM.diag_prob_t10]\n\t\t[DEFAULT DNA: dna_diag_prob_100_exp_550000]\n");
	printf("  -v\tadd to each score (to prevent negative values) [DEFAULT 0]\n");
	printf("  -t\t\"even\" threshold for low score for sequences alignment \n \t\t[DEFAULT PROTEIN: 4]\n\t\t[DEFAULT DNA: 0]\n");
	printf("  -n\tmaximum number of consecutive positions for window containing\n     \tlow scoring positions \n \t\t[DEFAULT PROTEIN: 4]\n\t\t[DEFAULT DNA: 4]\n");
	printf("  -g\tglobal minimum fragment length for stop criterion \n \t\t[DEFAULT PROTEIN: 40] \n\t\t[DEFAULT DNA: 40]\n");
	printf("  -m\tminimal allowed average score in frag window containing low \n     \tscoring positions \n \t\t[DEFAULT PROTEIN: 4.0]\n\t\t[DEFAULT DNA: 0.25]\n");
	printf("  -o\twhether overlap weights are calculated or not [DEFAULT 0]\n");
	printf("  -f\tminimum fragment length [DEFAULT 1]\n");
	printf("  -r\tthreshold weight to consider the fragment at all [DEFAULT 0.0]\n");
	printf("  -u\t[DEFAULT 0]\n");
	printf("    \t\t1: only use a sqrt(amount_of_seqs) stripe of neighbour\n     \t\t   sequences to calculate pairwise alignments (increase performance)\n");
	printf("    \t\t0: all pairwise alignments will be calculated\n");
	printf("  -A\toptional anchor file [DEFAULT none]\n");

	printf("  -D\tinput is DNA-sequence\n");
	printf("  -T\ttranslate DNA into aminoacids from begin to end (length will be cut to mod 3 = 0)\n\tWARNING: Do not use -D with this option \n\t(Default values for PROTEIN input will be loaded)\n");
 	printf("  -L\tcompare only longest Open Reading Frame\n\tWARNING: Do not use -D with this option \n\t(Default values for PROTEIN input will be loaded)\n");
 	printf("  -O\ttranslate DNA to aminoacids, reading frame for each sequence calculated due to its longest ORF\n\tWARNING: Do not use -D with this option \n\t(Default values for PROTEIN input will be loaded)\n");
	printf("  -P\toutput in aminoacids, no retranslation of DNA sequences\n\t[DEFAULT: input = output]\n");
 	printf("  -F\tfast mode (implies -l0, since it already significantly reduces sensitivity)\n");
 	printf("  -C\tgenerate probability table saved in <config_dir>/prob_table and exit\n");
	printf("  -H -h\tprint this message\n\n");
	exit(1);
}

void parameters(int argc, char** argv)
{
	init_parameters();	
	int opt_flag = 0;
	int opt;
	opterr = 0;
	while((opt = getopt(argc, argv, "D") )!=-1){
		switch(opt){
			case 'D': opt_flag = 1; break;
			case '?': break;
		}
	}
					
	if(opt_flag)
		set_parameters_dna();
	
	check_input(argc, argv);
}


