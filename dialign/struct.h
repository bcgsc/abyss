/**
 *
 * struct.h: Basic data structures
 *
 * Author:  A.R.Subramanian
 */


/**
 * score matrix (e.g. BLOSUM62)
 */
struct scr_matrix
{
  int length;   // number of amino acids
  int max_score; // maximum among all entries in the data array
  int *char2num; // resolves the character of an amino acid to its number
  int *num2char; // resolves the number of an amino acid to its character
  int *data;     // contains the matrix indexed by the number of the particular amino acid
  int *dist;   // number of pairs of amino acids (i,j) having score at equal to              // the  value of the index
  long double **raw_dist;
  double avg_sim_score;

};

/**
 * raw sequence
 */
struct seq 
{
  char *data;  // sequence data  
  char *name;  // name/description of the sequence
  int num;     // number of the sequence
  int length;  // length of sequence
  int max_seen;
  char *dna_num; // Numbers fo retranslation from Protein to DNA
  int orf_frame; // reading frame of the longest orf
  char crick_strand; // orf translation or sequence translation on crickstrand
};


/**
 * sequence collection
 */
struct seq_col 
{
  struct seq *seqs; // array of the sequences
  int avg_length;   // average length of sequences
  int length;       // number of sequences
};

/**
 * probability distribution of scores in diagonals
 */
struct prob_dist {
  struct scr_matrix *smatrix; // pointer to the associated score matrix
  long double **data;  // distribution of scores dist[i][j] contains the     
                       // probability of a diags of length i having score >=j
  double **log_data; // distribution of scores dist[i][j] contains the     
                  // -log(probability of a diags of length i having score >=j)
  //  long double *expect;      // score expectancy for each diaglen
  unsigned int max_dlen;     // maximum diaglength in dist
};


/**
 * part of a sequence (auxiliary data structure)
 */
struct seq_part {
  int num;         // a number that indicates a position in an array
  struct seq* sq;  // the pointer to the sequence
  int startpos;   // startpos in the sequence
  //int leftmargin;
  //int rightmargin;
};


/**
 * diagonal in the dotmatrix
 */
struct diag {
  struct seq_part seq_p1; // first sequence part
  struct seq_part seq_p2; // seconde sequence part
  unsigned int length;     // length of the diag
  long score;     // score of the diag
  long orig_score;     // orig score of the diag

  struct diag *pred_diag;  // predecessor diag for dynamic programming
  struct diag *col_pred_diag;  // col predecessor diag for dynamic programming
  int pool_pos;            // position in diag pool

  char meetsThreshold;	   // whether diag meets threshold

  // for vertex cover
  int degree;
  int max_degree;
  struct diag **neighbours;

  char anchor;              // if this is an anchor diag
  char marked;              // marking flag for arbitrary use
  char multi_dg;            // is >0 if this is a multi dg
  struct diag **multi_cont; // the contained dgs of this is a multi dg
  int multi_length;         // size of multi_cont

  //  char onlyOverThres;

  double weight;           // weight of the diag = -log(prob)
  double weight_sum;       // weight sum for dialign
  double weight_fac;       // weight factor
  double ov_weight;       // overlap weight
  double total_weight;       // total_weight = weight+o_weight
};

/**
 * collection of diag 
 */
struct simple_diag_col {
  unsigned int length;       // number of diags
  double total_weight;       // total weight
  double weight_fac;         // weight factor
  struct diag** data;        // the array of diags
};

/**
 * guide tree node
 */
struct gt_node {
  char isLeaf;    // whether it is leaf
  int *seq_num;    // the sequence numbers 
  int seq_num_length; // length of sequence numbers array
  struct gt_node *succ1;  // successor nodes
  struct gt_node *succ2;
};


/**
 * vertex cover node
struct vc_node {
  double weight;
  struct diag *dg;
  int degree;
  struct vc_node *adjacents;
}
 */


/**
 * collection of all diagonals sorted by the sequences 
 */ 
struct diag_col {
  int seq_amount;            // number of sequences involved
  struct simple_diag_col** diag_matrix; // diag_matrix[i +seq_amount *j] contains
                             // all diags found involving the sequences i and j

  double total_weight;       // total weight
  double average_weight;     // average_weight
  struct diag** diags; // all diags unordered
  unsigned int diag_amount; // number of diags found 
  struct gt_node *gt_root;
};


/**
 * diag container
 */
struct diag_cont {
  struct diag* dg;
  struct diag_cont *next;
};

/**
 * alignment position
 */
struct algn_pos {
  //  int seq_num;               // sequence number
  //  unsigned long pos_in_seq;  // position in the sequence
  char state;         // orphane: not aligned to any pos, 

  struct diag_cont *dg_cont; // diags that are aligned with that position

  int row;
  int col;
  //  unsigned int succFPos;  // if orphane, the position holding the succF
  // unsigned int predFPos;  // analogous to succFPos
  //  char* isAli;            // if alignemnt at the positions exist
  int* predF;    // predecessor frontier, only filled if non-orphane
  int* succF;    // successor frontier, only filled if non-orphane

  char *proceed;    // for output
  //char isInherited; // whether the pointers are inherited
  int predFPos;  // in case of orphane, where to find the predF or succF
  int succFPos;

  int *eqcAlgnPos; // equivalence class minimum alignment position (for output)
  struct algn_pos *eqcParent; // equivalence class parent
  int eqcRank; // equivalence class rank (>= maximum number of children)
  //  unsigned int *maxpos; // needed for output of the alignment 
};

/**
 * alignment
 */
struct alignment {
  //int seq_amount;            // number of sequences involved
  //char *redo_seqs;           // which pairs of sequences are to be aligned again
  char *seq_is_orphane;      // boolean array indicating for each sequence 
                             // whether it is orphane or not
  int max_pos;               // the greatest position in the alignment (including all -'s)
                             // if <0: the alignment has not yet been prepared
  struct seq_col *scol;      // all the sequences involved
  struct algn_pos **algn;    // the alignment
  double total_weight;       // the total weight of the alignment
  struct alignment *next;    // pointer to next alignment in the sorted linked list
  //struct alignment *prev;    // pointer to previous alignment in the sorted linked list
  //unsigned long pos;         // position in the sorted linked list

  //struct diag** aligned_diags; // all aligned diags
  //int aligned_diags_amount;
  //int max_aligned_diags_amount;
  //int orig_max_aligned_diags_amount;
  
  //struct diag_cont* backlog_diags; // all backlog diags
};
