#ifndef DIALIGN_H
#define DIALIGN_H 1

extern "C" {
#include "dialign/io.h"
#include "dialign/parameters.h"
#include "dialign/struct.h"

extern struct scr_matrix *smatrix;
extern struct prob_dist *pdist;
extern const double dna_diag_prob_100_exp_550000[5151];

struct alignment* create_empty_alignment(struct seq_col *scol);
struct diag_col *find_all_diags(struct scr_matrix *smatrix,
	struct prob_dist *pdist,
	struct seq_col *in_seq_col, struct alignment *algn, int round);
struct alignment* guided_aligner(struct alignment *palgn,
	struct seq_col *scol, struct diag_col *dcol,
	struct scr_matrix* smatrix,
	struct prob_dist *pdist,
	struct gt_node *gtn,
	int round);
char simple_aligner(struct seq_col *scol, struct diag_col *dcol,
	struct scr_matrix* smatrix,
	struct prob_dist *pdist,
	struct alignment *algn, int round);
void prepare_alignment(struct alignment *algn);
struct algn_pos *find_eqc(struct algn_pos **ap, int seqnum, int pos);
void free_alignment(struct alignment* algn);
void free_diag_col(struct diag_col* dcol);
}

#include <string>
#include <vector>

void initDialign();
void free_prob_dist(struct prob_dist* pdist);
std::string dialign(const std::vector<std::string>& amb_seqs,
		std::string& alignment, unsigned& matches);

#endif
