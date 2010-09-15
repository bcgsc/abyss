#ifndef DIALIGN_H
#define DIALIGN_H 1

extern "C" {
#include "dialign/io.h"
#include "dialign/parameters.h"
#include "dialign/struct.h"

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
void free_alignment(struct alignment* algn);
void free_diag_col(struct diag_col* dcol);
}

#include <string>
#include <vector>

void free_prob_dist(struct prob_dist* pdist);
void free_seq_col(struct seq_col* scol);
void get_alignment_consensus(struct alignment* algn, std::string& consensus);
struct seq_col* read_seqs(const std::vector<std::string>& amb_seqs);

#endif
