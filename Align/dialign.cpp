#include "dialign.h"
#include "Common/Options.h"
#include <cassert>
#include <cmath> // for log
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime> // for clock
#include <iostream>

using namespace std;

/** Score matrix. */
scr_matrix* smatrix;

/** Diagonal length probability distribution. */
prob_dist* pdist;

static scr_matrix* newDefaultScoreMatrix()
{
	string s("ACGT?#$");
	struct scr_matrix* p = (scr_matrix*)calloc(1, sizeof *smatrix);
	p->length = s.size();
	p->num2char = (int*)calloc(256, sizeof(int));
	p->char2num = (int*)calloc(256, sizeof(int));
	for (unsigned i = 0; i < s.size(); ++i) {
		unsigned c = s[i];
		p->num2char[i] = c;
		p->char2num[c] = i;
	}

	p->data = (int*)calloc(s.size() * s.size(), sizeof(int));
	unsigned n = s.size() - 3; // ignore ?#$
	// Set the diagonal to 1.
	for (unsigned i = 0; i < n; i++)
		p->data[s.size() * i + i] = 1;
	p->max_score = 1;
	p->avg_sim_score = para->PROT_SIM_SCORE_THRESHOLD;
	p->dist = (int*)calloc(2, sizeof(int));
	p->dist[0] = n * n - n;
	p->dist[1] = n;
	return p;
}

/** Initialize dialign. */
void initDialign()
{
	smatrix = strlen(para->SCR_MATRIX_FILE_NAME) > 0
		? read_scr_matrix(para->SCR_MATRIX_FILE_NAME)
		: newDefaultScoreMatrix();

	// print the score matrix
	if (para->DEBUG >5)
		print_scr_matrix(smatrix);

	// read the probability distribution for diagonals
	pdist = read_diag_prob_dist(smatrix, para->DIAG_PROB_FILE_NAME);
}

static void free_scr_matrix(struct scr_matrix* smatrix)
{
	free(smatrix->dist);
	free(smatrix->data);
	free(smatrix->char2num);
	free(smatrix->num2char);
	free(smatrix);
}

void free_prob_dist(struct prob_dist* pdist)
{
	unsigned int length = pdist->max_dlen;
	unsigned int i;
	for (i=1; i<=length; i++) {
		free(pdist->data[i]);
		free(pdist->log_data[i]);
	}
	free(pdist->data);
	free(pdist->log_data);
	free_scr_matrix(pdist->smatrix);
	free(pdist);
}

static void free_seq_col(struct seq_col* scol)
{
	unsigned int length = scol->length;
	unsigned int i;
	for (i=0; i<length; i++)
		free((scol->seqs[i]).data);
	free(scol->seqs);
	free(scol);
}

/** Print a dialign alignment. */
static ostream& print(ostream& out, const alignment& o,
		const string& consensus)
{
	const seq_col& scol = *o.scol;
	vector<int> proc(scol.length);
	algn_pos **ap = o.algn;
	for (int s = 0; s < scol.length; s++) {
		const seq& sq = scol.seqs[s];
		for (int j = 0; j < o.max_pos; j++) {
			if (proc[s] < sq.length) {
				const algn_pos& ap1 = *find_eqc(ap, s, proc[s]);
				assert(j <= *ap1.eqcAlgnPos);
				if (*ap1.eqcAlgnPos == j) {
					char c = sq.data[proc[s]];
					if (toupper(c) == toupper(consensus[j]))
						out << '.';
					else if (ap1.state & para->STATE_ORPHANE)
						out << (char)tolower(c);
					else
						out << c;
					proc[s]++;
				} else
					out << '*';
			} else
				out << '*';
		}
		out << '\n';
	}
	return out;
}

static struct seq_col* read_seqs(const vector<string>& amb_seqs)
{
	struct seq_col* scol = (struct seq_col*)calloc(1, sizeof(struct seq_col));
	struct seq* seqs = (scol->seqs = (struct seq*)calloc(amb_seqs.size(), sizeof(struct seq)));
	if(scol==NULL || seqs==NULL) {
		cerr << "read_seqs(): Out of memory !\n";
		exit(EXIT_FAILURE);
	}
	scol->length = amb_seqs.size();
	scol->avg_length = 0;

	seq* seq;
	for (size_t i=0; i<amb_seqs.size(); i++) {
		assert(!amb_seqs[i].empty());
		seq = &(scol->seqs[i]);
		seq->max_seen = 0;
		//seq->name = calloc(rlen, sizeof(char)); //do I need this?
		seq->num = i;
		seq->orf_frame=0;
		seq->crick_strand=0;
		//strncpy(seq->name, &(rline[1]), rlen-2);
		seq->data = (char*)calloc(amb_seqs[i].length()+1, sizeof(char));
		if (seq->data == NULL) {
			cerr << "seq->data out of memory !\n";
			exit(EXIT_FAILURE);
		}
		strcpy(seq->data, amb_seqs[i].c_str());
		seq->length = amb_seqs[i].length();
		scol->avg_length += amb_seqs[i].length();
		if(para->DEBUG >1) printf("DEBUG: seq:%s\n", seq->data);
	}
	scol->avg_length /= scol->length;
	if(para->DEBUG >1) printf("DEBUG: total # of amb_seqs: %i, avg_length: %i\n", scol->length, scol->avg_length);
	return scol;
}

/** Return the consensus base.
 * @param counts a count of the characters "-ACGTN"
 */
static char make_consensus(unsigned counts[6])
{
	static const char IUPAC[16] = {
		'N', //----
		'A', //---A
		'C', //--C-
		'M', //--CA
		'G', //-G--
		'R', //-G-A
		'S', //-GC-
		'V', //-GCA
		'T', //T---
		'W', //T--A
		'Y', //T-C-
		'H', //T-CA
		'K', //TG--
		'D', //TG-A
		'B', //TGC-
		'N', //TGCA
	};
	unsigned bases = 0;
	for (unsigned i = 1, mask = 1; i < 5; i++, mask <<= 1)
		if (counts[i] > 0)
			bases |= mask;
	if (bases == 0)
		assert(counts[5] > 0);
	assert(bases < 16);
	char c = IUPAC[bases];
	return counts[0] > 0 ? tolower(c) : c;
}

// assume initial sequences contain only a/c/g/t/n
static string get_alignment_consensus(struct alignment *algn)
{
  struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  int j;
  unsigned int s,max;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;

  prepare_alignment(algn);
  max = algn->max_pos;
  if (para->DEBUG > 5) printf("slen is %u, max pos is %u\n", slen, max);
  struct algn_pos *ap1;

	max = algn->max_pos;
	int* proc = new int[slen];
	for (j=0; j<(int)slen; j++)
		proc[j] = 0;
	string consensus;
	for (j=0; j<(int)max; j++) {
		unsigned int chars[6]; //store count for -,a,c,g,t,n
		for (s=0; s<6; s++)
			chars[s]=0;
		for(s=0;s<slen;s++) {
			sq = &(scol->seqs[s]);
			if(proc[s] < sq->length) {
				ap1 = find_eqc(ap,s,proc[s]);
				if(*ap1->eqcAlgnPos==j) {
					char cur_char = toupper(sq->data[proc[s]]);
					switch (cur_char) {
					case 'A':
						chars[1]++; break;
					case 'C':
						chars[2]++; break;
					case 'G':
						chars[3]++; break;
					case 'T':
						chars[4]++; break;
					case 'N':
						chars[5]++; break;
					default:
						cerr << "error: unexpected character: `"
							<< cur_char << "'\n";
						assert(false);
						exit(EXIT_FAILURE);
					}
					proc[s]++;
				} else {
					chars[0]++;
				}
			} else {
				chars[0]++;
			}
		}
		char c = make_consensus(chars);
		consensus += c;
		if (para->DEBUG > 5) {
			for (s=0; s<6; s++) printf("chars[%u]:%u; ", s, chars[s]);
			printf("\nconsensus: %c\n", c);
		}
	}
	delete[] proc;
	return consensus;
}

/** Align multiple sequences using DIALIGN-TX. */
string dialign(const vector<string>& amb_seqs)
{
	int i;
	struct seq_col *in_seq_col = NULL;
	double tim = clock();

	in_seq_col = read_seqs(amb_seqs);

	// fast mode has higher threshold weights
	struct parameters *dialign_para = para;
	if(dialign_para->FAST_MODE)
		dialign_para->PROT_SIM_SCORE_THRESHOLD += 0.25;

	// Consider Anchors -> default for DNA: DO_ANCHOR = 0;
	struct alignment *algn = NULL;
	if (!dialign_para->FAST_MODE)
		algn = create_empty_alignment(in_seq_col);
	struct alignment *salgn = create_empty_alignment(in_seq_col);
	if (dialign_para->DEBUG > 1)
		printf("empty alignments created\n");

	// Compute pairwise diagonals
	struct diag_col *all_diags = find_all_diags(smatrix, pdist,
		in_seq_col, salgn, 1);
	double duration = (clock()-tim)/CLOCKS_PER_SEC;
	if (dialign_para->DEBUG > 1)
		printf("Found %i diags in %f secs\n",
			all_diags->diag_amount, duration);
	int diag_amount = all_diags->diag_amount;

	// Compute alignment
	double tim2 = clock();
	if (!dialign_para->FAST_MODE) {
		struct diag *cp_diags[all_diags->diag_amount];
		for(i = 0; i < diag_amount; i++) {
			cp_diags[i] = (diag*)malloc(sizeof(struct diag));
			*(cp_diags[i]) = *(all_diags->diags[i]);
		}
		guided_aligner(algn, in_seq_col, all_diags, smatrix,
			pdist, all_diags->gt_root, 1);

		for(i = 0; i < diag_amount; i++)
			all_diags->diags[i] = cp_diags[i];

		all_diags->diag_amount = diag_amount;
	}
	simple_aligner(in_seq_col, all_diags, smatrix, pdist,
		salgn, 1);
	duration = (clock()-tim2)/CLOCKS_PER_SEC;

	if (!dialign_para->FAST_MODE) {
		if (dialign_para->DEBUG > 1)
			printf("First alignment after %f secs. "
					"simple: %f guided: %f\n",
				duration, salgn->total_weight, algn->total_weight);
		else
			if (dialign_para->DEBUG > 1)
				printf("First alignment after %f secs. simple: %f \n",
					duration, salgn->total_weight);
	}

	free_diag_col(all_diags);

	dialign_para->DO_ANCHOR = 0; // anchors done

	// round 2+
	int round;
	char newFound = 0;
	int type;

	// consider sensitivity level
	if (!dialign_para->FAST_MODE) {
		if (dialign_para->SENS_MODE == 0) {
			dialign_para->DIAG_THRESHOLD_WEIGHT = 0.0;
		} else if (dialign_para->SENS_MODE == 1) {
			dialign_para->DIAG_THRESHOLD_WEIGHT
				= -log(0.75);//-log(.875+0.125/2.0);
		} else if (dialign_para->SENS_MODE == 2) {
			dialign_para->DIAG_THRESHOLD_WEIGHT
				= -log(0.5);//-log(0.875);
		}
	}

	int stype = (dialign_para->FAST_MODE ? 1 : 0);
	for (type = stype; type < 2; type++) {
		for (round = 2; round <= 20; round++) {
			tim2 = clock();
			all_diags = find_all_diags(smatrix, pdist,
				in_seq_col, (type ? salgn : algn), round);
			duration = (clock()-tim2)/CLOCKS_PER_SEC;
			if (dialign_para->DEBUG > 1)
				printf("Found %i diags after %f secs\n",
					all_diags->diag_amount, duration);
			if (all_diags->diag_amount == 0) {
				free_diag_col(all_diags);
				break;
			} else {
			// round 2 and further we use the simple aligner
				newFound = simple_aligner(in_seq_col,
					all_diags, smatrix, pdist,
					(type ? salgn : algn), round);
				free_diag_col(all_diags);
				if (!newFound)
					break;
			}
		}
	}
	if (dialign_para->DEBUG > 1)
		printf("Alignment ready!\n");

	if (!dialign_para->FAST_MODE) {
		if (dialign_para->DEBUG > 1)
			printf("Final alignment simple: %f guided: %f\n",
				salgn->total_weight, algn->total_weight);
	} else {
		if (dialign_para->DEBUG > 1)
			printf("Final alignment simple: %f \n",
				salgn->total_weight);
	}

	if (dialign_para->FAST_MODE
			|| salgn->total_weight > algn->total_weight) {
		if (!dialign_para->FAST_MODE)
			free_alignment(algn);
		algn = salgn;
	} else {
		free_alignment(salgn);
	}

	if (opt::verbose > 2)
		simple_print_alignment_default(algn);
	string consensus = get_alignment_consensus(algn);
	if (opt::verbose > 0)
		print(cerr, *algn, consensus);

	if (dialign_para->DEBUG > 0) {
		duration = (clock()-tim)/CLOCKS_PER_SEC;
		cerr << "Total time: " << duration << " s\n"
			"Total weight: " << algn->total_weight << '\n';
	}

	free_alignment(algn);
	free_seq_col(in_seq_col);
	return consensus;
}
