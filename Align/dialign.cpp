#include "dialign.h"
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

using namespace std;

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

void free_seq_col(struct seq_col* scol)
{
	unsigned int length = scol->length;
	unsigned int i;
	for (i=0; i<length; i++)
		free((scol->seqs[i]).data);
	free(scol->seqs);
	free(scol);
}

struct seq_col* read_seqs(const vector<string>& amb_seqs)
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
void get_alignment_consensus(struct alignment *algn, string& consensus)
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
	unsigned int chars[6]; //store count for -,a,c,g,t,n
	for (j=0; j<(int)max; j++) {
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
}
