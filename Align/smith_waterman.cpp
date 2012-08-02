//////////////////////////////////////////////////////////////////////////////
// Smith_waterman implementation, with modification, to find sequence overlap
//
// Written by Rong She (rshe@bcgsc.ca)
// Last modified: Jul 7, 2010
//////////////////////////////////////////////////////////////////////////////

#include "smith_waterman.h"
#include "Sequence.h"
#include "Align/Options.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat> // for DBL_MAX
#include <iostream>

using namespace std;

namespace opt {
	/** The score of a match. */
	int match = 5;

	/** The score of a mismatch. */
	int mismatch = -4;

	/** gap open penalty */
	int gap_open = -12;

	/** gap extend penalty */
	int gap_extend = -4;
}
/** Print the specified alignment. */
static ostream& printAlignment(ostream& out,
		const string& aseq, const string& bseq,
		unsigned alignPos[], SMAlignment align)
{
	unsigned astart = alignPos[0], aend = alignPos[1] + 1;
	unsigned bstart = alignPos[2], bend = alignPos[3] + 1;
	assert(aend == aseq.size());
	assert(bstart == 0);
	(void)aend; (void)bstart;
	out << aseq.substr(0, astart) << align.query_align << '\n'
		<< string(astart, ' ');
	for (unsigned i = 0; i < align.match_align.size(); i++) {
		char a = align.query_align[i],
			b = align.target_align[i],
			c = align.match_align[i];
		out << (c == a || c == b ? '.' : c);
	}
	return out << '\n' << string(astart, ' ')
		<< align.target_align << bseq.substr(bend) << '\n';
}

/** Index comparison functor. */
template<class T>
struct index_cmp {
	const T arr;
	index_cmp(const T arr) : arr(arr) {}
	bool operator()(const int a, const int b) const { return arr[a] > arr[b]; }
};

/** Return whether the characters a and b match.
 * @param c [out] the consensus character
 */
static bool isMatch(char a, char b, char& c)
{
	if (a == b) {
		c = a;
	} else if (toupper(a) == toupper(b)) {
		c = islower(a) || islower(b) ? tolower(a) : a;
	} else if (a == 'N' || a == 'n') {
		c = b;
	} else if (b == 'N' || b == 'n') {
		c = a;
	} else {
		c = ambiguityOr(a, b);
		return ambiguityIsSubset(a, b);
	}
	return true;
}

/** Return the score of the alignment of a and b. */
static int matchScore(const char a, const char b)
{
	char consensus;
	return isMatch(a, b, consensus) ? opt::match : opt::mismatch;
}

/** Return the score of a gap, either newly opened or extended. */
static int gapScore(bool prev_is_gap)
{
	return prev_is_gap ? opt::gap_extend : opt::gap_open;
}

//the backtrack step in smith_waterman
unsigned Backtrack(const int i_max, const int j_max, int** I_i, int** I_j,
		const string& seq_a, const string& seq_b, SMAlignment& align, unsigned* align_pos)
{
	// Backtracking from H_max
	int current_i=i_max,current_j=j_max;
	int next_i=I_i[current_i][current_j];
	int next_j=I_j[current_i][current_j];
	string consensus_a(""), consensus_b(""), match("");
	unsigned num_of_match = 0;
	while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){
		if(next_i==current_i) {
			consensus_a += '-'; //deletion in A
			match += tolower(seq_b[current_j-1]);
			consensus_b += seq_b[current_j-1]; //b must be some actual char, cannot be '-' aligns with '-'!
		}
		else {
			consensus_a += seq_a[current_i-1]; // match/mismatch in A
			if(next_j==current_j) {
				consensus_b += '-'; // deletion in B
				match += tolower(seq_a[current_i-1]);
			}
			else {
				consensus_b += seq_b[current_j-1]; // match/mismatch in B
				char consensus_char;
				if (isMatch(seq_a[current_i-1], seq_b[current_j-1],
							consensus_char)) {
					match += consensus_char;
					num_of_match++;
				}
				else {
					match += ambiguityOr(
							seq_a[current_i-1], seq_b[current_j-1]);
				}
			}
		}

		current_i = next_i;
		current_j = next_j;
		next_i=I_i[current_i][current_j];
		next_j=I_j[current_i][current_j];
	}

	//check whether the alignment is what we want (pinned at the ends), modified version of SW (i_max is already fixed)
	if (current_j > 1)
		return 0;

	//record the last one
	consensus_a += seq_a[current_i-1];
	consensus_b += seq_b[current_j-1];
	char consensus_char;
	if (isMatch(seq_a[current_i-1], seq_b[current_j-1],
				consensus_char)) {
		match += consensus_char;
		num_of_match++;
	}
	else {
		match += ambiguityOr(seq_a[current_i-1], seq_b[current_j-1]);
	}
	align_pos[0] = current_i-1;
	align_pos[2] = current_j-1;

	align_pos[1] = i_max-1;
	align_pos[3] = j_max-1;

	reverse(consensus_a.begin(), consensus_a.end());
	reverse(consensus_b.begin(), consensus_b.end());
	reverse(match.begin(), match.end());
	align.query_align = consensus_a;
	align.match_align = match;
	align.target_align = consensus_b;

	return num_of_match;
}

/* This is the Smith-Waterman algorithm (a variation of
 * Needleman-Wunsch algorithm), finds one optimal local alignment
 * Modified to find overlap of seq_a and seq_b (alignment that is
 * pinned at the end of seq_a and beginning of seq_b).
 * Actually, this should be a variation of needleman algorithm, that
 * looks for a global alignment, but without penalizing overhangs...
 * and make sure the alignment is end-to-end (end of seqA to beginning
 * of seqB).
 */
void alignOverlap(const string& seq_a, const string& seq_b, unsigned seq_a_start_pos,
	vector<overlap_align>& overlaps, bool multi_align, bool verbose)
{
	// get the actual lengths of the sequences
	int N_a = seq_a.length();
	int N_b = seq_b.length();

	// initialize H
	int i, j;
	double** H;
	int **I_i, **I_j;
	H = new double*[N_a+1];
	I_i = new int*[N_a+1];
	I_j = new int*[N_a+1];
	bool** V = new bool*[N_a+1];

	for(i=0;i<=N_a;i++){
		H[i] = new double[N_b+1];
		I_i[i] = new int[N_b+1];
		I_j[i] = new int[N_b+1];
		H[i][0]=0; //only need to initialize first row and first column
		I_i[i][0] = i-1;
		V[i] = new bool[N_b+1];
		V[i][0] = true; //valid start
	}

	for (j = 0; j <= N_b; j++) {
		H[0][j] = 0; //initialize first column
		I_j[0][j] = j-1;
		V[0][j] = false; //wrong start, not overlap
	}
	V[0][0] = true;

	for(i=1;i<=N_a;i++){
		for(j=1;j<=N_b;j++){
			char a = seq_a[i-1], b = seq_b[j-1];
			double scores[3] = {
				V[i-1][j-1] ? H[i-1][j-1] + matchScore(a, b)
					: -DBL_MAX, // match or mismatch
				V[i-1][j] ? H[i-1][j] + gapScore(I_j[i-1][j] == j)
					: -DBL_MAX, // deletion in sequence A
				V[i][j-1] ? H[i][j-1] + gapScore(I_i[i][j-1] == i)
					: -DBL_MAX // deletion in sequence B
			};
			double* pMax = max_element(scores, scores + 3);
			H[i][j] = *pMax;
			switch (pMax - scores) {
			  case 0: // match or mismatch
				I_i[i][j] = i-1;
				I_j[i][j] = j-1;
				break;
			  case 1: // deletion in sequence A
				I_i[i][j] = i-1;
				I_j[i][j] = j;
				break;
			  case 2: // deletion in sequence B
				I_i[i][j] = i;
				I_j[i][j] = j-1;
				break;
			}
			V[i][j] = H[i][j] == -DBL_MAX ? false : true;
		}
	}

	// search H for the maximal score
	unsigned num_of_match = 0;
	double H_max = 0.;
	int i_max=N_a, j_max;
	int* j_max_indexes=new int[N_b]; //this array holds the index of j_max in H[N_a]
	for (j=0; j<N_b; j++)
		j_max_indexes[j]=j+1;

	//sort H[N_a], store the sorted index in j_max_indexes
	sort(j_max_indexes, j_max_indexes+N_b, index_cmp<double*>(H[N_a]));

	//find ALL overlap alignments, starting from the highest score j_max
	j = 0;
	bool found = false;
	while (j < N_b) {
		j_max = j_max_indexes[j];
		H_max = H[N_a][j_max];
		if (H_max == 0)
			break;

		SMAlignment align;
		unsigned align_pos[4];
		num_of_match = Backtrack(i_max, j_max, I_i, I_j, seq_a, seq_b, align, align_pos);
		if (num_of_match) {
			overlaps.push_back(overlap_align(seq_a_start_pos+align_pos[0], align_pos[3], align.match_align, num_of_match));
			if (!found) {
				if (verbose)
					printAlignment(cerr, seq_a, seq_b,
							align_pos, align);
				found = true;
				if (!multi_align
						|| (j+1 < N_b
							&& H[N_a][j_max_indexes[j+1]] < H_max))
					break;
			}
		}
		j++;
	}
	delete [] j_max_indexes;
	for(i=0;i<=N_a;i++){
		delete [] H[i];
		delete [] I_i[i];
		delete [] I_j[i];
		delete [] V[i];
	}
	delete [] H;
	delete [] I_i;
	delete [] I_j;
	delete [] V;
}
