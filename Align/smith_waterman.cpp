//////////////////////////////////////////////////////////////////////////////
// Smith_waterman implementation, with modification, to find sequence overlap
//
// Written by Rong She (rshe@bcgsc.ca)
// Last modified: Jul 7, 2010
//////////////////////////////////////////////////////////////////////////////

#include "smith_waterman.h"
#include "Sequence.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <iostream>

using namespace std;

/** Print the specified alignment. */
static ostream& printAlignment(ostream& out,
		const string& aseq, const string& bseq,
		unsigned alignPos[], SMAlignment align, unsigned matches)
{
	unsigned astart = alignPos[0], aend = alignPos[1] + 1;
	unsigned bstart = alignPos[2], bend = alignPos[3] + 1;
	assert(aend == aseq.size());
	assert(bstart == 0);
	out << aseq.substr(0, astart) << align.query_align << '\n'
		<< string(astart, ' ');
	for (unsigned i = 0; i < align.match_align.size(); i++) {
		char a = align.query_align[i],
			b = align.target_align[i],
			c = align.match_align[i];
		out << (c == a || c == b ? '.' : c);
	}
	const string& consensus = align.match_align;
	return out << '\n'
		<< string(astart, ' ')
			<< align.target_align << bseq.substr(bend) << '\n'
		<< matches << " / " << consensus.size() << " = "
			<< (float)matches / consensus.size() << '\n';
}

#ifndef __GNUC__
#ifndef NGNUC_INDEX
#define NGNUC_INDEX
#define INDEX_ONE_DIM(num_rows, row, col) (row * num_rows + col)
#endif
#endif

template <class T>
T find_array_max(T array[],int length, int& ind);

template <class T>
void AssignScores_std(T* temp, T H_im1_jm1, T H_im1_j, T H_i_jm1, const char seq_a_im1, const char seq_b_jm1);

//index comparison functor
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

//the backtrack step in smith_waterman
#ifdef __GNUC__
unsigned Backtrack(const int i_max, const int j_max, int** I_i, int** I_j,
		const string& seq_a, const string& seq_b, SMAlignment& align, unsigned* align_pos)
#else
unsigned Backtrack(const int i_max, const int j_max, int* I_i, int* I_j, const int I_i_rows, const int I_j_rows,
		const string& seq_a, const string& seq_b, SMAlignment& align, unsigned* align_pos)
#endif
{
	// Backtracking from H_max
	int current_i=i_max,current_j=j_max;
#ifdef __GNUC__
	int next_i=I_i[current_i][current_j];
	int next_j=I_j[current_i][current_j];
#else
	int next_i=I_i[INDEX_ONE_DIM(I_i_rows, current_i, current_j)];
	int next_j=I_j[INDEX_ONE_DIM(I_j_rows, current_i, current_j)];
#endif
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
#ifdef __GNUC__
		next_i=I_i[current_i][current_j];
		next_j=I_j[current_i][current_j];
#else
		next_i = I_i[INDEX_ONE_DIM(I_i_rows, current_i, current_j)];
		next_j = I_j[INDEX_ONE_DIM(I_j_rows, current_i, current_j)];
#endif
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

//This is the Smith-Waterman algirthm (a variation of Needleman-Wunsch algorithm), finds one optimal local alignment
//Modified to find overlap of seq_a and seq_b (alignment that is pinned at the end of seq_a and beginning of seq_b)
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

	for(i=0;i<=N_a;i++){
		H[i] = new double[N_b+1];
		I_i[i] = new int[N_b+1];
		I_j[i] = new int[N_b+1];
		H[i][0]=0; //only need to initialize first row and first column
	}

	for (j=0;j<=N_b;j++)
		H[0][j]=0; //initialize first column

	double temp[4];
	int ind;
	for(i=1;i<=N_a;i++){
		for(j=1;j<=N_b;j++){
			AssignScores_std<double>(temp, H[i-1][j-1], H[i-1][j], H[i][j-1], seq_a[i-1], seq_b[j-1]);
			temp[3] = 0.;
			H[i][j] = find_array_max<double>(temp,4, ind);
			switch(ind){
			case 0: // score in (i,j) stems from a match/mismatch
				I_i[i][j] = i-1;
				I_j[i][j] = j-1;
				break;
			case 1: // score in (i,j) stems from a deletion in sequence A
				I_i[i][j] = i-1;
				I_j[i][j] = j;
				break;
			case 2: // score in (i,j) stems from a deletion in sequence B
				I_i[i][j] = i;
				I_j[i][j] = j-1;
				break;
			case 3: // (i,j) is the beginning of a subsequence
				I_i[i][j] = i;
				I_j[i][j] = j;
				break;
			}
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

		//if H_max is 0, report failure
		if (H_max == 0) {
			if (verbose)
				cout //<< seq_a << endl << seq_b << endl
					<< "H_max is 0, no more alignment!" << endl;
			break;
		}

		SMAlignment align;
		unsigned align_pos[4];
#ifdef __GNUC__
		num_of_match = Backtrack(i_max, j_max, I_i, I_j, seq_a, seq_b, align, align_pos);
#else
		num_of_match = Backtrack(i_max, j_max, I_i, I_j, N_a+1, N_a+1,
					seq_a, seq_b, align, align_pos);
#endif
		if (num_of_match) {
			overlaps.push_back(overlap_align(seq_a_start_pos+align_pos[0], align_pos[3], align.match_align, num_of_match));
			if (!found) {
				if (verbose)
					printAlignment(cout, seq_a, seq_b,
							align_pos, align, num_of_match) << '\n';
				found = true;
				if (!multi_align)
					break;
			}
		}
		j++;
	}
	delete [] j_max_indexes;

	if (verbose && !found)
		cout <<"no alignment found" << endl;

	//clean up memory
	for(i=0;i<=N_a;i++){
		delete [] H[i];
		delete [] I_i[i];
		delete [] I_j[i];
	}
	delete [] H;
	delete [] I_i;
	delete [] I_j;
}

/** Return the score of the alignment of a and b. */
static int sim_score(const char a, const char b)
{
	char consensus;
	return a == '-' || b == '-' ? -2 // gap penalty
		: !isMatch(a, b, consensus) ? -1 // mismatch penalty
		: 1; // match score
}

// Find the max value in array[], use 'ind' to hold the index of max
template <class T>
T find_array_max(T array[],int length, int& ind)
{

  T max = array[0];            // start with max = first element
  ind = 0;

  for(int i = 1; i<length; i++){
      if(array[i] > max){
		max = array[i];
		ind = i;
      }
  }
  return max;                    // return highest value in array
}

/** Return the alignment score. */
template <class T>
void AssignScores_std(T* temp, T H_im1_jm1, T H_im1_j, T H_i_jm1, const char seq_a_im1, const char seq_b_jm1)
{
  temp[0] = H_im1_jm1+sim_score(seq_a_im1,seq_b_jm1);
  temp[1] = H_im1_j+sim_score(seq_b_jm1, '-');
  temp[2] = H_i_jm1+sim_score(seq_a_im1, '-');
}
