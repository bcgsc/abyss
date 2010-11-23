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

static inline ostream& printAlignment(ostream& out,
		const string& aseq, const string& bseq,
		unsigned alignPos[], SMAlignment align, unsigned matches)
{
	unsigned astart = alignPos[0], aend = alignPos[1] + 1;
	unsigned bstart = alignPos[2], bend = alignPos[3] + 1;
	assert(aend == aseq.size());
	assert(bstart == 0);
	const string& consensus = align.match_align;
	out << aseq.substr(0, astart) << align.query_align << '\n'
		<< string(astart, ' ');
	for (string::const_iterator it = align.match_align.begin();
			it != align.match_align.end(); ++it) {
		switch (*it) {
		  case 'A': case 'C': case 'G': case 'T':
			out << '.';
			break;
		  default:
			out << *it;
		}
	}
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

template <class T>
void AssignScores_scorematrix(T* temp, T H_im1_jm1, T H_im1_j, T H_i_jm1, char seq_a_im1, char seq_b_jm1,
	bool& opengap_im1_jm1, bool& opengap_im1_j, bool& opengap_i_jm1, bool* temp_opengap);


//index comparison functor
template<class T>
struct index_cmp {
	const T arr;

	index_cmp(const T arr) : arr(arr) {}
	bool operator()(const int a, const int b) const { return arr[a] > arr[b]; }
};

static inline int convert_to_index(char c)
{
	if ('A' <= c && c <= 'Z')
		return c - 'A';
	else
		return c - 'a';
}

static inline bool isMatch(char a, char b, char& consensus)
{
	if (a == b || convert_to_index(a) == convert_to_index(b)) {
		consensus = a;
		return true;
	} else if (a == '-' || b == '-') {
		consensus = '-';
		return false;
	} else if (a == 'N') {
		consensus = b;
		return true;
	} else if (b == 'N') {
		consensus = a;
		return true;
	} else {
		consensus = 'x';
		return false;
	}
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
//mismatch or gap positions are filled with 'N's in match_align string
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

	//bool **opengap; //use this to keep track whether a gap is opengap
	//opengap = new bool*[N_a+1];

	for(i=0;i<=N_a;i++){
		H[i] = new double[N_b+1];
		I_i[i] = new int[N_b+1];
		I_j[i] = new int[N_b+1];
		//opengap[i] = new bool[N_b+1];
		H[i][0]=0; //only need to initialize first row and first column
		//opengap[i][0] = false;
	}

	for (j=0;j<=N_b;j++){
		H[0][j]=0; //initialize first column
		//opengap[0][j] = false;
	}

	double temp[4];
	//bool temp_opengap[3];
	int ind;

	// here comes the actual algorithm
	for(i=1;i<=N_a;i++){
		for(j=1;j<=N_b;j++){
			AssignScores_std<double>(temp, H[i-1][j-1], H[i-1][j], H[i][j-1], seq_a[i-1], seq_b[j-1]);
			/*AssignScores_scorematrix<double>(temp, H[i-1][j-1], H[i-1][j], H[i][j-1], seq_a[i-1], seq_b[j-1],
				opengap[i-1][j-1], opengap[i-1][j], opengap[i][j-1], temp_opengap);*/

			temp[3] = 0.;
			H[i][j] = find_array_max<double>(temp,4, ind);
			switch(ind){
			case 0: // score in (i,j) stems from a match/mismatch
				I_i[i][j] = i-1;
				I_j[i][j] = j-1;
				//opengap[i][j] = temp_opengap[0];
				break;
			case 1: // score in (i,j) stems from a deletion in sequence A
				I_i[i][j] = i-1;
				I_j[i][j] = j;
				//opengap[i][j] = temp_opengap[1];
				break;
			case 2: // score in (i,j) stems from a deletion in sequence B
				I_i[i][j] = i;
				I_j[i][j] = j-1;
				//opengap[i][j] = temp_opengap[2];
				break;
			case 3: // (i,j) is the beginning of a subsequence
				I_i[i][j] = i;
				I_j[i][j] = j;
				//opengap[i][j] = false;
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
		//delete [] opengap[i];
	}
	delete [] H;
	delete [] I_i;
	delete [] I_j;
	//delete [] opengap;

} // END of smith_waterman

// The simple score matrix, match:1, mismatch/gap:-1
int sim_score(const char a, const char b)
{
	char consensus;
	if (isMatch(a, b, consensus))
		return 1;
	else
		if (consensus == '-')
			return -2; //gap penalty
		else
			return -1; //mismatch penalty
}

/*
// find the index of 'a' in the amino-acid substitution matrix, a is (a-z)|(A-Z)|-|*
int score_index(char a)
{
  if (a <= 'Z' && a >= 'A')
	  return a - 'A';
  else
	  if (a <= 'z' && a >= 'a')
		return a - 'a';
	  else
		if (a == '-')
			return 26;
		else
		{
			if (a == '*')
				return 27;
			else
			{
				cout << "cannot find score when running SW local alignment: " << "char:" << a << endl;
				exit(-1);
			}
		}
}

// find the score of (a,b) using amino acid substitution matrix, also update opengap status
int similarity_score(char a,char b, bool opengap, bool& opengap_cur)
{
  //ALIGN_SCORE_MATRIX[28][28] stores the align scores, 26 Uppercase letters (a.a.), '-', '*'
  int ai = score_index(a);
  int bi = score_index(b);
  if (ai != 26 && bi != 26)
  {
	  opengap_cur = false;
	  return ALIGN_SCORE_MATRIX[ai][bi];
  }
  else
  {
	  opengap_cur = true;
	  if (opengap) //previous is gap
		  return ALIGN_SCORE_MATRIX[26][1]; //gap extension
	  else
		  return ALIGN_SCORE_MATRIX[26][0]; //gap opening
  }
}
*/

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

//compute alignment score based on simple score matrix: match is 1, mismatch/deletion/insertion is -1
template <class T>
void AssignScores_std(T* temp, T H_im1_jm1, T H_im1_j, T H_i_jm1, const char seq_a_im1, const char seq_b_jm1)
{
  temp[0] = H_im1_jm1+sim_score(seq_a_im1,seq_b_jm1);
  temp[1] = H_im1_j+sim_score(seq_b_jm1, '-');
  temp[2] = H_i_jm1+sim_score(seq_a_im1, '-');
}

/*
//compute alignment score based on a substitution score matrix
template <class T>
void AssignScores_scorematrix(T* temp, T H_im1_jm1, T H_im1_j, T H_i_jm1, char seq_a_im1, char seq_b_jm1,
	bool& opengap_im1_jm1, bool& opengap_im1_j, bool& opengap_i_jm1, bool* temp_opengap)
{
	temp[0] = H_im1_jm1+similarity_score(seq_a_im1,seq_b_jm1, opengap_im1_jm1, temp_opengap[0]);
	temp[1] = H_im1_j+similarity_score(seq_a_im1, '-', opengap_im1_j, temp_opengap[1]);
	temp[2] = H_i_jm1+similarity_score(seq_b_jm1, '-', opengap_i_jm1, temp_opengap[2]);
	//temp[3] = 0.;
}
*/

//////////////////////////////////////////////////////////////////////////////

