/** Global sequence alignment with an affine gap penalty using the
 * Needleman-Wunsch algorithm and the improvement by Gotoh.
 * @author Shaun Jackman <sjackman@bcgsc.ca>
 */

#include "alignGlobal.h"
#include "Sequence.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib> // for abort

using namespace std;

/** A character representing a gap. */
static const char GAP = '*';

/** The score of a match. */
static const int MATCH = 5;

/** The penalty of a mismatch. */
static const int MISMATCH = -4;

/** The penalty of opening a gap. */
static const int GAP_OPEN = -12;

/** The penalty of extending a gap. */
static const int GAP_EXTEND = -4;

/** Return the score of the alignment of a and b.
 * @param [out] consensus the consensus of a and b
 * @return the score
 */
static int score(char a, char b, char& c)
{
	if (a == b) {
		c = a;
		return MATCH;
	} else {
		c = ambiguityOr(a, b);
		return c == a || c == b ? MATCH : MISMATCH;
	}
}

/** Return the score of the alignment of a and b. */
static int score(char a, char b)
{
	char c;
	return score(a, b, c);
}

/** Find the optimal alignment from the score matrices.
 * @param[out] align the alignment
 * @return the number of matches
 */
static unsigned backtrack(int** f, int** g, int** h,
		const string& seqA, const string& seqB, NWAlignment& align)
{
	string alignmentA, alignmentB, consensus;
	unsigned matches = 0;
	unsigned i = seqA.size(), j = seqB.size();
	while (i > 0 && j > 0) {
		int fij = f[i][j];
		char a = seqA[i-1], b = seqB[j-1], c;
		int s = score(a, b, c);
		if (fij == f[i-1][j-1] + s) {
			alignmentA += a;
			alignmentB += b;
			consensus += c;
			if (s == MATCH)
				matches++;
			i--;
			j--;
		} else if (fij == f[i-1][j] + GAP_OPEN
				|| fij == g[i-1][j] + GAP_EXTEND) {
			while (g[i][j] == g[i-1][j] + GAP_EXTEND) {
				char a = seqA[i-1];
				alignmentA += a;
				alignmentB += GAP;
				consensus += tolower(a);
				i--;
				assert(i > 0);
			}
			assert(g[i][j] == f[i-1][j] + GAP_OPEN);
			char a = seqA[i-1];
			alignmentA += a;
			alignmentB += GAP;
			consensus += tolower(a);
			i--;
		} else if (fij == f[i][j-1] + GAP_OPEN
				|| fij == h[i][j-1] + GAP_EXTEND) {
			while (h[i][j] == h[i][j-1] + GAP_EXTEND) {
				char b = seqB[j-1];
				alignmentA += GAP;
				alignmentB += b;
				consensus += tolower(b);
				j--;
				assert(j > 0);
			}
			assert(h[i][j] == f[i][j-1] + GAP_OPEN);
			char b = seqB[j-1];
			alignmentA += GAP;
			alignmentB += b;
			consensus += tolower(b);
			j--;
		} else {
			assert(false);
			abort();
		}
	}

	while (i > 0) {
		char a = seqA[i-1];
		alignmentA += a;
		alignmentB += GAP;
		consensus += tolower(a);
		i--;
	}

	while (j > 0) {
		char b = seqB[j-1];
		alignmentA += GAP;
		alignmentB += b;
		consensus += tolower(b);
		j--;
	}

	reverse(alignmentA.begin(), alignmentA.end());
	reverse(alignmentB.begin(), alignmentB.end());
	reverse(consensus.begin(), consensus.end());
	align.query_align = alignmentA;
	align.target_align = alignmentB;
	align.match_align = consensus;
	return matches;
}

/** Find the optimal global alignment of the two sequences using the
 * Needleman-Wunsch algorithm and the improvement by Gotoh to use an
 * affine gap penalty rather than a linear gap penalty.
 * @param[out] align the alignment
 * @return the number of matches
 */
unsigned alignGlobal(const string& seqA, const string& seqB,
		NWAlignment& align)
{
	unsigned lenA = seqA.size();
	unsigned lenB = seqB.size();
	int** f = new int*[lenA + 1];
	int** g = new int*[lenA + 1];
	int** h = new int*[lenA + 1];
	for (unsigned i = 0; i <= lenA; i++) {
		f[i] = new int[lenB + 1];
		g[i] = new int[lenB + 1];
		h[i] = new int[lenB + 1];
	}

	// Initialize the score matrix.
	for (unsigned i = 0; i <= lenA; i++) {
		f[i][0] = g[i][0] = i == 0 ? 0
			: GAP_OPEN + GAP_EXTEND * ((int)i - 1);
		h[i][0] = INT_MIN/2;
	}
	for (unsigned j = 0; j <= lenB; j++) {
		f[0][j] = h[0][j] = j == 0 ? 0
			: GAP_OPEN + GAP_EXTEND * ((int)j - 1);
		g[0][j] = INT_MIN/2;
	}

	// Calculate the score matrix.
	for (unsigned i = 1; i <= lenA; i++) {
		for (unsigned j = 1; j <= lenB; j++) {
			g[i][j] = max(
					f[i-1][j] + GAP_OPEN,
					g[i-1][j] + GAP_EXTEND);
			h[i][j] = max(
					f[i][j-1] + GAP_OPEN,
					h[i][j-1] + GAP_EXTEND);
			f[i][j] = max(
					f[i-1][j-1] + score(seqA[i-1], seqB[j-1]),
					max(g[i][j], h[i][j]));
		}
	}

	// Find the best alignment.
	unsigned matches = backtrack(f, g, h, seqA, seqB, align);

	for (unsigned i = 0; i <= lenA; i++) {
		delete[] f[i];
		delete[] g[i];
		delete[] h[i];
	}
	delete[] f;
	delete[] g;
	delete[] h;

	return matches;
}
