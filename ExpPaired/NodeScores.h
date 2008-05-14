#ifndef NODESCORES_H
#define NODESCORES_H

#include <map>
#include "CommonDefs.h"


struct TreeScore
{
	double indScore;
	double subtreeScore;
};

typedef std::map<PackedSeq, TreeScore> ScoreRecord;

class NodeScores
{
	public:
		NodeScores();
		void createNode(const PackedSeq& seq);
		void addScore(const PackedSeq& seq, bool indscore, double weight);
		TreeScore getScore(const PackedSeq& seq);
		
	private:
		ScoreRecord m_scoreRecord;
};

#endif
