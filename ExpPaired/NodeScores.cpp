#include "NodeScores.h"	
#include "PackedSeq.h"

NodeScores::NodeScores()
{
	
}

void NodeScores::createNode(const PackedSeq& seq)
{
	TreeScore temp;
	temp.indScore = 0.f;
	temp.subtreeScore = 0.f;
	m_scoreRecord[seq] = temp;
}
		
void NodeScores::addScore(const PackedSeq& seq, bool indscore, double weight)
{
	ScoreRecord::iterator iter = m_scoreRecord.find(seq);
	if(iter != m_scoreRecord.end())
	{
		if(indscore)
		{
			iter->second.indScore += weight;
		}
		
		iter->second.subtreeScore += weight;
	}
}

TreeScore NodeScores::getScore(const PackedSeq& seq)
{
	ScoreRecord::iterator iter = m_scoreRecord.find(seq);
	assert(iter != m_scoreRecord.end());
	
	return iter->second;
}

void NodeScores::reset()
{
	for(ScoreRecord::iterator iter = m_scoreRecord.begin(); iter != m_scoreRecord.end(); iter++)
	{
		iter->second.clear();
	}
}

void NodeScores::clearScore(const PackedSeq& seq)
{
	ScoreRecord::iterator iter = m_scoreRecord.find(seq);
	if(iter != m_scoreRecord.end())
	{
		iter->second.clear();
	}	
}

