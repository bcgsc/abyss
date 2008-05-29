#include "VisitAlgorithms.h"
#include "Timer.h"

void ContigDataFunctions::merge(const ContigID& targetKey, ContigData& targetData, const ContigID& /*slaveKey*/, const ContigData& slaveData, extDirection dir, bool reverse)
{
	// should the slave be reversed?
	Sequence slaveSeq = slaveData.seq;
	if(reverse)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	Sequence* leftSeq;
	Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &targetData.seq;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &targetData.seq;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - m_overlap, m_overlap);
	PackedSeq rightBegin = rightSeq->substr(0, m_overlap);
	//printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
	// ensure that there is a legitimate k-1 overlap between these sequences
	
	if(leftEnd != rightBegin)
	{
		printf("merge called data1: %s %s (%d, %d)\n", targetData.seq.c_str(), slaveData.seq.c_str(), dir, reverse);	
		printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
	}
	assert(leftEnd == rightBegin);
	
	// TODO: make this in-place?
	// generate the merged sequence
	Sequence merged = *leftSeq;
	merged.append(rightSeq->substr(m_overlap));
	
	targetData.seq = merged;
	
	// merge the seq sets
	targetData.m_seqSet.insert(slaveData.m_seqSet.begin(), slaveData.m_seqSet.end());
	
	// Now, update the lookup table
	m_pDatabase->addKeys(slaveData.m_seqSet, targetKey);
	
	//printf("merge called data1: %s %s (%d, %d)\n", target.seq.c_str(), slave.seq.c_str(), dir, reverse);	
	//printf("new target seq: %s\n", merged.c_str());
}

void ContigDataFunctions::deleteCallback(const ContigID& slaveKey, const ContigData& slaveData)
{
	m_pDatabase->deleteKeys(slaveData.m_seqSet, slaveKey);
}

bool ContigDataFunctions::check(ContigData& target, const ContigData& slave, extDirection dir, bool reverse)
{
	// should the slave be reversed?
	Sequence slaveSeq = slave.seq;
	if(reverse)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	Sequence* leftSeq;
	Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &target.seq;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &target.seq;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - m_overlap, m_overlap);
	PackedSeq rightBegin = rightSeq->substr(0, m_overlap);
	
	return leftEnd == rightBegin;
}

void DBGenerator::visit(const ContigID& id, const ContigData& data)
{
	// add a record for every sequence in the data
	m_pDatabase->addKeys(data.m_seqSet, id);
}

void ContigDataOutputter::visit(const ContigID& id, const ContigData& data)
{
	// super hacky
	int iID = atoi(id.c_str());
	m_pWriter->WriteSequence(data.seq, iID, 0);	
}

int PairedResolver::resolve(const ContigData& data, extDirection dir, ContigDataComponents& components)
{
	printf("Resolving with %zu pairs\n", m_pPairs->getNumReads());
	
	// First, generate the scoring sets, one for each component
	ScoreSets scoreSets;
	{
		Timer sTimer("ScoreSets");
		for(ContigDataComponents::iterator compIter = components.begin(); compIter != components.end(); ++compIter)
		{
			
			PSeqSet newSet;
			scoreSets.push_back(newSet);
			for(ContigDataCollection::iterator dataIter = compIter->begin(); dataIter != compIter->end(); ++dataIter)
			{
				// union the current set in
				scoreSets.back().insert((*dataIter)->m_seqSet.begin(), (*dataIter)->m_seqSet.end());
			}
			
			printf("Component has %zu kmers\n", scoreSets.back().size());
		}
	}
	
	
	// Break the reference contig into kmers
	PSequenceVector refKmers;
	for(size_t idx = 0; idx < data.seq.length() - m_kmer; idx++)
	{
		PackedSeq kmer(data.seq.substr(idx, m_kmer));
		refKmers.push_back(kmer);
	}
	
	// Compute the indices to use
	size_t numNodesInBranch = refKmers.size();
	size_t firstIndex;
	size_t lastIndex;
	size_t numNodesToUse = m_maxlength;
	
	if(dir == SENSE)
	{
		// use the last n sequences
		lastIndex = (int)numNodesInBranch - 1;						
		if((int)lastIndex - (int)numNodesToUse < 0)
		{
			firstIndex = 0;
		}
		else
		{
			firstIndex = lastIndex - numNodesToUse;
		}
	}
	else
	{
		// use the first n sequences	
		firstIndex = 0;
		if(firstIndex + numNodesToUse > numNodesInBranch)
		{
			lastIndex = numNodesInBranch - 1;	
		}
		else
		{
			lastIndex = firstIndex + numNodesToUse - 1;
		}
	}
	
	// Filter the possible pairs by only accepting pairs of sequences that only hit one of the possible branches
	PSequenceVector validPairs;
	int firstSupported = -1;
	int lastSupported = -1;
	
	printf("Support string: \n");
	
	for(size_t idx = firstIndex; idx <= lastIndex; ++idx)
	{
		// Get the pairs for this kmer
		
		PSequenceVector seqPairs;
		if(dir == SENSE)
		{
			m_pPairs->getPairsWithComp(refKmers[idx], false, seqPairs);
		}
		else
		{
			m_pPairs->getPairsWithComp(refKmers[idx], true, seqPairs);
		}
		

		ContigIDSet idCollection;
		for(PSequenceVector::const_iterator seqIter = seqPairs.begin(); seqIter != seqPairs.end(); ++seqIter)
		{
			// append into the id collection
			m_pDatabase->getSet(*seqIter, idCollection);
		}
		
		// score these pairs against the sequence sets
		int numSupported = 0;
		for(ScoreSets::iterator iter = scoreSets.begin(); iter != scoreSets.end(); ++iter)
		{
			size_t currScore = scorePairsToSet(seqPairs, *iter);
			if(currScore > 0)
			{
				numSupported++;
			}
		}
		
		// Get the set of contigs supported
		printf("%d ", numSupported);
		m_pDatabase->printSet(idCollection);
		printf("\n");
		
		//printf("Num supported: %d\n", numSupported);
		
		// if only one branch is supported, add its pairs to the valid set
		if(numSupported == 1)
		{
			validPairs.insert(validPairs.end(), seqPairs.begin(), seqPairs.end());
			
			if(firstSupported == -1)
			{
				firstSupported = idx;
			}
			if((int)idx > lastSupported)
			{
				lastSupported = idx;
			}
		}
	}
	printf("\n");
	
	printf("length: %zu valid index range: [%d %d] valid length: %d, num valid pairs: %zu\n", refKmers.size(), firstSupported, lastSupported, lastSupported - firstSupported, validPairs.size());
	
	const size_t minSupport = 20;
	size_t numSupported = 0;
	int supportedIndex = -1;
	for(size_t idx = 0; idx < scoreSets.size(); idx++)
	{
		size_t currScore = scorePairsToSet(validPairs, scoreSets[idx]);
		printf("idx: %zu score: %zu\n", idx, currScore);
		
		if(currScore > minSupport)
		{
			numSupported++;
			supportedIndex = (int)idx;
		}
	}
	
	if(numSupported == 1)
	{
		return supportedIndex;
	}
	else
	{
		return -1;
	}
}

size_t PairedResolver::scorePairsToSet(PSequenceVector& pairs, PSeqSet& seqSet)
{
	size_t score = 0;	
	for(PSequenceVector::iterator iter = pairs.begin(); iter != pairs.end(); ++iter)
	{
		//printf("Scoring pair %s\n", iter->decode().c_str());
		// look up the seq in the set
		// as the reverse complements are present in the set we dont have to look them up
		if(seqSet.find(*iter) != seqSet.end())
		{
			score++;
		}
	}
	
	return score;
}

