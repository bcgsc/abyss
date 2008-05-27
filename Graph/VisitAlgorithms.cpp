#include "VisitAlgorithms.h"

void ContigDataFunctions::merge(ContigData& target, const ContigData& slave, extDirection dir, bool reverse)
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
	//printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
	// ensure that there is a legitimate k-1 overlap between these sequences
	
	if(leftEnd != rightBegin)
	{
		printf("merge called data1: %s %s (%d, %d)\n", target.seq.c_str(), slave.seq.c_str(), dir, reverse);	
		printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
	}
	assert(leftEnd == rightBegin);
	
	// TODO: make this in-place?
	// generate the merged sequence
	Sequence merged = *leftSeq;
	merged.append(rightSeq->substr(m_overlap));
	
	target.seq = merged;
	
	//printf("merge called data1: %s %s (%d, %d)\n", target.seq.c_str(), slave.seq.c_str(), dir, reverse);	
	//printf("new target seq: %s\n", merged.c_str());
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
	for(ContigDataComponents::iterator compIter = components.begin(); compIter != components.end(); ++compIter)
	{
		PSeqSet newSet;
		scoreSets.push_back(newSet);
		for(ContigDataCollection::iterator dataIter = compIter->begin(); dataIter != compIter->end(); ++dataIter)
		{
			// break the data into kmers
			Sequence& seq = dataIter->seq;
			//printf("Seq len %zu\n", seq.length());
			for(size_t idx = 0; idx < seq.length() - m_kmer; idx++)
			{
				PackedSeq kmer(seq.substr(idx, m_kmer));
				//printf("adding %s\n", kmer.decode().c_str());
				scoreSets.back().insert(kmer);
				scoreSets.back().insert(reverseComplement(kmer));
			}
		}
		
		printf("Component has %zu kmers\n", scoreSets.back().size());
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
	
	if(dir == SENSE)
	{
		// use the last n sequences
		lastIndex = (int)numNodesInBranch - 1;						
		if((int)lastIndex - (int)m_maxlength < 0)
		{
			firstIndex = 0;
		}
		else
		{
			firstIndex = lastIndex - m_maxlength;
		}
	}
	else
	{
		// use the first n sequences	
		firstIndex = 0;
		if(firstIndex + m_maxlength > numNodesInBranch)
		{
			lastIndex = numNodesInBranch - 1;	
		}
		else
		{
			lastIndex = firstIndex + m_maxlength - 1;
		}
	}
	
	// Filter the possible pairs by only accepting pairs of sequences that only hit one of the possible branches
	PSequenceVector validPairs;
	int firstSupported = -1;
	int lastSupported = -1;
	
	printf("Support string: ");
	for(size_t idx = firstIndex; idx <= lastIndex; ++idx)
	{
		// Get the pairs for this kmer
		
		// TODO: make this a reference
		PSequenceVector seqPairs = m_pPairs->getPairs(refKmers[idx]);
		
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
		printf("%d", numSupported);
		
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

