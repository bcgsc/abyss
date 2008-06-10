#include "VisitAlgorithms.h"
#include "Timer.h"
#include "SetOperations.h"

void ContigDataFunctions::merge(const ContigID& targetKey, ContigData& targetData, const ContigID& slaveKey, const ContigData& slaveData, extDirection dir, bool reverse, bool removeMerge)
{
	// should the slave be reversed?
	Sequence slaveSeq = slaveData.m_seq;
	if(reverse)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	Sequence* leftSeq;
	Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &targetData.m_seq;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &targetData.m_seq;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - m_overlap, m_overlap);
	PackedSeq rightBegin = rightSeq->substr(0, m_overlap);
	//printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
	// ensure that there is a legitimate k-1 overlap between these sequences
	
	if(leftEnd != rightBegin)
	{
		printf("merge called data1: %s %s (%d, %d)\n", targetData.m_seq.c_str(), slaveData.m_seq.c_str(), dir, reverse);	
		printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
	}
	assert(leftEnd == rightBegin);
	
	// TODO: make this in-place?
	// generate the merged sequence
	Sequence merged = *leftSeq;
	merged.append(rightSeq->substr(m_overlap));
	
	targetData.m_seq = merged;
	
	// merge the seq sets
	
	// Add the sequence data for the slave into the target, correctly setting the usable flag
	// usable = true if the sequences are equivalent (the child is being deleted), else false
	
	CKDataVec::iterator insertPos;
	if(dir == SENSE)
	{
		// add to end
		insertPos = targetData.m_kmerVec.end();
	}
	else
	{
		// add to beginning
		insertPos = targetData.m_kmerVec.begin();
	}
	
	// add the new sequences to the contigdata
	targetData.addSeqs(slaveSeq, insertPos, removeMerge);
	
	// add the id, indicating the contig is made up of multiple joined contigs
	targetData.addID(slaveKey);
	
	/*
	// Now, update the lookup table
	PSeqSet updateSet;
	slaveData.copySeqSet(updateSet);
	m_pDatabase->addKeys(updateSet, targetKey);
	*/
	
	(void)targetKey;
	(void)slaveKey;
	//printf("merging (%s, %s)\n", targetKey.c_str(), slaveKey.c_str());
	//printf("merge called data1: %s %s (%d, %d)\n", targetData.seq.c_str(), slaveData.seq.c_str(), dir, reverse);	
	//printf("new target seq: %s\n", merged.c_str());
	
	// validate the new target data
	targetData.validateSequences();
}

void ContigDataFunctions::deleteCallback(const ContigID& slaveKey, const ContigData& slaveData)
{
	PSeqSet seqSet;
	slaveData.copySeqSet(seqSet);
	m_pDatabase->deleteKeys(seqSet, slaveKey);
}

bool ContigDataFunctions::check(ContigData& target, const ContigData& slave, extDirection dir, bool reverse)
{
	// should the slave be reversed?
	Sequence slaveSeq = slave.m_seq;
	if(reverse)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	Sequence* leftSeq;
	Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &target.m_seq;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &target.m_seq;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - m_overlap, m_overlap);
	PackedSeq rightBegin = rightSeq->substr(0, m_overlap);
	
	target.validateSequences();
	return leftEnd == rightBegin;
}

void DBGenerator::visit(const ContigID& id, const ContigData& data)
{
	// add a record for every sequence in the data
	PSeqSet seqSet;
	data.copySeqSet(seqSet);
	m_pDatabase->addKeys(seqSet, id);
}

void ContigDataOutputter::visit(const ContigID& id, const ContigData& data)
{
	// super hacky
	int iID = atoi(id.c_str());
	m_pWriter->WriteSequence(data.m_seq, iID, 0);	
}

int PairedResolver::resolve(DirectedGraph<ContigID, ContigData>* pGraph, const ContigID id)
{
	int numMerged = 0;
	printf("\n\n***ATTEMPTING TO RESOLVE %s***\n\n", id.c_str());
	
	// Functors
	SequenceDataCost dataCost(m_kmer);
	LinkFinder linkFinder(m_pPairs, m_pDatabase);
	ContigDataFunctions dataFunctor(m_kmer, m_pDatabase);
	
	for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
	{
		ContigIDSet reachableSet;
		extDirection dir = (extDirection)dirIdx;
		const ContigData& data = pGraph->getDataForVertex(id);
		
		// get the reachable set using the paired info 
		linkFinder.getReachableSet(id, data, dir, reachableSet);
		
		printf("Nodes reachable by pairs from %s: ", id.c_str());
		SetOperations::printSet(reachableSet);
		printf("\n");

		// check if (the shortest?) path from the root node to any node in the set is a super set of all the nodes
		ContigIDVec superPath;	
		bool found = pGraph->findSuperpath(id, dir, reachableSet, superPath, dataCost);

		// Create the list of sequences that have been merged so far
		ContigIDSet mergedPathSet;
		
		// add the source node to the merged path set
		mergedPathSet.insert(id);
		
		
		if(found)
		{
			for(ContigIDVec::const_iterator iter = superPath.begin(); iter != superPath.end(); ++iter)
			{
				// check if the node should be appended (target node stays) or merged (target node removed)
				ContigIDSet reachableSet;
				
				// Get the edge from the parent to the child
				
				EdgeDescription parentEdgeDesc = pGraph->getUniqueEdgeDesc(id, *iter, dir);
				EdgeDescription relativeDesc = parentEdgeDesc.makeRelativeDesc();
								
				printf("Parent dir: %d relative dir: %d\n", parentEdgeDesc.dir, relativeDesc.dir);
				
				// Get the reachable set for the target node in the direction towards the parents
				ContigIDSet targetReachSet;
				ContigID currID = *iter;
				const ContigData& currData = pGraph->getDataForVertex(currID);
				
				
				linkFinder.getReachableSet(currID, currData, !relativeDesc.dir, targetReachSet);
				printf("Nodes reachable by pairs from %s: ", currID.c_str());
				SetOperations::printSet(targetReachSet);
				printf("\n");
				
				// Get the difference between the merged so far this set
				
				// If the resulting set is empty, the defined path is a superset of this node and the node can be removed
				ContigIDSet pathDiff;
				SetOperations::difference(targetReachSet, mergedPathSet, pathDiff);
				
				printf("Difference from merged path: ");
				SetOperations::printSet(pathDiff);
				printf("\n");
				
				// if the path child's reachable set is a strict subset of the parent's path, it can be removed
				bool removeChild = pathDiff.empty() && !targetReachSet.empty();
				
				if(currID == "2313")
					removeChild = false;
				pGraph->mergePath(id, *iter, dir, removeChild, dataFunctor);
				
				mergedPathSet.insert(currID);
				numMerged++;
			}
		}
	}
	
	return numMerged;
}


void LinkFinder::visit(const ContigID& /*id*/, ContigData& /*data*/)
{
	assert(false);
}

void LinkFinder::getReachableSet(const ContigID& /*id*/, const ContigData& data, extDirection dir, ContigIDSet& idSet)
{
	// Copy out the usable sequences
	PSeqSet usableSeqs;
	data.copyUsableSeqSet(usableSeqs);
	
	for(PSeqSet::iterator iter = usableSeqs.begin(); iter != usableSeqs.end(); ++iter)
	{		
		PSequenceVector seqPairs;
		if(dir == SENSE)
		{
			m_pPairs->getPairsWithComp(*iter, false, seqPairs);
		}
		else
		{
			m_pPairs->getPairsWithComp(*iter, true, seqPairs);
		}
		

		for(PSequenceVector::const_iterator seqIter = seqPairs.begin(); seqIter != seqPairs.end(); ++seqIter)
		{
			// append into the id collection
			m_pDatabase->getSet(*seqIter, idSet);
		}
	}
	
	// Remove the IDs that are in the exclusion list (they have been merged already)
	SetOperations::removeElements(idSet, data.getExclusionSet());

	return;
}

