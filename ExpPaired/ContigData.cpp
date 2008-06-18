#include "ContigData.h"
#include "SetOperations.h"
#include "VisitAlgorithms.h"

ContigData::ContigData(const ContigID& id, const Sequence& s, size_t kmer, int copyNumber, AlignmentCache* pDB) : m_seq(s), m_id(id), m_kmer(kmer), m_copyNumber(copyNumber), m_pDatabase(pDB)
{
	// Create the seq set
	size_t numNewSeqs = getNumKmers(m_seq);	
	m_kmerVec.reserve(numNewSeqs);
	
	// Create the kmers for this sequence
	CKDataVec tempVec;
	for(size_t idx = 0; idx < numNewSeqs; ++idx)
	{
		// Create the new sequence for insertion
		PackedSeq insSeq(m_seq.substr(idx, m_kmer));
		
		// Create the data structure
		// Note the pairs are generated seperately (for now)
		ContigKmerData ndata;
		ndata.seq = insSeq;
		ndata.usable = true;
		
		// 
		tempVec.push_back(ndata);
	}	
	
	// Add the initial sequences, this will populate the database with the intial kmer alignments
	addSeqs(m_kmerVec.begin(), tempVec);
	
	addID(id, SENSE);
	addID(id, ANTISENSE);
}

//
//
//
void ContigData::merge(const ContigData& other, extDirection dir, bool isReversed, bool isUsable)
{
	// First decide the insertion place, the beginning or end of the sequence
	CKDataVec::iterator insertPos;
	if(dir == SENSE)
	{
		insertPos = m_kmerVec.end();
	}
	else
	{
		insertPos = m_kmerVec.begin();
	}
	
	// Create the vector of sequences to insert, reversing them if necessary
	// Should the pairs be marked as resolved at this point?
	CKDataVec tempVec;
	for(CKDataVec::const_iterator otherIter = other.m_kmerVec.begin(); otherIter != other.m_kmerVec.end(); ++otherIter)
	{
		ContigKmerData nData;
		if(isReversed)
		{
			nData.seq = reverseComplement(otherIter->seq);
			
			// swap the pairs
			nData.pairs[0] = otherIter->pairs[1];
			nData.pairs[1] = otherIter->pairs[0];
			
		}
		else
		{
			nData.seq = otherIter->seq;
			nData.pairs[0] = otherIter->pairs[0];
			nData.pairs[1] = otherIter->pairs[1];
		}
		
		nData.usable = isUsable;
		tempVec.push_back(nData);
	}
	
	// Finally, reverse the vector if the other contig is reverse-comp
	if(isReversed)
	{
		std::reverse(tempVec.begin(), tempVec.end());
	}
	
	// Perform the insertion
	addSeqs(insertPos, tempVec);
	
	// Now, append the sequence string into m_seq
	appendSeqString(other.m_seq, dir, isReversed);
	
	// Add the id of the other contig
	addID(other.m_id, dir);
}

//
// add sequences to the kmer record and update the db
//
void ContigData::addSeqs(CKDataVec::iterator position, CKDataVec& newSeqs)
{	
	// Remove all the alignments to this contig for the current sequences
	size_t currPosition = 0;
	for(CKDataVec::iterator removeIter = m_kmerVec.begin(); removeIter != m_kmerVec.end(); ++removeIter)
	{
		m_pDatabase->removeAlignment(removeIter->seq, m_id, currPosition);
		currPosition++;
	}

	// Perform the actual insertion
	m_kmerVec.insert(position, newSeqs.begin(), newSeqs.end());
	
	// Now readd the alignments for the new contig
	currPosition = 0;
	for(CKDataVec::iterator addIter = m_kmerVec.begin(); addIter != m_kmerVec.end(); ++addIter)
	{
		m_pDatabase->addAlignment(addIter->seq, m_id, currPosition);
		currPosition++;
	}
	
	validateDatabase(m_pDatabase);
}

//
//
//
void ContigData::validateSequences() const
{
	assert((m_seq.length() - m_kmer + 1) == m_kmerVec.size());
	
	for(size_t idx = 0; idx < m_seq.length() - m_kmer + 1; ++idx)
	{
		PackedSeq refSeq(m_seq.substr(idx, m_kmer));
		const PackedSeq& kmer = m_kmerVec[idx].seq;
		
		assert(refSeq == kmer);
	}
}

//
//
//
void ContigData::copySeqSet(PSeqSet& outSeqs) const
{
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		outSeqs.insert(iter->seq);
	}
}

//
//
//
void ContigData::copyUsableSeqSet(PSeqSet& outSeqs) const
{
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		if(iter->usable)
		{
			outSeqs.insert(iter->seq);
		}
	}
}

//
// Add the pairs of the sequences of this contig to the record
//
void ContigData::addPairs(PairRecord* pairRecord)
{
	for(CKDataVec::iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		for(size_t idx = 0; idx <= 1; ++idx)
		{
			bool rc = (idx == 0) ? false : true;
			// Look up the pairs for this sequence
			
			PSequenceVector seqPairs;
			pairRecord->getPairsWithComp(iter->seq, rc, seqPairs);
			
			// add the pairs to the record for this sequence
			for(PSequenceVector::iterator pairIter = seqPairs.begin(); pairIter != seqPairs.end(); ++pairIter)
			{
				PairData data;
				data.seq = *pairIter;
				data.resolved = false;
				
				iter->pairs[idx].push_back(data);
			}
		}
	}
}

//
//
//
void ContigData::getUnresolvedUsableSet(extDirection dir, ContigSupportMap& idSupport) const
{
	size_t pairIdx = dir2Idx(dir);
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		if(iter->usable)
		{
			const PairVector& pairVec = iter->pairs[pairIdx];	
			for(PairVector::const_iterator pairIter = pairVec.begin(); pairIter != pairVec.end(); ++pairIter)
			{
				if(!pairIter->resolved)
				{
					AlignmentPairVec alignPairs;
					getAlignmentPairs(iter->seq, pairIter->seq, alignPairs);
					
					if(alignPairs.empty())
					{
						continue;
					}
					
					// Disallow alignments of the same sequence to different contigs
					// Maybe disallow multiple in general?
					bool uniqueAlign = true;
					ContigID uniqueContig = alignPairs.front().pairSeqAlign.alignment.contigID;
					
					for(AlignmentPairVec::iterator apIter = alignPairs.begin(); apIter != alignPairs.end(); ++apIter)
					{
						if(uniqueContig != apIter->pairSeqAlign.alignment.contigID)
						{
							uniqueAlign = false;
							break;
						}
					}
					
					if(uniqueAlign)
					{
						idSupport[uniqueContig]++;
					}
				}
			}
		}		
	}
}

//
//
//
void ContigData::getAllAlignmentPairs(const ContigKmerData& refKmerData, extDirection dir, AlignmentPairVec& outAligns) const
{
	size_t pairIdx = dir2Idx(dir);
	
	// Now get the alignments for its pairs
	// Iterate over the pairs for this sequence
	for(PairVector::const_iterator pairIter = refKmerData.pairs[pairIdx].begin(); pairIter != refKmerData.pairs[pairIdx].end(); ++pairIter)
	{
		getAlignmentPairs(refKmerData.seq, pairIter->seq, outAligns);
	}	
}


//
// This could be more efficient
//
void ContigData::getAlignmentPairs(const PackedSeq& refSeq, const PackedSeq& pairSeq, AlignmentPairVec& outAligns) const
{
	// Get the alignments for each
	AlignSet seqAlignments;
	m_pDatabase->getAlignments(refSeq, seqAlignments);
	
	AlignSet pairAlignments;
	m_pDatabase->getAlignments(pairSeq, pairAlignments);
	
	// Make up all the pairs
	for(AlignSet::const_iterator saIter = seqAlignments.begin(); saIter != seqAlignments.end(); ++saIter)
	{
		for(AlignSet::const_iterator paIter = pairAlignments.begin(); paIter != pairAlignments.end(); ++paIter)
		{
			SeqAlignment ref = {refSeq, *saIter};
			SeqAlignment pair = {pairSeq, *paIter};
			
			AlignmentPair alignPair = {ref, pair};
			outAligns.push_back(alignPair);
		}
	}
}

//
// Mark pairs as resolved if their alignment is consistent with 
//
void ContigData::resolvePairs(PairedResolvePolicy* pResolvePolicy, extDirection dir)
{
	size_t pairIdx = dir2Idx(dir);	
	for(CKDataVec::iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		// Only resolve usable pairs
		if(iter->usable)
		{
			// Get the pairs for this read, in the specified direction
			PairVector& pairVec = iter->pairs[pairIdx];
	
			for(PairVector::iterator pairIter = pairVec.begin(); pairIter != pairVec.end(); ++pairIter)
			{
				// pairs can only be resolved once
				if(!pairIter->resolved)
				{
					// Check if any alignment between these two pairs are resolved, if so, mark them as such
					bool resolved = false;
					AlignmentPairVec alignPairs;
					getAlignmentPairs(iter->seq, pairIter->seq, alignPairs);
					
					// Check if these two reads are aligned within the expected distance distribution
					for(AlignmentPairVec::iterator apIter = alignPairs.begin(); apIter != alignPairs.end(); ++apIter)
					{
						if(apIter->refSeqAlign.alignment.contigID == apIter->pairSeqAlign.alignment.contigID)
						{
							// Calculate the distance between the pairs
							size_t pos1 = apIter->refSeqAlign.alignment.position;
							size_t pos2 = apIter->pairSeqAlign.alignment.position;
							size_t distance = (pos1 > pos2) ? pos1 - pos2 : pos2 - pos1;
							
							// Check if the pair is resolved
							if(pResolvePolicy->isResolved(apIter->refSeqAlign.alignment.contigID, apIter->pairSeqAlign.alignment.contigID, distance))
							{
								resolved = true;
								break;
							}
						}
					}
					
					// Mark the pair as resolved
					if(resolved)
					{
						pairIter->resolved = true;
					}
				}
			}
		}
	}
}

//
// Add the pairs where both sequences are on this contig
//
void ContigData::addSelfPairsToHist(Histogram* pHist)
{
	// We only check one direction as they are equivalent and checking both would be redundant
	extDirection dir = SENSE;
	for(CKDataVec::iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		AlignmentPairVec alignPairs;
		getAllAlignmentPairs(*iter, dir, alignPairs);
		
		for(AlignmentPairVec::iterator apIter = alignPairs.begin(); apIter != alignPairs.end(); ++apIter)
		{
			if(apIter->refSeqAlign.alignment.contigID == apIter->pairSeqAlign.alignment.contigID)
			{
				int difference = abs(apIter->refSeqAlign.alignment.position - apIter->pairSeqAlign.alignment.position);
				pHist->addDataPoint(difference);
			}
		}
	}
}

//
//
//
void ContigData::printPairAlignments(extDirection dir, bool showResolved) const
{
	size_t pairIdx = dir2Idx(dir);
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		PairVector pairVec = iter->pairs[pairIdx];
		for(PairVector::const_iterator pairIter = pairVec.begin(); pairIter != pairVec.end(); ++pairIter)
		{
			AlignmentPairVec alignPairs;
			getAlignmentPairs(iter->seq, pairIter->seq, alignPairs);
			
			for(AlignmentPairVec::iterator apIter = alignPairs.begin(); apIter != alignPairs.end(); ++apIter)
			{		
				if(!pairIter->resolved || showResolved)
				{
					std::cout << *apIter << " Resolved: " << pairIter->resolved << " " << iter->seq.decode() << " " << pairIter->seq.decode() << std::endl;
				}
			}
		}
	}
}

void ContigData::computeTestStat(const PDF& empDist)
{
	Histogram testHist;
	addSelfPairsToHist(&testHist);
	ChiSquare(empDist, testHist);
}

//
//
//
void ContigData::writeToAlignDB(AlignmentCache* pDB) const
{
	int position = 0;
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		pDB->addAlignment(iter->seq, m_id, position);
		++position;
	}
}

//
//
//
void ContigData::removeFromAlignDB(AlignmentCache* pDB) const
{
	int position = 0;
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		pDB->removeAlignment(iter->seq, m_id, position);
		++position;
	}
}

//
//
//
void ContigData::validateDatabase(AlignmentCache* pDB) const
{
	AlignSet alignments;
	
	int position = 0;
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		pDB->getAlignments(iter->seq, alignments);

		// Make sure that the alignment set contains the alignment to this sequence
		AlignData lookupAlign;
		lookupAlign.contigID = m_id;
		lookupAlign.position = position;
		
		// Make sure an alignment to this position exists
		if(alignments.find(lookupAlign) == alignments.end())
		{
			std::cout << "Expected alignment of " << lookupAlign << " failed!" << std::endl;
			pDB->printAlignmentsForSeq(iter->seq);
			assert(false);
		}
		++position;
	}
}

void ContigData::appendSeqString(const Sequence& seq, extDirection dir, bool isReversed)
{
	size_t overlap = m_kmer - 1;
	
	// should the slave be reversed?
	Sequence slaveSeq = seq;
	if(isReversed)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	Sequence* leftSeq;
	Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &m_seq;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &m_seq;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - overlap, overlap);
	PackedSeq rightBegin = rightSeq->substr(0, overlap);

	// ensure that there is a legitimate k-1 overlap between these sequences	
	if(leftEnd != rightBegin)
	{
		printf("merge called data1: %s %s (%d, %d)\n", m_seq.c_str(), seq.c_str(), dir, isReversed);	
		printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
		assert(leftEnd == rightBegin);
	}
	
	// TODO: make this in-place?
	// generate the merged sequence
	Sequence merged = *leftSeq;
	merged.append(rightSeq->substr(overlap));
	
	m_seq = merged;	
}
