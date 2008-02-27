#include "NetworkSequenceCollection.h"
#include "FastaReader.h"
#include "SeqRecord.h"

//
//
//
NetworkSequenceCollection::NetworkSequenceCollection(int myID, int numDataNodes, int kmerSize) : m_id(myID), m_numDataNodes(numDataNodes), m_state(NAS_LOADING)
{
	// Load the phase space
	m_pLocalSpace = new SequenceCollection();
	
	// Create the comm layer
	m_pComm = new CommLayer(myID, kmerSize);
}

//
//
//
NetworkSequenceCollection::~NetworkSequenceCollection()
{
	// Delete the objects created in the constructor
	delete m_pLocalSpace;
	m_pLocalSpace = 0;
	
	delete m_pComm;
	m_pComm = 0;
}

//
//
//
void NetworkSequenceCollection::run(int readLength, int kmerSize)
{
	bool stop = false;
	while(!stop)
	{
		switch(m_state)
		{
			case NAS_LOADING:
			{
				// Spin in the message loop, waiting for sequences
				pumpNetwork();
				break;
			}
			case NAS_FINALIZE:
			{
				finalize();
				SetState(NAS_GEN_ADJ);
				break;
			}
			case NAS_GEN_ADJ:
			{
				generateAdjacency();
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();
				//SendControlMessageTo 
				break;
			}
			case NAS_TRIM:
			{
				performTrim(readLength, kmerSize);
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();				
				break;
			}
			case NAS_ASSEMBLE:
			{
				printf("doing noncontrol assemble\n");
				assemble();
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();	
				break;
			}
			case NAS_WAITING:
			{
				pumpNetwork();
				break;
			}
			case NAS_DONE:
			{
				stop = true;
				break;
			}						
		}
	}
}

//
// The main 
//
void NetworkSequenceCollection::runControl(std::string fastaFile, int readLength, int kmerSize)
{
	bool stop = false;
	while(!stop)
	{
		switch(m_state)
		{
			case NAS_LOADING:
			{
				readSequences(fastaFile, readLength, kmerSize);
				SetState(NAS_FINALIZE);
				
				// Tell the rest of the loaders that the load is finished
				m_pComm->SendControlMessage(m_numDataNodes, APC_DONELOAD);
				break;
			}
			case NAS_FINALIZE:
			{
				finalize();
				SetState(NAS_GEN_ADJ);
				break;
			}
			case NAS_GEN_ADJ:
			{
				generateAdjacency();
				m_numReachedCheckpoint++;
				while(!checkpointReached())
				{
					pumpNetwork();
				}
				
				SetState(NAS_TRIM);
				m_pComm->SendControlMessage(m_numDataNodes, APC_TRIM);
				break;
			}
			case NAS_TRIM:
			{
				performTrim(readLength, kmerSize);
				m_numReachedCheckpoint++;
				while(!checkpointReached())
				{
					pumpNetwork();
				}
				
				SetState(NAS_ASSEMBLE);
				m_pComm->SendControlMessage(m_numDataNodes, APC_ASSEMBLE);							
				break;
			}
			case NAS_ASSEMBLE:
			{
				printf("doing control assemble\n");
				assemble();
				m_numReachedCheckpoint++;
				while(!checkpointReached())
				{
					pumpNetwork();
				}
				
				SetState(NAS_DONE);
				m_pComm->SendControlMessage(m_numDataNodes, APC_FINISHED);				
				break;
			}
			case NAS_DONE:
			{
				printf("in done state\n");
				stop = true;
				break;
			}											
		}
	}
}

//
// Set the state
//
void NetworkSequenceCollection::SetState(NetworkAssemblyState newState)
{
	m_state = newState;
	
	// Reset the checkpoint counter
	m_numReachedCheckpoint = 0;
}

APResult NetworkSequenceCollection::pumpNetwork()
{
	int senderID;
	
	APMessage msg = m_pComm->CheckMessage(senderID);
	if(msg != APM_NONE)
	{
		switch(msg)
		{
			case APM_SEQ:
				{
					parseSeqMessage(senderID);
					return APR_NONE;
				}
			case APM_SEQ_FLAG:
				{
					parseSeqFlagMessage(senderID);
					return APR_NONE;
				}
			case APM_SEQ_EXT:
				{
					parseSeqExtMessage(senderID);
					return APR_NONE;	
				}
			case APM_CONTROL:
				{
					parseControlMessage(senderID);
					return APR_NONE;
				}
			case APM_RESULT:
				{
					//PrintDebug(0, "Result got\n");
					ResultMessage msg = m_pComm->ReceiveResultMessage();
					return msg.result;
				}				
			default:
				printf("Unhandled message code: %d\n", msg);
		}
	}
	else
	{
		//printf("flushing network\n");
		m_pComm->flush();	
	}
	return APR_NONE;
}

//
//
//
bool NetworkSequenceCollection::pumpUntilResult()
{
	while(true)
	{
		APResult result = pumpNetwork();
		if(result == APR_TRUE)
		{
			return true;
		}
		else if(result == APR_FALSE)
		{
			return false;
		}
	}	
}

void NetworkSequenceCollection::parseSeqMessage(int senderID)
{
	SeqMessage seqMsg = m_pComm->ReceiveSeqMessage();
	switch(seqMsg.operation)
	{
		case APO_ADD:
		{
			add(seqMsg.seq);
			break;
		}
		case APO_REMOVE:
		{
			remove(seqMsg.seq);
			break;
		}
		case APO_CHECKSEQ:
		{
			//PrintDebug(0, "Check seq got\n");
			ResultPair result = m_pLocalSpace->exists(seqMsg.seq);
			m_pComm->SendResultMessage(senderID, result.forward || result.reverse);
			break;
		}
		case APO_HAS_PARENT:
		{
			bool result = m_pLocalSpace->hasParent(seqMsg.seq);
			m_pComm->SendResultMessage(senderID, result);			
			break;	
		}
		case APO_HAS_CHILD:
		{
			bool result = m_pLocalSpace->hasChild(seqMsg.seq);
			m_pComm->SendResultMessage(senderID, result);				
			break;	
		}
	}	
}

//
//
//
void NetworkSequenceCollection::parseSeqFlagMessage(int senderID)
{
	SeqFlagMessage flagMsg = m_pComm->ReceiveSeqFlagMessage();
	switch(flagMsg.operation)
	{
		case APSFO_CHECK:
		{
			ResultPair result = m_pLocalSpace->checkFlag(flagMsg.seq, flagMsg.flag);
			m_pComm->SendResultMessage(senderID, (result.forward || result.reverse));
			break;
		}
		case APSFO_SET:
		{
			m_pLocalSpace->setFlag(flagMsg.seq, flagMsg.flag);
			break;
		}
	}	
}

//
//
//
void NetworkSequenceCollection::parseSeqExtMessage(int senderID)
{
	SeqExtMessage extMsg = m_pComm->ReceiveSeqExtMessage();
	switch(extMsg.operation)
	{
		case APSEO_CHECK:
		{
			ResultPair result = m_pLocalSpace->checkExtension(extMsg.seq, extMsg.dir, extMsg.base);
			m_pComm->SendResultMessage(senderID, (result.forward || result.reverse));			
			break;
		}
		case APSEO_SET:
		{
			m_pLocalSpace->setExtension(extMsg.seq, extMsg.dir, extMsg.ext);
			break;
		}
		case APSEO_REMOVE:
		{
			m_pLocalSpace->removeExtension(extMsg.seq, extMsg.dir, extMsg.base);
			break;
		}		
	}	
}

//
//
//
void NetworkSequenceCollection::parseControlMessage(int senderID)
{
	ControlMessage controlMsg = m_pComm->ReceiveControlMessage();
	switch(controlMsg.msgType)
	{
		case APC_LOAD:
		{
			break;	
		}
		case APC_DONELOAD:
		{
			PrintDebug(0, "%d: got done load message\n");
			SetState(NAS_FINALIZE);
			break;	
		}
		case APC_CHECKPOINT:
		{
			m_numReachedCheckpoint++;
			break;	
		}
		case APC_FINISHED:
		{
			SetState(NAS_DONE);
			PrintDebug(0, "finished message received, exiting\n");							
			break;	
		}
		case APC_TRIM:
		{
			SetState(NAS_TRIM);
			break;				
		}
		case APC_ASSEMBLE:
		{
			SetState(NAS_ASSEMBLE);
			break;	
		}
	}
}

/*
//
//
//
NetworkResult NetworkSequenceCollection::pumpNetwork()
{
	int sendID;
	APMessage msg = m_pComm->CheckMessage(sendID);
	bool stop = false;
	
	// a message has been found, process it
	if(msg != APM_NONE)
	{
		switch(msg)
		{
			case APM_SEQADD:
				{
					PackedSeq seq = m_pComm->ReceiveSequence();
					add(seq);
					return NR_NONE;
				}
			case APM_SEQDEL:
				{
					PackedSeq seq = m_pComm->ReceiveSequence();
					remove(seq);
					return NR_NONE;
				}					
			case APM_DONELOAD:
				{
					m_pComm->ClearControlMessage();
					SetState(NAS_FINALIZE);
					return NR_NONE;
				}
			case APM_SEQCHECK:
				{
					PackedSeq seq = m_pComm->ReceiveSequence();
					ResultPair result = m_pLocalSpace->exists(seq);
					m_pComm->SendBool(sendID, result.forward || result.reverse);
					return NR_NONE;
				}				
			case APM_CHKRESPONSE:
				{
					bool found = m_pComm->ReceiveBool(sendID);
					if(found)
					{
						return NR_SEQFOUND;
					}
					else
					{
						return NR_SEQNOTFOUND;
					}
				}
			case APM_CHECKPOINT:
				{
					m_numReachedCheckpoint++;
					m_pComm->ClearControlMessage();
					return NR_NONE;
				}
			case APM_FINISHED:
				{
					m_pComm->ClearControlMessage();
					SetState(NAS_DONE);
					PrintDebug(0, "finished message received, exiting\n");
					return NR_NONE;
				}				
			default:
				printf("Unhandled message code: %d\n", msg);
		}
	}
	return NR_NONE;
}
*/

//
//
//
void NetworkSequenceCollection::removeSequence(const PackedSeq& seq)
{
	// This removes the reverse complement as well
	remove(seq);
	// Remove this sequence as an extension to the adjacent sequences
	for(int i = 0; i <= 1; i++)
	{
		extDirection dir = (i == 0) ? SENSE : ANTISENSE;
		extDirection oppDir = oppositeDirection(dir);	
			
		for(int i = 0; i < NUM_BASES; i++)
		{	
			char currBase = BASES[i];
			// does this sequence have an extension to the deleted seq?
			
			bool hasExt  = checkExtension(seq, dir, currBase);
			if(hasExt)
			{
				PackedSeq tempSeq(seq);	
				// generate the sequence that the extension is to
				char extBase = tempSeq.rotate(dir, currBase);				
				// remove the extension, this removes the reverse complement as well
				removeExtension(tempSeq, oppDir, extBase);
			}
		}
	}	
	
}


//
//
//
void NetworkSequenceCollection::add(const PackedSeq& seq)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->add(seq);
	}
	else
	{

		int nodeID = computeNodeID(seq);		
		m_pComm->SendSeqMessage(nodeID, seq, APO_ADD);
	}
}

//
//
//
void NetworkSequenceCollection::remove(const PackedSeq& seq)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		m_pLocalSpace->remove(seq);
	}
	else
	{
		
		int nodeID = computeNodeID(seq);	
		m_pComm->SendSeqMessage(nodeID, seq, APO_REMOVE);
	}	
}

//
//
//
void NetworkSequenceCollection::finalize()
{
	// this command is broadcast from the controller so we only perform a local finalize
	PrintDebug(1, "loaded: %d sequences\n", m_pLocalSpace->count());	
	m_pLocalSpace->finalize();
}

//
//
//
bool NetworkSequenceCollection::checkForSequence(const PackedSeq& seq)
{
	// Check if this sequence is local
	bool result;
	if(isLocal(seq))
	{
		ResultPair rp = m_pLocalSpace->exists(seq);
		result = (rp.forward || rp.reverse);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		m_pComm->SendSeqMessage(nodeID, seq, APO_CHECKSEQ);
		result = pumpUntilResult();
	}
	return result;
}

//
//
//
bool NetworkSequenceCollection::checkpointReached() const
{
	if(m_numReachedCheckpoint == m_numDataNodes)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//
//
//
void NetworkSequenceCollection::setFlag(const PackedSeq& seq, SeqFlag flag)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->setFlag(seq, flag);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_pComm->SendSeqFlagMessage(nodeID, seq, APSFO_SET, flag);
	}
}

//
//
//
bool NetworkSequenceCollection::checkFlag(const PackedSeq& seq, SeqFlag flag)
{
	// Check if this sequence is local
	bool result;
	if(isLocal(seq))
	{
		ResultPair rp = m_pLocalSpace->checkFlag(seq, flag);
		result = (rp.forward || rp.reverse);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		m_pComm->SendSeqFlagMessage(nodeID, seq, APSFO_CHECK, flag);
		result = pumpUntilResult();
	}
	return result;
}

//
//
// COMBINE THESE TWO FUNCTIONS
//
//

//
//
//
bool NetworkSequenceCollection::hasParent(const PackedSeq& seq)
{
	// Check if this sequence is local
	bool result;
	if(isLocal(seq))
	{
		result = m_pLocalSpace->hasParent(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		m_pComm->SendSeqMessage(nodeID, seq, APO_HAS_PARENT);
		result = pumpUntilResult();
	}
	return result;
}

//
//
//
bool NetworkSequenceCollection::hasChild(const PackedSeq& seq)
{
	// Check if this sequence is local
	bool result;
	if(isLocal(seq))
	{
		result = m_pLocalSpace->hasChild(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		m_pComm->SendSeqMessage(nodeID, seq, APO_HAS_CHILD);
		result = pumpUntilResult();
	}

	return result;
}

//
//
//
void NetworkSequenceCollection::setExtension(const PackedSeq& seq, extDirection dir, SeqExt extension)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->add(seq);
	}
	else
	{

		int nodeID = computeNodeID(seq);		
		// base is ignored for now
		m_pComm->SendSeqExtMessage(nodeID, seq, APSEO_SET, dir, extension);
	}
}

//
//
//
int NetworkSequenceCollection::count() const
{
	return m_pLocalSpace->count();
}

//
//
//
void NetworkSequenceCollection::removeExtension(const PackedSeq& seq, extDirection dir, char base)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->removeExtension(seq, dir, base);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		SeqExt dummy;
		m_pComm->SendSeqExtMessage(nodeID, seq, APSEO_REMOVE, dir, dummy, base);
	}
}

//
//
//
bool NetworkSequenceCollection::checkExtension(const PackedSeq& seq, extDirection dir, char base)
{
	// Check if this sequence is local
	bool result;
	if(isLocal(seq))
	{
		ResultPair rp = m_pLocalSpace->checkExtension(seq, dir, base);
		result = (rp.forward || rp.reverse);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		SeqExt dummy;
		m_pComm->SendSeqExtMessage(nodeID, seq, APSEO_CHECK, dir, dummy, base);
		result = pumpUntilResult();
	}

	return result;
}

//
//
//
void  NetworkSequenceCollection::performTrim(int readLen, int kmerSize)
{
	int start = 2;
	int step = 2;
	int maxBranch =  4 * (readLen - kmerSize + 1);
	
	
	while(start <= maxBranch)
	{
		//trimSequences(pSS, minCoords, maxCoords);
		trimSequences(start);
		//trimSequences3(pSS, minCoords, maxCoords, i+1, 3 * (readLen - kmerSize + 1));
		start <<= 1;
	}
	
	// Now trim at the max branch length
	for(int i = 0; i < 2; i++)
	{
		trimSequences(maxBranch);
	}
	
	// finally, trim at the min contig length
	bool stop = false;
	while(!stop)
	{
		int numRemoved = trimSequences(100);
		if(numRemoved <= 0)
		{
			stop = true;
		}
	}	
}

//
//
//

int NetworkSequenceCollection::trimSequences(int maxBranchCull)
{
	const int MAX_DEAD_LENGTH = maxBranchCull;
	PrintDebug(0, "trimming max branch: %d\n", maxBranchCull);	
	int numBranchesRemoved = 0;
	int count = 0;
	SequenceCollectionIter endIter  = m_pLocalSpace->getEndIter();
	for(SequenceCollectionIter iter = m_pLocalSpace->getStartIter(); iter != endIter; ++iter)
	{
		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
		
		if(count % 10000 == 0)
		{
			PrintDebug(0, "trimmed: %d\n", count);
		}
		count++;
				
		bool child = hasChild(*iter);
		bool parent = hasParent(*iter);
		
		extDirection dir;
		if(!child && !parent)
		{
			// remove this sequence, it has no extensions
			removeSequence(*iter);
			continue;
		}
		else if(!child)
		{
			dir = ANTISENSE;
		}
		else if(!parent)
		{
			dir = SENSE;
		}
		else
		{
			continue;	
		}

		count++;
		PSequenceVector branchElements;
		
		extDirection oppositeDir = oppositeDirection(dir);
					
		PackedSeq currSeq = *iter;
		SeqRecord loopCheck;
		bool stop = false;
		
		while(!stop)
		{		
				
			HitRecord hr = calculateExtension(currSeq, dir);
			HitRecord oppHr = calculateExtension(currSeq, oppositeDir);
			
			if(branchElements.size() == MAX_DEAD_LENGTH + 1 )
			{
				// no ext
				stop = true;
				//printf("stopped because of too long: %d\n", branchElements.size());
			}
			else if(oppHr.getNumHits() > 1)
			{
				//printf("stopped because of reverse branch\n");
				stop = true;	
			}
			else if(hr.getNumHits() == 0 || hr.getNumHits() > 1)
			{
				branchElements.push_back(currSeq);
				//printf("stopped because of noext/ambi branch\n");
				stop = true;
			}
			else
			{
				branchElements.push_back(currSeq);
					
				// good ext
				currSeq = hr.getFirstHit().seq;
			}
		}
		
		//printf("	branch has size: %d\n", branchElements.size());
		if(branchElements.size() <= MAX_DEAD_LENGTH && branchElements.size() > 0)
		{
			//printf("		removing\n");
			numBranchesRemoved++;
			for(PSequenceVectorIterator bIter = branchElements.begin(); bIter != branchElements.end(); bIter++)
			{
				removeSequence(*bIter);
			}
		}
		
		pumpNetwork();
	}
	
	
	
	printf("seqs after trimming: %d\n", this->count());
	printf("num branches removed: %d\n", numBranchesRemoved);
	return numBranchesRemoved;
}

//
//
//
void NetworkSequenceCollection::assemble()
{
	// create file writer
	//FastaWriter writer("contigs.fa");
	int noext = 0;
	int ambiext = 0;
	MPI_File handle;
	MPI_Info info;
	//MPI_File_open(MPI_COMM_WORLD, "pcontigs.fa", MPI_MODE_WRONLY, info, &handle);

	printf("starting assembly\n");
	int count = 0;
	SequenceCollectionIter endIter  = m_pLocalSpace->getEndIter();
	for(SequenceCollectionIter iter = m_pLocalSpace->getStartIter(); iter != endIter; ++iter)
	{
		if(!iter->isFlagSet(SF_SEEN) && !iter->isFlagSet(SF_DELETE))
		{
			// the record of extensions			
			PSequenceVector extensions[2];
							
			for(int i = 0; i <= 1; i++)
			{
				bool stop = false;
				extDirection dir = (i == 0) ? SENSE : ANTISENSE;
				PackedSeq currSeq = *iter;
				SeqRecord loopCheck;			
				
				while(!stop)
				{
					
					// Mark the flag for the selected sequence
					setFlag(currSeq, SF_SEEN);
									
					HitRecord hr = calculateExtension(currSeq, dir);
					if(hr.getNumHits() == 0)
					{
						// no ext
						stop = true;
						noext++;
					}
					else if(hr.getNumHits() == 1)
					{
						// good ext
						currSeq = hr.getFirstHit().seq;
						
						if(loopCheck.contains(currSeq))
						{
							stop = true;
						}
						else
						{
							//printf("good ext (%s)\n", currSeq.decode().c_str());
							extensions[i].push_back(currSeq);
						}
						
					}
					else
					{
						// ambi ext
						stop = true;
						ambiext++;
					}
					
					loopCheck.addSequence(currSeq);
				}
			}
			
			Sequence contig = BuildContig(extensions, *iter);
			
			// is this contig worth outputting?
			if(contig.length() >= 100)
			{
				count++;
				const int bufferSize = 10*1024;
				char buffer[bufferSize];
				int numChars = sprintf(buffer, ">%d\n%s\n", count, contig.c_str());
				printf("%s", buffer);
				//MPI_Status status;
				//MPI_File_write(handle, buffer, numChars, MPI::CHAR, &status);
				//writer.WriteSequence(contig);
			}
		}
		pumpNetwork();	
	}
	//MPI_File_close(&handle);
	printf("noext: %d, ambi: %d\n", noext, ambiext);	
	
}

Sequence NetworkSequenceCollection::BuildContig(PSequenceVector* extensions, const PackedSeq& originalSeq)
{
	Sequence contig;
	contig.reserve(originalSeq.getSequenceLength() + extensions[0].size() + extensions[1].size());
	
	// output the contig
	// output all the antisense extensions
	for(PSequenceVector::reverse_iterator asIter = extensions[1].rbegin(); asIter != extensions[1].rend(); asIter++)
	{
		contig.append(1, asIter->getFirstBase());
	}
	
	// output the current sequence itself
	contig.append(originalSeq.decode());
	
	// output the sense extensions
	for(PSequenceVector::iterator sIter = extensions[0].begin(); sIter != extensions[0].end(); sIter++)
	{
		contig.append(1, sIter->getLastBase());
	}	
	return contig;
}


//
//
//
void NetworkSequenceCollection::generateAdjacency()
{
	PrintDebug(0, "generating adjacency info\n");
	int count = 0;
	SequenceCollectionIter endIter  = m_pLocalSpace->getEndIter();
	for(SequenceCollectionIter iter = m_pLocalSpace->getStartIter(); iter != endIter; ++iter)
	{
		if(count % 10000 == 0)
		{
			printf("generated for %d\n", count);
		}
		count++;
		
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			SeqExt extension;
			for(int j = 0; j < NUM_BASES; j++)
			{
				char currBase = BASES[j];
				PackedSeq testSeq(*iter);
				testSeq.rotate(dir, currBase);
				
				if(checkForSequence(testSeq))
				{
					extension.SetBase(currBase);
				}
			}
			m_pLocalSpace->setExtension(*iter, dir, extension);			
		}
		
		
		//iter->printExtension();
		pumpNetwork();
	}
}

//
//
//
HitRecord NetworkSequenceCollection::calculateExtension(const PackedSeq& currSeq, extDirection dir)
{	
	
	// Create the return structure
	HitRecord hitRecord;
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		
		// SLOW
		bool hasExt = checkExtension(currSeq, dir, currBase);
			
		// Does this sequence have an extension?
		if(hasExt)
		{
			PackedSeq extSeq(currSeq);
			extSeq.rotate(dir, currBase);
			
			// is there a forward extension?
			//if(hasExt.forward)
			if(true)
			{
				hitRecord.addHit(extSeq, false);
			}
			else
			{
				// extension is of the reverse complement
				hitRecord.addHit(extSeq, true);	
			}
		}
	}
		
	return hitRecord;
}

//
//
//
void NetworkSequenceCollection::readSequences(std::string fastaFile, int readLength, int kmerSize)
{
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	int count = 0;
	
	// Read the sequences and add them to the network sequence space
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();		
		assert(kmerSize <= seq.getSequenceLength());
		
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);
			
			// Add the sequence to the network space
			add(sub);
			
			if(count % 10000 == 0)
			{
				printf("sent %d sequences\n", count);
			}
			count++;
		}
	}
}

//
// Check if this sequence belongs in the local space
//
bool NetworkSequenceCollection::isLocal(const PackedSeq& seq) const
{
	int id = computeNodeID(seq);
	return id == m_id;
}

//
// 
//
int NetworkSequenceCollection::computeNodeID(const PackedSeq& seq) const
{
	unsigned int code = seq.getCode();
	unsigned int id = code % m_numDataNodes;
	return id;
}

//
//
//
int NetworkSequenceCollection::PrintDebug(int level,char* fmt, ...) const
{	
	int retval=0;
	if(level <= 3)
	{
		printf("%d: ", m_id);
		va_list ap;
		va_start(ap, fmt); /* Initialize the va_list */
		retval = vprintf(fmt, ap); /* Call vprintf */
		
		va_end(ap); /* Cleanup the va_list */
	}

	return retval;
}

