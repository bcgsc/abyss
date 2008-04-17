#include "NetworkSequenceCollection.h"

//
//
//
NetworkSequenceCollection::NetworkSequenceCollection(int myID, int numDataNodes, int kmerSize, int readLen) : m_id(myID), m_numDataNodes(numDataNodes), m_state(NAS_LOADING), m_kmer(kmerSize), m_readLen(readLen)
{
	// Load the phase space
	m_pLocalSpace = new SequenceCollectionHash();
	
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
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage();				
				break;
			}
			case NAS_GEN_ADJ:
			{
				generateAdjacency(this);
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();
				//SendControlMessageTo 
				break;
			}
			case NAS_TRIM:
			{
				performTrim(this, readLength, kmerSize);
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();				
				break;
			}
			case NAS_ASSEMBLE:
			{
				printf("doing noncontrol assemble\n");
				assemble(this, m_readLen, m_kmer);
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
// The main loop for the controller (rank = 0 process)
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
				loadSequences(this, fastaFile, readLength, kmerSize);
				SetState(NAS_FINALIZE);
				
				// Tell the rest of the loaders that the load is finished
				m_pComm->SendControlMessage(m_numDataNodes, APC_DONELOAD);
				break;
			}
			case NAS_FINALIZE:
			{
				finalize();
				m_numReachedCheckpoint++;
				while(!checkpointReached())
				{
					pumpNetwork();
				}
				
				SetState(NAS_GEN_ADJ);
				m_pComm->SendControlMessage(m_numDataNodes, APC_GEN_ADJ);				
				break;
			}
			case NAS_GEN_ADJ:
			{
				generateAdjacency(this);
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
				performTrim(this, readLength, kmerSize);
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
				assemble(this, m_readLen, m_kmer);
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
			case NAS_WAITING:
			{
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
	m_pComm->flush();	
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
				
			// APM_RESULT is handled in pumpUntilResult()
			case APM_RESULT:
			default:
				{
					return APR_NONE;	
				}
		}
	}
	else
	{
		//printf("flushing network\n");

	}

	return APR_NONE;
}

//
//
//
ResultPair NetworkSequenceCollection::pumpUntilResult()
{
	while(true)
	{
		// Service incoming requests
		pumpNetwork();
		
		// Check if the result has returned
		int senderID;
		m_pComm->flush();	
		APMessage msg = m_pComm->CheckMessage(senderID);
		
		// The result has arrived
		if(msg == APM_RESULT)
		{		
	
			//PrintDebug(0, "Result got\n");
			ResultMessage msg = m_pComm->ReceiveResultMessage();
			
			// Convert the result to a result pair
			ResultPair rp;
			rp.forward = (msg.result[0] == APR_TRUE) ? true : false;
			rp.reverse = (msg.result[1] == APR_TRUE) ? true : false;
			
			// Return the result
			return rp;
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
			bool result = m_pLocalSpace->exists(seqMsg.seq);
			m_pComm->SendResultMessage(senderID, result);
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
			bool result = m_pLocalSpace->checkFlag(flagMsg.seq, flagMsg.flag);
			m_pComm->SendResultMessage(senderID, result);
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
// Parse a sequence extension message
//
void NetworkSequenceCollection::parseSeqExtMessage(int senderID)
{
	SeqExtMessage extMsg = m_pComm->ReceiveSeqExtMessage();
	switch(extMsg.operation)
	{
		case APSEO_CHECK:
		{
			ResultPair result = m_pLocalSpace->checkExtension(extMsg.seq, extMsg.dir, extMsg.base);
			m_pComm->SendResultMessage(senderID, result);			
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
		case APSEO_CLEAR_ALL:
		{
			printf("Clearing extension\n");
			m_pLocalSpace->clearExtensions(extMsg.seq, extMsg.dir);
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
			PrintDebug(0, "got done load message\n");
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
		case APC_GEN_ADJ:
		{
			SetState(NAS_GEN_ADJ);
			break;	
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
		//PrintDebug(3, "done add\n");
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
bool NetworkSequenceCollection::exists(const PackedSeq& seq)
{
	// Check if this sequence is local
	bool result;
	if(isLocal(seq))
	{
		result = m_pLocalSpace->exists(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		m_pComm->SendSeqMessage(nodeID, seq, APO_CHECKSEQ);
		//PrintDebug(1, "after send\n");
		ResultPair rp = pumpUntilResult();
		result = rp.forward || rp.reverse;
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
		result = m_pLocalSpace->checkFlag(seq, flag);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		//PrintDebug(1, "checking for non-local seq %s\n", seq.decode().c_str());
		m_pComm->SendSeqFlagMessage(nodeID, seq, APSFO_CHECK, flag);
		ResultPair rp = pumpUntilResult();
		result = rp.forward || rp.reverse;
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
		ResultPair rp = pumpUntilResult();
		result = rp.forward || rp.reverse;
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
		ResultPair rp = pumpUntilResult();
		result = rp.forward || rp.reverse;
	}

	return result;
}

//
// set the extension for the sequence
//
void NetworkSequenceCollection::setExtension(const PackedSeq& seq, extDirection dir, SeqExt extension)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->setExtension(seq, dir, extension);
	}
	else
	{

		int nodeID = computeNodeID(seq);		
		// base is ignored for now
		m_pComm->SendSeqExtMessage(nodeID, seq, APSEO_SET, dir, extension);
	}
}

//
// clear the extensions for the sequence
//
void NetworkSequenceCollection::clearExtensions(const PackedSeq& seq, extDirection dir)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->clearExtensions(seq, dir);
	}
	else
	{

		int nodeID = computeNodeID(seq);		
		// base is ignored for now
		SeqExt dummyExt;
		char dummyBase = 'A';
		m_pComm->SendSeqExtMessage(nodeID, seq, APSEO_CLEAR_ALL, dir, dummyExt, dummyBase);
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
ResultPair NetworkSequenceCollection::checkExtension(const PackedSeq& seq, extDirection dir, char base)
{
	// Check if this sequence is local
	ResultPair result;
	if(isLocal(seq))
	{
		result = m_pLocalSpace->checkExtension(seq, dir, base);
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
// GetMultiplicity
//
int NetworkSequenceCollection::getMultiplicity(const PackedSeq& seq)
{
	// does nothing for now
	return 0;	
}

//
// The iterator functions return the local space's begin/end
//
SequenceCollectionIterator NetworkSequenceCollection::getStartIter()
{
	return m_pLocalSpace->getStartIter();	
}

//
// The iterator functions return the local space's begin/end
//
SequenceCollectionIterator NetworkSequenceCollection::getEndIter()
{
	return m_pLocalSpace->getEndIter();	
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

