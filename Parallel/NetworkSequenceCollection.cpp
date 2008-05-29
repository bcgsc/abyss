#include "NetworkSequenceCollection.h"
#include "Options.h"
#include <sstream>

//
//
//
NetworkSequenceCollection::NetworkSequenceCollection(int myID, 
											 int numDataNodes, 
											 int kmerSize, 
											 int readLen) : m_id(myID), m_numDataNodes(numDataNodes), 
											 m_kmer(kmerSize), m_readLen(readLen), 
											 m_numBasesAdjSet(0), m_startTrimLen(-1), m_trimStep(0), m_numAssembled(0), m_numOutstandingRequests(0), m_timer("Total")
{
	// Load the phase space
	m_pLocalSpace = new SequenceCollectionHash();
	
	// Create the comm layer
	m_pComm = new CommLayer(myID, kmerSize);
	
	// Create the message buffer
	m_pMsgBuffer = new MessageBuffer(numDataNodes, m_pComm);
	
	stringstream strStrm;
	strStrm << "log_" << myID << ".txt";
	m_pLog = new Log(strStrm.str());
	
}

//
//
//
NetworkSequenceCollection::~NetworkSequenceCollection()
{
	// write the last message to the log
	m_pLog->write(m_timer.toString().c_str());
	
	// Delete the objects created in the constructor
	delete m_pLocalSpace;
	m_pLocalSpace = 0;
	
	delete m_pComm;
	m_pComm = 0;
	
	delete m_pMsgBuffer;
	m_pMsgBuffer = 0;
	
	delete m_pLog;
	m_pLog = 0;
}

void NetworkSequenceCollection::loadSequences()
{
	Timer timer("LoadSequences");
	for (unsigned i = m_id;
			i < opt::inFiles.size();
			i += m_numDataNodes)
		AssemblyAlgorithms::loadSequences(this, opt::inFiles[i]);
	m_pLog->write(timer.toString().c_str());
}

int adjSet = 0;
int numAdjMessageSent = 0;
int numAdjMessageParsed = 0;

//
//
//
void NetworkSequenceCollection::run()
{
	SetState(NAS_LOADING);

	bool stop = false;
	while(!stop)
	{
		switch(m_state)
		{
			case NAS_LOADING:
				loadSequences();
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage();
				break;
			case NAS_FINALIZE:
			{
				finalize();
				
				// Cleanup any messages that are pending
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage();				
				break;
			}
			case NAS_GEN_ADJ:
			{
				networkGenerateAdjacency(this);
				// Cleanup any messages that are pending
				EndState();
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();
				break;
			}
			case NAS_ERODE:
			{
				AssemblyAlgorithms::erodeEnds(this);
				// Cleanup any messages that are pending
				EndState();
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage();
				break;
			}
			case NAS_TRIM:
			case NAS_TRIM2:
			{					
				assert(m_trimStep != 0);
				
				// Perform the trim at the branch length that the control node sent
				int numRemoved = performNetworkTrim(this, m_trimStep);
				// Cleanup any messages that are pending
				EndState();				
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage(numRemoved);
				break;
			}
			case NAS_POPBUBBLE:
			{
				unsigned numPopped
					= performNetworkBubblePop(this, opt::kmerSize);
				m_pComm->SendCheckPointMessage(numPopped);
				// Cleanup any messages that are pending
				EndState();				
				SetState(NAS_WAITING);	
				break;	
			}
			case NAS_SPLIT:
			{
				AssemblyAlgorithms::splitAmbiguous(this);	
				m_pComm->SendCheckPointMessage();
				// Cleanup any messages that are pending
				EndState();				
				SetState(NAS_WAITING);
				break;
			}
			case NAS_ASSEMBLE:
			{
				// The slave node opens the file in append mode
				FastaWriter* writer = new FastaWriter("pcontigs.fa", true);
				unsigned numAssembled = performNetworkAssembly(this, writer);
				
				// Close the writer
				delete writer;
				
				// Cleanup any messages that are pending
				EndState();				
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage(numAssembled);
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
			assert(false);			
		}
	}
}

//
// The main loop for the controller (rank = 0 process)
//
void NetworkSequenceCollection::runControl()
{
	SetState(NAS_LOADING);

	bool stop = false;
	while(!stop)
	{
		switch(m_state)
		{
			case NAS_LOADING:
				loadSequences();
				EndState();
				m_numReachedCheckpoint++;
				while (!checkpointReached(m_numDataNodes))
					pumpNetwork();
				SetState(NAS_FINALIZE);
				m_pComm->SendControlMessage(m_numDataNodes, APC_DONELOAD);
				break;
			case NAS_FINALIZE:
			{
				finalize();
				
				// Cleanup any messages that are pending
				EndState();
								
				m_numReachedCheckpoint++;
				while(!checkpointReached(m_numDataNodes))
				{
					pumpNetwork();
				}
						
				SetState(NAS_GEN_ADJ);
				m_pComm->SendControlMessage(m_numDataNodes, APC_GEN_ADJ);		
				
				//SetState(NAS_DONE);
				//m_pComm->SendControlMessage(m_numDataNodes, APC_FINISHED);		
				break;
			}
			case NAS_GEN_ADJ:
			{
				puts("Generating adjacency");
				networkGenerateAdjacency(this);
				
				// Cleanup any messages that are pending
				EndState();
								
				m_numReachedCheckpoint++;
				while(!checkpointReached(m_numDataNodes))
				{
					pumpNetwork();
				}
				
				// should erosion be performed?
				if(!opt::disableErosion)
				{
					SetState(NAS_ERODE);
				}
				else
				{
					m_startTrimLen = 2;	
					SetState(NAS_TRIM);
				}
				//SetState(NAS_DONE);
				//m_pComm->SendControlMessage(m_numDataNodes, APC_FINISHED);		
				break;
			}
			case NAS_ERODE:
			{
				puts("Eroding");
				int numErodes = opt::readLen - opt::kmerSize + 1;
				for(int i = 0; i < numErodes; i++)
				{
					m_pComm->SendControlMessage(m_numDataNodes, APC_ERODE);
					AssemblyAlgorithms::erodeEnds(this);
					
					// Cleanup any messages that are pending
					EndState();							

					m_numReachedCheckpoint++;
					while(!checkpointReached(m_numDataNodes))
					{
						pumpNetwork();
					}
					
					// All checkpoints are reached, reset the state
					SetState(NAS_ERODE);					
				}
				
				// Cleanup any messages that are pending
				EndState();				
				
				// erosion has been completed
				m_startTrimLen = numErodes + 1;
				SetState(NAS_TRIM);				
				
			}
			case NAS_TRIM:
			{		
				// The control node drives the trimming and passes the value to trim at to the other nodes
				int start = m_startTrimLen;
				
				assert(start > 1);
				
				bool stopTrimming = false;
				while(!stopTrimming)
				{
					printf("Trimming short branches: %d\n", start);	
					m_pComm->SendControlMessage(m_numDataNodes, APC_TRIM, start);
					
					// perform the trim
					int numRemoved = performNetworkTrim(this, start);
					
					if (start < opt::trimLen)
						start <<= 1;
					if (start > opt::trimLen)
						start = opt::trimLen;

					// Wait for all the nodes to hit the checkpoint
					if(numRemoved == 0)
					{
						stopTrimming = true;
					}	
					
					// Cleanup any messages that are pending
					EndState();							

					m_numReachedCheckpoint++;
					while(!checkpointReached(m_numDataNodes))
					{
						pumpNetwork();
					}
					
					// All checkpoints are reached, reset the state
					SetState(NAS_TRIM);
				}
				
				// Cleanup any messages that are pending
				EndState();				
				
				// Trimming has been completed
				SetState(NAS_POPBUBBLE);					
				break;
			}
			case NAS_POPBUBBLE:
				puts("Popping bubbles");	
				while (controlPopBubbles() > 0);
				SetState(NAS_TRIM2);
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_TRIM, opt::trimLen);
				break;
			case NAS_TRIM2:
			{
				printf("Trimming short branches: %d\n", opt::trimLen);
				performNetworkTrim(this, opt::trimLen);
				
				// Cleanup any messages that are pending
				EndState();
								
				m_numReachedCheckpoint++;
				while(!checkpointReached(m_numDataNodes))
				{
					pumpNetwork();
				}
								
				SetState(NAS_SPLIT);
				m_pComm->SendControlMessage(m_numDataNodes, APC_SPLIT);						
			}
			case NAS_SPLIT:
			{
				puts("Splitting ambiguous branches");
				AssemblyAlgorithms::splitAmbiguous(this);

				// Cleanup any messages that are pending
				EndState();
								
				m_numReachedCheckpoint++;
				while(!checkpointReached(m_numDataNodes))
				{
					pumpNetwork();
				}
				
				SetState(NAS_ASSEMBLE);			
				break;				
			}	
			case NAS_ASSEMBLE:
			{
				puts("Assembling");
				// Perform a round-robin assembly
				// The assembly operations cannot be concurrent since a contig can have a start in two different
				// nodes and would therefore result in a collision (and it would be output twice)
				// Using a round-robin assembly like this will have a minimal impact on performance since most of the heavy
				// computation is already done and the output of the contigs will mostly be bounded by file io

				// First, assemble the local sequences
				
				// Note: all other nodes will be in a waiting state so they will service network requests
				
				// The master opens the file in truncate mode
				FastaWriter* writer = new FastaWriter("pcontigs.fa");
				unsigned numAssembled = performNetworkAssembly(this, writer);
				
				// Close the writer
				delete writer;

				// Cleanup any messages that are pending
				EndState();
								
				// Now tell all the slave nodes to perform the assemble one by one
				for(unsigned int i = 1; i < m_numDataNodes; ++i)
				{
					m_pComm->SendControlMessageToNode(i,
							APC_ASSEMBLE, numAssembled);
					
					// Wait for this node to return
					int slaveNumAssembled = 0;
					while(!checkpointReached(1))
					{
						pumpNetwork(&slaveNumAssembled);
					}
					numAssembled += slaveNumAssembled;
					
					// Cleanup any messages that are pending
					EndState();
					
					//Reset the state and loop
					SetState(NAS_ASSEMBLE);
				}
				
				// Cleanup any messages that are pending
				EndState();
				
				SetState(NAS_DONE);
				m_pComm->SendControlMessage(m_numDataNodes, APC_FINISHED);				
				break;
			}
			case NAS_DONE:
			{
				puts("Done.");
				stop = true;
				break;
			}
			case NAS_WAITING:
			{
				break;
			}
			assert(false);							
		}
	}
}

void NetworkSequenceCollection::EndState()
{
	// Flush the message buffer
	m_pMsgBuffer->flush();
}

//
// Set the state
//
void NetworkSequenceCollection::SetState(NetworkAssemblyState newState)
{
	// Ensure there are no pending messages
	assert(m_pMsgBuffer->empty());
	
	m_state = newState;
	
	// Reset the checkpoint counter
	m_numReachedCheckpoint = 0;
}

APResult NetworkSequenceCollection::pumpNetwork(int* pArg)
{
	int senderID;
	m_pComm->flush();	
	//int numReceived = 0;
	bool stop = false;
	
	// Loop until either a) a APM_RESULT message is found in which case return and allow pump_until_result to handle it or b) no more messages are waiting
	while(!stop)
	{
		
		APMessage msg = m_pComm->CheckMessage(senderID);
		
		switch(msg)
		{
			case APM_CONTROL:
				{
					int arg = parseControlMessage();
					if (pArg != NULL)
						*pArg = arg;
					break;
				}
			case APM_BUFFERED:
				{
					MessagePtrVector msgs;
					m_pComm->ReceiveBufferedMessage(msgs);
					for(MessagePtrVector::iterator iter = msgs.begin(); iter != msgs.end(); iter++)
					{
						// Handle each message based on its type
						(*iter)->handle(senderID, *this);
						
						// Delete the message
						delete (*iter);
						*iter = 0;
					}
								
					break;	
				}
			case APM_NONE:
				{
					// either there is no message waiting or it is an AP_RESULT in which case pump_until_result will handle it
					stop = true;
					break;
				}
		}
	}
	return APR_NONE;
}


//
//
//
void NetworkSequenceCollection::handleSeqOpMessage(int /*senderID*/, const SeqOpMessage& seqMsg)
{
	switch(seqMsg.m_operation)
	{
		case MO_ADD:
		{
			add(seqMsg.m_seq);
			break;
		}
		case MO_REMOVE:
		{
			remove(seqMsg.m_seq);
			break;
		}
		case MO_EXIST:
		{
			assert(false);
			break;
		}
		default:
		{
			assert(false);
			break;
		}	
	}
}

//
//
//
void NetworkSequenceCollection::handleSetFlagMessage(int /*senderID*/, const SetFlagMessage& message)
{
	assert(isLocal(message.m_seq));
	m_pLocalSpace->setFlag(message.m_seq, message.m_flag);
}

void NetworkSequenceCollection::handleSetBaseMessage(int /*senderID*/, const SetBaseMessage& message)
{
	numAdjMessageParsed++;
	assert(isLocal(message.m_seq));
	setBaseExtension(message.m_seq, message.m_dir, message.m_base);
}

//
// Parse a sequence extension message
//
void NetworkSequenceCollection::handleRemoveExtensionMessage(int /*senderID*/, const RemoveExtensionMessage& message)
{
	assert(isLocal(message.m_seq));	
	m_pLocalSpace->removeExtension(message.m_seq, message.m_dir, message.m_base);
}

//
//
//
int NetworkSequenceCollection::parseControlMessage()
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
			break;	
		}
		case APC_TRIM:
		{
			// This message came from the control node along with an argument indicating the maximum branch to trim at
			m_trimStep = controlMsg.argument;
			SetState(NAS_TRIM);
			break;				
		}
		case APC_ERODE:
		{
			SetState(NAS_ERODE);
			break;	
		}
		case APC_POPBUBBLE:
		{		
			SetState(NAS_POPBUBBLE);
			break;	
		}
		case APC_SPLIT:
		{
			SetState(NAS_SPLIT);
			break;	
		}
		case APC_ASSEMBLE:
		{
			// This message came from the control node along with an argument indicating the number of sequences assembled so far
			m_numAssembled = controlMsg.argument;			
			SetState(NAS_ASSEMBLE);
			break;	
		}
		case APC_GEN_ADJ:
		{
			SetState(NAS_GEN_ADJ);
			break;	
		}		
	}
	return controlMsg.argument;
}

//
// Parse an extension request
//
void NetworkSequenceCollection::handleSequenceDataRequest(int senderID, SeqDataRequest& message)
{	
	// Get the extension for this sequence
	assert(isLocal(message.m_seq));
	
	ExtensionRecord extRec;
	int multiplicity = -1;
	bool found = m_pLocalSpace->getSeqData(message.m_seq, extRec, multiplicity);
	assert(found);
	(void)found;

	// Return the extension to the sender
	m_pMsgBuffer->sendSeqDataResponse(senderID, message.m_group, message.m_id, message.m_seq, extRec, multiplicity);
}

//
// Parse an extension request
//
void NetworkSequenceCollection::handleSequenceDataResponse(int /*senderID*/, SeqDataResponse& message)
{
	processSequenceExtension(message.m_group, message.m_id, message.m_seq, message.m_extRecord, message.m_multiplicity);
}

//
// Generate adjacency - Network version
// This function operates in the same manner as AssemblyAlgorithms::GenerateAdjacency but has been rewritten to hide latency between nodes
//
void NetworkSequenceCollection::networkGenerateAdjacency(ISequenceCollection* seqCollection)
{
	Timer timer("NetworkGenerateAdjacency");
	
	int count = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		count++;
		if(count % 100000 == 0) {
			PrintDebug(0, "Generating adjacency: %d sequences\n",
					count);
			if (m_numOutstandingRequests > 0)
				PrintDebug(0, "Back log of %zu requests\n",
						m_numOutstandingRequests);
		}
		
		const PackedSeq& currSeq = *iter;
		//printf("gen for: %s\n", iter->decode().c_str());
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			extDirection oppDir = !dir;
			SeqExt extension;
			
			PackedSeq testSeq(currSeq);
			char adjBase = testSeq.rotate(dir, 'A');
			
			for(int j = 0; j < NUM_BASES; j++)
			{
				char currBase = BASES[j];
				testSeq.setLastBase(dir, currBase);
				
				// Here is the divergence from the common adjacency generation function
				// We only generate a request for the existance of the sequence at this moment and then carry on
				// When the data meanders over the network and eventually returns to use, THEN the adjacency is set
				// See:: SendAdjancencyRequest/ParseAdjancencyResponse
				//computeAdjacency(currSeq, testSeq, dir, currBase);
				
				// Optimistically send a message over the network that there is an extension from testSeq to adjBase
				setBaseExtension(testSeq, oppDir, adjBase);
				
				pumpNetwork();
			}	
		}
	}
	
	// Wait for all the requests to be filled
	while(m_numOutstandingRequests != 0)
	{
		pumpNetwork();
	}

	m_pLog->write(timer.toString().c_str());
}

//
// Distributed trimming function
// 
int NetworkSequenceCollection::performNetworkTrim(ISequenceCollection* seqCollection, int maxBranchCull)
{
	Timer timer("NetworkTrim");
	int numBranchesRemoved = 0;

	// The branch ids
	uint64_t branchGroupID = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{	
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// remove this sequence, it has no extensions
			AssemblyAlgorithms::removeSequenceAndExtensions(seqCollection, *iter);
		}
		
		// Sequence is trimmable, create a new branch for it
		BranchGroup newGroup(branchGroupID, dir, 1);
		BranchRecord newBranch(dir, maxBranchCull);
		newGroup.addBranch(0, newBranch);
		m_activeBranchGroups[branchGroupID] = newGroup;

		// Generate the first extension request
		generateExtensionRequest(branchGroupID, 0, *iter);
		branchGroupID++;
		
		// Process the active branches
		numBranchesRemoved += processBranchesTrim();
		
		// Service any waiting network events
		seqCollection->pumpNetwork();
		
		// Primitive load balancing
		if(m_activeBranchGroups.size() > MAX_ACTIVE)
		{
			while(m_activeBranchGroups.size() > LOW_ACTIVE)
			{
				seqCollection->pumpNetwork();
				numBranchesRemoved += processBranchesTrim();
			}
		}
	}
	
	// Clear out the remaining branches
	while(!m_activeBranchGroups.empty())
	{
		numBranchesRemoved += processBranchesTrim();
		seqCollection->pumpNetwork();
	}		
	
	m_pLog->write(timer.toString().c_str());
	PrintDebug(0, "Trimmed %d branches\n", numBranchesRemoved);
	return numBranchesRemoved;	
}

//
// Process current branches, removing those that are finished
// returns true if the branch list has branches remaining
//
int NetworkSequenceCollection::processBranchesTrim()
{
	int numBranchesRemoved = 0;
	std::vector<BranchGroupMap::iterator> removeBranches;
	// Check if any of the current branches have gone inactive
	for(BranchGroupMap::iterator iter = m_activeBranchGroups.begin(); iter != m_activeBranchGroups.end(); iter++)
	{
		if(!iter->second.isActive())
		{
			// In the trimming context, the group should have 1 and only 1 branch
			assert(iter->second.getNumBranches() == 1);
			
			// Trim the branch if possible
			
			// Get lastBranch returns the only branch in the group (see assert above)
			if(AssemblyAlgorithms::processTerminatedBranchTrim(this, iter->second.getBranch(0)))
			{
				numBranchesRemoved++;
			}
			
			// Mark the group for removal
			removeBranches.push_back(iter);
		}	
	}
	
	// Remove all the finished branches
	for(std::vector<BranchGroupMap::iterator>::iterator rmIter = removeBranches.begin(); rmIter != removeBranches.end(); rmIter++)
	{
		//printf("erased branch %llu\n", (*rmIter)->first);
		m_activeBranchGroups.erase(*rmIter);	
	}	
	
	return numBranchesRemoved;
}

//
// pop bubbles
//
int NetworkSequenceCollection::performNetworkBubblePop(ISequenceCollection* seqCollection, int kmerSize)
{
	Timer timer("NetworkPopBubbles");
	
	// The branch ids
	uint64_t branchGroupID = 0;
	
	// make sure the branch group structure is initially empty
	assert(m_activeBranchGroups.empty());
	
	int numPopped = 0;
	
	int count = 0;

	// Set the cutoffs
	const unsigned int expectedBubbleSize = 2*(kmerSize + 1);
	const unsigned int maxNumBranches = 3;
	
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		// Skip sequences that have already been deleted	
		if(iter->isFlagSet(SF_DELETE))
		{		
			continue;
		}

		count++;
		if (count % 100000 == 0)
			PrintDebug(0, "Popping bubbles: %d sequences\n", count);
		
		// Get the extensions for this sequence, this function populates the extRecord structure
		ExtensionRecord extRec;
		int multiplicity = -1;
		// THIS CALL TO GET EXTENSIONS IS GUARENTEED TO BE LOCAL SO WE DO NOT HAVE TO WAIT FOR THE RETURN
		bool success = seqCollection->getSeqData(*iter, extRec, multiplicity);
		assert(success);
		(void)success;
		
		// Check for ambiguity
		for(int i = 0; i <= 1; ++i)
		{	
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			
			if(extRec.dir[dir].IsAmbiguous())
			{
				// Found a potential bubble, examine each branch
				
				// Create the branch group
				BranchGroup branchGroup(branchGroupID, dir, maxNumBranches);
				
				// insert the new group into the active group map
				BranchGroupMap::iterator groupIter = m_activeBranchGroups.insert(std::pair<uint64_t, BranchGroup>(branchGroupID,branchGroup)).first;
				
				// initiate the new group
				AssemblyAlgorithms::initiateBranchGroup(groupIter->second, *iter, extRec.dir[dir], multiplicity, expectedBubbleSize);
				
				// generate a sequence extension request for each sequence in the group
				// this will be handled in process sequence extension
				size_t maxID = groupIter->second.getNumBranches();
				for(size_t id = 0; id < maxID; ++id)
				{
					generateExtensionRequest(groupIter->first, id, groupIter->second.getBranch(id).getLastSeq());
				}
				
				// increment the group id
				branchGroupID++;
			}
		}

		// Process groups that may be finished	
		numPopped += processBranchesPop();
		
		seqCollection->pumpNetwork();
	}
	
	// Clear out the remaining branches
	while(!m_activeBranchGroups.empty())
	{
		numPopped += processBranchesPop();
		seqCollection->pumpNetwork();
	}
	
	m_pLog->write(timer.toString().c_str());
	PrintDebug(1, "Removed %d bubbles\n", numPopped);
	return numPopped;	
}

//
// Process groups that are finished searching for bubbles
//
int NetworkSequenceCollection::processBranchesPop()
{
	int numPopped = 0;
	std::vector<BranchGroupMap::iterator> removeBranches;
	// Check if any of the current branches have gone inactive
	for(BranchGroupMap::iterator iter = m_activeBranchGroups.begin(); iter != m_activeBranchGroups.end(); iter++)
	{
		bool remove = false;
		// All branches have been extended one sequence, check the stop conditions
		
		// First, check if the group hit a no-extension in one of the sequences
		if(iter->second.isNoExt())
		{
			remove = true;
		}
		else
		{
			// Update status is called in processSequenceExtensionPop(), check the status here
			BranchGroupStatus status = iter->second.getStatus();
			
			// Check if a stop condition was met
			if(status == BGS_TOOLONG || status == BGS_LOOPFOUND || status == BGS_TOOMANYBRANCHES || status == BGS_NOEXT)
			{
				remove = true;
			}
			else if(status == BGS_JOINED)
			{
				AssemblyAlgorithms::collapseJoinedBranches(this, iter->second);
				numPopped++;
				remove = true;
			}
			else
			{										
				// the branch is still active, continue
				assert(status == BGS_ACTIVE);
			}
		}
		
		if(remove)
		{
			// Mark the group for removal
			removeBranches.push_back(iter);
		}	
	}
	
	// Remove all the finished branches
	for(std::vector<BranchGroupMap::iterator>::iterator rmIter = removeBranches.begin(); rmIter != removeBranches.end(); rmIter++)
	{
		//printf("erased branch %llu\n", (*rmIter)->first);
		m_finishedGroups.insert((*rmIter)->first);
		m_activeBranchGroups.erase(*rmIter);	
	}	
	
	return numPopped;		
}

int NetworkSequenceCollection::controlPopBubbles()
{
	// Perform a round-robin bubble pop to avoid concurrency issues
	unsigned numPopped = performNetworkBubblePop(this, opt::kmerSize);

	// Now tell all the slave nodes to perform the pop one by one
	for(unsigned i = 1; i < m_numDataNodes; ++i) {
		m_pComm->SendControlMessageToNode(i, APC_POPBUBBLE);

		// Cleanup any messages that are pending
		EndState();

		// Wait for this node to return
		int slaveNumPopped;
		while (!checkpointReached(1))
			pumpNetwork(&slaveNumPopped);
		numPopped += slaveNumPopped;

		//Reset the state and loop
		SetState(NAS_POPBUBBLE);
	}

	// Cleanup any messages that are pending
	EndState();

	printf("Removed %d bubbles\n", numPopped);
	return numPopped;
}

//
// Perform a network assembly
//

//
// Distributed trimming function
// 
unsigned NetworkSequenceCollection::performNetworkAssembly(ISequenceCollection* seqCollection, IFileWriter* fileWriter)
{
	Timer timer("NetworkAssembly");
	
	unsigned numAssembled = 0;
	
	// The branch ids
	uint64_t branchGroupID = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{	
		extDirection dir;
		// dir will be set to the assembly direction if the sequence can be assembled
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// singleton, ignore for now
			continue;
		}
		
		// Sequence is trimmable, create a new branch for it
		BranchGroup newGroup(branchGroupID, dir, 1);
		BranchRecord newBranch(dir, -1);
		newGroup.addBranch(0, newBranch);
		m_activeBranchGroups[branchGroupID] = newGroup;

		// Generate the first extension request
		generateExtensionRequest(branchGroupID, 0, *iter);
		branchGroupID++;
		
		// Process the active branches
		numAssembled += processBranchesAssembly(seqCollection, fileWriter, numAssembled);
		
		// Service any waiting network events
		seqCollection->pumpNetwork();
		
		// Primitive load balancing
		if(m_activeBranchGroups.size() > MAX_ACTIVE)
		{
			while(m_activeBranchGroups.size() > LOW_ACTIVE)
			{
				seqCollection->pumpNetwork();
				numAssembled += processBranchesAssembly(seqCollection, fileWriter, numAssembled);
			}
		}
	}
	
	// Clear out the remaining branches
	while(!m_activeBranchGroups.empty())
	{
		numAssembled += processBranchesAssembly(seqCollection, fileWriter, numAssembled);
		seqCollection->pumpNetwork();
	}		
	m_pLog->write(timer.toString().c_str());
	PrintDebug(0, "Assembled %d contigs\n", numAssembled);
	return numAssembled;
}

//
// Process current branches, removing those that are finished
// returns true if the branch list has branches remaining
//
int NetworkSequenceCollection::processBranchesAssembly(ISequenceCollection* seqCollection, IFileWriter* fileWriter, int currContigID)
{
	int numAssembled = 0;
	std::vector<BranchGroupMap::iterator> removeBranches;
	// Check if any of the current branches have gone inactive
	for(BranchGroupMap::iterator iter = m_activeBranchGroups.begin(); iter != m_activeBranchGroups.end(); iter++)
	{
		if(!iter->second.isActive())
		{
			// In this context, the group should have 1 and only 1 branch
			assert(iter->second.getNumBranches() == 1);
			
			// check if the branch is redundant, assemble if so, else it will simply be removed
			if(!isBranchRedundant(iter->second.getBranch(0)))
			{
				// Assemble the contig
				Sequence contig;
				AssemblyAlgorithms::processTerminatedBranchAssemble(seqCollection, iter->second.getBranch(0), contig);
				numAssembled++;
				
				// Output the contig
				fileWriter->WriteSequence(contig,
						m_numAssembled + currContigID++, 0);
			}
			
			// Mark the group for removal
			removeBranches.push_back(iter);
		}	
	}
	
	// Remove all the finished branches
	for(std::vector<BranchGroupMap::iterator>::iterator rmIter = removeBranches.begin(); rmIter != removeBranches.end(); rmIter++)
	{
		//printf("erased branch %llu\n", (*rmIter)->first);
		m_activeBranchGroups.erase(*rmIter);	
	}	
	
	return numAssembled;
}

// Check if a branch is redundant with a previously output branch
bool NetworkSequenceCollection::isBranchRedundant(BranchRecord& branch)
{
	// Since branches are assembled simulatenously it is possibly to start a branch from both ends at the same time
	// This can only happen if both ends of the branch are in the same node (since we assemble one node at a time)
	// To get around that, check 1) if both ends are local and 2) if either end has been seen
	// If so, the branch will be discarded
	
	// Since flag sets are atomic within a node, this logic stands up
	
	if(branch.empty())
	{
		// empty branches are trivally non-redundant (.....or are they?)
		return false;
	}
	
	// Check if the first and last sequences are in the local node
	// note: the first node always should be!
	if(isLocal(branch.getFirstSeq()) && isLocal(branch.getLastSeq()))
	{
		// Check if either sequences have the seen flag set
		// note: this flag returns immediately since both sequenes are local by definition
		if(checkFlag(branch.getFirstSeq(), SF_SEEN) || checkFlag(branch.getLastSeq(), SF_SEEN))
		{
			return true;
		}
		else
		{
			return false;	
		}
	}
	else
	{
		// Unless both sequences are local, the branch is not redundant
		return false;	
	}
}


//
// Generate a request for a sequence's extension, it will be handled in parseSequenceExtensionResponse
//
void NetworkSequenceCollection::generateExtensionRequest(uint64_t groupID, uint64_t branchID, const PackedSeq& seq)
{
	// Check if the test sequence is local
	if(isLocal(seq))
	{
		// simply look up the sequence in the local space
		ExtensionRecord extRec;
		int multiplicity = -1;
		bool success = m_pLocalSpace->getSeqData(seq, extRec, multiplicity);
		assert(success);
		(void)success;
		
		// process the message
		processSequenceExtension(groupID, branchID, seq, extRec, multiplicity);
	}
	else
	{
		// Send the request
		int nodeID = computeNodeID(seq);
		assert(nodeID != m_id);
		
		// Send the request, it will be processed in the callback
		m_pMsgBuffer->sendSeqDataRequest(nodeID, groupID, branchID, seq);
	}
}

//
//
//
void NetworkSequenceCollection::processSequenceExtension(uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec, int multiplicity)
{
	switch(m_state)
	{
		case NAS_TRIM:
		case NAS_TRIM2:
		case NAS_ASSEMBLE:
			return processLinearSequenceExtension(groupID, branchID, seq, extRec, multiplicity);
			break;
		case NAS_POPBUBBLE:
			return processSequenceExtensionPop(groupID, branchID, seq, extRec, multiplicity);
			break;
		default:
			assert(false);
			break;
	}	
}

//
// Process a sequence extension for trimming
//
void NetworkSequenceCollection::processLinearSequenceExtension(uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec, int multiplicity)
{
	//printf("processing %llu\n", id);
	// Find the branch by its ID	
	BranchGroupMap::iterator iter = m_activeBranchGroups.find(groupID);
	
	// should always exist
	assert(iter != m_activeBranchGroups.end());

	PackedSeq currSeq = seq;
	bool active = AssemblyAlgorithms::processLinearExtensionForBranch(iter->second.getBranch(branchID), currSeq, extRec, multiplicity);
	
	// if the branch is still active generate a new request
	if(active)
	{
		return generateExtensionRequest(groupID, branchID, currSeq);
	}
	else
	{
		return;	
	}
}

//
// Process a sequence extension for popping
//
void NetworkSequenceCollection::processSequenceExtensionPop(uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec, int multiplicity)
{
	//printf("processing %llu\n", id);
	// Find the branch by its ID	
	BranchGroupMap::iterator iter = m_activeBranchGroups.find(groupID);
		
	// If the iterator was not found we finished with that branch already, ensure this is so
	if(iter == m_activeBranchGroups.end())
	{
		assert(m_finishedGroups.find(groupID) != m_finishedGroups.end());
		// do nothing
		return;
	}
	
	PackedSeq currSeq = seq;
	bool extendable = AssemblyAlgorithms::processBranchGroupExtension(iter->second, branchID, currSeq, extRec, multiplicity);

	// The extendable flag indicates that one round of extension has happened and each branch is equal length
	if(extendable)
	{
		// Update the status of the branch
		BranchGroupStatus status = iter->second.updateStatus();
			
		// if the group is still active, generate new requests for all the branches in the group
		if(status == BGS_ACTIVE)
		{
			size_t numBranches = iter->second.getNumBranches();
			for(size_t i = 0; i < numBranches; ++i)
			{
				generateExtensionRequest(groupID, i, iter->second.getBranch(i).getLastSeq());
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
		//PrintDebug(3, "done add\n");
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_pMsgBuffer->sendSeqOpMessage(nodeID, seq, MO_ADD);			
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
		m_pMsgBuffer->sendSeqOpMessage(nodeID, seq, MO_REMOVE);
	}	
}

//
//
//
void NetworkSequenceCollection::finalize()
{
	// this command is broadcast from the controller so we only perform a local finalize
	PrintDebug(1, "Loaded %d sequences\n", m_pLocalSpace->count());	
	m_pLocalSpace->finalize();
}

//
//
//
bool NetworkSequenceCollection::exists(const PackedSeq& seq)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		return m_pLocalSpace->exists(seq);
	}
	else
	{
		assert(false);
		return false;
		/*
		//PrintDebug(1, "after send\n");
		ResultPair rp = pumpUntilResult();
		return rp.forward || rp.reverse;
		*/
	}
}

//
//
//
bool NetworkSequenceCollection::checkpointReached(int numRequired) const
{
	if(m_numReachedCheckpoint == numRequired)
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
		m_pMsgBuffer->sendSetFlagMessage(nodeID, seq, flag);
	}
}

//
//
//
bool NetworkSequenceCollection::checkFlag(const PackedSeq& seq, SeqFlag flag)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		return m_pLocalSpace->checkFlag(seq, flag);
	}
	else
	{
		// Check flag should be for local sequences only
		assert(false);
		return false;
	}
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
	if(isLocal(seq))
	{
		return m_pLocalSpace->hasParent(seq);
	}
	else
	{
		// Never should be called for non-local sequences
		assert(false);
		return false;
	}
}

//
//
//
bool NetworkSequenceCollection::hasChild(const PackedSeq& seq)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		return m_pLocalSpace->hasChild(seq);
	}
	else
	{
		// Never should be called for non-local sequences
		assert(false);
		return false;
	}
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
		assert(false);
	}
}

//
// 
//
bool NetworkSequenceCollection::setBaseExtension(const PackedSeq& seq, extDirection dir, char base)
{
	if(isLocal(seq))
	{
		if(m_pLocalSpace->setBaseExtension(seq, dir, base))
		{
			m_numBasesAdjSet++;
		}
	}
	else
	{	
		numAdjMessageSent++;
		int nodeID = computeNodeID(seq);
		m_pMsgBuffer->sendSetBaseExtension(nodeID, seq, dir, base);
	}
	
	// As this call delegates, the return value is meaningless so return false
	return false;
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
		assert(false);
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
		m_pMsgBuffer->sendRemoveExtension(nodeID, seq, dir, base);
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
		//PrintDebug(1, "checking for local seq %s\n", seq.decode().c_str());
		result = m_pLocalSpace->checkExtension(seq, dir, base);
	}
	else
	{
		assert(false);
	}

	return result;
}

//
// get the extensions of the sequence
//
bool NetworkSequenceCollection::getSeqData(const PackedSeq& seq, ExtensionRecord& extRecord, int& multiplicity)
{
	// This function can only be called locally, the distributed version is through generateSequenceExtensionMessage
	if(isLocal(seq))
	{
		m_pLocalSpace->getSeqData(seq, extRecord, multiplicity);
		return true;
	}
	else
	{
		assert(false);
		return false;
	}
}


//
// GetMultiplicity
//
int NetworkSequenceCollection::getMultiplicity(const PackedSeq& /*seq*/)
{
	// does nothing for now
	return 0;	
}

//
//
//
void NetworkSequenceCollection::wipeFlag(SeqFlag flag)
{
	m_pLocalSpace->wipeFlag(flag);
}

//
// The iterator functions return the local space's begin/end
//
SequenceCollectionIterator NetworkSequenceCollection::getStartIter() const
{
	return m_pLocalSpace->getStartIter();	
}

//
// The iterator functions return the local space's begin/end
//
SequenceCollectionIterator NetworkSequenceCollection::getEndIter() const
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
