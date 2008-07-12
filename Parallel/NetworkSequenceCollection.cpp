#include "NetworkSequenceCollection.h"
#include "Options.h"
#include <sstream>
#include <iostream>

//
//
//
NetworkSequenceCollection::NetworkSequenceCollection(int myID, 
		int numDataNodes, int kmerSize, int readLen) :
	m_id(myID),
	m_numDataNodes(numDataNodes),
	m_kmer(kmerSize),
	m_readLen(readLen),
	m_numBasesAdjSet(0),
	m_startTrimLen(-1),
	m_trimStep(0),
	m_numPopped(0),
	m_numAssembled(0),
	m_numOutstandingRequests(0),
	m_timer("Total")
{
	// Load the phase space
	m_pLocalSpace = new SequenceCollectionHash();
	
	// Create the comm layer
	m_pComm = new CommLayer(myID, kmerSize);
	
	// Create the message buffer
	m_pMsgBuffer = new MessageBuffer(numDataNodes, m_pComm);
	m_pComm->setMsgBuffer(m_pMsgBuffer);
	
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
				unsigned numEroded
					= AssemblyAlgorithms::erodeEnds(this);
				// Cleanup any messages that are pending
				EndState();
				SetState(NAS_WAITING);
				
				// Tell the control process this checkpoint has been reached
				m_pComm->SendCheckPointMessage(numEroded);
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
			case NAS_DISCOVER_BUBBLES:
			{
				unsigned numDiscovered
					= performNetworkDiscoverBubbles(this,
							opt::kmerSize);
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage(numDiscovered);
				break;
			}
			case NAS_POPBUBBLE:
			{
				unsigned numPopped
					= performNetworkPopBubbles(this);
				EndState();				
				SetState(NAS_WAITING);	
				m_pComm->SendCheckPointMessage(numPopped);
				break;	
			}
			case NAS_SPLIT:
			{
				AssemblyAlgorithms::splitAmbiguous(this);	
				EndState();				
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage();
				break;
			}
			case NAS_ASSEMBLE:
			{
				FastaWriter* writer = new FastaWriter(
						opt::contigsPath.c_str());
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
				
				if (opt::erode > 0) {
					SetState(NAS_ERODE);
				} else {
					m_startTrimLen = 2;	
					SetState(NAS_TRIM);
				}
				break;
			}
			case NAS_ERODE:
			{
				puts("Eroding");
				unsigned totalEroded = 0;
				for (int i = 0; i < opt::erode; i++) {
					m_pComm->SendControlMessage(m_numDataNodes, APC_ERODE);
					unsigned numEroded 
						= AssemblyAlgorithms::erodeEnds(this);
					
					// Cleanup any messages that are pending
					EndState();							

					m_numReachedCheckpoint++;
					while(!checkpointReached(m_numDataNodes))
					{
						pumpNetwork();
					}
					numEroded += m_checkpointSum;
					printf("Eroded %d tips\n", numEroded);
					totalEroded += numEroded;
					
					// All checkpoints are reached, reset the state
					SetState(NAS_ERODE);					
				}
				printf("Eroded %d tips in total\n",
						totalEroded);
				
				// Cleanup any messages that are pending
				EndState();				
				
				// erosion has been completed
				m_startTrimLen = opt::erode + 1;
				SetState(NAS_TRIM);				
				
			}
			case NAS_TRIM:
			{		
				// The control node drives the trimming and passes the value to trim at to the other nodes
				int start = m_startTrimLen;
				
				assert(start > 1);
				
				unsigned totalRemoved = 0;
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

					// Cleanup any messages that are pending
					EndState();							

					// Wait for all the nodes to hit the checkpoint
					m_numReachedCheckpoint++;
					while(!checkpointReached(m_numDataNodes))
					{
						pumpNetwork();
					}
					numRemoved += m_checkpointSum;
					
					if (numRemoved == 0)
						stopTrimming = true;
					else
						printf("Trimmed %d branches\n", numRemoved);
					totalRemoved += numRemoved;
					
					// All checkpoints are reached, reset the state
					SetState(NAS_TRIM);
				}
				printf("Trimmed %d branches in total\n",
						totalRemoved);
				
				// Cleanup any messages that are pending
				EndState();				
				
				// Trimming has been completed
				SetState(NAS_POPBUBBLE);					
				break;
			}
			case NAS_DISCOVER_BUBBLES:
				// This state is only used by the slaves.
				assert(false);
				exit(EXIT_FAILURE);
			case NAS_POPBUBBLE:
			{
				puts("Popping bubbles");
				unsigned totalPopped = 0;
				int i;
				for (i = 0; i < opt::bubbles; i++) {
					unsigned numPopped = controlPopBubbles();
					if (numPopped == 0)
						break;
					totalPopped += numPopped;
				}
				assert(totalPopped == m_numPopped);
				printf("Removed %d bubbles in %d rounds\n",
						totalPopped, i);
				SetState(NAS_TRIM2);
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_TRIM, opt::trimLen);
				break;
			}
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
				FastaWriter* writer = new FastaWriter(
						opt::contigsPath.c_str());
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_ASSEMBLE);
				unsigned numAssembled = performNetworkAssembly(this,
						writer);
				delete writer;
				EndState();

				m_numReachedCheckpoint++;
				while (!checkpointReached(m_numDataNodes))
					pumpNetwork();
				numAssembled += m_checkpointSum;
				printf("Assembled %u contigs\n", numAssembled);

				SetState(NAS_DONE);
				m_pComm->SendControlMessage(m_numDataNodes, APC_FINISHED);				
				break;
			}
			case NAS_DONE:
			{
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
	m_checkpointSum = 0;
}

APResult NetworkSequenceCollection::pumpNetwork()
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
					parseControlMessage();
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
void NetworkSequenceCollection::parseControlMessage()
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
			m_checkpointSum += controlMsg.argument;
			break;	
		}
		case APC_BARRIER:
		{
			assert(m_state == NAS_WAITING);
			m_pComm->barrier();
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
		case APC_DISCOVER_BUBBLES:
		{
			SetState(NAS_DISCOVER_BUBBLES);
			break;
		}
		case APC_POPBUBBLE:
		{		
			m_numPopped = controlMsg.argument;			
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
			PrintDebug(1, "Generating adjacency: %d sequences\n",
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
			continue;
		}
		
		// Sequence is trimmable, create a new branch for it
		BranchGroup newGroup(branchGroupID, dir, 1, *iter);
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

/** Write a packed sequence checkpoint. */
static void writePackedSeqCheckpoint(
		ISequenceCollection* seqCollection)
{
	Timer timer("PackedSequenceCheckpoint");
	ostringstream s;
	s << "trimmed-" << opt::rank << ".psq";
	AssemblyAlgorithms::outputPackedSequences(
			s.str().c_str(), seqCollection);
}

/** Discover bubbles to pop. */
int NetworkSequenceCollection::performNetworkDiscoverBubbles(ISequenceCollection* seqCollection, int kmerSize)
{
	writePackedSeqCheckpoint(seqCollection);

	Timer timer("NetworkDiscoverBubbles");
	
	// The branch ids
	uint64_t branchGroupID = 0;
	m_finishedGroups.clear();
	
	// make sure the branch group structure is initially empty
	assert(m_activeBranchGroups.empty());
	
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
			PrintDebug(1, "Popping bubbles: %d sequences\n", count);
		
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
				BranchGroup branchGroup(branchGroupID, dir, maxNumBranches, *iter);
				
				// insert the new group into the active group map
				BranchGroupMap::iterator groupIter = m_activeBranchGroups.insert(std::pair<uint64_t, BranchGroup>(branchGroupID,branchGroup)).first;
				
				// initiate the new group
				AssemblyAlgorithms::initiateBranchGroup(groupIter->second, *iter, extRec.dir[dir], multiplicity, expectedBubbleSize);
				
				// Disallow any further branching.
				unsigned numInitialBranches
					= groupIter->second.getNumBranches();
				if (numInitialBranches <= maxNumBranches)
					groupIter->second.setMaxNumBranches(
							numInitialBranches);

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
		
		// Primitive load balancing
		if (m_activeBranchGroups.size() > MAX_ACTIVE) {
			while (m_activeBranchGroups.size() > LOW_ACTIVE) {
				seqCollection->pumpNetwork();
				processBranchesDiscoverBubbles();
			}
		}

		// Process groups that may be finished	
		processBranchesDiscoverBubbles();
		seqCollection->pumpNetwork();
	}
	
	// Wait until the groups finish extending.
	while (processBranchesDiscoverBubbles())
		seqCollection->pumpNetwork();
	assert(m_activeBranchGroups.empty());

	unsigned numDiscovered = m_bubbles.size();
	PrintDebug(1, "Discovered %d bubbles\n", numDiscovered);
	return numDiscovered;
}

/** Pop bubbles discovered previously. */
int NetworkSequenceCollection::performNetworkPopBubbles(ISequenceCollection* /*seqCollection*/)
{
	Timer timer("NetworkPopBubbles");

	unsigned numPopped = 0;
	for (BranchGroupMap::iterator iter = m_bubbles.begin();
			iter != m_bubbles.end(); iter++) {
		assert(iter->second.getStatus() == BGS_JOINED);
		// Check whether this bubble has already been popped.
		if (!iter->second.isAmbiguous(this))
			continue;
		numPopped++;
		AssemblyAlgorithms::writeSNP(
				iter->second, m_numPopped + numPopped);
		AssemblyAlgorithms::collapseJoinedBranches(
				this, iter->second);
		assert(m_pComm->empty());
	}
	m_bubbles.clear();

	if (opt::snpFile != NULL)
		fflush(opt::snpFile);

	m_pLog->write(timer.toString().c_str());
	PrintDebug(0, "Removed %d bubbles\n", numPopped);
	return numPopped;
}

//
// Process groups that are finished searching for bubbles
//
bool NetworkSequenceCollection::processBranchesDiscoverBubbles()
{
	bool active = false;
	// Check if any of the current branches have gone inactive
	BranchGroupMap::iterator iter = m_activeBranchGroups.begin();
	while (iter != m_activeBranchGroups.end()) {
		// All branches have been extended one sequence. Check the
		// stop conditions. updateStatus() is called in
		// processSequenceExtensionPop().
		BranchGroupStatus status = iter->second.isNoExt() ? BGS_NOEXT
			: iter->second.getStatus();
		bool finished = false;
		switch (status) {
			case BGS_TOOLONG:
			case BGS_LOOPFOUND:
			case BGS_TOOMANYBRANCHES:
			case BGS_NOEXT:
				finished = true;
				break;
			case BGS_JOINED:
				m_bubbles.insert(*iter);
				finished = true;
				break;
			case BGS_ACTIVE:
				active = true;
				break;
			default:
				assert(false);
		}
		if (finished) {
			m_finishedGroups.insert(iter->first);
			m_activeBranchGroups.erase(iter++);
		} else
			iter++;
	}
	return active;
}

/** Discover bubbles to pop. */
unsigned NetworkSequenceCollection::controlDiscoverBubbles()
{
	SetState(NAS_DISCOVER_BUBBLES);
	m_pComm->SendControlMessage(m_numDataNodes,
			APC_DISCOVER_BUBBLES);

	unsigned numDiscovered = performNetworkDiscoverBubbles(this,
			opt::kmerSize);
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached(m_numDataNodes))
		pumpNetwork();
	numDiscovered += m_checkpointSum;
	SetState(NAS_POPBUBBLE);
	if (numDiscovered > 0 && opt::verbose > 0)
		printf("Discovered %d bubbles\n", numDiscovered);
	return numDiscovered;
}

/** Pop the bubbles discovered previously. */
int NetworkSequenceCollection::controlPopBubbles()
{
	unsigned numDiscovered = controlDiscoverBubbles();
	if (numDiscovered == 0)
		return 0;

	// Perform a round-robin bubble pop to avoid concurrency issues
	unsigned numPopped = performNetworkPopBubbles(this);
	EndState();

	// Now tell all the slave nodes to perform the pop one by one
	for(unsigned i = 1; i < m_numDataNodes; ++i) {
		m_pComm->SendControlMessage(m_numDataNodes, APC_BARRIER);
		m_pComm->barrier();
		SetState(NAS_POPBUBBLE);
		m_pComm->SendControlMessageToNode(i, APC_POPBUBBLE,
				m_numPopped + numPopped);
		while (!checkpointReached(1))
			pumpNetwork();
		numPopped += m_checkpointSum;
		EndState();
	}

	m_numPopped += numPopped;
	if (numPopped > 0)
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
			// Output the singleton contig.
			BranchRecord currBranch(SENSE, -1);
			currBranch.addSequence(*iter, iter->getMultiplicity());
			currBranch.terminate(BS_NOEXT);
			Sequence contig;
			AssemblyAlgorithms::processTerminatedBranchAssemble( 
					seqCollection, currBranch, contig);
			fileWriter->WriteSequence(contig,
					m_numAssembled + numAssembled++,
					currBranch.calculateBranchMultiplicity());
			continue;
		}
		
		// Sequence is trimmable, create a new branch for it
		BranchGroup newGroup(branchGroupID, dir, 1, *iter);
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
			BranchRecord& currBranch = iter->second.getBranch(0);
			if (currBranch.isCanonical()) {
				// Assemble the contig
				Sequence contig;
				AssemblyAlgorithms::processTerminatedBranchAssemble(
						seqCollection, currBranch, contig);
				numAssembled++;

				// Output the contig
				fileWriter->WriteSequence(contig,
						m_numAssembled + currContigID++,
						currBranch.calculateBranchMultiplicity());
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
		case NAS_DISCOVER_BUBBLES:
			return processSequenceExtensionPop(groupID, branchID, seq, extRec, multiplicity);
			break;
		case NAS_WAITING:
			if(m_finishedGroups.find(groupID) == m_finishedGroups.end())
			{
				// The extension message is not in the finished groups list therefore it is unexpected
				// Print some debug info and return
				std::cerr << "Unexpected sequence extension message! gid: " << groupID << " bid: " << branchID << " seq: " << seq.decode() << " Aborting...\n";
				assert(false);
			}
			break;
		default:
			std::cerr << "Unexpected sequence extension message! State: " << m_state << " gid: " << groupID << " bid: " << branchID << " seq: " << seq.decode() << " Aborting...\n";
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
	PrintDebug(0, "Loaded %d sequences\n", m_pLocalSpace->count());	
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
