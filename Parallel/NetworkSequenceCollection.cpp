#include "NetworkSequenceCollection.h"
#include "AssemblyAlgorithms.h"
#include "FastaWriter.h"
#include "Log.h"
#include "Options.h"
#include <cstdlib>
#include <iostream>

using namespace std;

//
//
//
NetworkSequenceCollection::NetworkSequenceCollection(
		int myID, int numDataNodes) :
	m_id(myID),
	m_numDataNodes(numDataNodes),
	m_state(NAS_WAITING),
	m_numBasesAdjSet(0),
	m_trimStep(0),
	m_numPopped(0),
	m_numAssembled(0)
{
	// Load the phase space
	m_pLocalSpace = new SequenceCollectionHash();

	// Create the comm layer
	m_pComm = new CommLayer(myID);

	// Create the message buffer
	m_pMsgBuffer = new MessageBuffer(numDataNodes, m_pComm);
	m_pComm->setMsgBuffer(m_pMsgBuffer);
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
	
	delete m_pMsgBuffer;
	m_pMsgBuffer = 0;
}

void NetworkSequenceCollection::loadSequences()
{
	Timer timer("LoadSequences");
	for (unsigned i = m_id;
			i < opt::inFiles.size();
			i += m_numDataNodes)
		AssemblyAlgorithms::loadSequences(this, opt::inFiles[i]);
}

/** Receive, process, send, and synchronize.
 * @return the number of packets received
 */
unsigned NetworkSequenceCollection::pumpFlushReduce()
{
	m_pMsgBuffer->flush(); // Send.
	m_pComm->barrier(); // Synchronize.
	unsigned count = pumpNetwork(); // Receive and process.
	if (count == 0)
		assert(m_pMsgBuffer->empty());
	return m_pComm->reduce(count); // Reduce.
}

/** Receive packets and process them until no more work exists for any
 * slave processor.
 */
void NetworkSequenceCollection::completeOperation()
{
	Timer timer("completeOperation");

	while (pumpFlushReduce() > 0)
		;

	assert(m_pMsgBuffer->empty()); // Nothing to send.
	m_pComm->barrier(); // Synchronize.
	assert(m_pComm->empty()); // Nothing to receive.
}

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
			case NAS_LOAD_COMPLETE:
				m_pComm->barrier();
				pumpNetwork();
				PrintDebug(0, "Loaded %zu sequences\n",
					m_pLocalSpace->count());
				m_pLocalSpace->printLoad();
				m_pComm->reduce(m_pLocalSpace->count());
				EndState();
				SetState(NAS_WAITING);
				break;
			case NAS_GEN_ADJ:
				AssemblyAlgorithms::generateAdjacency(this);
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage();
				break;
			case NAS_ADJ_COMPLETE:
				m_pComm->barrier();
				pumpNetwork();
				PrintDebug(0, "Generated %u edges\n",
						m_numBasesAdjSet);
				m_pComm->reduce(m_numBasesAdjSet);
				EndState();
				SetState(NAS_WAITING);
				break;
			case NAS_ERODE:
			{
				m_pComm->barrier();
				unsigned numEroded
					= AssemblyAlgorithms::erodeEnds(this);
				EndState();
				SetState(NAS_ERODE_WAITING);
				m_pComm->SendCheckPointMessage(numEroded);
				break;
			}
			case NAS_ERODE_WAITING:
				pumpNetwork();
				break;
			case NAS_ERODE_COMPLETE:
				completeOperation();
				m_pComm->reduce(AssemblyAlgorithms::getNumEroded());

				m_pComm->reduce(m_pLocalSpace->cleanup());
				m_pLocalSpace->printLoad();
				m_pComm->barrier();

				SetState(NAS_WAITING);
				break;
			case NAS_TRIM:
			{
				assert(m_trimStep != 0);
				m_pComm->barrier();
				int numRemoved = performNetworkTrim(this, m_trimStep);
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage(numRemoved);
				break;
			}
			case NAS_REMOVE_MARKED: {
				m_pComm->barrier();
				unsigned count
					= AssemblyAlgorithms::removeMarked(this);
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage(count);
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
				m_pComm->barrier();
				assert(m_pComm->empty());
				m_pComm->reduce(
						AssemblyAlgorithms::markAmbiguous(this));
				assert(m_pMsgBuffer->empty());
				assert(m_pComm->empty());
				m_pComm->barrier();
				unsigned count
					= AssemblyAlgorithms::splitAmbiguous(this);
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage(count);
				break;
			}
			case NAS_ASSEMBLE:
			{
				FastaWriter* writer = new FastaWriter(
						opt::contigsTempPath.c_str());
				unsigned numAssembled = performNetworkAssembly(this, writer);
				delete writer;
				EndState();
				SetState(NAS_WAITING);
				m_pComm->SendCheckPointMessage(numAssembled);
				break;
			}
			case NAS_WAITING:
				pumpNetwork();
				break;
			case NAS_DONE:
				stop = true;
				break;
			assert(false);
		}
	}
}

unsigned NetworkSequenceCollection::controlErode()
{
	SetState(NAS_ERODE);
	m_pComm->SendControlMessage(m_numDataNodes, APC_ERODE);
	m_pComm->barrier();
	unsigned numEroded = AssemblyAlgorithms::erodeEnds(this);
	EndState();

	// Do not call SetState, because it clears the
	// checkpoint information.
	//SetState(NAS_ERODE_WAITING);
	m_state = NAS_ERODE_WAITING;

	m_numReachedCheckpoint++;
	while (!checkpointReached(m_numDataNodes))
		pumpNetwork();
	numEroded += m_checkpointSum;
	EndState();

	if (numEroded == 0) {
		SetState(NAS_WAITING);
		m_pComm->SendControlMessage(m_numDataNodes, APC_WAIT);
		m_pComm->barrier();
		return 0;
	}

	SetState(NAS_ERODE_COMPLETE);
	m_pComm->SendControlMessage(m_numDataNodes,
			APC_ERODE_COMPLETE);
	completeOperation();
	numEroded += m_pComm->reduce(
			AssemblyAlgorithms::getNumEroded());
	printf("Eroded %u tips\n", numEroded);

	unsigned removed = m_pComm->reduce(m_pLocalSpace->cleanup());
	m_pLocalSpace->printLoad();
	m_pComm->barrier();
	assert(removed == numEroded);

	SetState(NAS_WAITING);
	return numEroded;
}

/** Remove marked sequences.
 * @return the number of sequences removed
 */
unsigned NetworkSequenceCollection::controlRemoveMarked()
{
	if (opt::verbose > 0)
		puts("Sweeping");
	SetState(NAS_REMOVE_MARKED);
	m_pComm->SendControlMessage(m_numDataNodes, APC_REMOVE_MARKED);
	m_pComm->barrier();
	unsigned count = AssemblyAlgorithms::removeMarked(this);
	m_checkpointSum += count;
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached(m_numDataNodes))
		pumpNetwork();
	return m_checkpointSum;
}

/** Perform a single round of trimming at the specified length. */
unsigned NetworkSequenceCollection::controlTrimRound(unsigned trimLen)
{
	printf("Trimming short branches: %d\n", trimLen);
	SetState(NAS_TRIM);
	m_pComm->SendControlMessage(m_numDataNodes, APC_TRIM, trimLen);
	m_pComm->barrier();
	unsigned numRemoved = performNetworkTrim(this, trimLen);
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached(m_numDataNodes))
		pumpNetwork();
	numRemoved += m_checkpointSum;

	unsigned numSweeped = controlRemoveMarked();

	if (numRemoved > 0)
		printf("Trimmed %u sequences in %u branches\n",
				numSweeped, numRemoved);
	return numRemoved;
}

/** Perform multiple rounds of trimming until complete. */
void NetworkSequenceCollection::controlTrim(unsigned start)
{
	unsigned rounds = 0, total = 0;
	for (int trim = start; trim < opt::trimLen; trim *= 2) {
		rounds++;
		total += controlTrimRound(trim);
	}
	unsigned count;
	while ((count = controlTrimRound(opt::trimLen)) > 0) {
		rounds++;
		total += count;
	}
	printf("Trimmed %u branches in %u rounds\n", total, rounds);
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

				SetState(NAS_LOAD_COMPLETE);
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_LOAD_COMPLETE);
				m_pComm->barrier();
				pumpNetwork();
				PrintDebug(0, "Loaded %zu sequences\n",
					m_pLocalSpace->count());
				m_pLocalSpace->printLoad();
				printf("Loaded %u sequences\n",
						m_pComm->reduce(m_pLocalSpace->count()));
				EndState();

				SetState(m_pLocalSpace->isAdjacencyLoaded()
						? NAS_ERODE : NAS_GEN_ADJ);
				break;
			case NAS_GEN_ADJ:
				puts("Generating adjacency");
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_GEN_ADJ);
				AssemblyAlgorithms::generateAdjacency(this);
				
				// Cleanup any messages that are pending
				EndState();
								
				m_numReachedCheckpoint++;
				while(!checkpointReached(m_numDataNodes))
				{
					pumpNetwork();
				}

				SetState(NAS_ADJ_COMPLETE);
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_ADJ_COMPLETE);
				m_pComm->barrier();
				pumpNetwork();
				PrintDebug(0, "Generated %u edges\n",
						m_numBasesAdjSet);
				printf("Generated %u edges\n",
						m_pComm->reduce(m_numBasesAdjSet));
				EndState();

				SetState(opt::erode > 0 ? NAS_ERODE : NAS_TRIM);
				break;
			case NAS_ERODE:
				assert(opt::erode > 0);
				puts("Eroding tips");
				controlErode();
				assert(controlErode() == 0);
				SetState(NAS_TRIM);
				break;
			case NAS_LOAD_COMPLETE:
			case NAS_ADJ_COMPLETE:
			case NAS_REMOVE_MARKED:
			case NAS_ERODE_WAITING:
			case NAS_ERODE_COMPLETE:
			case NAS_DISCOVER_BUBBLES:
				// These states are used only by the slaves.
				assert(false);
				exit(EXIT_FAILURE);

			case NAS_TRIM:
				controlTrim();
				SetState(NAS_POPBUBBLE);
				break;

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

				// Another round of trimming to remove tips created by
				// popping complex bubbles.
				if (totalPopped > 0)
					controlTrim(opt::trimLen);

				SetState(NAS_SPLIT);
				break;
			}
			case NAS_SPLIT:
			{
				m_pComm->SendControlMessage(m_numDataNodes,
						APC_SPLIT);
				m_pComm->barrier();
				unsigned marked = controlMarkAmbiguous();
				unsigned split = controlSplitAmbiguous();
				assert(marked == split);
				SetState(NAS_ASSEMBLE);
				break;
			}
			case NAS_ASSEMBLE:
			{
				puts("Assembling");
				FastaWriter* writer = new FastaWriter(
						opt::contigsTempPath.c_str());
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
				stop = true;
				break;
			case NAS_WAITING:
				break;
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
	PrintDebug(2, "SetState %u (was %u)\n", newState, m_state);

	// Ensure there are no pending messages
	assert(m_pMsgBuffer->empty());
	
	m_state = newState;
	
	// Reset the checkpoint counter
	m_numReachedCheckpoint = 0;
	m_checkpointSum = 0;
}

/** Receive and dispatch packets.
 * @return the number of packets received
 */
unsigned NetworkSequenceCollection::pumpNetwork()
{
	for (unsigned count = 0; ; count++) {
		int senderID;
		APMessage msg = m_pComm->CheckMessage(senderID);
		switch(msg)
		{
			case APM_CONTROL:
				parseControlMessage();
				// Deal with the control packet before we continue
				// processing further packets.
				return ++count;
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
				return count;
		}
	}
}

/** Call the observers of the specified sequence. */
void NetworkSequenceCollection::notify(
		const PackedSeq& key)
{
	const PackedSeq& seq = m_pLocalSpace->getSeqAndData(key);
	switch (m_state) {
		case NAS_ERODE:
		case NAS_ERODE_WAITING:
		case NAS_ERODE_COMPLETE:
			AssemblyAlgorithms::erode(this, seq);
			break;
		default:
			// Nothing to do.
			break;
	}
}

//
//
//
void NetworkSequenceCollection::handleSeqOpMessage(int /*senderID*/, const SeqOpMessage& seqMsg)
{
	switch(seqMsg.m_operation)
	{
		case MO_ADD:
			add(seqMsg.m_seq);
			break;
		case MO_REMOVE:
			remove(seqMsg.m_seq);
			break;
		default:
			assert(false);
			break;
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
	assert(isLocal(message.m_seq));
	setBaseExtension(message.m_seq, message.m_dir, message.m_base);
}

//
// Parse a sequence extension message
//
void NetworkSequenceCollection::handleRemoveExtensionMessage(int /*senderID*/, const RemoveExtensionMessage& message)
{
	assert(isLocal(message.m_seq));	
	bool found = m_pLocalSpace->removeExtension(
			message.m_seq, message.m_dir, message.m_base);
	if (found)
		notify(message.m_seq);
}

//
//
//
void NetworkSequenceCollection::parseControlMessage()
{
	ControlMessage controlMsg = m_pComm->ReceiveControlMessage();
	switch(controlMsg.msgType)
	{
		case APC_LOAD_COMPLETE:
			SetState(NAS_LOAD_COMPLETE);
			break;
		case APC_CHECKPOINT:
			m_numReachedCheckpoint++;
			m_checkpointSum += controlMsg.argument;
			break;	
		case APC_WAIT:
			SetState(NAS_WAITING);
			m_pComm->barrier();
			break;
		case APC_BARRIER:
			assert(m_state == NAS_WAITING);
			m_pComm->barrier();
			break;
		case APC_FINISHED:
			SetState(NAS_DONE);
			break;	
		case APC_TRIM:
			m_trimStep = controlMsg.argument;
			SetState(NAS_TRIM);
			break;
		case APC_REMOVE_MARKED:
			SetState(NAS_REMOVE_MARKED);
			break;
		case APC_ERODE:
			SetState(NAS_ERODE);
			break;
		case APC_ERODE_COMPLETE:
			assert(m_state == NAS_ERODE_WAITING);
			EndState();
			SetState(NAS_ERODE_COMPLETE);
			break;
		case APC_DISCOVER_BUBBLES:
			SetState(NAS_DISCOVER_BUBBLES);
			break;
		case APC_POPBUBBLE:
			m_numPopped = controlMsg.argument;			
			SetState(NAS_POPBUBBLE);
			break;	
		case APC_SPLIT:
			SetState(NAS_SPLIT);
			break;	
		case APC_ASSEMBLE:
			// This message came from the control node along with an argument indicating the number of sequences assembled so far
			m_numAssembled = controlMsg.argument;			
			SetState(NAS_ASSEMBLE);
			break;	
		case APC_GEN_ADJ:
			SetState(NAS_GEN_ADJ);
			break;	
		case APC_ADJ_COMPLETE:
			SetState(NAS_ADJ_COMPLETE);
			break;
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

/** Discover bubbles to pop. */
int NetworkSequenceCollection::performNetworkDiscoverBubbles(ISequenceCollection* seqCollection, int kmerSize)
{
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
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			if (extRec.dir[dir].isAmbiguous()) {
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

	// Deal with any packets still in the queue. The barrier
	// synchronization guarantees that the packets have been
	// delivered, but we may not have dealt with them yet.
	pumpNetwork();
	assert(m_pComm->empty());

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

		if(!m_pComm->empty())
		{
			int sendID;
			std::cerr << " COMM NOT EMPTY MESSAGE WAITING ID: " << (int)m_pComm->CheckMessage(sendID) << " sender " << sendID << std::endl;
			std::cerr << " Attempting pump " << std::endl;
			pumpNetwork();
			std::cerr << " Pump returned ok " << std::endl;
			assert(false);
		}
		assert(m_pComm->empty());
	}
	m_bubbles.clear();

	if (opt::snpFile != NULL)
		fflush(opt::snpFile);

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
	m_checkpointSum = performNetworkPopBubbles(this);
	EndState();

	// Now tell all the slave nodes to perform the pop one by one
	for(unsigned i = 1; i < m_numDataNodes; ++i) {
		m_pComm->SendControlMessage(m_numDataNodes, APC_BARRIER);
		m_pComm->barrier();
		m_numReachedCheckpoint = 0;
		m_pComm->SendControlMessageToNode(i, APC_POPBUBBLE,
				m_numPopped + m_checkpointSum);
		while (!checkpointReached(1))
			pumpNetwork();
	}

	unsigned numPopped = m_checkpointSum;
	m_numPopped += numPopped;
	if (numPopped > 0)
		printf("Removed %d bubbles\n", numPopped);
	return numPopped;
}

/** Mark ambiguous branches. */
unsigned NetworkSequenceCollection::controlMarkAmbiguous()
{
	puts("Marking ambiguous branches");
	assert(m_pComm->empty());
	unsigned count = m_pComm->reduce(
			AssemblyAlgorithms::markAmbiguous(this));
	assert(m_pMsgBuffer->empty());
	assert(m_pComm->empty());
	m_pComm->barrier();
	printf("Marked %u ambiguous branches\n", count);
	return count;
}

/** Remove ambiguous branches. */
unsigned NetworkSequenceCollection::controlSplitAmbiguous()
{
	puts("Splitting ambiguous branches");
	unsigned count = AssemblyAlgorithms::splitAmbiguous(this);
	m_checkpointSum += count;
	EndState();
	m_numReachedCheckpoint++;
	while (!checkpointReached(m_numDataNodes))
		pumpNetwork();
	printf("Split %u ambiguous branches\n",
			m_checkpointSum);
	return m_checkpointSum;
}

/** Assemble contigs. */
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
		case NAS_ASSEMBLE:
			return processLinearSequenceExtension(groupID, branchID, seq, extRec, multiplicity);
		case NAS_DISCOVER_BUBBLES:
			return processSequenceExtensionPop(groupID, branchID, seq, extRec, multiplicity);
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
bool NetworkSequenceCollection::exists(const PackedSeq& seq)
{
	assert(isLocal(seq));
	return m_pLocalSpace->exists(seq);
}

//
//
//
bool NetworkSequenceCollection::checkpointReached(int numRequired) const
{
	return m_numReachedCheckpoint == numRequired;
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

bool NetworkSequenceCollection::checkFlag(const PackedSeq& seq,
		SeqFlag flag) const
{
	assert(isLocal(seq));
	return m_pLocalSpace->checkFlag(seq, flag);
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
	assert(isLocal(seq));
	return m_pLocalSpace->hasParent(seq);
}

//
//
//
bool NetworkSequenceCollection::hasChild(const PackedSeq& seq)
{
	assert(isLocal(seq));
	return m_pLocalSpace->hasChild(seq);
}

bool NetworkSequenceCollection::setBaseExtension(
		const PackedSeq& seq, extDirection dir, uint8_t base)
{
	if (isLocal(seq)) {
		if (m_pLocalSpace->setBaseExtension(seq, dir, base))
			m_numBasesAdjSet++;
	} else {
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
	assert(isLocal(seq));
	m_pLocalSpace->clearExtensions(seq, dir);
}

//
//
//
size_t NetworkSequenceCollection::count() const
{
	return m_pLocalSpace->count();
}

/** Remove the specified extension from the specified sequence if it
 * exists in this collection.
 * @return true if the specified sequence is local and exists and
 * false otherwise
 */
bool NetworkSequenceCollection::removeExtension(
		const PackedSeq& seq, extDirection dir, uint8_t base)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		bool found = m_pLocalSpace->removeExtension(seq, dir, base);
		if (found)
			notify(seq);
		return found;
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_pMsgBuffer->sendRemoveExtension(nodeID, seq, dir, base);
		return false;
	}
}

// get the extensions of the sequence
bool NetworkSequenceCollection::getSeqData(const PackedSeq& seq,
		ExtensionRecord& extRecord, int& multiplicity) const
{
	// This function can only be called locally, the distributed
	// version is through generateSequenceExtensionMessage.
	assert(isLocal(seq));
	m_pLocalSpace->getSeqData(seq, extRecord, multiplicity);
	return true;
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
