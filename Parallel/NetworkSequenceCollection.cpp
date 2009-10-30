#include "NetworkSequenceCollection.h"
#include "Assembly/Options.h"
#include "AssemblyAlgorithms.h"
#include "Common/Options.h"
#include "FastaWriter.h"
#include "Histogram.h"
#include "Log.h"
#include <cmath> // for roundf
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <utility>

using namespace std;

NetworkSequenceCollection::NetworkSequenceCollection()
	: m_pLocalSpace(new SequenceCollectionHash()),
	m_state(NAS_WAITING),
	m_trimStep(0),
	m_numPopped(0),
	m_numAssembled(0)
{
}

NetworkSequenceCollection::~NetworkSequenceCollection()
{
	// Delete the objects created in the constructor
	delete m_pLocalSpace;
	m_pLocalSpace = 0;
}

void NetworkSequenceCollection::loadSequences()
{
	Timer timer("LoadSequences");
	for (unsigned i = opt::rank; i < opt::inFiles.size();
			i += opt::numProc)
		AssemblyAlgorithms::loadSequences(this, opt::inFiles[i]);
}

/** Receive, process, send, and synchronize.
 * @return the number of inflight messages
 */
unsigned NetworkSequenceCollection::pumpFlushReduce()
{
	pumpNetwork();
	m_comm.flush();
	return m_comm.reduceInflight();
}

/** Receive packets and process them until no more work exists for any
 * slave processor.
 */
void NetworkSequenceCollection::completeOperation()
{
	Timer timer("completeOperation");

	while (pumpFlushReduce() > 0)
		;

	assert(m_comm.transmitBufferEmpty()); // Nothing to send.
	m_comm.barrier(); // Synchronize.
	assert(m_comm.receiveEmpty()); // Nothing to receive.
	assert(m_comm.reduceInflight() == 0);
}

/** Run the assembly state machine. */
void NetworkSequenceCollection::run()
{
	/** The number of contigs and k-mer assembled. */
	pair<unsigned, unsigned> numAssembled;

	ofstream bubbleFile;

	SetState(NAS_LOADING);
	while (m_state != NAS_DONE) {
		switch (m_state) {
			case NAS_LOADING:
				m_pLocalSpace->setColourSpace(
						m_comm.receiveBroadcast());
				loadSequences();
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage();
				break;
			case NAS_LOAD_COMPLETE:
			{
				m_comm.barrier();
				pumpNetwork();
				PrintDebug(0, "Loaded %zu k-mer\n",
					m_pLocalSpace->count());
				m_pLocalSpace->printLoad();
				m_comm.reduce(m_pLocalSpace->count());

				Histogram h = m_comm.reduce(
						AssemblyAlgorithms::coverageHistogram(
							*m_pLocalSpace));
				AssemblyAlgorithms::determineMinimumCoverage(h);
				EndState();
				SetState(NAS_WAITING);
				break;
			}
			case NAS_GEN_ADJ:
				m_comm.barrier();
				m_numBasesAdjSet = 0;
				AssemblyAlgorithms::generateAdjacency(this);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage();
				break;
			case NAS_ADJ_COMPLETE:
				m_comm.barrier();
				pumpNetwork();
				PrintDebug(0, "Generated %u edges\n",
						m_numBasesAdjSet);
				m_comm.reduce(m_numBasesAdjSet);
				EndState();
				SetState(NAS_WAITING);
				break;
			case NAS_ERODE:
			{
				m_comm.barrier();
				unsigned numEroded
					= AssemblyAlgorithms::erodeEnds(this);
				EndState();
				SetState(NAS_ERODE_WAITING);
				m_comm.sendCheckPointMessage(numEroded);
				break;
			}
			case NAS_ERODE_WAITING:
				pumpNetwork();
				break;
			case NAS_ERODE_COMPLETE:
				completeOperation();
				m_comm.reduce(AssemblyAlgorithms::getNumEroded());

				m_comm.reduce(m_pLocalSpace->cleanup());
				m_pLocalSpace->printLoad();
				m_comm.barrier();

				SetState(NAS_WAITING);
				break;
			case NAS_TRIM:
			{
				assert(m_trimStep != 0);
				m_comm.barrier();
				int numRemoved = performNetworkTrim(this, m_trimStep);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage(numRemoved);
				break;
			}
			case NAS_REMOVE_MARKED: {
				m_comm.barrier();
				unsigned count
					= AssemblyAlgorithms::removeMarked(this);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage(count);
				break;
			}

			case NAS_COVERAGE:
			{
				m_comm.barrier();
				m_lowCoverageContigs = 0;
				m_lowCoverageKmer = 0;
				numAssembled = performNetworkAssembly(this);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage();
				break;
			}
			case NAS_COVERAGE_COMPLETE:
				m_comm.barrier();
				pumpNetwork();
				m_pLocalSpace->wipeFlag(
						SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));
				m_comm.reduce(numAssembled.first);
				m_comm.reduce(numAssembled.second);
				m_comm.reduce(m_lowCoverageContigs);
				m_comm.reduce(m_lowCoverageKmer);
				m_comm.reduce(m_pLocalSpace->cleanup());
				opt::coverage = 0;
				EndState();
				SetState(NAS_WAITING);
				break;

			case NAS_DISCOVER_BUBBLES:
			{
				unsigned numDiscovered
					= performNetworkDiscoverBubbles(this,
							opt::kmerSize);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage(numDiscovered);
				break;
			}
			case NAS_POPBUBBLE:
			{
				if (!bubbleFile.is_open())
					AssemblyAlgorithms::openBubbleFile(bubbleFile);
				unsigned numPopped
					= performNetworkPopBubbles(bubbleFile);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage(numPopped);
				break;
			}
			case NAS_SPLIT:
			{
				m_comm.barrier();
				assert(m_comm.receiveEmpty());
				m_comm.reduce(
						AssemblyAlgorithms::markAmbiguous(this));
				assert(m_comm.transmitBufferEmpty());
				assert(m_comm.receiveEmpty());
				m_comm.barrier();
				unsigned count
					= AssemblyAlgorithms::splitAmbiguous(this);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage(count);
				break;
			}
			case NAS_ASSEMBLE:
			{
				m_comm.barrier();
				FastaWriter writer(opt::contigsTempPath.c_str());
				numAssembled = performNetworkAssembly(this, &writer);
				EndState();
				SetState(NAS_WAITING);
				m_comm.sendCheckPointMessage();
				break;
			}
			case NAS_ASSEMBLE_COMPLETE:
				m_comm.reduce(numAssembled.first);
				m_comm.reduce(numAssembled.second);
				EndState();
				SetState(NAS_DONE);
				break;
			case NAS_WAITING:
				pumpNetwork();
				break;
			case NAS_DONE:
				break;
		}
	}
}

unsigned NetworkSequenceCollection::controlErode()
{
	SetState(NAS_ERODE);
	m_comm.sendControlMessage(APC_ERODE);
	m_comm.barrier();
	unsigned numEroded = AssemblyAlgorithms::erodeEnds(this);
	EndState();

	// Do not call SetState, because it clears the
	// checkpoint information.
	//SetState(NAS_ERODE_WAITING);
	m_state = NAS_ERODE_WAITING;

	m_numReachedCheckpoint++;
	while (!checkpointReached())
		pumpNetwork();
	numEroded += m_checkpointSum;
	EndState();

	if (numEroded == 0) {
		SetState(NAS_WAITING);
		m_comm.sendControlMessage(APC_WAIT);
		m_comm.barrier();
		return 0;
	}

	SetState(NAS_ERODE_COMPLETE);
	m_comm.sendControlMessage(APC_ERODE_COMPLETE);
	completeOperation();
	numEroded += m_comm.reduce(
			AssemblyAlgorithms::getNumEroded());
	printf("Eroded %u tips\n", numEroded);

	unsigned removed = m_comm.reduce(m_pLocalSpace->cleanup());
	m_pLocalSpace->printLoad();
	m_comm.barrier();
	assert(removed == numEroded);

	SetState(NAS_WAITING);
	return numEroded;
}

/** Remove marked k-mer.
 * @return the number of k-mer removed
 */
unsigned NetworkSequenceCollection::controlRemoveMarked()
{
	if (opt::verbose > 0)
		puts("Sweeping");
	SetState(NAS_REMOVE_MARKED);
	m_comm.sendControlMessage(APC_REMOVE_MARKED);
	m_comm.barrier();
	unsigned count = AssemblyAlgorithms::removeMarked(this);
	m_checkpointSum += count;
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached())
		pumpNetwork();
	return m_checkpointSum;
}

/** Perform a single round of trimming at the specified length. */
unsigned NetworkSequenceCollection::controlTrimRound(unsigned trimLen)
{
	assert(trimLen > 0);
	printf("Trimming short branches: %u\n", trimLen);
	SetState(NAS_TRIM);
	m_comm.sendControlMessage(APC_TRIM, trimLen);
	m_comm.barrier();
	unsigned numRemoved = performNetworkTrim(this, trimLen);
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached())
		pumpNetwork();
	numRemoved += m_checkpointSum;

	unsigned numSweeped = controlRemoveMarked();

	if (numRemoved > 0)
		printf("Trimmed %u k-mer in %u branches\n",
				numSweeped, numRemoved);
	return numRemoved;
}

/** Perform multiple rounds of trimming until complete. */
void NetworkSequenceCollection::controlTrim(unsigned start)
{
	if (opt::trimLen == 0)
		return;
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

/** Remove low-coverage contigs. */
void NetworkSequenceCollection::controlCoverage()
{
	assert(opt::coverage > 0);

	// Split ambiguous branches.
	SetState(NAS_SPLIT);
	controlSplit();

	// Remove low-coverage contigs.
	printf("Removing low-coverage contigs "
			"(mean k-mer coverage < %f)\n", opt::coverage);
	SetState(NAS_COVERAGE);
	m_comm.sendControlMessage(APC_COVERAGE);
	m_comm.barrier();
	m_lowCoverageContigs = 0;
	m_lowCoverageKmer = 0;
	pair<unsigned, unsigned> numAssembled
		= performNetworkAssembly(this);
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached())
		pumpNetwork();

	// Count the number of low-coverage contigs.
	SetState(NAS_COVERAGE_COMPLETE);
	m_comm.sendControlMessage(APC_COVERAGE_COMPLETE);
	m_comm.barrier();
	pumpNetwork();
	m_pLocalSpace->wipeFlag(
			SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));

	numAssembled.first = m_comm.reduce(numAssembled.first);
	numAssembled.second = m_comm.reduce(numAssembled.second);
	printf("Found %u k-mer in %u contigs "
			"before removing low-coverage contigs\n",
			numAssembled.second, numAssembled.first);

	unsigned lowCoverageContigs = m_comm.reduce(m_lowCoverageContigs);
	unsigned lowCoverageKmer = m_comm.reduce(m_lowCoverageKmer);
	unsigned removed = m_comm.reduce(m_pLocalSpace->cleanup());
	printf("Removed %u k-mer in %u low-coverage contigs\n",
			lowCoverageKmer, lowCoverageContigs);
	printf("Removed %u marked k-mer\n", removed);
	opt::coverage = 0;
	EndState();
}

/** Run the assembly state machine for the controller (rank = 0). */
void NetworkSequenceCollection::runControl()
{
	SetState(NAS_LOADING);
	while (m_state != NAS_DONE) {
		switch (m_state) {
			case NAS_LOADING:
			{
				loadSequences();
				assert(m_pLocalSpace->count() > 0);
				EndState();

				m_numReachedCheckpoint++;
				while (!checkpointReached())
					pumpNetwork();

				SetState(NAS_LOAD_COMPLETE);
				m_comm.sendControlMessage(APC_LOAD_COMPLETE);
				m_comm.barrier();
				pumpNetwork();
				PrintDebug(0, "Loaded %zu k-mer\n",
					m_pLocalSpace->count());
				m_pLocalSpace->printLoad();
				printf("Loaded %lu k-mer\n",
						m_comm.reduce(m_pLocalSpace->count()));

				Histogram h = m_comm.reduce(
						AssemblyAlgorithms::coverageHistogram(
							*m_pLocalSpace));
				AssemblyAlgorithms::determineMinimumCoverage(h);
				EndState();

				SetState(m_pLocalSpace->isAdjacencyLoaded()
						? NAS_ERODE : NAS_GEN_ADJ);
				break;
			}
			case NAS_GEN_ADJ:
				puts("Generating adjacency");
				m_comm.sendControlMessage(APC_GEN_ADJ);
				m_comm.barrier();
				m_numBasesAdjSet = 0;
				AssemblyAlgorithms::generateAdjacency(this);
				EndState();

				m_numReachedCheckpoint++;
				while (!checkpointReached())
					pumpNetwork();

				SetState(NAS_ADJ_COMPLETE);
				m_comm.sendControlMessage(APC_ADJ_COMPLETE);
				m_comm.barrier();
				pumpNetwork();
				PrintDebug(0, "Generated %u edges\n",
						m_numBasesAdjSet);
				printf("Generated %lu edges\n",
						m_comm.reduce(m_numBasesAdjSet));
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
			case NAS_COVERAGE_COMPLETE:
			case NAS_DISCOVER_BUBBLES:
			case NAS_ASSEMBLE_COMPLETE:
			case NAS_WAITING:
				// These states are used only by the slaves.
				assert(false);
				exit(EXIT_FAILURE);

			case NAS_TRIM:
				controlTrim();
				SetState(opt::coverage > 0 ? NAS_COVERAGE
						: opt::bubbles > 0 ? NAS_POPBUBBLE
						: NAS_SPLIT);
				break;

			case NAS_COVERAGE:
				controlCoverage();
				SetState(NAS_GEN_ADJ);
				break;

			case NAS_POPBUBBLE:
			{
				ofstream out;
				AssemblyAlgorithms::openBubbleFile(out);

				puts("Popping bubbles");
				unsigned totalPopped = 0;
				int i;
				for (i = 0; i < opt::bubbles; i++) {
					unsigned numPopped = controlPopBubbles(out);
					if (numPopped == 0)
						break;
					totalPopped += numPopped;
				}
				assert(totalPopped == m_numPopped);
				assert(out.good());
				out.close();
				printf("Removed %u bubbles in %u rounds\n",
						totalPopped, i);

				// Another round of trimming to remove tips created by
				// popping complex bubbles.
				if (totalPopped > 0)
					controlTrim(opt::trimLen);

				SetState(NAS_SPLIT);
				break;
			}
			case NAS_SPLIT:
				controlSplit();
				SetState(NAS_ASSEMBLE);
				break;
			case NAS_ASSEMBLE:
			{
				puts("Assembling");
				m_comm.sendControlMessage(APC_ASSEMBLE);
				m_comm.barrier();
				FastaWriter* writer = new FastaWriter(
						opt::contigsTempPath.c_str());
				pair<unsigned, unsigned> numAssembled
					= performNetworkAssembly(this, writer);
				delete writer;
				EndState();

				m_numReachedCheckpoint++;
				while (!checkpointReached())
					pumpNetwork();

				SetState(NAS_ASSEMBLE_COMPLETE);
				m_comm.sendControlMessage(APC_ASSEMBLE_COMPLETE);

				numAssembled.first = m_comm.reduce(
						numAssembled.first);
				numAssembled.second = m_comm.reduce(
						numAssembled.second);
				printf("Assembled %u k-mer in %u contigs\n",
						numAssembled.second, numAssembled.first);

				SetState(NAS_DONE);
				break;
			}
			case NAS_DONE:
				break;
		}
	}
}

void NetworkSequenceCollection::EndState()
{
	// Flush the message buffer
	m_comm.flush();
}

//
// Set the state
//
void NetworkSequenceCollection::SetState(NetworkAssemblyState newState)
{
	PrintDebug(2, "SetState %u (was %u)\n", newState, m_state);

	// Ensure there are no pending messages
	assert(m_comm.transmitBufferEmpty());

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
		APMessage msg = m_comm.checkMessage(senderID);
		switch(msg)
		{
			case APM_CONTROL:
				parseControlMessage(senderID);
				// Deal with the control packet before we continue
				// processing further packets.
				return ++count;
			case APM_BUFFERED:
				{
					MessagePtrVector msgs;
					m_comm.receiveBufferedMessage(msgs);
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
void NetworkSequenceCollection::notify(const Kmer& key)
{
	switch (m_state) {
		case NAS_ERODE:
		case NAS_ERODE_WAITING:
		case NAS_ERODE_COMPLETE:
			AssemblyAlgorithms::erode(this,
					m_pLocalSpace->getSeqAndData(key));
			break;
		default:
			// Nothing to do.
			break;
	}
}

void NetworkSequenceCollection::handleSeqOpMessage(int /*senderID*/, const SeqOpMessage& seqMsg)
{
	assert(isLocal(seqMsg.m_seq));
	switch(seqMsg.m_operation)
	{
		case MO_ADD:
			m_pLocalSpace->add(seqMsg.m_seq);
			break;
		case MO_REMOVE:
			m_pLocalSpace->remove(seqMsg.m_seq);
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
	m_pLocalSpace->removeExtension(
			message.m_seq, message.m_dir, message.m_ext);
	notify(message.m_seq);
}

//
//
//
void NetworkSequenceCollection::parseControlMessage(int source)
{
	ControlMessage controlMsg = m_comm.receiveControlMessage();
	switch(controlMsg.msgType)
	{
		case APC_LOAD_COMPLETE:
			SetState(NAS_LOAD_COMPLETE);
			break;
		case APC_CHECKPOINT:
			PrintDebug(4, "checkpoint from %u: %u\n",
					source, controlMsg.argument);
			m_numReachedCheckpoint++;
			m_checkpointSum += controlMsg.argument;
			break;	
		case APC_WAIT:
			SetState(NAS_WAITING);
			m_comm.barrier();
			break;
		case APC_BARRIER:
			assert(m_state == NAS_WAITING);
			m_comm.barrier();
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
		case APC_COVERAGE:
			SetState(NAS_COVERAGE);
			break;
		case APC_COVERAGE_COMPLETE:
			SetState(NAS_COVERAGE_COMPLETE);
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
			m_numAssembled = controlMsg.argument;			
			SetState(NAS_ASSEMBLE);
			break;
		case APC_ASSEMBLE_COMPLETE:
			SetState(NAS_ASSEMBLE_COMPLETE);
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
	if (!found)
		cout << "error: from " << senderID << ' ' << message << endl;
	assert(found);
	(void)found;

	// Return the extension to the sender
	m_comm.sendSeqDataResponse(senderID, message.m_group, message.m_id, message.m_seq, extRec, multiplicity);
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

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(
				*iter, dir);
		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// remove this sequence, it has no extensions
			AssemblyAlgorithms::removeSequenceAndExtensions(seqCollection, *iter);
			numBranchesRemoved++;
			continue;
		}

		// Sequence is trimmable, create a new branch for it
		BranchGroup newGroup(branchGroupID, dir, 1, iter->first);
		BranchRecord newBranch(dir, maxBranchCull);
		newGroup.addBranch(0, newBranch);
		m_activeBranchGroups[branchGroupID] = newGroup;

		// Generate the first extension request
		generateExtensionRequest(branchGroupID, 0, iter->first);
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

	PrintDebug(0, "Trimmed %u branches\n", numBranchesRemoved);
	return numBranchesRemoved;	
}

//
// Process current branches, removing those that are finished
// returns true if the branch list has branches remaining
//
int NetworkSequenceCollection::processBranchesTrim()
{
	int numBranchesRemoved = 0;
	vector<BranchGroupMap::iterator> removeBranches;
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
	for (vector<BranchGroupMap::iterator>::iterator rmIter
				= removeBranches.begin();
			rmIter != removeBranches.end(); rmIter++)
		m_activeBranchGroups.erase(*rmIter);	
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

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		if (++count % 100000 == 0)
			PrintDebug(1, "Popping bubbles: %u k-mer\n", count);

		ExtensionRecord extRec = iter->second.extension();
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			if (extRec.dir[dir].isAmbiguous()) {
				BranchGroup branchGroup(branchGroupID, dir,
						maxNumBranches, iter->first);

				// insert the new group into the active group map
				BranchGroupMap::iterator groupIter = m_activeBranchGroups.insert(pair<uint64_t, BranchGroup>(branchGroupID,branchGroup)).first;

				AssemblyAlgorithms::initiateBranchGroup(
						groupIter->second, iter->first,
						extRec.dir[dir], expectedBubbleSize);

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
	PrintDebug(1, "Discovered %u bubbles\n", numDiscovered);
	return numDiscovered;
}

/** Pop bubbles discovered previously. */
int NetworkSequenceCollection::performNetworkPopBubbles(ostream& out)
{
	Timer timer("NetworkPopBubbles");

	// Deal with any packets still in the queue. The barrier
	// synchronization guarantees that the packets have been
	// delivered, but we may not have dealt with them yet.
	pumpNetwork();
	assert(m_comm.receiveEmpty());

	unsigned numPopped = 0;
	for (BranchGroupMap::iterator iter = m_bubbles.begin();
			iter != m_bubbles.end(); iter++) {
		assert(iter->second.getStatus() == BGS_JOINED);
		// Check whether this bubble has already been popped.
		if (!iter->second.isAmbiguous(this))
			continue;
		numPopped++;
		AssemblyAlgorithms::writeBubble(out,
				iter->second, m_numPopped + numPopped);
		AssemblyAlgorithms::collapseJoinedBranches(
				this, iter->second);

		if (!m_comm.receiveEmpty()) {
			int sendID;
			cerr << " COMM NOT EMPTY MESSAGE WAITING ID: "
				<< (int)m_comm.checkMessage(sendID)
				<< " sender " << sendID << endl;
			cerr << " Attempting pump " << endl;
			pumpNetwork();
			cerr << " Pump returned ok " << endl;
		}
		assert(m_comm.receiveEmpty());
	}
	m_bubbles.clear();
	out.flush();
	assert(out.good());

	PrintDebug(0, "Removed %u bubbles\n", numPopped);
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
	m_comm.sendControlMessage(APC_DISCOVER_BUBBLES);

	unsigned numDiscovered = performNetworkDiscoverBubbles(this,
			opt::kmerSize);
	EndState();

	m_numReachedCheckpoint++;
	while (!checkpointReached())
		pumpNetwork();
	numDiscovered += m_checkpointSum;
	SetState(NAS_POPBUBBLE);
	if (numDiscovered > 0 && opt::verbose > 0)
		printf("Discovered %u bubbles\n", numDiscovered);
	return numDiscovered;
}

/** Pop the bubbles discovered previously. */
int NetworkSequenceCollection::controlPopBubbles(ostream& out)
{
	controlDiscoverBubbles();

	// Perform a round-robin bubble pop to avoid concurrency issues
	m_checkpointSum = performNetworkPopBubbles(out);
	EndState();

	// Now tell all the slave nodes to perform the pop one by one
	for(int i = 1; i < opt::numProc; i++) {
		m_comm.sendControlMessage(APC_BARRIER);
		m_comm.barrier();
		m_numReachedCheckpoint = 0;
		m_comm.sendControlMessageToNode(i, APC_POPBUBBLE,
				m_numPopped + m_checkpointSum);
		while (!checkpointReached(1))
			pumpNetwork();
	}

	unsigned numPopped = m_checkpointSum;
	m_numPopped += numPopped;
	if (numPopped > 0)
		printf("Removed %u bubbles\n", numPopped);
	return numPopped;
}

/** Mark ambiguous branches. */
unsigned NetworkSequenceCollection::controlMarkAmbiguous()
{
	puts("Marking ambiguous branches");
	assert(m_comm.receiveEmpty());
	unsigned count = m_comm.reduce(
			AssemblyAlgorithms::markAmbiguous(this));
	assert(m_comm.transmitBufferEmpty());
	assert(m_comm.receiveEmpty());
	m_comm.barrier();
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
	while (!checkpointReached())
		pumpNetwork();
	printf("Split %u ambiguous branches\n",
			m_checkpointSum);
	return m_checkpointSum;
}

/** Split ambiguous edges. */
unsigned NetworkSequenceCollection::controlSplit()
{
	m_comm.sendControlMessage(APC_SPLIT);
	m_comm.barrier();
	unsigned marked = controlMarkAmbiguous();
	unsigned split = controlSplitAmbiguous();
	assert(marked == split);
	return split;
}

/** Assemble a contig. */
void NetworkSequenceCollection::assembleContig(
		ISequenceCollection* seqCollection, FastaWriter* writer,
		BranchRecord& branch, unsigned id)
{
	unsigned removed = AssemblyAlgorithms::assembleContig(
			seqCollection, writer, branch, id);
	if (removed > 0) {
		m_lowCoverageContigs++;
		m_lowCoverageKmer += removed;
	}
}

namespace std {
	pair<unsigned, unsigned>& operator +=(pair<unsigned, unsigned>& a,
			pair<unsigned, unsigned> b)
	{
		a.first += b.first;
		a.second += b.second;
		return a;
	}
};

/** Assemble contigs.
 * @return the number of contigs and k-mer assembled
 */
pair<unsigned, unsigned> NetworkSequenceCollection::
performNetworkAssembly(ISequenceCollection* seqCollection,
		FastaWriter* fileWriter)
{
	Timer timer("NetworkAssembly");
	pair<unsigned, unsigned> numAssembled(0, 0);
	uint64_t branchGroupID = 0;
	assert(m_activeBranchGroups.empty());

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		extDirection dir;
		// dir will be set to the assembly direction if the sequence can be assembled
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(
				*iter, dir);
		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// Output the singleton contig.
			BranchRecord currBranch(SENSE, -1);
			currBranch.addSequence(*iter);
			currBranch.terminate(BS_NOEXT);
			assembleContig(seqCollection, fileWriter, currBranch,
					m_numAssembled + numAssembled.first);
			numAssembled.first++;
			numAssembled.second += currBranch.getLength();
			continue;
		}
		
		// Sequence is trimmable, create a new branch for it
		BranchGroup newGroup(branchGroupID, dir, 1, iter->first);
		BranchRecord newBranch(dir, -1);
		newGroup.addBranch(0, newBranch);
		m_activeBranchGroups[branchGroupID] = newGroup;

		// Generate the first extension request
		generateExtensionRequest(branchGroupID, 0, iter->first);
		branchGroupID++;

		numAssembled += processBranchesAssembly(seqCollection,
				fileWriter, numAssembled.first);
		seqCollection->pumpNetwork();

		if(m_activeBranchGroups.size() > MAX_ACTIVE)
		{
			while(m_activeBranchGroups.size() > LOW_ACTIVE)
			{
				seqCollection->pumpNetwork();
				numAssembled += processBranchesAssembly(seqCollection,
						fileWriter, numAssembled.first);
			}
		}
	}
	
	// Clear out the remaining branches
	while(!m_activeBranchGroups.empty())
	{
		numAssembled += processBranchesAssembly(seqCollection,
				fileWriter, numAssembled.first);
		seqCollection->pumpNetwork();
	}

	if (opt::coverage > 0) {
		PrintDebug(0, "Found %u k-mer in %u contigs before removing "
				"low-coverage contigs\n",
				numAssembled.second, numAssembled.first);
		PrintDebug(0, "Removed %u k-mer in %u low-coverage contigs\n",
				m_lowCoverageKmer, m_lowCoverageContigs);
	} else
		PrintDebug(0, "Assembled %u k-mer in %u contigs\n",
				numAssembled.second, numAssembled.first);
	return numAssembled;
}

/** Processes branches that are in progress, removing those that have
 * completed.
 * @return the number of contigs and k-mer assembled
 */
pair<unsigned, unsigned> NetworkSequenceCollection::
processBranchesAssembly(ISequenceCollection* seqCollection,
		FastaWriter* fileWriter, int currContigID)
{
	unsigned assembledContigs = 0, assembledKmer = 0;
	vector<BranchGroupMap::iterator> removeBranches;
	// Check if any of the current branches have gone inactive
	for(BranchGroupMap::iterator iter = m_activeBranchGroups.begin(); iter != m_activeBranchGroups.end(); iter++)
	{
		if(!iter->second.isActive())
		{
			// In this context, the group should have 1 and only 1 branch
			assert(iter->second.getNumBranches() == 1);
			
			// check if the branch is redundant, assemble if so, else it will simply be removed
			BranchRecord& currBranch = iter->second.getBranch(0);
			assert(currBranch.getState() == BS_NOEXT);
			if (currBranch.isCanonical()) {
				assembledContigs++;
				assembledKmer += currBranch.getLength();
				assembleContig(seqCollection, fileWriter, currBranch,
						m_numAssembled + currContigID++);
			}

			// Mark the group for removal
			removeBranches.push_back(iter);
		}	
	}

	// Remove all the finished branches
	for (vector<BranchGroupMap::iterator>::iterator rmIter
				= removeBranches.begin();
			rmIter != removeBranches.end(); rmIter++)
		m_activeBranchGroups.erase(*rmIter);	
	return make_pair(assembledContigs, assembledKmer);
}

//
// Generate a request for a sequence's extension, it will be handled in parseSequenceExtensionResponse
//
void NetworkSequenceCollection::generateExtensionRequest(
		uint64_t groupID, uint64_t branchID, const Kmer& seq)
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
		int nodeID = computeNodeID(seq);
		m_comm.sendSeqDataRequest(nodeID, groupID, branchID, seq);
	}
}

void NetworkSequenceCollection::processSequenceExtension(
		uint64_t groupID, uint64_t branchID, const Kmer& seq,
		const ExtensionRecord& extRec, int multiplicity)
{
	switch(m_state)
	{
		case NAS_TRIM:
		case NAS_ASSEMBLE:
		case NAS_COVERAGE:
			return processLinearSequenceExtension(groupID, branchID, seq, extRec, multiplicity);
		case NAS_DISCOVER_BUBBLES:
			return processSequenceExtensionPop(groupID, branchID, seq, extRec, multiplicity);
		case NAS_WAITING:
			if(m_finishedGroups.find(groupID) == m_finishedGroups.end())
			{
				cerr << "Unexpected sequence extension message! gid: " << groupID << " bid: " << branchID << " seq: " << seq.decode() << " Aborting...\n";
				assert(false);
			}
			break;
		default:
			cerr << "Unexpected sequence extension message! State: " << m_state << " gid: " << groupID << " bid: " << branchID << " seq: " << seq.decode() << " Aborting...\n";
			assert(false);
			break;
	}	
}

/** Process a sequence extension for trimming. */
void NetworkSequenceCollection::processLinearSequenceExtension(
		uint64_t groupID, uint64_t branchID, const Kmer& seq,
		const ExtensionRecord& extRec, int multiplicity)
{
	BranchGroupMap::iterator iter = m_activeBranchGroups.find(groupID);
	assert(iter != m_activeBranchGroups.end());
	Kmer currSeq = seq;
	bool active = AssemblyAlgorithms::processLinearExtensionForBranch(iter->second.getBranch(branchID), currSeq, extRec, multiplicity);
	if (active)
		generateExtensionRequest(groupID, branchID, currSeq);
}

/** Process a sequence extension for popping. */
void NetworkSequenceCollection::processSequenceExtensionPop(
		uint64_t groupID, uint64_t branchID, const Kmer& seq,
		const ExtensionRecord& extRec, int multiplicity)
{
	BranchGroupMap::iterator iter = m_activeBranchGroups.find(groupID);
	// If the iterator was not found we finished with that branch already, ensure this is so
	if(iter == m_activeBranchGroups.end())
	{
		assert(m_finishedGroups.find(groupID) != m_finishedGroups.end());
		// do nothing
		return;
	}

	Kmer currSeq = seq;
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
				generateExtensionRequest(groupID, i,
						iter->second.getBranch(i).getLastSeq());
			}
		}
	}
}


/** Add a k-mer to this collection. */
void NetworkSequenceCollection::add(const Kmer& seq)
{
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->add(seq);
		//PrintDebug(3, "done add\n");
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_comm.sendSeqOpMessage(nodeID, seq, MO_ADD);
	}
}

/** Remove a k-mer from this collection. */
void NetworkSequenceCollection::remove(const Kmer& seq)
{
	if(isLocal(seq))
	{
		m_pLocalSpace->remove(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);	
		m_comm.sendSeqOpMessage(nodeID, seq, MO_REMOVE);
	}
}

bool NetworkSequenceCollection::checkpointReached() const
{
	return m_numReachedCheckpoint == opt::numProc;
}

bool NetworkSequenceCollection::checkpointReached(int numRequired) const
{
	return m_numReachedCheckpoint == numRequired;
}

void NetworkSequenceCollection::setFlag(const Kmer& seq, SeqFlag flag)
{
	if(isLocal(seq))
	{
		//PrintDebug(3, "received local seq: %s\n", seq.decode().c_str());
		m_pLocalSpace->setFlag(seq, flag);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_comm.sendSetFlagMessage(nodeID, seq, flag);
	}
}

bool NetworkSequenceCollection::setBaseExtension(
		const Kmer& seq, extDirection dir, uint8_t base)
{
	if (isLocal(seq)) {
		if (m_pLocalSpace->setBaseExtension(seq, dir, base))
			m_numBasesAdjSet++;
	} else {
		int nodeID = computeNodeID(seq);
		m_comm.sendSetBaseExtension(nodeID, seq, dir, base);
	}

	// As this call delegates, the return value is meaningless so return false
	return false;
}

size_t NetworkSequenceCollection::count() const
{
	return m_pLocalSpace->count();
}

/** Remove the specified extensions from this k-mer. */
void NetworkSequenceCollection::removeExtension(
		const Kmer& seq, extDirection dir, SeqExt ext)
{
	if (isLocal(seq)) {
		m_pLocalSpace->removeExtension(seq, dir, ext);
		notify(seq);
	} else {
		int nodeID = computeNodeID(seq);
		m_comm.sendRemoveExtension(nodeID, seq, dir, ext);
	}
}

/** Return the data associated with the specified k-mer. */
bool NetworkSequenceCollection::getSeqData(const Kmer& seq,
		ExtensionRecord& extRecord, int& multiplicity) const
{
	assert(isLocal(seq));
	return m_pLocalSpace->getSeqData(seq, extRecord, multiplicity);
}

void NetworkSequenceCollection::wipeFlag(SeqFlag flag)
{
	m_pLocalSpace->wipeFlag(flag);
}

/** Return whether this sequence belongs to this process. */
bool NetworkSequenceCollection::isLocal(const Kmer& seq) const
{
	return computeNodeID(seq) == opt::rank;
}

/** Return the process ID to which the specified kmer belongs. */
int NetworkSequenceCollection::computeNodeID(const Kmer& seq) const
{
	return seq.getCode() % (unsigned)opt::numProc;
}
