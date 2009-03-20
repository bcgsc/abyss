#ifndef NETWORKSEQUENCECOLLECTION_H
#define NETWORKSEQUENCECOLLECTION_H 1

#include "PackedSeq.h"
#include "SequenceCollectionHash.h"
#include "BranchGroup.h"
#include "BranchRecord.h"
#include "CommLayer.h"
#include "IFileWriter.h"
#include "MessageBuffer.h"
#include "Timer.h"
#include <set>

enum NetworkAssemblyState
{
	NAS_LOADING, // loading sequences
	NAS_LOAD_COMPLETE, // loading is complete
	NAS_GEN_ADJ, // generating the sequence data
	NAS_ADJ_COMPLETE, // adjacency generation is complete
	NAS_ERODE, // erode the branch ends one sequence at a time
	NAS_ERODE_WAITING,
	NAS_ERODE_COMPLETE,
	NAS_TRIM, // trimming the data
	NAS_DISCOVER_BUBBLES, // discover read errors/SNPs
	NAS_POPBUBBLE, // remove read errors/SNPs
	NAS_TRIM2, // second trimming step after the bubble removal
	NAS_SPLIT, // remove ambiguous links (just before assembling into contigs)
	NAS_ASSEMBLE, // assembling the data
	NAS_WAITING, // non-control process is waiting, this just loops over the network function
	NAS_DONE // finished, clean up and exit
};

typedef std::map<uint64_t, BranchGroup> BranchGroupMap;

class NetworkSequenceCollection : public ISequenceCollection
{
	public:
	
		// Constructor/destructor
		NetworkSequenceCollection(int myID, int numDataNodes);
		~NetworkSequenceCollection();
		
		// This function is similar to AssemblyAlgorithms::performNetworkTrim but is optimized to hide latency
		int performNetworkTrim(ISequenceCollection* seqCollection, int maxBranchCull);
		
		// This function is similar to AssemblyAlgorithms::popBubbles but is optimized to hide latency
		int performNetworkDiscoverBubbles(ISequenceCollection* seqCollection, int kmerSize);
		int performNetworkPopBubbles(ISequenceCollection* seqCollection);

		unsigned controlErode();
		unsigned controlDiscoverBubbles();
		int controlPopBubbles();

		// Perform a network assembly
		unsigned performNetworkAssembly(ISequenceCollection* seqCollection, IFileWriter* fileWriter);

		// add a sequence to the collection
		void add(const PackedSeq& seq);
		
		// remove a sequence from the collection
		void remove(const PackedSeq& seq);
		
		// check if a sequence exists
		bool exists(const PackedSeq& seq);
		
		// Set flag for sequence seq
		void setFlag(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set
		bool checkFlag(const PackedSeq& seq, SeqFlag flag);
		
		// Clear the specified flag from every sequence in the collection
		void wipeFlag(SeqFlag flag);
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq);

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq);
		
		// Return the number of sequences in the collection
		int count() const;
		
		// remove the extension to the sequence
		bool removeExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// check if the extension exists
		ResultPair checkExtension(const PackedSeq& seq, extDirection dir, char base);

		// set a single base extension
		bool setBaseExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// remove all the extensions of this sequence
		void clearExtensions(const PackedSeq& seq, extDirection dir);		
		
		// get the multiplicity of the sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// get the extensions of the sequence
		bool getSeqData(const PackedSeq& seq, ExtensionRecord& extRecord, int& multiplicity);

		// Receive and dispatch packets.
		unsigned pumpNetwork();
		unsigned pumpFlushReduce();

		void completeOperation();

		// run the assembly
		void run();
		
		// run the assembly from the controller's point of view (rank 0 node)
		void runControl();
		
		// test if the checkpoint has been reached
		bool checkpointReached(int numRequired) const;
		
		// get an iterator to the first sequence
		SequenceCollectionIterator getStartIter() const;
		
		// get an iterator to the last sequence
		SequenceCollectionIterator getEndIter() const;
		
		// Message handlers, polymorphically called by the message types
		void handleSeqOpMessage(int senderID, const SeqOpMessage& seqMsg);
		void handleSetBaseMessage(int senderID, const SetBaseMessage& message);
		void handleSetFlagMessage(int senderID, const SetFlagMessage& message);
		void handleRemoveExtensionMessage(int senderID, const RemoveExtensionMessage& message);
		void handleSequenceDataRequest(int senderID, SeqDataRequest& message);
		void handleSequenceDataResponse(int senderID, SeqDataResponse& message);		
		
		// Observer pattern, not implemented.
		virtual void attach(SeqObserver f) { (void)f; }
		virtual void detach(SeqObserver f) { (void)f; }
	private:
		// Observer pattern
		void notify(const PackedSeq& seq);

		void loadSequences();
	
		// Branch processing
		int processBranchesTrim();
		bool processBranchesDiscoverBubbles();
		int processBranchesAssembly(ISequenceCollection* seqCollection, IFileWriter* fileWriter, int currContigID);
		
		void generateExtensionRequest(uint64_t groupID, uint64_t branchID, const PackedSeq& seq);
		void processSequenceExtension(uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec, int multiplicity);
		void processLinearSequenceExtension(uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec, int multiplicity);
		void processSequenceExtensionPop(uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec, int multiplicity);
		
		// Check if a branch is redundant with a previously output branch
		bool isBranchRedundant(const BranchRecord& branch);
		
		// Network message parsers
		void parseControlMessage();
		
		// Check if this sequence belongs in the local phase space
		bool isLocal(const PackedSeq& seq) const;
		
		// Get the node id that this sequence should reside in
		int computeNodeID(const PackedSeq& seq) const;
		
		// Clean up the state
		void EndState();
		
		// Set the state of the network assembly
		void SetState(NetworkAssemblyState newState);
		
		// Pointer to the local sequence space
		// These sequences are held in memory on this process
		SequenceCollectionHash* m_pLocalSpace;
		
		// The communications layer implements the functions over the network
		CommLayer* m_pComm;
		
		// The ID for this node, assigned by MPI
		int m_id;
		
		// The number of nodes in the network
		unsigned int m_numDataNodes;
		
		// the state of the assembly
		NetworkAssemblyState m_state;
		
		// The number of processes that have sent a checkpoint reached message, this is used by the control process to determine the state flow
		int m_numReachedCheckpoint;

		/** The sum of the values returned by the slave nodes in their
		 * checkpoint messages.
		 */
		int m_checkpointSum;
		
		// the number of bases of adjacency set
		int m_numBasesAdjSet;

		// the current length to trim on (comes from the control node)
		int m_trimStep;
		
		/** The number of bubbles popped so far. */
		unsigned m_numPopped;

		// the number of sequences assembled so far
		int m_numAssembled;
		
		// The current branches that are active
		BranchGroupMap m_activeBranchGroups;

		/** Bubbles, which are branch groups that have joined. */
		BranchGroupMap m_bubbles;
		
		// List of IDs of finished groups, used for sanity checking during bubble popping
		std::set<uint64_t> m_finishedGroups;
		
		// Message buffer
		MessageBuffer* m_pMsgBuffer;

		static const size_t MAX_ACTIVE = 50;
		static const size_t LOW_ACTIVE = 10;		
};

#endif
