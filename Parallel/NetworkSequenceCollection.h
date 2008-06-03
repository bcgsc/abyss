#ifndef NETWORKSEQUENCECOLLECTION_H
#define NETWORKSEQUENCECOLLECTION_H

#include "PackedSeq.h"
#include "SequenceCollectionHash.h"
#include "CommLayer.h"
#include "SeqRecord.h"
#include "AssemblyAlgorithms.h"
#include "MessageBuffer.h"
#include "Log.h"
#include "Timer.h"

enum NetworkAssemblyState
{
	NAS_LOADING, // loading sequences
	NAS_FINALIZE, // finalizing the sequence data and getting ready for processing
	NAS_GEN_ADJ, // generating the sequence data
	NAS_ERODE, // erode the branch ends one sequence at a time
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
		NetworkSequenceCollection(int myID, int numDataNodes, int kmerSize, int readLen);
		~NetworkSequenceCollection();
		
		// This function operates in the same manner as AssemblyAlgorithms::GenerateAdjacency 
		// but has been rewritten to hide latency between nodes
		void networkGenerateAdjacency(ISequenceCollection* seqCollection);		
		
		// This function is similar to AssemblyAlgorithms::performNetworkTrim but is optimized to hide latency
		int performNetworkTrim(ISequenceCollection* seqCollection, int maxBranchCull);
		
		// This function is similar to AssemblyAlgorithms::popBubbles but is optimized to hide latency
		int performNetworkDiscoverBubbles(ISequenceCollection* seqCollection, int kmerSize);
		int performNetworkPopBubbles(ISequenceCollection* seqCollection);
		unsigned controlDiscoverBubbles();
		int controlPopBubbles();
		
		// Perform a network assembly
		unsigned performNetworkAssembly(ISequenceCollection* seqCollection, IFileWriter* fileWriter);

		// add a sequence to the collection
		void add(const PackedSeq& seq);
		
		// remove a sequence from the collection
		void remove(const PackedSeq& seq);
		
		// end the data load and make the sequence space ready for data read
		void finalize();
		
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
		void removeExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// check if the extension exists
		ResultPair checkExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// Set the extension of this sequence
		void setExtension(const PackedSeq& seq, extDirection dir, SeqExt extension);
		
		// set a single base extension
		bool setBaseExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// remove all the extensions of this sequence
		void clearExtensions(const PackedSeq& seq, extDirection dir);		
		
		// get the multiplicity of the sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// get the extensions of the sequence
		bool getSeqData(const PackedSeq& seq, ExtensionRecord& extRecord, int& multiplicity);
		
		// The loop to run the network code
		APResult pumpNetwork(int* pArg = NULL);
		
		// Loop over the pumping function while waiting for a result from the network
		ResultPair pumpUntilResult(); 
		
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
		
	private:
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
		bool isBranchRedundant(BranchRecord& branch);
		
		// Network message parsers
		int parseControlMessage();
		
		// Read a fasta file and distribute the sequences
		void readSequences(std::string fastaFile, int readLength, int kmerSize);
	
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
		
		// the size of k used in the assembly
		int m_kmer;
		
		// the original read length
		int m_readLen;
		
		// the number of bases of adjacency set
		int m_numBasesAdjSet;
		
		// the starting trim value
		int m_startTrimLen;
		
		// the current length to trim on (comes from the control node)
		int m_trimStep;
		
		// the number of sequences assembled so far
		int m_numAssembled;
		
		// The number of requests that are pending
		size_t m_numOutstandingRequests;
		
		// The current branches that are active
		BranchGroupMap m_activeBranchGroups;

		/** Bubbles, which are branch groups that have joined. */
		BranchGroupMap m_bubbles;
		
		// List of IDs of finished groups, used for sanity checking during bubble popping
		std::set<uint64_t> m_finishedGroups;
		
		// Message buffer
		MessageBuffer* m_pMsgBuffer;
		
		// Log file
		Log* m_pLog;
		
		// Timer for the entire lifetime of the object
		Timer m_timer;
		
		static const size_t MAX_ACTIVE = 50;
		static const size_t LOW_ACTIVE = 10;		
};

#endif
