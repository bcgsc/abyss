#ifndef NETWORKSEQUENCECOLLECTION_H
#define NETWORKSEQUENCECOLLECTION_H 1

#include "SequenceCollection.h"
#include "BranchGroup.h"
#include "BranchRecord.h"
#include "CommLayer.h"
#include "FastaWriter.h"
#include "MessageBuffer.h"
#include "Timer.h"
#include <ostream>
#include <set>
#include <utility>

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
	NAS_REMOVE_MARKED, // remove marked sequences
	NAS_COVERAGE, // remove low-coverage contigs
	NAS_COVERAGE_COMPLETE,
	NAS_DISCOVER_BUBBLES, // discover read errors/SNPs
	NAS_POPBUBBLE, // remove read errors/SNPs
	NAS_MARK_AMBIGUOUS, // mark ambiguous branches
	NAS_SPLIT_AMBIGUOUS, // split ambiguous branches
	NAS_CLEAR_FLAGS, // clear the flags
	NAS_ASSEMBLE, // assembling the data
	NAS_ASSEMBLE_COMPLETE, // assembling is complete
	NAS_WAITING, // non-control process is waiting, this just loops over the network function
	NAS_DONE // finished, clean up and exit
};

typedef std::map<uint64_t, BranchGroup> BranchGroupMap;

class NetworkSequenceCollection : public ISequenceCollection
{
	public:
		NetworkSequenceCollection();
		~NetworkSequenceCollection();

		// This function is similar to AssemblyAlgorithms::performNetworkTrim but is optimized to hide latency
		int performNetworkTrim(ISequenceCollection* seqCollection, int maxBranchCull);

		int performNetworkDiscoverBubbles(ISequenceCollection* c);
		int performNetworkPopBubbles(std::ostream& out);

		unsigned controlErode();
		unsigned controlTrimRound(unsigned trimLen);
		void controlTrim(unsigned start = 1);
		unsigned controlRemoveMarked();
		void controlCoverage();
		unsigned controlDiscoverBubbles();
		int controlPopBubbles(std::ostream& out);
		unsigned controlMarkAmbiguous();
		unsigned controlSplitAmbiguous();
		unsigned controlSplit();

		// Perform a network assembly
		std::pair<unsigned, unsigned> performNetworkAssembly(
				ISequenceCollection* seqCollection,
				FastaWriter* fileWriter = NULL);

		void add(const Kmer& seq);
		void remove(const Kmer& seq);
		void setFlag(const Kmer& seq, SeqFlag flag);
		void wipeFlag(SeqFlag flag);

		// Return the number of sequences in the collection
		size_t count() const;

		void printLoad() const { m_pLocalSpace->printLoad(); }

		void removeExtension(const Kmer& seq, extDirection dir,
				SeqExt ext);
		bool setBaseExtension(const Kmer& seq, extDirection dir,
				uint8_t base);

		bool getSeqData(const Kmer& seq,
				ExtensionRecord& extRecord, int& multiplicity) const;
		const value_type& getSeqAndData(const Kmer& key) const
		{
			assert(isLocal(key));
			return m_pLocalSpace->getSeqAndData(key);
		}

		// Receive and dispatch packets.
		unsigned pumpNetwork();
		unsigned pumpFlushReduce();

		void completeOperation();

		// run the assembly
		void run();
		
		// run the assembly from the controller's point of view (rank 0 node)
		void runControl();
		
		// test if the checkpoint has been reached
		bool checkpointReached() const;
		bool checkpointReached(int numRequired) const;

		void handle(int senderID, const SeqAddMessage& message);
		void handle(int senderID, const SeqRemoveMessage& message);
		void handle(int senderID, const SetBaseMessage& message);
		void handle(int senderID, const SetFlagMessage& message);
		void handle(int senderID, const RemoveExtensionMessage& message);
		void handle(int senderID, const SeqDataRequest& message);
		void handle(int senderID, const SeqDataResponse& message);

		// Observer pattern, not implemented.
		void attach(SeqObserver f) { (void)f; }
		void detach(SeqObserver f) { (void)f; }

		/** Load this collection from disk. */
		void load(const char *path)
		{
			m_pLocalSpace->load(path);
		}

		/** Indicate that this is a colour-space collection. */
		void setColourSpace(bool flag)
		{
			m_pLocalSpace->setColourSpace(flag);
			m_comm.broadcast(flag);
		}

		iterator begin() { return m_pLocalSpace->begin(); }
		const_iterator begin() const { return m_pLocalSpace->begin(); }
		iterator end() { return m_pLocalSpace->end(); }
		const_iterator end() const { return m_pLocalSpace->end(); }

	private:
		// Observer pattern
		void notify(const Kmer& seq);

		void loadSequences();

		std::pair<unsigned, unsigned> processBranchesAssembly(
				ISequenceCollection* seqCollection,
				FastaWriter* fileWriter, int currContigID);
		int processBranchesTrim();
		bool processBranchesDiscoverBubbles();

		void generateExtensionRequest(
				uint64_t groupID, uint64_t branchID, const Kmer& seq);
		void generateExtensionRequests(uint64_t groupID,
				BranchGroup::const_iterator first,
				BranchGroup::const_iterator last);
		void processSequenceExtension(
				uint64_t groupID, uint64_t branchID, const Kmer& seq,
				const ExtensionRecord& extRec, int multiplicity);
		void processLinearSequenceExtension(
				uint64_t groupID, uint64_t branchID, const Kmer& seq,
				const ExtensionRecord& extRec, int multiplicity);
		void processSequenceExtensionPop(
				uint64_t groupID, uint64_t branchID, const Kmer& seq,
				const ExtensionRecord& extRec, int multiplicity);

		void assembleContig(ISequenceCollection* seqCollection,
				FastaWriter* fileWriter,
				BranchRecord& branch, unsigned id);

		// Check if a branch is redundant with a previously output branch
		bool isBranchRedundant(const BranchRecord& branch);

		void parseControlMessage(int source);

		bool isLocal(const Kmer& seq) const;
		int computeNodeID(const Kmer& seq) const;

		void EndState();
		
		// Set the state of the network assembly
		void SetState(NetworkAssemblyState newState);
		
		// Pointer to the local sequence space
		// These sequences are held in memory on this process
		SequenceCollectionHash* m_pLocalSpace;

		// The communications layer implements the functions over the network
		MessageBuffer m_comm;

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

		/** The number of low-coverage contigs removed. */
		unsigned m_lowCoverageContigs;

		/** The number of low-coverage k-mer removed. */
		unsigned m_lowCoverageKmer;

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

		static const size_t MAX_ACTIVE = 50;
		static const size_t LOW_ACTIVE = 10;		
};

#endif
