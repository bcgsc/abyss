#ifndef NETWORKSEQUENCECOLLECTION_H
#define NETWORKSEQUENCECOLLECTION_H 1

#include "SequenceCollection.h"
#include "BranchRecord.h"
#include "BranchGroup.h"
#include "CommLayer.h"
#include "FastaWriter.h"
#include "MessageBuffer.h"
#include "Timer.h"
#include <ostream>
#include <set>
#include <utility>
#if _SQL
#include "Common/InsOrderedMap.h"

namespace NSC
{
	typedef InsOrderedMap<std::string, int> dbMap;
	dbMap moveFromAaStatMap();
}
#endif

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
	NAS_WAITING, // non-control process is waiting
	NAS_DONE // finished, clean up and exit
};

typedef std::map<uint64_t, BranchGroup> BranchGroupMap;

/** A distributed map of vertices to vertex properties. */
class NetworkSequenceCollection
{
	public:
		typedef SequenceDataHash::key_type V;
		typedef SequenceDataHash::mapped_type VP;

		typedef SequenceDataHash::key_type key_type;
		typedef SequenceDataHash::mapped_type mapped_type;
		typedef SequenceDataHash::value_type value_type;
		typedef SequenceDataHash::iterator iterator;
		typedef SequenceDataHash::const_iterator const_iterator;

		typedef mapped_type::Symbol Symbol;
		typedef mapped_type::SymbolSet SymbolSet;
		typedef mapped_type::SymbolSetPair SymbolSetPair;

		typedef key_type vertex_descriptor;
		typedef mapped_type vertex_bundled;
		typedef std::pair<V, V> edge_descriptor;

		typedef boost::directed_tag directed_category;
		typedef boost::disallow_parallel_edge_tag edge_parallel_category;
		struct traversal_category
			: boost::adjacency_graph_tag, boost::vertex_list_graph_tag
			{ };

		NetworkSequenceCollection()
			: m_state(NAS_WAITING), m_trimStep(0),
			m_numPopped(0), m_numAssembled(0) { }

		size_t performNetworkTrim();

		size_t performNetworkDiscoverBubbles();
		size_t performNetworkPopBubbles(std::ostream& out);

		size_t controlErode();
		size_t controlTrimRound(unsigned trimLen);
		void controlTrim(unsigned&, unsigned start = 1);
		size_t controlRemoveMarked();
		void controlCoverage();
		size_t controlDiscoverBubbles();
		size_t controlPopBubbles(std::ostream& out);
		size_t controlMarkAmbiguous();
		size_t controlSplitAmbiguous();
		size_t controlSplit();

		// Perform a network assembly
		std::pair<size_t, size_t> performNetworkAssembly(
				FastaWriter* fileWriter = NULL);

		void add(const V& seq, unsigned coverage = 1);
		void remove(const V& seq);
		void setFlag(const V& seq, SeqFlag flag);

		/** Mark the specified sequence in both directions. */
		void mark(const V& seq)
		{
			setFlag(seq, SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));
		}

		/** Mark the specified sequence. */
		void mark(const V& seq, extDirection sense)
		{
			setFlag(seq, sense == SENSE
					? SF_MARK_SENSE : SF_MARK_ANTISENSE);
		}

		/** Return true if this container is empty. */
		bool empty() const { return m_data.empty(); }

		void printLoad() const { m_data.printLoad(); }

		void removeExtension(const V& seq, extDirection dir,
				SymbolSet ext);

		/** Remove the specified edge of this vertex. */
		void removeExtension(const V& seq, extDirection dir, Symbol base)
		{
			removeExtension(seq, dir, SymbolSet(base));
		}

		bool setBaseExtension(const V& seq, extDirection dir,
				uint8_t base);

		// Receive and dispatch packets.
		size_t pumpNetwork();
		size_t pumpFlushReduce();

		void completeOperation();

		// run the assembly
		void run();

		// run the assembly from the controller's point of view
		void runControl();

		// test if the checkpoint has been reached
		bool checkpointReached() const;
		bool checkpointReached(unsigned numRequired) const;

		void handle(int senderID, const SeqAddMessage& message);
		void handle(int senderID, const SeqRemoveMessage& message);
		void handle(int senderID, const SetBaseMessage& message);
		void handle(int senderID, const SetFlagMessage& message);
		void handle(int senderID, const RemoveExtensionMessage& m);
		void handle(int senderID, const SeqDataRequest& message);
		void handle(int senderID, const SeqDataResponse& message);

		/** The observer callback function. */
		typedef void (*SeqObserver)(SequenceCollectionHash* c,
				const value_type& seq);

		// Observer pattern, not implemented.
		void attach(SeqObserver) { }
		void detach(SeqObserver) { }

		/** Load this collection from disk. */
		void load(const char *path)
		{
			m_data.load(path);
		}

		/** Indicate that this is a colour-space collection. */
		void setColourSpace(bool flag)
		{
			m_data.setColourSpace(flag);
			m_comm.broadcast(flag);
		}

		iterator begin() { return m_data.begin(); }
		const_iterator begin() const { return m_data.begin(); }
		iterator end() { return m_data.end(); }
		const_iterator end() const { return m_data.end(); }

	private:
		// Observer pattern
		void notify(const V& seq);

		void loadSequences();

		std::pair<size_t, size_t> processBranchesAssembly(
				FastaWriter* fileWriter, unsigned currContigID);
		size_t processBranchesTrim();
		bool processBranchesDiscoverBubbles();

		void generateExtensionRequest(
				uint64_t groupID, uint64_t branchID, const V& seq);
		void generateExtensionRequests(uint64_t groupID,
				BranchGroup::const_iterator first,
				BranchGroup::const_iterator last);
		void processSequenceExtension(
				uint64_t groupID, uint64_t branchID, const V& seq,
				const SymbolSetPair& extRec, int multiplicity);
		void processLinearSequenceExtension(
				uint64_t groupID, uint64_t branchID, const V& seq,
				const SymbolSetPair& extRec, int multiplicity,
				unsigned maxLength);
		void processSequenceExtensionPop(
				uint64_t groupID, uint64_t branchID, const V& seq,
				const SymbolSetPair& extRec, int multiplicity,
				unsigned maxLength);

		void assembleContig(
				FastaWriter* fileWriter,
				BranchRecord& branch, unsigned id);

		// Check if a branch is redundant with a previously output
		// branch.
		bool isBranchRedundant(const BranchRecord& branch);

		void parseControlMessage(int source);

		bool isLocal(const V& seq) const;
		int computeNodeID(const V& seq) const;

		void EndState();

		// Set the state of the network assembly
		void SetState(NetworkAssemblyState newState);

		SequenceCollectionHash m_data;

		// The communications layer implements the functions over the
		// network.
		MessageBuffer m_comm;

		// The number of nodes in the network
		unsigned m_numDataNodes;

		// the state of the assembly
		NetworkAssemblyState m_state;

		// The number of processes that have sent a checkpoint reached
		// message, this is used by the control process to determine
		// the state flow.
		unsigned m_numReachedCheckpoint;

		/** The sum of the values returned by the slave nodes in their
		 * checkpoint messages.
		 */
		size_t m_checkpointSum;

		// the number of bases of adjacency set
		size_t m_numBasesAdjSet;

		// the current length to trim on (comes from the control node)
		unsigned m_trimStep;

		/** The number of low-coverage contigs removed. */
		size_t m_lowCoverageContigs;

		/** The number of low-coverage k-mer removed. */
		size_t m_lowCoverageKmer;

		/** The number of bubbles popped so far. */
		size_t m_numPopped;

		// the number of sequences assembled so far
		size_t m_numAssembled;

		// The current branches that are active
		BranchGroupMap m_activeBranchGroups;

		/** Bubbles, which are branch groups that have joined. */
		BranchGroupMap m_bubbles;

		// List of IDs of finished groups, used for sanity checking
		// during bubble popping.
		std::set<uint64_t> m_finishedGroups;

		static const size_t MAX_ACTIVE = 50;
		static const size_t LOW_ACTIVE = 10;
};

#endif
