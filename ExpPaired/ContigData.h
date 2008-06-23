#ifndef CONTIGDATA_H
#define CONTIGDATA_H

#include "CommonDefs.h"
#include "PackedSeq.h"
#include "AlignmentCache.h"
#include "PairRecord.h"
#include "Stats.h"

// Flags to indicate what sequences are desired
// By default there are no filters and every pair is returned
// Flags can by OR'd to filter on multiple groups
enum PairSelectFilter
{
	PSF_ALL = 0x0,
	PSF_USABLE_ONLY = 0x1,
	PSF_UNRESOLVED_ONLY = 0x2,
	PSF_UNIQUE_ALIGN_ONLY = 0x4
};

// Forward declare
struct PairedResolvePolicy;

struct SeqAlignment
{
	PackedSeq seq;
	AlignData alignment;	
};

struct AlignmentPair
{
	SeqAlignment refSeqAlign;
	SeqAlignment pairSeqAlign;

	friend std::ostream& operator<<(std::ostream& out, const AlignmentPair& object)
	{
		out << "Seq Align: " << object.refSeqAlign.alignment << " Pair Align: " << object.pairSeqAlign.alignment;
		return out;
	}	
};

typedef std::vector<AlignmentPair> PairAlignments;
typedef std::vector<PairAlignments> PairAlignmentsCollection;

// Information about whether a pair was resolved or not
struct ResolvedData
{
	public:
		ResolvedData() : m_isResolved(false), m_resolvedID("") {}
		void setResolved(ContigID id) { m_isResolved = true; m_resolvedID = id; }
		bool isResolved() const { return m_isResolved; }
		ContigID getResolvedID() const { return m_resolvedID; }
	
	private:
		bool m_isResolved;
		ContigID m_resolvedID;
};

struct PairData
{
	PackedSeq seq;
	ResolvedData resolvedData;
};

typedef std::vector<PairData> PairVector;
struct ContigKmerData
{
	PackedSeq seq;
	bool usable;
	
	// Keep a record of all the pairs of each sequence, one for the sequence and one for its reverse complement
	PairVector pairs[2];
};

typedef std::map<ContigID, int> ContigSupportMap;
typedef std::vector<ContigKmerData> CKDataVec;

class ContigData
{
	public:
		
		ContigData(const ContigID& id, const Sequence& s, size_t kmer, int copyNumber, AlignmentCache* pDB);
		
		// Merge the data in other into this data, taking into consideration directionality and complement
		// The isUsable flag is used to detemine whether the other contig's pairs will be usable for resolution (the contigs are equivalent
		// in copy number) or not (the child contig is a repeat element that is only used to bridge sequences)
		void merge(const ContigData& other, extDirection dir, bool isReversed, bool isUsable);
		
		// Add a contig id to the id set, indicating this contig was concataneted into this sequence
		void addID(const ContigID& id, extDirection dir) { m_idParts.insert(id); m_mergeRecord[dir].push_back(id); }
		
		// add a group of sequences
		void addSeqs(CKDataVec::iterator position, CKDataVec& newSeqs);
		
		// update the alignment db with this contig's data
		void writeToAlignDB(AlignmentCache* pDB) const;
		void removeFromAlignDB(AlignmentCache* pDB) const;
		
		// Copy out the sequence set or usable set
		void copySeqSet(PSeqSet& outSeqs) const;
		void copyUsableSeqSet(PSeqSet& outSeqs) const;
		
		// Get the ids of all the sequences that have been added
		ContigIDSet getExclusionSet() const { return m_idParts; }
		
		// get the sequence of this contig
		Sequence getSeq() const { return m_seq; }
		
		// get the contig length
		size_t getLength() const { return m_seq.length(); }
		
		// validate the kmer seqs against the reference
		void validateSequences() const;
		
		// validate that the database has correct sequences for this contig
		void validateDatabase(AlignmentCache* pDB) const;
		
		// add the pairs to this contig
		void addPairs(PairRecord* pairRecord);
		
		// add the distance between pairs that both reside on this contig to the histogram
		void addSelfPairsToHist(Histogram* pHist);
		
		// Get the number of kmers that make up the given sequence
		size_t getNumKmers(const Sequence& seq) { return (seq.length() - m_kmer + 1); } 
		
		// print all the pair alignments from this contig
		void printPairAlignments(extDirection dir, unsigned int filter) const;
		
		// Convert a direction to a pair index
		static size_t dir2Idx(extDirection dir) { return (dir == SENSE) ? 0 : 1; }
		
		// Resolve the pairs on this contig using the resolve policy passed in
		void resolvePairs(PairedResolvePolicy* pResolvePolicy, extDirection dir);
		
		// Extract the alignments of all the unresolved pairs in this direction
		void extractAlignments(extDirection dir, unsigned int filter, PairAlignmentsCollection& allAlignments) const;
		
		// Get the list of contigs that are supported by pairs in the specified direction
		void getSupportMap(extDirection dir, unsigned int filter, ContigSupportMap& outMap) const;
		
		// Convert a support map to an ID set
		static void supportMap2IDSet(const ContigSupportMap& inMap, ContigIDSet& outSet);
		
		// Compute a test stat about this contig, this is a debug function
		void computeTestStat(const PDF& empDist);
		
		// Get the id of this contig
		ContigID getID() const { return m_id; }
		
		
	private:
		
		// Get the alignment pairs for the two specified sequences
		void getAlignmentPairs(const PackedSeq& refSeq, const PackedSeq& pairSeq, unsigned int filter, PairAlignments& outAligns) const;
		
		// append a sequence into this sequence
		void appendSeqString(const Sequence& seq, extDirection dir, bool isReversed);
		
		// The sequence for the contig
		Sequence m_seq;
		
		// The vector of kmers that make up this contig
		CKDataVec m_kmerVec;
		
		// The record of initial contigs that made up this contig
		ContigID m_id;
		ContigIDSet m_idParts;
		
		// the kmer size
		size_t m_kmer;
		
		// the infered copy number for this sequence, 
		int m_copyNumber;
		
		// Merge record, for debug
		ContigIDVec m_mergeRecord[NUM_DIRECTIONS];
		
		// The alignment cache database
		AlignmentCache* m_pDatabase;
};

#endif
