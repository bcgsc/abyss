#ifndef BRANCHRECORD_H
#define BRANCHRECORD_H

#include "CommonDefs.h"
#include "PackedSeq.h"

enum BranchState
{
	
	BS_ACTIVE, // can be extended
	BS_NOEXT, // the branch has ended because of a lack of sequence to extend to
	BS_AMBI_SAME, // the branch has ended because the extension from this branch is ambigious
	BS_AMBI_OPP, // the branch has ended because the extension to this branch is ambigiuous
	BS_TOO_LONG, // the branch is too long
	BS_LOOP // the branch has a loop
	
};

typedef std::pair<PackedSeq, int> MultMapPair;
typedef std::vector<PackedSeq> BranchData;
typedef std::map<PackedSeq, int> BranchMultMap;
typedef BranchData::iterator BranchDataIter;


class BranchRecord
{
	public:
		
		// Constructor
		BranchRecord();
		BranchRecord(extDirection dir, int maxLength);
		BranchRecord(const BranchRecord&);
		
		// Assignment Operator
		BranchRecord& operator=(const BranchRecord& other);		
		
		// Add a single sequence to the branch
		void addSequence(const PackedSeq& seq, int multiplicity = -1);
		
		// Terminate the branch and indicate why
		void terminate(BranchState reason);
		
		// Get the branch length
		size_t getLength() const;
		
		// Is the branch active?
		bool isActive() const;
		
		// Get the state of the branch
		BranchState getState() const;
		
		// get the multiplicity of a sequence
		int getMultiplicity(const PackedSeq& seq) const;
		
		// Set the multiplicity of a sequence in the branch
		void setMultiplicity(const PackedSeq& seq, int multiplicity);
		
		/** Forget the multiplicity information. */
		void BranchRecord::clearMultiplicity();

		// Set the state of the branch
		void setState(BranchState state) { m_state = state; }
		
		// Check if the branch is empty
		bool empty() const { return m_data.empty(); }
		
		// Get the first sequence added to the branch
		const PackedSeq& getFirstSeq() const;
		
		// Get the last sequence added to this branch
		const PackedSeq& getLastSeq() const;
		
		const PackedSeq& getSeqByIndex(size_t index) const;
		
		// Get the iterators
		BranchDataIter getStartIter();
		BranchDataIter getEndIter();
		
		// Get the direction of extension
		extDirection getDirection() const;
		
		// Return the maximum branch length
		size_t getMaxLength() const { return m_maxLength; }
		
		// build a contig from the branch
		void buildContig(Sequence& outseq) const;
		
		// check if a sequence exists in the branch record
		bool exists(const PackedSeq& seq) const;
		
		// should the length of the branch be checked?
		bool doLengthCheck() const;
		
		/** Return whether the branch is too long. */
		bool isTooLong() const;
		
		// Does this branch have a loop?
		bool hasLoop() const { return m_loopDetected; }
		
		// Calculate the total multiplicity for this branch.
		int calculateBranchMultiplicity(bool ignorelast = false);

		// Return the precalculated multiplicity for this branch.
		int getBranchMultiplicity() const
		{
			assert(m_multiplicity > 0);
			return m_multiplicity;
		}

		bool isCanonical() const;
		
		// Print out the branch
		void printBranch(std::ostream& ostr);

	private:
				
		// BranchData is used for the ordering/length of the branch, BranchSet is used for the existance of sequences in the branch.
		// They are populate simulataneously
		BranchData m_data;
		BranchMultMap m_seqMap;
		extDirection m_dir;
		BranchState m_state;
		int m_maxLength;
		bool m_loopDetected;
		int m_multiplicity;
};

#endif
