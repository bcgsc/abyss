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
	BS_TOO_LONG // the branch is too long
	
};

typedef std::vector<PackedSeq> BranchData;
typedef std::set<PackedSeq> BranchSet;
typedef BranchData::iterator BranchDataIter;


class BranchRecord
{
	public:
		
		// Constructor
		BranchRecord();
		BranchRecord(extDirection dir, size_t maxLength);
		BranchRecord(const BranchRecord&);
		
		// Assignment Operator
		BranchRecord& operator=(const BranchRecord& other);		
		
		// Add a single sequence to the branch
		void addSequence(const PackedSeq& seq);
		
		// Terminate the branch and indicate why
		void terminate(BranchState reason);
		
		// Get the branch length
		size_t getLength() const;
		
		// Is the branch active?
		bool isActive() const;
		
		// Get the state of the branch
		BranchState getState() const;
		
		// Get the last sequence added to this branch
		const PackedSeq& getLastSeq() const;
		
		// Get the iterators
		BranchDataIter getStartIter();
		BranchDataIter getEndIter();
		
		// Get the direction of extension
		extDirection getDirection() const;
		
		// Return the maximum branch length
		size_t getMaxLength() const { return m_maxLength; }
		
		// check if a sequence exists in the branch record
		bool exists(const PackedSeq& seq) const;
		
		// Does this branch have a loop?
		bool hasLoop() const { return m_loopDetected; }
		
		
		
	private:
				
		// BranchData is used for the ordering/length of the branch, BranchSet is used for the existance of sequences in the branch.
		// They are populate simulataneously
		BranchData m_data;
		BranchSet m_set;
		extDirection m_dir;
		BranchState m_state;
		size_t m_maxLength;
		bool m_loopDetected;
};

#endif
