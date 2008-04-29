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

typedef std::list<PackedSeq> BranchData;
typedef BranchData::iterator BranchDataIter;


class BranchRecord
{
	public:
		
		// Constructor
		BranchRecord();
		BranchRecord(extDirection dir, size_t maxLength);
		BranchRecord(const BranchRecord&);
		
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
		
		// Get the iterators
		BranchDataIter getStartIter();
		BranchDataIter getEndIter();
		
		// Get the direction of extension
		extDirection getDirection() const;
		
		// Return the maximum branch length
		size_t getMaxLength() const { return m_maxLength; }
		
		
		
	private:
		BranchData m_data;
		extDirection m_dir;
		BranchState m_state;
		size_t m_maxLength;
};

#endif
