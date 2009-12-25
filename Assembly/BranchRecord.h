#ifndef BRANCHRECORD_H
#define BRANCHRECORD_H 1

#include "PackedSeq.h"
#include <set>
#include <vector>

enum BranchState
{
	
	BS_ACTIVE, // can be extended
	BS_NOEXT, // the branch has ended because of a lack of sequence to extend to
	BS_AMBI_SAME, // the branch has ended because the extension from this branch is ambigious
	BS_AMBI_OPP, // the branch has ended because the extension to this branch is ambigiuous
	BS_TOO_LONG, // the branch is too long
	BS_LOOP // the branch has a loop
	
};

class BranchRecord
{
	public:
		typedef std::vector<PackedSeq> BranchData;
		typedef std::set<Kmer> BranchSet;
		typedef BranchData::iterator iterator;
		typedef BranchData::const_iterator const_iterator;

		BranchRecord()
			: m_dir(SENSE), m_state(BS_ACTIVE), m_maxLength(-1),
			m_loopDetected(false), m_multiplicity(-1) { }
		BranchRecord(extDirection dir, int maxLength)
			: m_dir(dir), m_state(BS_ACTIVE), m_maxLength(maxLength),
			m_loopDetected(false), m_multiplicity(-1) { }

		operator Sequence() const;

		void addSequence(const PackedSeq& seq);

		/** Add a k-mer to the branch without data. */
		void addSequence(const Kmer& kmer)
		{
			addSequence(std::make_pair(kmer, KmerData()));
		}

		// Remove all the sequences including and following the
		// specified iterator.
		void truncate(iterator position);

		// Terminate the branch and indicate why
		void terminate(BranchState reason);
		
		// Get the branch length
		size_t getLength() const;
		
		// Is the branch active?
		bool isActive() const;
		
		// Get the state of the branch
		BranchState getState() const;

		// Set the data of a sequence in the branch.
		void setData(const PackedSeq& seq);

		/** Forget the multiplicity information. */
		void clearMultiplicity();

		// Get the first sequence added to the branch
		const Kmer& getFirstSeq() const;

		// Get the last sequence added to this branch
		const Kmer& getLastSeq() const;

		iterator begin() { return m_data.begin(); }
		iterator end() { return m_data.end(); }
		const_iterator begin() const { return m_data.begin(); }
		const_iterator end() const { return m_data.end(); }

		// Get the direction of extension
		extDirection getDirection() const;
		
		// Return the maximum branch length
		size_t getMaxLength() const { return m_maxLength; }

		// check if a sequence exists in the branch record
		bool exists(const Kmer& seq) const;

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

	private:
				
		// BranchData is used for the ordering/length of the branch, BranchSet is used for the existance of sequences in the branch.
		// They are populate simulataneously
		BranchData m_data;
		BranchSet m_seqMap;
		extDirection m_dir;
		BranchState m_state;
		int m_maxLength;
		bool m_loopDetected;
		int m_multiplicity;
};

#endif
