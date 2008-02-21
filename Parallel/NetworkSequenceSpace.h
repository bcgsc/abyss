#ifndef NETWORKSEQUENCESPACE_H
#define NETWORKSEQUENCESPACE_H

#include "PackedSeq.h"
#include "SimpleSequenceSpace.h"
#include "CommLayer.h"

class NetworkSequenceSpace
{
	public:
	
		// Constructor/destructor
		NetworkSequenceSpace(int myID, int numDataNodes, int kmerSize);
		~NetworkSequenceSpace();
		
		// add a sequence to the collection
		void addSequence(const PackedSeq& seq);
		
		// remove a sequence from the collection
		void removeSequence(const PackedSeq& seq);
		
		// end the data load and make the sequence space ready for data read
		void finalize();
		
		// generate the adjacency info for every sequence in the collection
		void generateAdjacency();
		
		// check if a sequence exists
		bool checkForSequence(const PackedSeq& seq) const;
		
		// Set flag for sequence seq
		void markSequence(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set
		bool checkSequenceFlag(const PackedSeq& seq, SeqFlag flag);

		// calculate whether this sequence has an extension in the phase space
		HitRecord calculateExtension(const PackedSeq& currSeq, extDirection dir) const;
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq) const;

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq) const;
		
		// Return the number of sequences in the collection
		int countAll() const;
		
	private:
	
		// Check if this sequence belongs in the local phase space
		bool isLocal(const PackedSeq& seq) const;
		
		// Get the node id that this sequence should reside in
		int computeNodeID(const PackedSeq& seq) const;
		
		// Pointer to the local sequence space
		// These sequences are held in memory on this process
		SimpleSequenceSpace* m_pLocalSpace;
		
		// The communications layer implements the functions over the network
		CommLayer* m_pComm;
		
		// The ID for this node, assigned by MPI
		int m_id;
		
		int m_numDataNodes;
};

#endif
