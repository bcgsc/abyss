#ifndef NETWORKSEQUENCECOLLECTION_H
#define NETWORKSEQUENCECOLLECTION_H

#include "PackedSeq.h"
#include "SequenceCollectionHash.h"
#include "CommLayer.h"
#include "SeqRecord.h"
#include "AssemblyAlgorithms.h"

enum NetworkAssemblyState
{
	NAS_LOADING, // loading sequences
	NAS_FINALIZE, // finalizing the sequence data and getting ready for processing
	NAS_GEN_ADJ, // generating the sequence data
	NAS_TRIM, // trimming the data
	NAS_POPBUBBLE, // remove read errors/SNPs
	NAS_SPLIT, // remove ambiguous links (just before assembling into contigs)
	NAS_ASSEMBLE, // assembling the data
	NAS_WAITING, // non-control process is waiting, this just loops over the network function
	NAS_DONE // finished, clean up and exit
};

class NetworkSequenceCollection : public ISequenceCollection
{
	public:
	
		// Constructor/destructor
		NetworkSequenceCollection(int myID, int numDataNodes, int kmerSize, int readLen);
		~NetworkSequenceCollection();
		
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
		
		// remove all the extensions of this sequence
		void clearExtensions(const PackedSeq& seq, extDirection dir);		
		
		// get the multiplicity of the sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// The loop to run the network code
		APResult pumpNetwork();
		
		// Loop over the pumping function while waiting for a result from the network
		ResultPair pumpUntilResult(); 
		
		// run the assembly
		void run(int readLength, int kmerSize);
		
		// run the assembly from the controller's point of view (rank 0 node)
		void runControl(std::string fastaFile, int readLength, int kmerSize);
		
		// test if the checkpoint has been reached
		bool checkpointReached() const;
		
		// get an iterator to the first sequence
		SequenceCollectionIterator getStartIter() const;
		
		// get an iterator to the last sequence
		SequenceCollectionIterator getEndIter() const;		
		
	private:
	
	
		// Network message parsers
		void parseControlMessage(int senderID);
		void parseSeqMessage(int senderID);
		void parseSeqFlagMessage(int senderID);
		void parseSeqExtMessage(int senderID);	
		
		// Read a fasta file and distribute the sequences
		void readSequences(std::string fastaFile, int readLength, int kmerSize);
	
		// Check if this sequence belongs in the local phase space
		bool isLocal(const PackedSeq& seq) const;
		
		// Get the node id that this sequence should reside in
		int computeNodeID(const PackedSeq& seq) const;
		
		// Print function
		int PrintDebug(int level,char* fmt, ...) const;
		
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
		int m_numDataNodes;
		
		// the state of the assembly
		NetworkAssemblyState m_state;
		
		// The number of processes that have sent a checkpoint reached message, this is used by the control process to determine the state flow
		int m_numReachedCheckpoint;
		
		// the size of k used in the assembly
		int m_kmer;
		
		// the original read length
		int m_readLen;
		
		// the current length to trim on (comes from the control node)
		int m_trimStep;
		
		// the number of sequences assembled so far
		int m_numAssembled;
};

#endif
