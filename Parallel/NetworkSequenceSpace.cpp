#include "NetworkSequenceSpace.h"

//
//
//
NetworkSequenceSpace::NetworkSequenceSpace(int myID, int numDataNodes, int kmerSize) : m_id(myID), m_numDataNodes(numDataNodes)
{
	// Load the phase space
	m_pLocalSpace = new SimpleSequenceSpace();
	
	// Create the comm layer
	m_pComm = new CommLayer(myID, kmerSize);
}

//
//
//
NetworkSequenceSpace::~NetworkSequenceSpace()
{
	// Delete the objects created in the constructor
	delete m_pLocalSpace;
	m_pLocalSpace = 0;
	
	delete m_pComm;
	m_pComm = 0;
}

//
//
//
void NetworkSequenceSpace::addSequence(const PackedSeq& seq)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		m_pLocalSpace->addSequence(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_pComm->SendSequence(nodeID, seq, APM_SEQADD);
	}
}

//
//
//
void NetworkSequenceSpace::removeSequence(const PackedSeq& seq)
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		m_pLocalSpace->removeSequence(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		m_pComm->SendSequence(nodeID, seq, APM_SEQDEL);
		
		// Perform network sequence extension removal
		assert(false);
	}	
}

//
//
//
void NetworkSequenceSpace::finalize()
{
	// this command is broadcast from the controller so we only perform a local finalize
	m_pLocalSpace->finalize();
}

//
//
//
void NetworkSequenceSpace::generateAdjacency()
{
	// generate the adjacency information for all sequences in the local sequence space
	for(SequenceCollectionIter iter = m_pLocalSpace->getStartIter(); iter != m_pLocalSpace->getEndIter(); iter++)
	{
		// do stuff
	}
}

//
//
//
bool NetworkSequenceSpace::checkForSequence(const PackedSeq& seq) const
{
	// Check if this sequence is local
	if(isLocal(seq))
	{
		return m_pLocalSpace->checkForSequence(seq);
	}
	else
	{
		int nodeID = computeNodeID(seq);
		// Queue request, wait for return?
		
	}
}

//
//
//
void NetworkSequenceSpace::markSequence(const PackedSeq& seq, SeqFlag flag)
{
	
}

//
//
//
bool NetworkSequenceSpace::checkSequenceFlag(const PackedSeq& seq, SeqFlag flag)
{
	
}

//
//
//
HitRecord NetworkSequenceSpace::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
{
	
}

//
//
//
bool NetworkSequenceSpace::hasParent(const PackedSeq& seq) const
{
	
}

//
//
//
bool NetworkSequenceSpace::hasChild(const PackedSeq& seq) const
{
	
}

//
//
//
int NetworkSequenceSpace::countAll() const
{
	
}

//
// Check if this sequence belongs in the local space
//
bool NetworkSequenceSpace::isLocal(const PackedSeq& seq) const
{
	int id = computeNodeID(seq);
	return id == m_id;
}

//
// 
//
int NetworkSequenceSpace::computeNodeID(const PackedSeq& seq) const
{
	int code = seq.getCode();
	int id = code % m_numDataNodes;
	return id;
}

