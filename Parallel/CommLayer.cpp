#include "CommLayer.h"

//
//
//
CommLayer::CommLayer(int id, int kmerSize) : m_id(id), m_kmerSize(kmerSize)
{
	// Cache the number of bytes per sequence
	m_numBytesPerSeq = PackedSeq::getNumCodingBytes(m_kmerSize);
	
	// Create the receive buffer
	m_pRecvBuffer = new char[m_numBytesPerSeq];
}

//
//
//
CommLayer::~CommLayer()
{
	// Destroy the buffer
	delete [] m_pRecvBuffer;
	m_pRecvBuffer = 0;
}

//
//
//
APMessage CommLayer::CheckMessage(int& sendID) const
{ 
	MPI::Status status;
	MPI::COMM_WORLD.Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, status);
	
	//printf("message found of type: %d\n", status.Get_tag());
	sendID = status.Get_source();
	return (APMessage)status.Get_tag();
}

//
//
//
void CommLayer::SendControlMessage(int destID, APMessage msg)
{
	// Is there a simpler way to do this?
	char dummy = 'A';
	MPI::COMM_WORLD.Send(&dummy, 1, MPI::CHAR, destID, msg);
}

//
//
//
void CommLayer::SendSequence(int destID, const PackedSeq& seq, APMessage msg)
{
	MPI::COMM_WORLD.Send(seq.getDataPtr(), m_numBytesPerSeq, MPI::BYTE, destID, msg);
	//printf("%d: sent %s\n", m_id, seq.decode().c_str());
}

//
//
//
PackedSeq CommLayer::ReceiveSequence()
{
	MPI::COMM_WORLD.Recv(m_pRecvBuffer, m_numBytesPerSeq, MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	
	// Return the received sequence
	PackedSeq seq(m_pRecvBuffer, m_kmerSize);
	//printf("%d: received %s\n", m_id, seq.decode().c_str());
	return seq;
}

//
//
//
bool CommLayer::CheckForSequence(int destID, PackedSeq& seq)
{
	SendSequence(destID, seq, APM_SEQCHECK);
	return ReceiveBool(destID);
}

//
//
//
bool CommLayer::ReceiveBool(int sendID)
{
	char b;
	MPI::COMM_WORLD.Recv(&b, 1, MPI::BYTE, sendID, APM_DATA);
	return (bool)b;
}

//
//
//
void CommLayer::SendBool(int destID, bool b)
{
	char m = (char)b;
	MPI::COMM_WORLD.Send(&m, 1, MPI::BYTE, destID, APM_DATA);
}


//
//
//
void CommLayer::ClearControlMessage()
{
	char c;
	MPI::COMM_WORLD.Recv(&c, 1, MPI::CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG);	
}
