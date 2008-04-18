#include "CommLayer.h"

//
//
//
CommLayer::CommLayer(int id, int kmerSize) : m_id(id), m_kmerSize(kmerSize)
{

	m_bufferSize = 200*1024*1024;
	// Create the  buffer
	
	m_buffer = new char[m_bufferSize];
	MPI_Buffer_attach( m_buffer, m_bufferSize); 
}

//
//
//
CommLayer::~CommLayer()
{
	
	
	// Destroy the buffer
	delete [] m_buffer;
}

void CommLayer::flush()
{

}

//
//
//
APMessage CommLayer::CheckMessage(int& sendID) const
{ 
	MPI::Status status;
	int flag = MPI::COMM_WORLD.Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, status);
	
	// check if a message waits
	if(flag)
	{
		//printf("message found of type: %d (flag: %d)\n", status.Get_tag(), flag);
		sendID = status.Get_source();
		return (APMessage)status.Get_tag();
	}
	else
	{
		// no messsage
		return APM_NONE;
	}
}


void CommLayer::SendCheckPointMessage(int argument) const
{
	ControlMessage msg;
	msg.msgType = APC_CHECKPOINT;
	msg.argument = argument;
	MPI::COMM_WORLD.Ssend(&msg, sizeof(ControlMessage), MPI::BYTE, CONTROL_ID, APM_CONTROL);
}

//
// Send a sequence to a specific id
//
void CommLayer::SendSeqMessage(int destID, const PackedSeq& seq, APSeqOperation operation) const
{
	SeqMessage msg;
	
	msg.seq = seq;
	msg.operation = operation;	
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SeqMessage), MPI::BYTE, destID, APM_SEQ);
	//PrintBufferAsHex((char*)&msg, sizeof(SeqMessage));
	
	req.Free();

}

//
// Receive a seq message
//
SeqMessage CommLayer::ReceiveSeqMessage()
{
	SeqMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(SeqMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	//PrintBufferAsHex((char*)&msg, sizeof(SeqMessage));	
	return msg;	
}

//
// Send a sequence extension message
//
void CommLayer::SendSeqExtMessage(int destID, const PackedSeq& seq, APSeqExtOperation operation, extDirection dir, SeqExt ext, char base) const
{
	SeqExtMessage msg;
	
	msg.seq = seq;
	msg.operation = operation;
	msg.ext = ext;
	msg.dir = dir;
	msg.base = base;

	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SeqExtMessage), MPI::BYTE, destID, APM_SEQ_EXT);
	req.Free();
	
}

//
// Receive a sequence extension message
//
SeqExtMessage CommLayer::ReceiveSeqExtMessage()
{
	SeqExtMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(SeqExtMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;	
}

//
// Send a sequence flag message
//
void CommLayer::SendSeqFlagMessage(int destID, const PackedSeq& seq, APSeqFlagOperation operation, SeqFlag flag) const
{
	SeqFlagMessage msg;
	
	msg.seq = seq;
	msg.operation = operation;
	msg.flag = flag;

	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SeqFlagMessage), MPI::BYTE, destID, APM_SEQ_FLAG);
	req.Free();
}

//
// Receive a seq message
//
SeqFlagMessage CommLayer::ReceiveSeqFlagMessage()
{
	SeqFlagMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(SeqFlagMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;	
}

//
// Send a control message
//
void CommLayer::SendControlMessage(int numNodes, APControl m, int argument) const
{
	ControlMessage msg;
	msg.msgType = m;
	msg.argument = argument;
	
	for(int i = 1; i < numNodes; i++)
	{
		// Control messages are synchronous
		MPI::COMM_WORLD.Ssend(&msg, sizeof(ControlMessage), MPI::BYTE, i, APM_CONTROL);
	}
}

//
// Receive a control message
//
ControlMessage CommLayer::ReceiveControlMessage()
{
	ControlMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(ControlMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;	
}

//
// Send a result
//
void CommLayer::SendResultMessage(int destID, ResultPair rp)
{
	ResultMessage msg;
	
	// Convert the result pair to a result message
	msg.result[0] = (rp.forward) ? APR_TRUE : APR_FALSE;
	msg.result[1] = (rp.reverse) ? APR_TRUE : APR_FALSE;
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(ResultMessage), MPI::BYTE, destID, APM_RESULT);
	req.Free();
}

//
// Send a bool result
//
void CommLayer::SendResultMessage(int destID, bool b)
{
	ResultMessage msg;
	
	// Convert the result pair to a result message
	msg.result[0] = (b) ? APR_TRUE : APR_FALSE;
	msg.result[1] = (b) ? APR_TRUE : APR_FALSE;
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(ResultMessage), MPI::BYTE, destID, APM_RESULT);
	req.Free();
}

//
// Receive a result message
//
ResultMessage CommLayer::ReceiveResultMessage()
{
	ResultMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(ResultMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
}

