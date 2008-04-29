#include "CommLayer.h"

//
//
//
CommLayer::CommLayer(int id, int kmerSize) : m_id(id), m_kmerSize(kmerSize), m_msgID(0)
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
	printf("%d: Comm layer performed %lu buffered sends\n", m_id, (unsigned long)m_msgID);
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


uint64_t CommLayer::SendCheckPointMessage(int argument)
{
	ControlMessage msg;
	msg.id = m_msgID++;
	msg.msgType = APC_CHECKPOINT;
	msg.argument = argument;

	MPI::COMM_WORLD.Ssend(&msg, sizeof(ControlMessage), MPI::BYTE, CONTROL_ID, APM_CONTROL);
	return msg.id;
}

//
// Send a sequence to a specific id
//
uint64_t CommLayer::SendSeqMessage(int destID, const PackedSeq& seq, APSeqOperation operation)
{
	SeqMessage msg;
	msg.id = m_msgID++;
	msg.seq = seq;
	msg.operation = operation;	
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SeqMessage), MPI::BYTE, destID, APM_SEQ);
	req.Free();
	return msg.id;	
}

//
// Send a sequence extension message
//
uint64_t CommLayer::SendSeqExtMessage(int destID, const PackedSeq& seq, APSeqExtOperation operation, extDirection dir, SeqExt ext, char base)
{
	SeqExtMessage msg;
	msg.id = m_msgID++;
	msg.seq = seq;
	msg.operation = operation;
	msg.ext = ext;
	msg.dir = dir;
	msg.base = base;

	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SeqExtMessage), MPI::BYTE, destID, APM_SEQ_EXT);
	req.Free();
	return msg.id;	
}


//
// Send a sequence flag message
//
uint64_t CommLayer::SendSeqFlagMessage(int destID, const PackedSeq& seq, APSeqFlagOperation operation, SeqFlag flag)
{
	SeqFlagMessage msg;
	msg.id = m_msgID++;
	msg.seq = seq;
	msg.operation = operation;
	msg.flag = flag;

	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SeqFlagMessage), MPI::BYTE, destID, APM_SEQ_FLAG);
	req.Free();
	return msg.id;	
}

//
// Send a control message
//
void CommLayer::SendControlMessage(int numNodes, APControl m, int argument)
{
	// i starts at 1 because the control node does not get the message
	for(int i = 1; i < numNodes; i++)
	{
		SendControlMessageToNode(i, m, argument);
	}
}

//
// Send a control message to a specific node
//
uint64_t CommLayer::SendControlMessageToNode(int nodeID, APControl m, int argument)
{
	ControlMessage msg;
	msg.id = m_msgID++;
	msg.msgType = m;
	msg.argument = argument;
	
	// Control messages are synchronous
	MPI::COMM_WORLD.Ssend(&msg, sizeof(ControlMessage), MPI::BYTE, nodeID, APM_CONTROL);
	return msg.id;	
}

//
// Send an adjacency request
//
uint64_t CommLayer::SendAdjacencyRequest(int destID, const PackedSeq& testSequence, const PackedSeq& originalSeq, extDirection dir, char b)
{
	AdjacencyMessage msg;
	msg.id = m_msgID++;
	msg.testSeq = testSequence;
	msg.originalSeq = originalSeq;
	msg.dir = dir;
	msg.base = b;
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(AdjacencyMessage), MPI::BYTE, destID, APM_ADJACENCY);
	req.Free();
	return msg.id;	
	
}

//
// Send an adjacency result message
//
void CommLayer::SendAdjacencyResult(int destID, uint64_t reqID, PackedSeq originalSeq, extDirection dir, char base, bool b)
{
	ResultPair rp;
	rp.forward = b;
	rp.reverse = b;	

	AdjacencyResultMessage msg;
	msg.originalSeq = originalSeq;
	msg.dir = dir;
	msg.base = base;
	msg.id = reqID;
	
	// Convert the result pair to a result message
	msg.result[0] = (rp.forward) ? APR_TRUE : APR_FALSE;
	msg.result[1] = (rp.reverse) ? APR_TRUE : APR_FALSE;
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(AdjacencyResultMessage), MPI::BYTE, destID, APM_RESULT_ADJ);
	req.Free();		
}

//
// Send a request for the extension of a sequence
//
uint64_t CommLayer::SendExtensionRequest(int destID, uint64_t groupID, uint64_t branchID, const PackedSeq& seq)
{
	SequenceExtensionRequestMessage msg;
	msg.id = m_msgID++;
	msg.groupID = groupID;
	msg.branchID = branchID;
	msg.seq = seq;

	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SequenceExtensionRequestMessage), MPI::BYTE, destID, APM_SEQ_EXT_REQUEST);
	req.Free();
	return msg.id;	
}

//
// Send the response for the extension of a sequence
//
void CommLayer::SendExtensionResponse(int destID, uint64_t reqID, uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec)
{
	SequenceExtensionResponseMessage msg;
	msg.id = reqID;
	msg.groupID = groupID;
	msg.branchID = branchID;	
	msg.seq = seq;
	msg.extRec = extRec;

	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(SequenceExtensionResponseMessage), MPI::BYTE, destID, APM_SEQ_EXT_RESPONSE);
	req.Free();
}

//
// Send a bool result
//
void CommLayer::SendResultMessage(int destID, uint64_t reqID, bool b)
{
	ResultPair rp;
	rp.forward = b;
	rp.reverse = b;
	SendResultMessage(destID, reqID, rp);
}

//
// Send a result
//
void CommLayer::SendResultMessage(int destID, uint64_t reqID, ResultPair rp)
{
	ResultPairMessage msg;
	msg.id = reqID;
		
	// Convert the result pair to a result message
	msg.result[0] = (rp.forward) ? APR_TRUE : APR_FALSE;
	msg.result[1] = (rp.reverse) ? APR_TRUE : APR_FALSE;
	MPI::Request req = MPI::COMM_WORLD.Ibsend(&msg, sizeof(ResultPairMessage), MPI::BYTE, destID, APM_RESULT);
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
// Receive a sequence extension message
//
SeqExtMessage CommLayer::ReceiveSeqExtMessage()
{
	SeqExtMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(SeqExtMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
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
// Receive a result message
//
ResultPairMessage CommLayer::ReceiveResultMessage()
{
	ResultPairMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(ResultPairMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
}

//
// Receive an adjacency message
//
AdjacencyMessage CommLayer::ReceiveAdjacencyMessage()
{
	AdjacencyMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(AdjacencyMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
}

//
// Receive an adjacency result message
//
AdjacencyResultMessage CommLayer::ReceiveAdjacencyResultMessage()
{
	AdjacencyResultMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(AdjacencyResultMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
}

//
// Receive a sequence extension message
//
SequenceExtensionRequestMessage CommLayer::ReceiveSequenceExtensionMessage()
{
	SequenceExtensionRequestMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(SequenceExtensionRequestMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
}
		
// Receive a sequence extension result message
SequenceExtensionResponseMessage CommLayer::ReceiveSequenceExtensionResponseMessage()
{
	SequenceExtensionResponseMessage msg;
	MPI::COMM_WORLD.Recv(&msg, sizeof(SequenceExtensionResponseMessage), MPI::BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG);
	return msg;
}
