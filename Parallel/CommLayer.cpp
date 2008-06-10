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
	printf("%d: Sent %llu messages\n", m_id,
			(unsigned long long)m_msgID);
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
	MPI_Status status;
	int flag;
	MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
	
	// check if a message waits
	if(flag)
	{
		//printf("message found of type: %d (flag: %d)\n", status.Get_tag(), flag);
		sendID = status.MPI_SOURCE;
		return (APMessage)status.MPI_TAG;
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

	MPI_Ssend(&msg, sizeof(ControlMessage), MPI_BYTE, CONTROL_ID, APM_CONTROL, MPI_COMM_WORLD);
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
	MPI_Ssend(&msg, sizeof(ControlMessage), MPI_BYTE, nodeID, APM_CONTROL, MPI_COMM_WORLD);
	return msg.id;	
}

//
// Receive a control message
//
ControlMessage CommLayer::ReceiveControlMessage()
{
	ControlMessage msg;
	MPI_Status status;	
	MPI_Recv(&msg, sizeof(ControlMessage), MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	return msg;	
}

//
// Send a buffered collection of messages
//
void CommLayer::SendBufferedMessage(int destID, char* msg, size_t size)
{
	MPI_Request req;
			
	MPI_Ibsend(msg, size, MPI_BYTE, destID, APM_BUFFERED, MPI_COMM_WORLD, &req);
	MPI_Request_free(&req);	
	
	//printf("buffered send: %zu bytes\n", size);
}

// Receive a buffered message
void CommLayer::ReceiveBufferedMessage(MessagePtrVector& outmessages)
{
	int size;
	MPI_Status status;
	
	MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Get_count(&status, MPI_BYTE, &size);

	//printf("buffered recv: %zu bytes\n", size);	
	
	// Allocate a buffer to hold the messages
	char* buffer = new char[size];
	
	// Copy the messages (blocking)
	MPI_Recv(buffer, size, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	//printf("Received %d bytes\n", size);	
	//printf("Receive buffer: \n");
	//PrintBufferAsHex(buffer, size);	
	// Interpret the sequences
	int offset = 0;
	
	//printf("Received Messages: \n");
	while(offset < size)
	{
		// Read the type of the current message
		MessageType type = Message::readMessageType(buffer + offset);
		//printf("pos: %d type: %d\n", offset, (int)type);
		Message* pNewMessage;
		
		switch(type)
		{
			case MT_SEQ_OP:
			{
				pNewMessage = new SeqOpMessage();
				break;
			}
			case MT_SET_FLAG:
			{
				pNewMessage = new SetFlagMessage();
				break;
			}	
			case MT_REMOVE_EXT:
			{
				pNewMessage = new RemoveExtensionMessage();
				break;
			}					
			case MT_SEQ_DATA_REQUEST:
			{
				pNewMessage = new SeqDataRequest();
				break;
			}	
			case MT_SEQ_DATA_RESPONSE:
			{
				pNewMessage = new SeqDataResponse();
				break;
			}	
			case MT_SET_BASE:
			{
				pNewMessage = new SetBaseMessage();
				break;
			}									
			default:
			{
				assert(false);
				break;
			}
		}
		
		// Unserialize the new message from the buffer
		offset += pNewMessage->unserialize(buffer + offset);
		//pNewMessage->print();
		
		// Constructed message will be deleted in the NetworkSequenceCollection calling function
		outmessages.push_back(pNewMessage);
	}
	
	assert(offset == size);
	
	delete [] buffer;
	
}


