#include "config.h"
#include "CommLayer.h"
#include "Log.h"
#include <mpi.h>
#include <stdint.h>

static const unsigned TX_BUFSIZE = 200*1024*1024;
static const unsigned RX_BUFSIZE = 16*1024;

//
//
//
CommLayer::CommLayer(int id)
	: m_id(id), m_msgID(0),
	  m_txBuffer(new uint8_t[TX_BUFSIZE]),
	  m_rxBuffer(new uint8_t[RX_BUFSIZE]),
	  m_request(MPI_REQUEST_NULL),
	  m_pMsgBuffer(NULL)
{
	MPI_Buffer_attach(m_txBuffer, TX_BUFSIZE);
	assert(m_request == MPI_REQUEST_NULL);
	MPI_Irecv(m_rxBuffer, RX_BUFSIZE,
			MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			&m_request);
}

//
//
//
CommLayer::~CommLayer()
{
	printf("%d: Sent %llu messages\n", m_id,
			(unsigned long long)m_msgID);
	delete[] m_txBuffer;
	delete[] m_rxBuffer;
}

/** Return the tag of the received message or APM_NONE if no message
 * has been received. If a message has been received, this call should
 * be followed by a call to either ReceiveControlMessage or
 * ReceiveBufferedMessage.
 */
APMessage CommLayer::CheckMessage(int& sendID)
{
	int flag;
	MPI_Status status;
	MPI_Request_get_status(m_request, &flag, &status);
	if (flag)
		sendID = status.MPI_SOURCE;
	return flag ? (APMessage)status.MPI_TAG : APM_NONE;
}

/** Return true if no message has been received. */
bool CommLayer::empty()
{
	int sendID;
	return CheckMessage(sendID) == APM_NONE;
}

/** Block until all processes have reached this routine. */
void CommLayer::barrier()
{
	PrintDebug(4, "entering barrier\n");
	MPI_Barrier(MPI_COMM_WORLD);
	PrintDebug(4, "left barrier\n");
}

/** Block until all processes have reached this routine.
 * @return the sum of count from all processors
 */
unsigned CommLayer::reduce(unsigned count)
{
	PrintDebug(4, "entering reduce: %u\n", count);
	unsigned sum;
	MPI_Allreduce(&count, &sum, 1, MPI_UNSIGNED, MPI_SUM,
			MPI_COMM_WORLD);
	PrintDebug(4, "left reduce: %u\n", sum);
	return sum;
}

uint64_t CommLayer::SendCheckPointMessage(int argument)
{
	assert(m_pMsgBuffer->empty());
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
	assert(m_pMsgBuffer->empty());
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
	assert(m_pMsgBuffer->empty());
	ControlMessage msg;
	msg.id = m_msgID++;
	msg.msgType = m;
	msg.argument = argument;
	
	// Control messages are synchronous
	MPI_Ssend(&msg, sizeof(ControlMessage), MPI_BYTE, nodeID, APM_CONTROL, MPI_COMM_WORLD);
	return msg.id;	
}

/** Receive a control message. */
ControlMessage CommLayer::ReceiveControlMessage()
{
	int flag;
	MPI_Status status;
	MPI_Test(&m_request, &flag, &status);
	assert(flag);
	assert((APMessage)status.MPI_TAG == APM_CONTROL);

	int count;
	MPI_Get_count(&status, MPI_BYTE, &count);
	ControlMessage msg;
	assert(count == sizeof msg);
	memcpy(&msg, m_rxBuffer, sizeof msg);
	assert(m_request == MPI_REQUEST_NULL);
	MPI_Irecv(m_rxBuffer, RX_BUFSIZE,
			MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			&m_request);
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

/** Receive a buffered message. */
void CommLayer::ReceiveBufferedMessage(MessagePtrVector& outmessages)
{
	int flag;
	MPI_Status status;
	MPI_Test(&m_request, &flag, &status);
	assert(flag);
	assert((APMessage)status.MPI_TAG == APM_BUFFERED);

	int size;
	MPI_Get_count(&status, MPI_BYTE, &size);

	int offset = 0;
	while (offset < size) {
		MessageType type = Message::readMessageType(
				(char*)m_rxBuffer + offset);

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
		offset += pNewMessage->unserialize(
				(char*)m_rxBuffer + offset);
		//pNewMessage->print();
		
		// Constructed message will be deleted in the NetworkSequenceCollection calling function
		outmessages.push_back(pNewMessage);
	}
	assert(offset == size);

	assert(m_request == MPI_REQUEST_NULL);
	MPI_Irecv(m_rxBuffer, RX_BUFSIZE,
			MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			&m_request);
}
