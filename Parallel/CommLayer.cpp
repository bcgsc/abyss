#include "config.h"
#include "CommLayer.h"
#include "Common/Options.h"
#include "Log.h"
#include "MessageBuffer.h"
#include <mpi.h>
#include <cstring>
#include <vector>

using namespace std;

static const unsigned RX_BUFSIZE = 16*1024;

CommLayer::CommLayer()
	: m_msgID(0),
	  m_rxBuffer(new uint8_t[RX_BUFSIZE]),
	  m_request(MPI_REQUEST_NULL),
	  m_rxPackets(0), m_rxMessages(0), m_rxBytes(0),
	  m_txPackets(0), m_txMessages(0), m_txBytes(0)
{
	assert(m_request == MPI_REQUEST_NULL);
	MPI_Irecv(m_rxBuffer, RX_BUFSIZE,
			MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			&m_request);
}

CommLayer::~CommLayer()
{
	delete[] m_rxBuffer;
	logger(1) << "Sent " << m_msgID << " control, "
		<< m_txPackets << " packets, "
		<< m_txMessages << " messages, "
		<< m_txBytes << " bytes. "
		<< "Received " << m_rxPackets << " packets, "
		<< m_rxMessages << " messages, "
		<< m_rxBytes << " bytes.\n";
}

/** Return the status of an MPI request.
 * Wraps MPI_Request_get_status.
 */
static bool request_get_status(const MPI_Request& req,
		MPI_Status& status)
{
	int flag;
	MPI_Request_get_status(req, &flag, &status);
	// Work around a bug present in Open MPI 1.3.3 and earlier.
	// MPI_Request_get_status may return false on the first call even
	// though a message is waiting. The second call should work.
	if (!flag)
		MPI_Request_get_status(req, &flag, &status);
	return flag;
}

/** Return the tag of the received message or APM_NONE if no message
 * has been received. If a message has been received, this call should
 * be followed by a call to either ReceiveControlMessage or
 * ReceiveBufferedMessage.
 */
APMessage CommLayer::checkMessage(int& sendID)
{
	MPI_Status status;
	bool flag = request_get_status(m_request, status);
	if (flag)
		sendID = status.MPI_SOURCE;
	return flag ? (APMessage)status.MPI_TAG : APM_NONE;
}

/** Return true if no message has been received. */
bool CommLayer::receiveEmpty()
{
	MPI_Status status;
	return !request_get_status(m_request, status);
}

/** Block until all processes have reached this routine. */
void CommLayer::barrier()
{
	logger(4) << "entering barrier\n";
	MPI_Barrier(MPI_COMM_WORLD);
	logger(4) << "left barrier\n";
}

/** Broadcast a message. */
void CommLayer::broadcast(int message)
{
	assert(opt::rank == 0);
	MPI_Bcast(&message, 1, MPI_INT, 0, MPI_COMM_WORLD);
	barrier();
}

/** Receive a broadcast message. */
int CommLayer::receiveBroadcast()
{
	assert(opt::rank != 0);
	int message;
	MPI_Bcast(&message, 1, MPI_INT, 0, MPI_COMM_WORLD);
	barrier();
	return message;
}

/** Block until all processes have reached this routine.
 * @return the sum of count from all processors
 */
long long unsigned CommLayer::reduce(long long unsigned count)
{
	logger(4) << "entering reduce: " << count << '\n';
	long long unsigned sum;
	MPI_Allreduce(&count, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
			MPI_COMM_WORLD);
	logger(4) << "left reduce: " << sum << '\n';
	return sum;
}

/** Reduce the specified vector. */
vector<unsigned> CommLayer::reduce(const vector<unsigned>& v)
{
	logger(4) << "entering reduce\n";
	vector<unsigned> sum(v.size());
	MPI_Allreduce(const_cast<unsigned*>(&v[0]),
			&sum[0], v.size(), MPI_UNSIGNED, MPI_SUM,
			MPI_COMM_WORLD);
	logger(4) << "left reduce\n";
	return sum;
}

/** Reduce the specified vector. */
vector<long unsigned> CommLayer::reduce(
		const vector<long unsigned>& v)
{
	logger(4) << "entering reduce\n";
	vector<long unsigned> sum(v.size());
	MPI_Allreduce(const_cast<long unsigned*>(&v[0]),
			&sum[0], v.size(), MPI_UNSIGNED_LONG, MPI_SUM,
			MPI_COMM_WORLD);
	logger(4) << "left reduce\n";
	return sum;
}

uint64_t CommLayer::sendCheckPointMessage(int argument)
{
	logger(4) << "checkpoint: " << argument << '\n';
	assert(opt::rank != 0);
	ControlMessage msg;
	msg.id = m_msgID++;
	msg.msgType = APC_CHECKPOINT;
	msg.argument = argument;

	MPI_Send(&msg, sizeof msg,
			MPI_BYTE, 0, APM_CONTROL, MPI_COMM_WORLD);
	return msg.id;
}

/** Send a control message to every other process. */
void CommLayer::sendControlMessage(APControl m, int argument)
{
	for (int i = 0; i < opt::numProc; i++)
		if (i != opt::rank) // Don't send the message to myself.
			sendControlMessageToNode(i, m, argument);
}

/** Send a control message to a specific node. */
uint64_t CommLayer::sendControlMessageToNode(int nodeID,
		APControl m, int argument)
{
	assert(opt::rank == 0);
	ControlMessage msg;
	msg.id = m_msgID++;
	msg.msgType = m;
	msg.argument = argument;

	MPI_Send(&msg, sizeof msg,
			MPI_BYTE, nodeID, APM_CONTROL, MPI_COMM_WORLD);
	return msg.id;
}

/** Receive a control message. */
ControlMessage CommLayer::receiveControlMessage()
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

/** Send a buffered collection of messages. */
void CommLayer::sendBufferedMessage(int destID,
		char* msg, size_t size)
{
	MPI_Send(msg, size, MPI_BYTE, destID, APM_BUFFERED,
			MPI_COMM_WORLD);
}

/** Receive a buffered message. */
void CommLayer::receiveBufferedMessage(MessagePtrVector& outmessages)
{
	assert(outmessages.empty());
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
			case MT_ADD:
				pNewMessage = new SeqAddMessage();
				break;
			case MT_REMOVE:
				pNewMessage = new SeqRemoveMessage();
				break;
			case MT_SET_FLAG:
				pNewMessage = new SetFlagMessage();
				break;
			case MT_REMOVE_EXT:
				pNewMessage = new RemoveExtensionMessage();
				break;
			case MT_SEQ_DATA_REQUEST:
				pNewMessage = new SeqDataRequest();
				break;
			case MT_SEQ_DATA_RESPONSE:
				pNewMessage = new SeqDataResponse();
				break;
			case MT_SET_BASE:
				pNewMessage = new SetBaseMessage();
				break;
			default:
				assert(false);
				break;
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

	m_rxPackets++;
	m_rxMessages += outmessages.size();
	m_rxBytes += size;
}
