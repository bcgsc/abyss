#ifndef COMMLAYER_H
#define COMMLAYER_H 1

#include "NetworkDefs.h"
#include "Messages.h"
#include <mpi.h>

typedef std::vector<Message*> MessagePtrVector;

struct ControlMessage
{
	int64_t id;
	APControl msgType;
	int argument;
};

// The comm layer wraps inter-process communication operations
class CommLayer
{
	public:
		CommLayer();
		~CommLayer();

		// Check if a message exists, if it does return the type
		APMessage CheckMessage(int &sendID);

		// Return whether a message has been received.
		bool receiveEmpty();

		// Block until all processes have reached this routine.
		void barrier();

		// Block until all processes have reached this routine.
		unsigned reduce(unsigned count);

		// Send a control message
		void sendControlMessage(APControl m, int argument = 0);

		// Send a control message to a specific node
		uint64_t SendControlMessageToNode(int nodeID, APControl m, int argument = 0);

		// Receive a control message
		ControlMessage ReceiveControlMessage();
		
		// Send a message that the checkpoint has been reached
		uint64_t SendCheckPointMessage(int argument = 0);

		// Send a buffered message
		void SendBufferedMessage(int destID, char* msg, size_t size);

		// Receive a buffered sequence of messages
		void ReceiveBufferedMessage(MessagePtrVector& outmessages);

	private:
		uint64_t m_msgID;
		uint8_t* m_rxBuffer;
		MPI_Request m_request;
};

#endif
