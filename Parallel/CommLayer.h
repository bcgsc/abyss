#ifndef COMMLAYER_H
#define COMMLAYER_H 1

#include "NetworkDefs.h"
#include "Messages.h"
#include <mpi.h>
#include <vector>

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
		APMessage checkMessage(int &sendID);

		// Return whether a message has been received.
		bool receiveEmpty();

		// Block until all processes have reached this routine.
		void barrier();

		void broadcast(int message);
		int receiveBroadcast();

		// Block until all processes have reached this routine.
		long unsigned reduce(long unsigned count);
		std::vector<unsigned> reduce(const std::vector<unsigned>& v);

		// Send a control message
		void sendControlMessage(APControl m, int argument = 0);

		// Send a control message to a specific node
		uint64_t sendControlMessageToNode(int nodeID,
				APControl command, int argument = 0);

		// Receive a control message
		ControlMessage receiveControlMessage();

		// Send a message that the checkpoint has been reached
		uint64_t sendCheckPointMessage(int argument = 0);

		// Send a buffered message
		void sendBufferedMessage(int destID, char* msg, size_t size);

		// Receive a buffered sequence of messages
		void receiveBufferedMessage(MessagePtrVector& outmessages);

		uint64_t reduceInflight()
		{
			return reduce(m_txPackets - m_rxPackets);
		}

	private:
		uint64_t m_msgID;
		uint8_t* m_rxBuffer;
		MPI_Request m_request;

	protected:
		// Counters
		uint64_t m_rxPackets;
		uint64_t m_rxMessages;
		uint64_t m_rxBytes;
		uint64_t m_txPackets;
		uint64_t m_txMessages;
		uint64_t m_txBytes;
};

#endif
