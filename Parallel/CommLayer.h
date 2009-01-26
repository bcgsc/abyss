#ifndef COMMLAYER
#define COMMLAYER

class CommLayer;

#include <list>
#include "PackedSeq.h"
#include "NetworkDefs.h"
#include "Messages.h"
#include "MessageBuffer.h"

typedef std::vector<Message*> MessagePtrVector;

struct ControlMessage
{
	int64_t id;
	APControl msgType;
	int argument;
};

const int CONTROL_ID = 0;

// The comm layer wraps inter-process communication operations
class CommLayer
{
	public:
	
		// Constructor/Destructor
		CommLayer(int id, int kmerSize);
		~CommLayer();
	
		// Check if a message exists, if it does return the type
		APMessage CheckMessage(int &sendID) const;

		// Return whether the queue of messages is empty.
		bool empty() const;
		
		// Block until all processes have reached this routine.
		void barrier();
		
		// Send a control message
		void SendControlMessage(int numNodes, APControl m, int argument = 0);
		
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
		
		// Flush the buffer
		void flush();
		
		void setMsgBuffer(const MessageBuffer *pMsgBuffer)
		{
			assert(pMsgBuffer != NULL);
			m_pMsgBuffer = pMsgBuffer;
		}

	private:
		int m_id;
		int m_kmerSize;
		int m_numBytesPerSeq;
		int m_bufferSize;
		uint64_t m_msgID;
		char* m_buffer;
		const MessageBuffer *m_pMsgBuffer;
};

#endif
