#ifndef MESSAGE_BUFFER_H
#define MESSAGE_BUFFER_H 1

class MessageBuffer;

#include "CommLayer.h"
#include "Messages.h"
#include "NetworkDefs.h"
#include <vector>

typedef std::vector<Message*> MsgBuffer;
typedef std::vector<MsgBuffer> MessageQueues;

enum SendMode
{
	SM_BUFFERED,
	SM_IMMEDIATE
};

class MessageBuffer
{
	public:
		// Constructor, create a message buffer for every process
		MessageBuffer(CommLayer* pComm);

		void sendCheckPointMessage(int argument = 0)
		{
			assert(empty());
			m_pCommLayer->sendCheckPointMessage(argument);
		}

		void sendControlMessage(APControl command, int argument = 0)
		{
			assert(empty());
			m_pCommLayer->sendControlMessage(command, argument);
		}

		void sendControlMessageToNode(int dest,
				APControl command, int argument = 0)
		{
			assert(empty());
			m_pCommLayer->sendControlMessageToNode(dest,
					command, argument);
		}

		// send a sequence operation message
		void sendSeqOpMessage(int nodeID, const PackedSeq& seq, MessageOp op);
		
		// Send a set flag message
		void sendSetFlagMessage(int nodeID, const PackedSeq& seq, SeqFlag flag);

		// Send a remove extension message
		void sendRemoveExtension(int nodeID,
				const PackedSeq& seq, extDirection dir, uint8_t base);

		// Send a sequence data request
		void sendSeqDataRequest(int nodeID, IDType group, IDType id, const PackedSeq& seq);
		
		// Send a sequence data response
		void sendSeqDataResponse(int nodeID, IDType group, IDType id, const PackedSeq& seq, ExtensionRecord extRec, int multiplicity);

		// Send an adjacency message
		void sendSetBaseExtension(int nodeID,
				const PackedSeq& seq, extDirection dir, uint8_t base);

		// Send an existance response message
		void sendExistResponse(int nodeID, const PackedSeq& seq, bool result);
		
		// send all remaining messages
		void flush();
		
		// Queue a message
		void queueMessage(int nodeID, Message* message, SendMode mode);
		
		// clear out a queue
		void clearQueue(int nodeID);
		
		// check if the message buffer is empty (no messages pending)
		bool empty() const;
		
		// check if a queue is full, if so, send the messages
		// if the immediate mode flag is set, send even if the queue is not full
		void checkQueueForSend(int nodeID, SendMode mode);
		
	private:
		static const size_t MAX_MESSAGES = 100;
		MessageQueues m_msgQueues;
		CommLayer* m_pCommLayer;
};

#endif
