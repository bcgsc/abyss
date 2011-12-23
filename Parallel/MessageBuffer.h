#ifndef MESSAGE_BUFFER_H
#define MESSAGE_BUFFER_H 1

class MessageBuffer;

#include "CommLayer.h"
#include "Messages.h"
#include <vector>

typedef std::vector<Message*> MsgBuffer;
typedef std::vector<MsgBuffer> MessageQueues;

enum SendMode
{
	SM_BUFFERED,
	SM_IMMEDIATE
};

/** A buffer of Message. */
class MessageBuffer : public CommLayer
{
	public:
		MessageBuffer();

		void sendCheckPointMessage(int argument = 0)
		{
			assert(transmitBufferEmpty());
			CommLayer::sendCheckPointMessage(argument);
		}

		void sendControlMessage(APControl command, int argument = 0)
		{
			assert(transmitBufferEmpty());
			CommLayer::sendControlMessage(command, argument);
		}

		void sendControlMessageToNode(int dest,
				APControl command, int argument = 0)
		{
			assert(transmitBufferEmpty());
			CommLayer::sendControlMessageToNode(dest,
					command, argument);
		}

		void sendSeqAddMessage(int nodeID, const Kmer& seq);
		void sendSeqRemoveMessage(int nodeID, const Kmer& seq);
		void sendSetFlagMessage(int nodeID,
				const Kmer& seq, SeqFlag flag);
		void sendRemoveExtension(int nodeID,
				const Kmer& seq, extDirection dir, SeqExt ext);
		void sendSeqDataRequest(int nodeID,
				IDType group, IDType id, const Kmer& seq);
		void sendSeqDataResponse(int nodeID,
				IDType group, IDType id, const Kmer& seq,
				ExtensionRecord extRec, int multiplicity);
		void sendSetBaseExtension(int nodeID,
				const Kmer& seq, extDirection dir, uint8_t base);

		void flush();
		void queueMessage
			(int nodeID, Message* message, SendMode mode);

		// clear out a queue
		void clearQueue(int nodeID);
		bool transmitBufferEmpty() const;

		// check if a queue is full, if so, send the messages if the
		// immediate mode flag is set, send even if the queue is not
		// full.
		void checkQueueForSend(int nodeID, SendMode mode);
		
	private:
		static const size_t MAX_MESSAGES = 100;
		MessageQueues m_msgQueues;
};

#endif
