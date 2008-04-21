#ifndef COMMLAYER
#define COMMLAYER

#include <mpi.h>
#include <list>
#include "PackedSeq.h"
#include "NetworkDefs.h"


struct SeqMessage
{
	int64_t id;
	APSeqOperation operation;
	PackedSeq seq;
};

struct SeqFlagMessage
{
	int64_t id;
	APSeqFlagOperation operation;
	PackedSeq seq;
	SeqFlag flag;
};

struct SeqExtMessage
{
	int64_t id;
	APSeqExtOperation operation;
	PackedSeq seq;
	extDirection dir;	
	SeqExt ext;
	char base;
	
};

struct ControlMessage
{
	int64_t id;
	APControl msgType;
	int argument;
};

struct ResultMessage
{
	int64_t id;
	APResultType resultType;
	APResult result[2];
};

struct RequestBuffer
{
	char* buffer;
	MPI::Request request;
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
		
		// Send a control message
		void SendControlMessage(int numNodes, APControl m, int argument = 0);
		
		// Send a control message to a specific node
		uint64_t SendControlMessageToNode(int nodeID, APControl m, int argument = 0);
		
		// Send a message that the checkpoint has been reached
		uint64_t SendCheckPointMessage(int argument = 0);
		
		// Send a sequence to a specific id
		uint64_t SendSeqMessage(int destID, const PackedSeq& seq, APSeqOperation operation);
		
		// Send a sequence extension message
		uint64_t SendSeqExtMessage(int destID, const PackedSeq& seq, APSeqExtOperation operation, extDirection dir, SeqExt ext, char base = 0);
		
		// Send a sequence flag message
		uint64_t SendSeqFlagMessage(int destID, const PackedSeq& seq, APSeqFlagOperation operation, SeqFlag flag);
				
		// Send an adjacency result message
		void SendAdjacencyResult(int destID, uint64_t reqID, bool b);
		
		// Send a bool result
		void SendResultMessage(int destID, uint64_t reqID, bool b);
		
		// Send a result
		void SendResultMessage(int destID, uint64_t reqID, ResultPair rp);
		
		// Receive a seq message
		SeqMessage ReceiveSeqMessage();
		
		// Receive a seq message
		SeqExtMessage ReceiveSeqExtMessage();
		
		// Receive a seq message
		SeqFlagMessage ReceiveSeqFlagMessage();		
			
		// Receive a control message
		ControlMessage ReceiveControlMessage();
		
		// Receive a result message
		ResultMessage ReceiveResultMessage();
		
		// Flush the buffer
		void flush();
		
	private:
		int m_id;
		int m_kmerSize;
		int m_numBytesPerSeq;
		int m_bufferSize;
		uint64_t m_msgID;
		char* m_buffer;
};

#endif
