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

struct ResultPairMessage
{
	int64_t id;
	APResult result[2];
};

struct SequenceExtensionRequestMessage
{
	int64_t id;
	uint64_t groupID;
	uint64_t branchID;
	PackedSeq seq;
};

struct SequenceExtensionResponseMessage
{
	int64_t id;
	uint64_t groupID;
	uint64_t branchID;	
	PackedSeq seq;
	extDirection dir;	
	ExtensionRecord extRec;
};

struct AdjacencyMessage
{
	int64_t id;
	PackedSeq originalSeq;
	PackedSeq testSeq;	
	extDirection dir;
	char base;	
	
};

struct AdjacencyResultMessage
{
	int64_t id;
	PackedSeq originalSeq;
	extDirection dir;
	char base;
	APResult result[2];
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
		
		// Send a request to a remote node to see if testSequence is present which means originalSeq has extension b in direction dir
		uint64_t SendAdjacencyRequest(int destID, const PackedSeq& testSequence, const PackedSeq& originalSeq, extDirection dir, char b);
		
		// Send a request for the extension of a sequence
		uint64_t SendExtensionRequest(int destID, uint64_t groupID, uint64_t branchID, const PackedSeq& seq);
		
		// Send a request for the extension of a sequence
		void SendExtensionResponse(int destID, uint64_t reqID, uint64_t groupID, uint64_t branchID, const PackedSeq& seq, const ExtensionRecord& extRec);
		
		// Send an adjacency result message, send back the original sequence and the base that was requested
		// This is done to help hide the latency. By not requiring the sender to save the state it can proceed ahead at will
		void SendAdjacencyResult(int destID, uint64_t reqID, PackedSeq originalSeq, extDirection dir, char base, bool b);
		
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
		ResultPairMessage ReceiveResultMessage();
		
		// Receieve an adjacency message
		AdjacencyMessage ReceiveAdjacencyMessage();
		
		// Receive an adjacency result message
		AdjacencyResultMessage ReceiveAdjacencyResultMessage();
		
		// Receive a sequence extension message
		SequenceExtensionRequestMessage ReceiveSequenceExtensionMessage();
		
		// Receive a sequence extension result message
		SequenceExtensionResponseMessage ReceiveSequenceExtensionResponseMessage();
		
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
