#ifndef COMMLAYER
#define COMMLAYER

#include <mpi.h>
#include "PackedSeq.h"

// The types of messages that we can send
enum APMessage
{
	APM_NONE,
	APM_SEQ,
	APM_SEQ_FLAG,
	APM_SEQ_EXT,
	APM_CONTROL,	
	APM_RESULT
};

// The type of operations on whole sequences that can be performed
enum APSeqOperation
{
	APO_ADD,
	APO_REMOVE,
	APO_CHECKSEQ,
	APO_HAS_CHILD,
	APO_HAS_PARENT
};

enum APSeqFlagOperation
{
	APSFO_CHECK,
	APSFO_SET
};

enum APSeqExtOperation
{
	APSEO_CHECK,
	APSEO_SET,
	APSEO_REMOVE
};


// The type of control messages that can be sent
enum APControl
{
	APC_LOAD,
	APC_DONELOAD,
	APC_GEN_ADJ,
	APC_TRIM,
	APC_ASSEMBLE,
	APC_CHECKPOINT,	
	APC_FINISHED	
};

// The type of results that can be sent
enum APResult
{
	APR_NONE,
	APR_TRUE,
	APR_FALSE
};


struct SeqMessage
{
	APSeqOperation operation;
	PackedSeq seq;
};

struct SeqFlagMessage
{
	APSeqFlagOperation operation;
	PackedSeq seq;
	SeqFlag flag;
};

struct SeqExtMessage
{
	APSeqExtOperation operation;
	PackedSeq seq;
	extDirection dir;	
	SeqExt ext;
	char base;
	
};

struct ControlMessage
{
	APControl msgType;
};

struct ResultMessage
{
	APResult result;
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
		void SendControlMessage(int numNodes, APControl m) const;
		
		void SendCheckPointMessage() const;
		
		// Send a sequence to a specific id
		void SendSeqMessage(int destID, const PackedSeq& seq, APSeqOperation operation) const;
		
		// Send a sequence extension message
		void SendSeqExtMessage(int destID, const PackedSeq& seq, APSeqExtOperation operation, extDirection dir, SeqExt ext, char base = 0) const;
		
		// Send a sequence flag message
		void SendSeqFlagMessage(int destID, const PackedSeq& seq, APSeqFlagOperation operation, SeqFlag flag) const;
		
		// Send a result
		void SendResultMessage(int destID, bool r);
		
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
		char* m_buffer;
};

#endif
