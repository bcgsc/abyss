#ifndef COMMLAYER
#define COMMLAYER

#include <mpi.h>
#include "PackedSeq.h"

enum APMessage
{
	APM_DATA,
	APM_SEQADD,
	APM_SEQDEL,
	APM_SEQCHECK,
	APM_DONELOAD,
	APM_FINISHED
};

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
		void SendControlMessage(int destID, APMessage msg);
		
		// Send a sequence to a specific id
		void SendSequence(int destID, const PackedSeq& seq, APMessage msg);
		
		// Check if this sequence exists in the phase space
		bool CheckForSequence(int destID, PackedSeq& seq);
		
		// Receive a sequence
		PackedSeq ReceiveSequence();
		
		// Send/Receive a bool
		void SendBool(int destID, bool b);
		bool ReceiveBool(int sendID);
		
		
		// Clear out a control message
		void ClearControlMessage();
		
	private:
		int m_id;
		int m_kmerSize;
		int m_numBytesPerSeq;
		char* m_pRecvBuffer;
};

#endif
