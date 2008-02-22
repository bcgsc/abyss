#include "DistributedPhaseSpace.h"
#if 0
DistributedPhaseSpace::DistributedPhaseSpace(int myID, int kmerSize)
{
	int maxC = kmerSize - 1;
	Coord4 minCoords = {0, 0, 0, 0};
	Coord4 maxCoords = {maxC, maxC, maxC, maxC};
	
	// Load the phase space
	m_pSS = new SimpleSequenceSpace();
	
	// Create the comm layer
	pComm = new CommLayer(myID, kmerSize);
}

DistributedPhaseSpace::~DistributedPhaseSpace()
{
	delete m_pSS;
	m_pSS = 0;
	
	delete pComm;
	pComm = 0;
}

void DistributedPhaseSpace::MessageLoop()
{
	bool stop = false;
	
	while(!stop)
	{
		int sendID;
		APMessage msg = pComm->CheckMessage(sendID);
		
		switch(msg)
		{
			case APM_SEQADD:
				{
					PackedSeq seq = pComm->ReceiveSequence();
					m_pSS->addSequence(seq);
					break;
				}
			case APM_DONELOAD:
				{
					pComm->ClearControlMessage();
					printf("loaded: %d sequences\n", m_pSS->countAll());
					printf("finalizing\n");
					m_pSS->finalize();
					break;
				}
			case APM_SEQCHECK:
				{
					PackedSeq seq = pComm->ReceiveSequence();
					pComm->SendBool(sendID, m_pSS->checkForSequence(seq));
					break;
				}				
			case APM_FINISHED:
				{
					pComm->ClearControlMessage();
					printf("finished message received, exiting\n");
					stop = true;
					break;
				}	
			default:
				printf("Unhandled message code: %d\n", msg);
		}
	}
}
#endif