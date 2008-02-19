#ifndef DISTRIBUTEDPHASESPACE_H
#define DISTRIBUTEDPHASESPACE_H

#include <mpi.h>
#include "SimpleSequenceSpace.h"
#include "CommLayer.h"

// The distributed phase space is a wrapped for the phase space that can operate between various processes
class DistributedPhaseSpace
{
	public:
		DistributedPhaseSpace(int myID, int kmerSize);
		~DistributedPhaseSpace();
		
		void MessageLoop();
			
	private:
		DistributedPhaseSpace();
			
		CommLayer* pComm;
		
		SimpleSequenceSpace* m_pSS;
};

#endif
