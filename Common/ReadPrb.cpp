#include "ReadPrb.h"

double calculatePError(const ReadPrb& readPrb)
{
	// calculate the Perror of this read based on the PRB values
	double PCorrect = 1.0f;
	for(std::vector<Prb>::const_iterator iter = readPrb.begin(); iter != readPrb.end(); iter++)
	{
			double baseCorrect = 1.0f - Prb::QsToP(iter->getOrderStat(4));
			PCorrect *= baseCorrect;
	}	
	return 1.0f - PCorrect;
}
