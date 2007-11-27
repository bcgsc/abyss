#ifndef PRB_H
#define PRB_H

#include <string>
#include <assert.h>
#include <vector>

class Prb
{
	public:
		Prb(int A, int C, int G, int T);
		Prb(std::string tuple);
		
		// get value by the base
		int getValue(const char base) const;
		
		// get ranking value 1 = smallest, 4 = biggest
		int getOrderStat(int rank) const;
		
		static double QsToP(int value);
		
	private:
	
		// get a value by index
		int getValue(const int index) const;

		
		static int baseToIndex(char base);
		
		// the actual values
		std::vector<short> m_values;
		
		// cached sorted list of sorted values
		// this is computed after construction; if the setters are called afterwards this will be invalid
		// it is used for order statistics
		std::vector<short> m_sortedValues;
	
};

#endif
