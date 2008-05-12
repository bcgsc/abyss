#include <math.h>
#include "Prb.h"


Prb::Prb(std::string tuple)
{
	bool stop = false;
	int currPos = 0;
	int index = 0;
	while(!stop)
	{
		// the prb values are delimited by commas
		int nextPos = tuple.find(',', currPos);
		int len;
		if(nextPos != -1)
		{
			len = nextPos - currPos;
		}
		else
		{
			len = tuple.length() - currPos;
			stop = true;
		}
		
		// read the next value and add it to the prb values
		std::string valueStr = tuple.substr(currPos, len);
		int val = atoi(valueStr.c_str());
		
		m_values.push_back(val);
		
		currPos = nextPos + 1;
		index++;
	}
	
	// sort and cache the result
	m_sortedValues = m_values;
	std::sort(m_sortedValues.begin(), m_sortedValues.end());
}

int Prb::getValue(const char base) const
{
	int index = baseToIndex(base);
	return getValue(index);	
}

int Prb::getValue(const int index) const
{
	assert(index < (int)m_values.size());
	return m_values[index];	
}

double Prb::QsToP(int value)
{
	double temp = static_cast<double>(1) + pow(10.0f, static_cast<double>(value) / 10.0f);
	return 1.0f / temp;	
}

int Prb::baseToIndex(char base)
{
	switch(base)
	{
		case 'A':
		case 'a':
			return 0;
		case 'C':
		case 'c':
			return 1;
		case 'G':
		case 'g':
			return 2;
		case 'T':
		case 't':
			return 3;
	}
	
	assert(false);	
	return 0;
}

int Prb::getOrderStat(int rank) const
{
	rank -= 1;
	assert(rank < (int)m_sortedValues.size());
	return m_sortedValues[rank];
}
