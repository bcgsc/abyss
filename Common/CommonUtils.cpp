#include <math.h>
#include "CommonUtils.h"

extDirection oppositeDirection(extDirection dir)
{
	if(dir == SENSE)
	{
		return ANTISENSE;
	}
	else
	{
		return SENSE;
	}	
}

double entropy(const Sequence& s)
{
	std::map<char, int> vals;
	
	int len = s.length();
	for(int i = 0; i < len; i++)
	{
		vals[s[i]]++;
	}
	
	std::map<char, int>::iterator iter;
	
	double sum = 0.0f;
	for(iter = vals.begin(); iter != vals.end(); iter++)
	{
		double f = static_cast<double>(iter->second) / static_cast<double>(len);
		double t = log(f) / log(2.0f);
		
		sum += (f*t);
	}
	
	return -sum;
}

//
//
//
int min(const int& n1, const int& n2)
{
	return (n1 < n2) ? n1 : n2;	
}

//
//
//
int max(const int& n1, const int& n2)
{
	return (n1 > n2) ? n1 : n2;	
}

//
// Return the complement of the specified base.
//
char complementBaseChar(char b)
{
	switch (b) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default:
			assert(false);
			return 0;
	}			
}

void PrintBufferAsHex(char* buffer, int length)
{
	int ptr = 0;
	for(;ptr < length;ptr++)
	{
		printf("%02X ",(unsigned char)*(buffer+ptr));
	}
	printf("\n");
	
}
