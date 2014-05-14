/**
 * A minifloat like datatype for probablistic log counts (PLC) of elements
 * Unsigned generic implementation
 * Mantissa = 1 bits
 * Exponent = 7 bits
 * Copyright 2014 bcgsc
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>

using namespace std;

static const unsigned mantissa = 1;

class plc {
public:
	plc()
	{
		m_val = 0;
	}

	void operator++()
	{
		//check if at max value
		if (m_val > 0xFF) {
			return;
		}
		//from 0-1
		if (m_val <= mantissa) {
			++m_val;
		} else {
			//this shifts the first bit off and creates the value
			//need to get the correct transition probability
			unsigned shiftVal = 1 << ((m_val >> mantissa) - 1);
			if (rand() % shiftVal == 0) {
				++m_val;
			}
		}
	}

	float toFloat()
	{
		if (m_val <= mantissa)
			return float(m_val);
		//the following needs to modified if mantissa is changed
		return float((m_val & 0x01) | 0x02)
				* pow(2.0, (m_val >> mantissa) - mantissa);
	}

	/*
	 * return raw value of byte use to store value
	 */
	uint8_t rawValue()
	{
		return m_val;
	}

private:
	uint8_t m_val;
};

//static float toFloat(plc n)
//{
//	return n;
//}
