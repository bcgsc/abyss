/**
 * A minifloat like datatype for probablistic log counts (PLC) of elements
 * Unsigned generic implementation
 * Mantissa = 3 bits
 * Exponent = 5 bits
 * Copyright 2014 bcgsc
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>

using namespace std;

static const unsigned mantissa = 3;
static const uint8_t mantiMask = 0xFF >> (8 - mantissa);
static const uint8_t addMask = 0x80 >> (7 - mantissa);

class plc {
public:
	plc()
	{
		m_val = 0;
	}

	void operator++()
	{
		//from 0-1
		if (m_val <= mantiMask) {
			++m_val;
		} else {
			//this shifts the first bit off and creates the value
			//need to get the correct transition probability
			size_t shiftVal = 1 << ((m_val >> mantissa) - 1);
			if (rand() % shiftVal == 0) {
				++m_val;
			}
		}
	}

	float toFloat()
	{
		if (m_val <= mantiMask)
			return float(m_val);
		return ldexp((m_val & mantiMask) | addMask, (m_val >> mantissa) - 1);
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
