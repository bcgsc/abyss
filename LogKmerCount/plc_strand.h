/**
 * A minifloat like datatype for probablistic log counts (PLC) of elements
 * Contains unique functionalities for minimizing collision
 * Mantissa = 3 bits
 * Exponent = 4 bits
 * ControlBit = 1 bits
 * Copyright 2014 bcgsc
 */

static const uint8_t numericBitMask = 0x7F;
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
		//extract numeric section of plc
		uint8_t maskedVal = m_val & numericBitMask;

		//check if at max value
		if (maskedVal == numericBitMask) {
			return;
		}

		if (maskedVal <= mantissa) {
			++m_val;
		} else {
			//this shifts the first bit off and creates the value
			//need to get the correct transition probability
			unsigned shiftVal = 1 << ((maskedVal >> mantissa) - 1);
			if (rand() % shiftVal == 0) {
				++maskedVal;
			}
		}
	}

	bool operator==(plc val)
	{
		return val.rawValue() & numericBitMask == m_val & numericBitMask;
	}

	operator bool()
	{
		return m_val & numericBitMask;
	}

	float toFloat()
	{
		//extract numeric section of plc
		uint8_t maskedVal = m_val & numericBitMask;

		if (maskedVal <= mantissa)
			return float(maskedVal);
		return ldexp((maskedVal & mantiMask) | addMask,
				(maskedVal >> mantissa) - 1);
	}

	/*
	 * return raw value of byte use to store value
	 */
	uint8_t rawValue()
	{
		return m_val;
	}

	void setStrandBit()
	{
		m_val &= 0x80;
	}

private:
	uint8_t m_val;
};
