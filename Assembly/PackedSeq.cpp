#include "config.h"
#include "Common/Options.h"
#include "PackedSeq.h"
#include "HashFunction.h"
#include "Sequence.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;

//
// Serialize this packed seq to the buffer
// TODO: make this machine independent (using mpi datatypes?)
// Return the number of bytes read
//
size_t PackedSeq::serialize(char* buffer) const
{
	size_t offset = 0;
	
	memcpy(buffer + offset, &m_seq, sizeof(m_seq));
	offset += sizeof(m_seq);
	
	memcpy(buffer + offset, &m_length, sizeof(m_length));
	offset += sizeof(m_length);
	
	memcpy(buffer + offset, &m_flags, sizeof(m_flags));
	offset += sizeof(m_flags);
	
	memcpy(buffer + offset, m_multiplicity, sizeof(m_multiplicity));
	offset += sizeof(m_multiplicity);	
	
	memcpy(buffer + offset, &m_extRecord, sizeof(m_extRecord));
	offset += sizeof(m_extRecord);
	
	assert(offset == serialSize());

	return offset;		
}

//
// Unserialize this packed seq from the buffer
// TODO: make this machine independent (using mpi datatypes?)
//
size_t PackedSeq::unserialize(const char* buffer)
{
	size_t offset = 0;
	
	memcpy(m_seq, buffer + offset, sizeof(m_seq));
	offset += sizeof(m_seq);
	
	memcpy(&m_length, buffer + offset, sizeof(m_length));
	offset += sizeof(m_length);
	
	memcpy(&m_flags, buffer + offset, sizeof(m_flags));
	offset += sizeof(m_flags);
	
	memcpy(m_multiplicity, buffer + offset, sizeof(m_multiplicity));
	offset += sizeof(m_multiplicity);	
	
	memcpy(&m_extRecord, buffer + offset, sizeof(m_extRecord));
	offset += sizeof(m_extRecord);
	
	assert(offset == serialSize());
	
	return offset;			
}
