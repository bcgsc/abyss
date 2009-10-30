#include "Messages.h"
#include "NetworkSequenceCollection.h"
#include <cstdio>
#include <cstring>

static size_t serializeData(const void* ptr, char* buffer,
		size_t size)
{
	memcpy(buffer, ptr, size);
	return size;
}

static size_t unserializeData(void* ptr, const char* buffer,
		size_t size)
{
	memcpy(ptr, buffer, size);
	return size;
}

MessageType Message::readMessageType(char* buffer)
{
	return (MessageType)*(uint8_t*)buffer;
}

size_t Message::unserialize(const char* buffer)
{
	size_t offset = 0;
	offset++; // MessageType
	offset += m_seq.unserialize(buffer + offset);
	return offset;
}

size_t SeqAddMessage::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	return offset;
}

void SeqAddMessage::handle(int senderID,
		NetworkSequenceCollection& handler)
{
	handler.handleSeqAddMessage(senderID, *this);
}

void SeqAddMessage::print() const
{
	printf("Message type: %d Sequence: %s\n",
			TYPE, m_seq.decode().c_str());
}

size_t SeqRemoveMessage::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	return offset;
}

void SeqRemoveMessage::handle(int senderID,
		NetworkSequenceCollection& handler)
{
	handler.handleSeqRemoveMessage(senderID, *this);
}

void SeqRemoveMessage::print() const
{
	printf("Message type: %d Sequence: %s\n",
			TYPE, m_seq.decode().c_str());
}

size_t SetFlagMessage::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	offset += serializeData(&m_flag, buffer + offset, sizeof(m_flag));
	return offset;
}

size_t SetFlagMessage::unserialize(const char* buffer)
{
	size_t offset = 0;
	offset += Message::unserialize(buffer);
	offset += unserializeData(&m_flag, buffer + offset, sizeof(m_flag));
	return offset;
}

void SetFlagMessage::handle(int senderID, NetworkSequenceCollection& handler)
{
	handler.handleSetFlagMessage(senderID, *this);
}

void SetFlagMessage::print() const
{
	printf("Message type: %d Sequence: %s Flag: %d\n",
			TYPE, m_seq.decode().c_str(), (int)m_flag);
}

size_t RemoveExtensionMessage::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	offset += serializeData(&m_dir, buffer + offset, sizeof m_dir);
	offset += serializeData(&m_ext, buffer + offset, sizeof m_ext);
	return offset;
}

size_t RemoveExtensionMessage::unserialize(const char* buffer)
{
	size_t offset = 0;
	offset += Message::unserialize(buffer);
	offset += unserializeData(&m_dir, buffer + offset, sizeof m_dir);
	offset += unserializeData(&m_ext, buffer + offset, sizeof m_ext);
	return offset;
}

void RemoveExtensionMessage::handle(int senderID, NetworkSequenceCollection& handler)
{
	handler.handleRemoveExtensionMessage(senderID, *this);
}

void RemoveExtensionMessage::print() const
{
	printf("Message type: %d Sequence: %s Dir: %d ",
			TYPE, m_seq.decode().c_str(), (int)m_dir);
	m_ext.print();
}

size_t SetBaseMessage::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	offset += serializeData(&m_dir, buffer + offset, sizeof(m_dir));
	offset += serializeData(&m_base, buffer + offset, sizeof(m_base));
	return offset;
}

size_t SetBaseMessage::unserialize(const char* buffer)
{
	size_t offset = 0;
	offset += Message::unserialize(buffer);
	offset += unserializeData(&m_dir, buffer + offset, sizeof(m_dir));
	offset += unserializeData(&m_base, buffer + offset, sizeof(m_base));
	return offset;
}

void SetBaseMessage::handle(int senderID, NetworkSequenceCollection& handler)
{
	handler.handleSetBaseMessage(senderID, *this);
}

void SetBaseMessage::print() const
{
	printf("Message type: %d Sequence: %s Dir: %d Base: %c\n",
			TYPE, m_seq.decode().c_str(), (int)m_dir, m_base);
}

size_t SeqDataRequest::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	offset += serializeData(&m_group, buffer + offset, sizeof(m_group));
	offset += serializeData(&m_id, buffer + offset, sizeof(m_id));
	return offset;
}

size_t SeqDataRequest::unserialize(const char* buffer)
{
	size_t offset = 0;
	offset += Message::unserialize(buffer);
	offset += unserializeData(&m_group, buffer + offset, sizeof(m_group));
	offset += unserializeData(&m_id, buffer + offset, sizeof(m_id));
	return offset;
}

void SeqDataRequest::handle(int senderID, NetworkSequenceCollection& handler)
{
	handler.handleSequenceDataRequest(senderID, *this);
}

void SeqDataRequest::print() const
{
	printf("Message type: %d Sequence: %s group: %d id: %d\n",
			TYPE, m_seq.decode().c_str(), m_group, m_id);
}

size_t SeqDataResponse::serialize(char* buffer)
{
	size_t offset = 0;
	buffer[offset++] = TYPE;
	offset += m_seq.serialize(buffer + offset);
	offset += serializeData(&m_group, buffer + offset, sizeof(m_group));
	offset += serializeData(&m_id, buffer + offset, sizeof(m_id));
	offset += serializeData(&m_extRecord, buffer + offset, sizeof(m_extRecord));
	offset += serializeData(&m_multiplicity, buffer + offset, sizeof(m_multiplicity));
	return offset;
}

size_t SeqDataResponse::unserialize(const char* buffer)
{
	size_t offset = 0;
	offset += Message::unserialize(buffer);
	offset += unserializeData(&m_group, buffer + offset, sizeof(m_group));
	offset += unserializeData(&m_id, buffer + offset, sizeof(m_id));
	offset += unserializeData(&m_extRecord, buffer + offset, sizeof(m_extRecord));
	offset += unserializeData(&m_multiplicity, buffer + offset, sizeof(m_multiplicity));
	return offset;
}

void SeqDataResponse::handle(int senderID, NetworkSequenceCollection& handler)
{
	handler.handleSequenceDataResponse(senderID, *this);
}

void SeqDataResponse::print() const
{
	printf("Message type: %d Sequence: %s Multiplicity: %d group: %d id: %d\n",
			TYPE, m_seq.decode().c_str(),
			m_multiplicity, m_group, m_id);
	m_extRecord.dir[SENSE].print();
	m_extRecord.dir[ANTISENSE].print();
}
