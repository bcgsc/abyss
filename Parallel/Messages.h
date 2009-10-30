#ifndef MESSAGES_H
#define MESSAGES_H 1

#include "Kmer.h"
#include "KmerData.h"
#include "NetworkDefs.h"
#include <ostream>

class NetworkSequenceCollection;

enum MessageType
{
	MT_VOID,
	MT_SEQ_OP,
	MT_SET_FLAG,
	MT_REMOVE_EXT,
	MT_SEQ_DATA_REQUEST,
	MT_SEQ_DATA_RESPONSE,
	MT_SET_BASE
};

enum MessageOp
{
	MO_VOID,
	MO_ADD,
	MO_REMOVE,
};

size_t serializeData(void* ptr, char* buffer, size_t size);


typedef uint32_t IDType;

class Message
{
	public:
		Message(MessageType m_type) : m_type(m_type) { }
		Message(const Kmer& seq, MessageType m_type)
			: m_type(m_type), m_seq(seq) { }
		virtual ~Message() { }

		virtual void handle(int senderID, NetworkSequenceCollection& handler) = 0;

		virtual size_t getNetworkSize() const
		{
			return sizeof m_type + Kmer::serialSize();
		}

		static MessageType readMessageType(char* buffer);
		virtual size_t serialize(char* buffer);
		virtual size_t unserialize(const char* buffer);
		virtual void print() const { }

		friend std::ostream& operator <<(std::ostream& o,
				const Message& m)
		{
			m.print();
			return o;
		}

		MessageType m_type;
		Kmer m_seq;
};

class SeqOpMessage : public Message
{
	public:
		SeqOpMessage() : Message(MT_SEQ_OP) { }
		SeqOpMessage(const Kmer& seq, MessageOp operation)
			: Message(seq, MT_SEQ_OP), m_operation(operation) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize() + sizeof m_operation;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);
		void print() const;

		uint8_t m_operation; // MessageOp
};

class SetFlagMessage : public Message
{
	public:
		SetFlagMessage() : Message(MT_SET_FLAG) { }
		SetFlagMessage(const Kmer& seq, SeqFlag flag)
			: Message(seq, MT_SET_FLAG), m_flag(flag) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize() + sizeof m_flag;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);
		void print() const;

		uint8_t m_flag; // SeqFlag
};

class RemoveExtensionMessage : public Message
{
	public:
		RemoveExtensionMessage() : Message(MT_REMOVE_EXT) { }
		RemoveExtensionMessage(const Kmer& seq,
				extDirection dir, SeqExt ext)
			: Message(seq, MT_REMOVE_EXT), m_dir(dir), m_ext(ext) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_dir + sizeof m_ext;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);
		void print() const;

		uint8_t m_dir; // extDirection
		SeqExt m_ext;
};

class SeqDataRequest : public Message
{
	public:
		SeqDataRequest() : Message(MT_SEQ_DATA_REQUEST) { }
		SeqDataRequest(const Kmer& seq, IDType group, IDType id)
			: Message(seq, MT_SEQ_DATA_REQUEST),
				m_group(group), m_id(id) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_group + sizeof m_id;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);
		void print() const;

		IDType m_group;
		IDType m_id;
};

class SeqDataResponse : public Message
{
	public:
		SeqDataResponse() : Message(MT_SEQ_DATA_RESPONSE) { }
		SeqDataResponse(const Kmer& seq, IDType group, IDType id,
				ExtensionRecord& extRecord, int multiplicity) :
			Message(seq, MT_SEQ_DATA_RESPONSE),
			m_group(group), m_id(id),
			m_extRecord(extRecord), m_multiplicity(multiplicity) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_group + sizeof m_id
				+ sizeof m_extRecord + sizeof m_multiplicity;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);
		void print() const;

		IDType m_group;
		IDType m_id;
		ExtensionRecord m_extRecord;
		uint16_t m_multiplicity;
};

class SetBaseMessage : public Message
{
	public:
		SetBaseMessage() : Message(MT_SET_BASE) { }
		SetBaseMessage(const Kmer& seq,
				extDirection dir, uint8_t base)
			: Message(seq, MT_SET_BASE), m_dir(dir), m_base(base) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_dir + sizeof m_base;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);
		void print() const;

		uint8_t m_dir; // extDirection
		uint8_t m_base;
};

#endif
