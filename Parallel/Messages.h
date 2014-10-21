#ifndef MESSAGES_H
#define MESSAGES_H 1

#include "SequenceCollection.h"
#include <ostream>

class NetworkSequenceCollection;

enum MessageType
{
	MT_VOID,
	MT_ADD,
	MT_REMOVE,
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

typedef uint32_t IDType;

/** The base class of all interprocess messages. */
class Message
{
	public:
		typedef SequenceCollectionHash Graph;
		typedef graph_traits<Graph>::vertex_descriptor V;
		typedef Graph::Symbol Symbol;
		typedef Graph::SymbolSet SymbolSet;
		typedef Graph::SymbolSetPair SymbolSetPair;

		Message() { }
		Message(const V& seq) : m_seq(seq) { }
		virtual ~Message() { }

		virtual void handle(
				int senderID, NetworkSequenceCollection& handler) = 0;

		virtual size_t getNetworkSize() const
		{
			return sizeof (uint8_t) // MessageType
				+ V::serialSize();
		}

		static MessageType readMessageType(char* buffer);
		virtual size_t serialize(char* buffer) = 0;
		virtual size_t unserialize(const char* buffer);

		friend std::ostream& operator <<(std::ostream& out,
				const Message& message)
		{
			return out << message.m_seq.str() << '\n';
		}

		V m_seq;
};

/** Add a vertex. */
class SeqAddMessage : public Message
{
	public:
		SeqAddMessage() { }
		SeqAddMessage(const V& seq) : Message(seq) { }

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);

		static const MessageType TYPE = MT_ADD;
};

/** Remove a vertex. */
class SeqRemoveMessage : public Message
{
	public:
		SeqRemoveMessage() { }
		SeqRemoveMessage(const V& seq) : Message(seq) { }

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);

		static const MessageType TYPE = MT_REMOVE;
};

/** Set a flag. */
class SetFlagMessage : public Message
{
	public:
		SetFlagMessage() { }
		SetFlagMessage(const V& seq, SeqFlag flag)
			: Message(seq), m_flag(flag) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize() + sizeof m_flag;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);

		static const MessageType TYPE = MT_SET_FLAG;
		uint8_t m_flag; // SeqFlag
};

/** Remove an edge. */
class RemoveExtensionMessage : public Message
{
	public:
		RemoveExtensionMessage() { }
		RemoveExtensionMessage(const V& seq,
				extDirection dir, SymbolSet ext)
			: Message(seq), m_dir(dir), m_ext(ext) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_dir + sizeof m_ext;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);

		static const MessageType TYPE = MT_REMOVE_EXT;
		uint8_t m_dir; // extDirection
		SymbolSet m_ext;
};

/** Request vertex properties. */
class SeqDataRequest : public Message
{
	public:
		SeqDataRequest() { }
		SeqDataRequest(const V& seq, IDType group, IDType id)
			: Message(seq), m_group(group), m_id(id) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_group + sizeof m_id;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);

		static const MessageType TYPE = MT_SEQ_DATA_REQUEST;
		IDType m_group;
		IDType m_id;
};

/** The response to a request for vertex properties. */
class SeqDataResponse : public Message
{
	public:
		SeqDataResponse() { }
		SeqDataResponse(const V& seq, IDType group, IDType id,
				SymbolSetPair& extRecord, int multiplicity) :
			Message(seq), m_group(group), m_id(id),
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

		static const MessageType TYPE = MT_SEQ_DATA_RESPONSE;
		IDType m_group;
		IDType m_id;
		SymbolSetPair m_extRecord;
		uint16_t m_multiplicity;
};

/** Add an edge. */
class SetBaseMessage : public Message
{
	public:
		SetBaseMessage() { }
		SetBaseMessage(const V& seq,
				extDirection dir, Symbol base)
			: Message(seq), m_dir(dir), m_base(base) { }

		size_t getNetworkSize() const
		{
			return Message::getNetworkSize()
				+ sizeof m_dir + sizeof m_base;
		}

		void handle(int senderID, NetworkSequenceCollection& handler);
		size_t serialize(char* buffer);
		size_t unserialize(const char* buffer);

		static const MessageType TYPE = MT_SET_BASE;
		uint8_t m_dir; // extDirection
		Symbol m_base;
};

#endif
