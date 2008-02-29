#ifndef NETWORKDEFS_H
#define NETWORKDEFS_H

// The types of messages that we can send
enum APMessage
{
	APM_NONE,
	APM_SEQ,
	APM_SEQ_FLAG,
	APM_SEQ_EXT,
	APM_CONTROL,	
	APM_RESULT
};

// The type of operations on whole sequences that can be performed
enum APSeqOperation
{
	APO_ADD,
	APO_REMOVE,
	APO_CHECKSEQ,
	APO_HAS_CHILD,
	APO_HAS_PARENT
};

enum APSeqFlagOperation
{
	APSFO_CHECK,
	APSFO_SET
};

enum APSeqExtOperation
{
	APSEO_CHECK,
	APSEO_SET,
	APSEO_REMOVE
};


// The type of control messages that can be sent
enum APControl
{
	APC_LOAD,
	APC_DONELOAD,
	APC_GEN_ADJ,
	APC_TRIM,
	APC_ASSEMBLE,
	APC_CHECKPOINT,	
	APC_FINISHED	
};

// The type of results that can be sent
enum APResult
{
	APR_NONE,
	APR_TRUE,
	APR_FALSE
};

#endif
