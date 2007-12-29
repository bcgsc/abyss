#include "Path.h"
#include "CommonUtils.h"

Path::Path(Sequence seedSeq, extDirection dir) : m_seed(seedSeq), m_growDir(dir), m_length(0)
{	
	// add the seed sequence to both paths
	addToPath(m_seed, false);
}

const Sequence& Path::getCurrentNode() const
{	
	if(m_growDir == SENSE)
	{
		return m_extensions.back();
	}
	else
	{
		return m_extensions.front();	
	}
}


const Sequence& Path::getNode(const int index) const
{
	const_seq_list_iter iter = m_extensions.begin();
	for(int i = 0; i < index; i++)
	{
		iter++;
	}
	return *iter;
}		

int Path::getNumNodes() const
{
	return m_extensions.size();
}

void Path::addToPath(const Sequence& seq, bool antiDir)
{
	// if antidir is true, add the node in the opposite direction of growth
	extDirection dir = antiDir ? oppositeDirection(m_growDir) : m_growDir;
	
	if(dir == SENSE)
	{
		m_extensions.push_back(seq);	
	}
	else
	{
		m_extensions.push_front(seq);	
	}
	m_length++;
	m_seqCache.addSequence(seq);
}

bool Path::contains(const Sequence& seq) const
{
	return m_seqCache.contains(seq);
}

void Path::mergePath(const Path& path2, bool reverseComp, bool antiDir, bool skipFirst)
{

	std::list<Sequence> path2Full;
	path2.getPath(path2Full);
	
	printf("merging\npath1: ");
	print();
	printf("path2: ");		
	path2.print();

		
	for(const_seq_list_iter iter = path2Full.begin(); iter != path2Full.end(); iter++)
	{
		if(skipFirst)
		{
			skipFirst = false;
			continue;
		}
		
		Sequence newSeq;
		
		if(!reverseComp)
		{
			newSeq = *iter;
		}
		else
		{
			newSeq = reverseComplement(*iter);	
		}
		
		addToPath(newSeq, antiDir);
	}
	
	printf("merged path: ");
	print();
}

void Path::getPath(std::list<Sequence>& outPath) const
{
	// push all sequences in this path onto the output vector in the antisense -> sense direction
	const_seq_list_iter iter;
	for(iter = m_extensions.begin(); iter != m_extensions.end(); iter++)
	{
		outPath.push_back(*iter);
	}
	
	if(m_growDir == ANTISENSE)
	{
		outPath.reverse();
	}
}

Sequence Path::getSequence() const
{
	Sequence fullSeq("");
	
	// append the first full sequence
	const_seq_list_iter iter = m_extensions.begin();
	seqAppend(fullSeq, *iter);
	iter++;
	
	Sequence prevSeq;
	for(; iter != m_extensions.end(); iter++)
	{
		int len = iter->length();
		seqAppendBase(fullSeq, (*iter)[len - 1]);
		
		if(prevSeq == *iter)
		{
			printf("SEQUENCE MATCH!! %s %s\n", prevSeq.c_str(), iter->c_str());
		}
		
		prevSeq = *iter;
	}
	
	return fullSeq;
	
}

extDirection Path::getGrowthDirection() const
{
	return m_growDir;	
}

int Path::getPathLength() const
{
	return m_length;	
}

void Path::print() const
{
	printf("path: %s\n", getSequence().c_str());
}
