#ifndef PATH_H
#define PATH_H

#include <list>
#include "Sequence.h"
#include "SeqRecord.h"


typedef std::list<Sequence>::iterator seq_list_iter;
typedef std::list<Sequence>::const_iterator const_seq_list_iter;

class Path
{

	public:
		Path(Sequence seedSeq, extDirection dir);
		
		// add a sequence to the path
		void addToPath(const Sequence& seq, bool antiDir);
		
		// merge two paths together
		void mergePath(const Path& path2, bool reverseComp, bool antiDir, bool skipFirst);
		
		// get the last sequence added to the vector
		const Sequence& getCurrentNode() const;
		
		// get an arbitrary node
		const Sequence& getNode(const int index) const;
		
		// get the number of nodes
		int getNumNodes() const;
		
		// get the extension
		void getPath(std::list<Sequence>& outPath) const;
		
		// get a string representation of the contig
		Sequence getSequence() const;
		
		// check if this path contains the requested sequence
		bool contains(const Sequence& seq) const;
		
		// return the direction this path is extending in
		extDirection getGrowthDirection() const;
		
		// get the length of the path
		int getPathLength() const;
		
		// print the full path
		void print() const;
		
	private:
	
		// extension vectors, one for each direction
		std::list<Sequence> m_extensions;
		
		// a cache of the sequences in this path for quick access
		SeqRecord m_seqCache;
		
		// The sequence this path started from
		Sequence m_seed;
		
		// The direction this path is extending
		extDirection m_growDir;
		
		// the length of the path
		int m_length;
};

#endif
