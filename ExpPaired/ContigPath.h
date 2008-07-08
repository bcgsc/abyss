#ifndef CONTIGPATH_H
#define CONTIGPATH_H

#include <iostream>
#include <sstream>
#include <iterator>
#include <list>
#include "CommonDefs.h"
#include "DirectedGraph.h"

struct MergeNode
{
	LinearNumKey id;
	bool isRC;
	
	void flip() { isRC = (isRC) ? 0 : 1; }
	friend std::ostream& operator<<(std::ostream& out, const MergeNode& object)
	{
		out << object.id << "," << object.isRC;
		return out;
	} 
  
	friend std::istream& operator>>(std::istream& in, MergeNode& object)
	{
		// Read 1 record from the stream
		std::string record;
		in >> record;
		
		// parse the record
		std::stringstream recss(record);
		std::stringstream convertor;
		std::string data;
	
		getline(recss, data, ',');
		convertor.str(data);
		convertor >> object.id;
	
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.isRC;
		return in;
	}  	
};

class ContigPath
{
	public:
		ContigPath();
		
		// Add a node to this path
		void appendNode(const MergeNode& mn);
		
		// prepend a path
		void prependPath(const ContigPath& other);
		
		// append a path
		void appendPath(const ContigPath& other);
		
		// Get the number of nodes in the path
		size_t getNumNodes() const { return m_path.size(); }
		
		// Get the node with a specified index
		MergeNode& getNode(size_t idx)
		{
			assert(idx < m_path.size());
			return m_path[idx];
		}
		
		// Get the node with a specified index
		const MergeNode& getNode(size_t idx) const
		{
			assert(idx < m_path.size());
			return m_path[idx];
		}
		
		// Returns the index of the first node that matches id
		size_t findFirstOf(LinearNumKey id);
		
		// Returns the index of the last node that matches id
		size_t findLastOf(LinearNumKey id);
		
		// reverse the path
		void reverse(bool flipNodes);
		
		// Extract a subset from the path
		ContigPath extractNodes(size_t start, size_t end);
	
		// Insertion/Extraction operators
		friend std::ostream& operator<<(std::ostream& out, const ContigPath& object);
		friend std::istream& operator>>(std::istream& in, ContigPath& object);

	private:
		std::vector<MergeNode> m_path;
};

#endif
