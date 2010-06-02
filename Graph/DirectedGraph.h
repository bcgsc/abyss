#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "ContigNode.h"
#include "ContigPath.h"
#include "Sense.h"
#include <ostream>
#include <map>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

typedef std::vector<ContigPath> ContigPaths;

typedef uint32_t LinearNumKey;

template<typename K, typename D>
struct Vertex
{
	typedef Vertex<K,D> VertexType;

	Vertex(const K& k, const D& d) : m_key(k), m_data(d) { }

	struct EdgeData
	{
		VertexType* pVertex;
		bool reverse;

		bool operator==(const EdgeData& e2) const
		{
			return (this->pVertex == e2.pVertex && this->reverse == e2.reverse);
		}
	};

	// Compare operators (needed?)
	inline bool operator==(const Vertex& other) const { return m_data == other.m_data; }
	inline bool operator<(const Vertex& other) const { return m_data < other.m_data; }

	/** Return the number of edges in the specified direction. */
	size_t numEdges(bool sense) const
	{
		return m_edges[sense].size();
	}

	// add an edge to the vertex in the specified direction
	void addEdge(VertexType* pNode, extDirection dir, bool reverse);

	K m_key;
	D m_data;

	typedef typename std::vector<EdgeData> EdgeCollection;
	EdgeCollection m_edges[NUM_DIRECTIONS];
};

typedef std::pair<ContigNode, unsigned> Constraint;
typedef std::vector<Constraint> Constraints;

template<typename D>
class DirectedGraph
{
	public:
		typedef typename std::vector<Vertex<LinearNumKey, D> > VertexTable;
		typedef Vertex<LinearNumKey, D> VertexType;

		/** Return the vertex specified by the given key. */
		const VertexType& operator[](const LinearNumKey& key) const
		{
			const VertexType* pVertex = findVertex(key);
			assert(pVertex != NULL);
			return *pVertex;
		}

		/** Return the vertex specified by the given key. */
		VertexType& operator[](const LinearNumKey& key)
		{
			VertexType* pVertex = findVertex(key);
			assert(pVertex != NULL);
			return *pVertex;
		}

		/** Return the vertex data specified by the given key. */
		const D& getDataForVertex(const LinearNumKey& key) const
		{
			return (*this)[key].m_data;
		}

		void addEdge(const LinearNumKey& parent, extDirection dir,
				const ContigNode& child);
		void addVertex(const LinearNumKey& key, const D& data);

		bool findSuperpaths(const LinearNumKey& sourceKey,
				extDirection dir, Constraints& constraints,
				ContigPaths& superPaths, unsigned& compCost) const;

		// return the number of vertices
		size_t getNumVertices() const { return m_vertexTable.size(); }

		// count the number of edges (SLOW)
		size_t countEdges() const;

		// Calculate the length of this path
		size_t calculatePathLength(const ContigPath& path) const;

		// Make a map of the distances to each node
		void makeDistanceMap(const ContigPath& path,
				std::map<ContigNode, int>& distanceMap) const;

	private:
		bool ConstrainedDFS(const VertexType* pCurrVertex,
				extDirection dir, Constraints& constraints,
				Constraints::const_iterator nextConstraint,
				unsigned satisfied,
				ContigPath& path, ContigPaths& solutions,
				size_t currLen, unsigned& visitedCount) const;

		VertexType* findVertex(const LinearNumKey& key);
		const VertexType* findVertex(const LinearNumKey& key) const;
		VertexTable m_vertexTable;
};

#endif
