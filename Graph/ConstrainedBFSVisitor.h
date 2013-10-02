#ifndef CONSTRAINED_BFS_VISITOR_H
#define CONSTRAINED_BFS_VISITOR_H

#include "Common/UnorderedMap.h"
#include "Common/IOUtil.h"
#include "Graph/DefaultColorMap.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <iostream>
#include <sstream>
#include <vector>

enum PathSearchResult { UNIQUE_PATH, MULTIPLE_PATHS, NO_PATH };

template <typename G>
class ConstrainedBFSVisitor : public boost::default_bfs_visitor
{
public:

	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;
	typedef std::vector<V> Path;
	typedef std::vector<Path> PathList;
	typedef unsigned short depth_t;

private:

	typedef std::vector<V> Predecessors;
	typedef unordered_map<V, Predecessors> PredecessorMap;
	typedef unordered_map<V, depth_t> DepthMap;

	static const int NO_MAX = -1;

	PredecessorMap m_predecessors;
	DepthMap m_depthMap;
	const V& m_start;
	const V& m_goal;
	depth_t m_minDepth;
	depth_t m_maxDepth;
	DefaultColorMap<G>& m_colorMap;
	bool m_bFoundGoal;

	depth_t m_maxDepthVisited;

public:

  ConstrainedBFSVisitor(
		  const V& start,
		  const V& goal,
		  depth_t minDepth,
		  depth_t maxDepth,
		  DefaultColorMap<G>& colorMap) :
			  m_start(start),
			  m_goal(goal),
			  m_minDepth(minDepth),
			  m_maxDepth(maxDepth),
			  m_colorMap(colorMap),
			  m_bFoundGoal(false)
  {
	  m_maxDepthVisited = 0;
  }

  void examine_edge(const E& e, const G& g)
  {
	  V u = source(e, g);
	  V v = target(e, g);

	  m_predecessors[v].push_back(u);

	  if (get(m_colorMap, v) == boost::white_color) // tree edge
		  m_depthMap[v] = m_depthMap[u] + 1;

	  if (m_depthMap[v] >= m_maxDepth) // limit depth of traversal
		  put(m_colorMap, v, boost::black_color);

	  if (m_depthMap[v] > m_maxDepthVisited)
	  	  m_maxDepthVisited = m_depthMap[v];
  }

  static std::string pathToString(const Path& path)
  {
		std::stringstream s;
		for (unsigned int i = 0; i < path.size(); i++) {
			if (i > 0)
				s << ",";
			s << path[i];
		}
		return s.str();
  }

  PathSearchResult uniquePathToGoal(Path& uniquePath)
  {
//	  std::cout << "num nodes in traversal: " << m_predecessors.size() << std::endl;
//	  std::cout << "minDepth: " << m_minDepth << std::endl;
//	  std::cout << "maxDepth: " << m_maxDepth << std::endl;

	  Path reversePath;
	  reversePath.push_back(m_goal);
	  PathList solutions;
	  bool exceededMaxPaths = false;
	  unsigned long fullDepthPaths = 0;
	  unsigned long deadEndPaths = 0;
	  pathsToGoal(reversePath, solutions, 1, exceededMaxPaths, fullDepthPaths, deadEndPaths);
//	  std::cout << "num full depth paths: " << fullDepthPaths << std::endl;
//	  std::cout << "num dead end paths: " << deadEndPaths << std::endl;

	  if (exceededMaxPaths) {
		  uniquePath.clear();
//		  std::cout << "result: MULTIPLE_PATHS" << std::endl;
		  return MULTIPLE_PATHS;
	  } else if (solutions.size() == 0) {
		  uniquePath.clear();
//		  std::cout << "result: NO_PATH" << std::endl;
		  return NO_PATH;
	  } else {
		  uniquePath = solutions[0];
//		  std::cout << "result: UNIQUE_PATH" << std::endl;
		  return UNIQUE_PATH;
	  }
  }

  /*
  PathList pathsToGoal()
  {
	  Path reversePath;
	  reversePath.push_back(m_goal);
	  PathList solutions;
	  bool exceededMaxPaths = false;
	  pathsToGoal(reversePath, solutions, NO_MAX, exceededMaxPaths);

	  return solutions;
  }
  */

  depth_t getMaxDepthVisited()
  {
	  return m_maxDepthVisited;
  }

private:

  void pathsToGoal(Path& pathToStart, PathList& pathList,
		  int maxPaths, bool& exceededMaxPaths, unsigned long& fullDepthPaths, unsigned long& deadEndPaths)
  {
	  if (pathToStart.size() > (depth_t)(m_maxDepth + 1)) {
		  fullDepthPaths++;
		  return;
	  }

	  if (exceededMaxPaths)
		  return;

	  V back = pathToStart.back();

	  if (back == m_start) {
		  if (pathToStart.size() > m_minDepth) {
			  if (maxPaths == NO_MAX || (int)(pathList.size()) < maxPaths) {
				  pathList.push_back(reversePath(pathToStart));
			  } else {
				  exceededMaxPaths = true;
			  }
		  }
	  }

	  if (back != m_start && m_predecessors[back].size() == 0) {
		  deadEndPaths++;
	  }

	  for (unsigned int i = 0; i < m_predecessors[back].size(); i++) {
		  pathToStart.push_back(m_predecessors[back][i]);
		  pathsToGoal(pathToStart, pathList, maxPaths, exceededMaxPaths, fullDepthPaths, deadEndPaths);
		  pathToStart.pop_back();
	  }
  }


  Path reversePath(Path path)
  {
	  Path reversed;
	  for (int i = path.size() - 1; i >= 0; i--)
		  reversed.push_back(path[i]);
	  return reversed;
  }

};


#endif /* CONSTRAINED_BFS_VISITOR_H */
