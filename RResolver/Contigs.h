#ifndef RRESOLVER_CONTIGS_H
#define RRESOLVER_CONTIGS_H 1

#include "Common/ConstString.h"
#include "Common/ContigNode.h"
#include "Common/ContigPath.h"
#include "Common/ContigProperties.h"
#include "Common/Histogram.h"
#include "Common/IOUtil.h"
#include "Common/Options.h"
#include "Common/Sequence.h"
#include "Common/StringUtil.h"
#include "DataLayer/FastaReader.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"

#include <set>
#include <string>
#include <vector>

typedef ContigGraph<DirectedGraph<ContigProperties, Distance>> Graph;
typedef Graph::vertex_descriptor vertex_descriptor;
typedef Graph::edge_descriptor edge_descriptor;
typedef Graph::adjacency_iterator adjacency_iterator;
typedef Graph::in_edge_iterator in_edge_iterator;
typedef Graph::out_edge_iterator out_edge_iterator;

class ContigSequence : public const_string
{

  public:
	ContigSequence(const std::string& s)
	  : const_string(s)
	{
		strSize = const_string::size();
	}

	size_t size() const { return strSize; }

  private:
	size_t strSize;
};

typedef std::vector<std::pair<std::string, std::string>> ContigSequencesInfo;
typedef std::vector<ContigSequence> ContigSequences;
typedef std::vector<ContigNode> ContigPath;
typedef std::vector<ContigPath> ContigPaths;
typedef std::vector<std::pair<ContigNode, int>>
    ImaginaryContigPath; // Each element is a pair of contig and distance to the next
typedef std::set<ImaginaryContigPath> ImaginaryContigPaths;

extern Graph g_contigGraph;
extern std::vector<std::string> g_contigComments;
extern ContigSequences g_contigSequences;

int
distanceBetween(const ContigNode& node1, const ContigNode& node2);
const ContigSequence&
getContigSequence(const ContigNode& node);
int
getContigSize(const ContigNode& node);
const std::string&
getContigComment(const ContigNode& node);
Sequence
getPathSequence(const ContigPath& path);
Sequence
getPathSequence(const ImaginaryContigPath& path);
long
getContigNumOfKmers(const ContigNode& node);
double
getContigBaseCoverage(const ContigNode& node);

void
loadContigGraph(const std::string& contigGraphPath);
void
loadContigs(const std::string& contigsPath);
void
storeContigGraph(
    const std::string& contigGraphPath,
    const std::string& program,
    const std::string& commandLine);
void
storeContigs(const std::string& contigsPath);

unsigned
num_vertices_removed(const Graph& contigGraph);

void
assembleContigs();

#endif
