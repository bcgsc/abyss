#include "Contigs.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

Graph g_contigGraph;
std::vector<std::string> g_contigComments;
ContigSequences g_contigSequences;

int distanceBetween(const ContigNode& node1, const ContigNode& node2) {
  return get(edge_bundle, g_contigGraph, node1, node2).distance;
}

const ContigSequence& getContigSequence(const ContigNode &node) {
  const long idx = node.index();
  assert(idx >= 0);
  assert(idx < long(g_contigSequences.size()));
  return g_contigSequences[idx];
}

int getContigSize(const ContigNode& node) {
  return getContigSequence(node).size();
}

const std::string& getContigComment(const ContigNode &node) {
  const long id = node.id();
  assert(id >= 0);
  assert(id < long(g_contigComments.size()));
  return g_contigComments[id];
}

Sequence getPathSequence(const ContigPath &path) {
  assert(path.size() >= 1);
  Sequence sequence = Sequence(getContigSequence(path[0]));
  for (size_t i = 1; i < path.size(); i++) {
    const auto &node = path[i];
    const auto distance = distanceBetween(path[i - 1], path[i]);
    const auto overlap = -distance;
    assert(overlap >= 0);
    Sequence newSequence = Sequence(getContigSequence(node));
    assert(int(sequence.size()) >= overlap);
    assert(int(newSequence.size()) >= overlap);
    assert(sequence.substr(sequence.size() - overlap) == newSequence.substr(0, overlap));
    sequence += newSequence.substr(overlap);
  }
  return sequence;
}

Sequence getPathSequence(const ImaginaryContigPath &path) {
  assert(path.size() >= 1);
  Sequence sequence = Sequence(getContigSequence(path[0].first));
  for (size_t i = 1; i < path.size(); i++) {
    const auto &node = path[i].first;
    const auto &distance = path[i].second;
    const auto overlap = -distance;
    assert(overlap >= 0);
    Sequence newSequence = Sequence(getContigSequence(node));
    assert(int(sequence.size()) >= overlap);
    assert(int(newSequence.size()) >= overlap);
    assert(sequence.substr(sequence.size() - overlap) == newSequence.substr(0, overlap));
    sequence += newSequence.substr(overlap);
  }
  return sequence;
}

long getContigNumOfKmers(const ContigNode& node) {
  return get(vertex_bundle, g_contigGraph, node).coverage;
}

double getContigBaseCoverage(const ContigNode& node) {
  return double(getContigNumOfKmers(node)) * opt::k / double(getContigSequence(node).size() - opt::k + 1);
}

void loadContigGraph(const std::string &contigGraphPath) {
  if (opt::verbose) std::cerr << "Loading contig graph from `" << contigGraphPath << "'...\n";
  std::ifstream fin(contigGraphPath);
  assert_good(fin, contigGraphPath);
  fin >> g_contigGraph;
  assert(fin.eof());
  g_contigNames.lock();
  if (opt::verbose) {
    std::cerr << "Contig graph loaded.\n";
    printGraphStats(std::cerr, g_contigGraph);
  }
}

void loadContigs(const std::string &contigsPath) {
  if (opt::verbose) std::cerr << "Loading contigs from `" << contigsPath << "'...\n";
  FastaReader in(contigsPath.c_str(), FastaReader::NO_FOLD_CASE);
  for (FastaRecord rec; in >> rec;) {
    if (g_contigNames.count(rec.id) == 0) continue;
    assert(g_contigSequences.size() / 2 == get(g_contigNames, rec.id));
    g_contigComments.push_back(rec.comment);
    g_contigSequences.push_back(rec.seq);
    g_contigSequences.push_back(reverseComplement(rec.seq));
  }
  assert(in.eof());
  assert(!g_contigSequences.empty());
  opt::colourSpace = isdigit(g_contigSequences.front()[0]);
  if (opt::verbose) {
    std::cerr << "Contigs loaded.\n";
  }
}

void storeContigGraph(const std::string &contigGraphPath, const std::string &program, const std::string& commandLine) {
    if (opt::verbose) { std::cerr << "Storing contig graph to `" << contigGraphPath << "'...\n"; }
		std::ofstream fout(contigGraphPath.c_str());
		assert_good(fout, contigGraphPath);
		write_graph(fout, g_contigGraph, program, commandLine);
		assert_good(fout, contigGraphPath);
    if (opt::verbose) { std::cerr << "Contig graph stored.\n"; }
}

void storeContigs(const std::string &contigsPath) {
    if (opt::verbose) { std::cerr << "Storing contigs to `" << contigsPath << "'...\n"; }
		std::ofstream fout(contigsPath.c_str());
		assert_good(fout, contigsPath);

    Graph::vertex_iterator vertexStart, vertexEnd;
    boost::tie(vertexStart, vertexEnd) = vertices(g_contigGraph);
    for (auto it = vertexStart; it != vertexEnd; ++it) {
      auto node = *it;
      if (!get(vertex_removed, g_contigGraph, node) && !node.sense()) {
        FastaRecord rec;

        std::stringstream ss;
        ss << get(vertex_name, g_contigGraph, node);
        std::string name = ss.str();
        
        rec.id = name.substr(0, name.size() - 1);
        rec.comment = getContigComment(node);
        rec.seq = getContigSequence(node);

        fout << rec;
      }
    }
		assert_good(fout, contigsPath);
    if (opt::verbose) { std::cerr << "Contigs stored.\n"; }
}

unsigned num_vertices_removed(const Graph& graph) {
  unsigned removed = 0;
  Graph::vertex_iterator vertexStart, vertexEnd;
  boost::tie(vertexStart, vertexEnd) = vertices(graph);
  for (auto it = vertexStart; it != vertexEnd; ++it) {
    auto node = *it;
    if (get(vertex_removed, graph, node)) {
      removed++;
    }
  }
  return removed;
}

void assembleContigs() {
  if (opt::verbose) {
    std::cerr << "Assembling contigs... ";
  }

  std::vector<std::string> pathIDs;
  ContigPaths paths;

  Graph newContigGraph(g_contigGraph);
  assemble(newContigGraph, back_inserter(paths));
  assert(num_vertices(newContigGraph) >= num_vertices(g_contigGraph));

  // Record all the contigs that are in a path.
	std::vector<bool> seen(g_contigGraph.num_vertices() / 2);

  g_contigNames.unlock();
  int counter = 0;
  for (const auto &path : paths) {
    assert(!path.empty());
    ContigNode pathNode(g_contigGraph.num_vertices() / 2 + counter, path[0].sense());
    std::string name = createContigName();
    pathIDs.push_back(name);
    put(vertex_name, newContigGraph, pathNode, name);
    for (const auto &node : path) {
      seen[node.id()] = true;
    }
    counter++;
  }
  g_contigNames.lock();

	for (const auto &path : paths) {
		assert(!path.empty());
    Sequence sequence = getPathSequence(path);
    unsigned coverage = 0;
    for (const auto &node : path) {
      if (!node.ambiguous()) {
        coverage += g_contigGraph[node].coverage;
      }
    }

    std::stringstream ss;
    ss << sequence.size() << ' ' << coverage << ' ';
    ss << get(vertex_name, g_contigGraph, path.front());
    if (path.size() == 1)
      return;
    else if (path.size() == 3)
      ss << ',' << get(vertex_name, g_contigGraph, path[1]);
    else if (path.size() > 3)
      ss << ",...";
    ss << ',' << get(vertex_name, g_contigGraph, path.back());

    g_contigSequences.push_back(sequence);
    g_contigSequences.push_back(reverseComplement(sequence));
    g_contigComments.push_back(ss.str());
	}
  
  g_contigGraph.swap(newContigGraph);

  if (opt::verbose) {
    std::cerr << "Done!\n";
  }
}
