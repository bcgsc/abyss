#include "SequenceTree.h"

#include <vector>
#include <queue>
#include <memory>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <functional>

class SequenceTreeNode: public ContigNode {

public:

  SequenceTreeNode(const ContigNode& contigNode,
                   const int overlap, const int maxLength, const bool forward);

  std::vector<SequenceTreeNode> getChildren(const int maxChildren = std::numeric_limits<int>::max()) const;
  const Sequence treeSequence() const;

  int treeSequenceStart;
  int treeSequenceLength;
  int maxLength;
  bool forward;

};

const int EXPECTED_BASES_PER_NODE = 4;

SequenceTreeNode::SequenceTreeNode(const ContigNode& contigNode,
                  const int overlap, const int maxLength, const bool forward):
  ContigNode(contigNode), maxLength(maxLength), forward(forward)
{
  assert(overlap >= 0);
  assert(maxLength > 0);
  treeSequenceStart = overlap;
  const auto &sequence = getContigSequence(*this);
  const int size = sequence.size();
  int treeSequenceEnd = std::min(overlap + maxLength, size);
  assert(treeSequenceStart >= 0);
  assert(treeSequenceEnd > 0);
  assert(treeSequenceEnd > treeSequenceStart);
  treeSequenceLength = treeSequenceEnd - treeSequenceStart;
  if (!forward) {
    treeSequenceStart = size - treeSequenceEnd;
    assert(treeSequenceStart >= 0);
  }
  assert(treeSequenceStart + treeSequenceLength <= size);
}

std::vector<SequenceTreeNode> SequenceTreeNode::getChildren(const int maxChildren) const {
  std::vector<SequenceTreeNode> children;
  if (maxLength > treeSequenceLength) {
    if (forward) {
      int step = out_degree(*this, g_contigGraph) / maxChildren;
      if (step <= 0) { step = 1; }
      out_edge_iterator it, itLast;
      for (std::tie(it, itLast) = out_edges(*this, g_contigGraph);
        it != itLast; std::advance(it, step))
      {
        ContigNode outig = target(*it, g_contigGraph);

        children.push_back(SequenceTreeNode(
                            outig, -distanceBetween(*this, outig),
                            maxLength - treeSequenceLength, forward));
        if (int(children.size()) >= maxChildren) {
          break;
        }
      }
    } else {
      int step = in_degree(*this, g_contigGraph) / maxChildren;
      if (step <= 0) { step = 1; }
      in_edge_iterator it, itLast;
      for (std::tie(it, itLast) = in_edges(*this, g_contigGraph);
        it != itLast; std::advance(it, step))
      {
        ContigNode intig = source(*it, g_contigGraph);

        children.push_back(SequenceTreeNode(
                            intig, -distanceBetween(intig, *this),
                            maxLength - treeSequenceLength, forward));
        if (int(children.size()) >= maxChildren) {
          break;
        }
      }
    }
  }
  return children;
}

const Sequence SequenceTreeNode::treeSequence() const {
  assert(treeSequenceStart >= 0);
  assert(treeSequenceLength > 0);
  return Sequence(getContigSequence(*this).c_str(), treeSequenceStart, treeSequenceLength);
}

typedef std::list<SequenceTreeNode> Trace;
typedef std::list<Trace> Traces;

static Traces
getTreeTraces(const ContigNode& start,
              const int overlap, const int maxLength,
              const bool forward, const int maxPaths)
{
  assert(maxLength > 0);

  Traces traces;
  std::queue<std::reference_wrapper<Trace>> queue;

  SequenceTreeNode root(start, overlap, maxLength, forward);

  int level = 1;
  int leaves = 1;
  traces.push_back({ root });
  queue.push(traces.back());
  while (!queue.empty()) {
    Trace &trace = queue.front();
    queue.pop();

    const SequenceTreeNode &node = trace.back();
    assert(node.maxLength > 0 && node.maxLength <= maxLength);

    assert(int(trace.size()) >= level);
    level = trace.size();

    auto children = node.getChildren();
    if ((children.size() > 0) && (leaves + int(children.size()) - 1 <= maxPaths)) {
      for (size_t i = 0; i < children.size(); i++) {
        const auto &child = children[i];
        assert(child.maxLength > 0);
        assert(child.maxLength < node.maxLength);
        if (i < children.size() - 1) {
          Trace childTrace = trace;
          childTrace.push_back(child);
          assert(childTrace.size() == trace.size() + 1);
          traces.push_back(childTrace);
          queue.push(traces.back());
        } else {
          trace.push_back(child);
          queue.push(trace);
        }
      }
      leaves += children.size() - 1;
    }
  }

  return traces;
}

std::list<Sequence>
getTreeSequences(const ContigNode& start,
                 const int overlap, const int maxLength,
                 const bool forward, const int maxPaths)
{
  assert(overlap >= 0);
  assert(maxLength > 0);
  assert(maxPaths >= 1);

  auto traces = getTreeTraces(start, overlap, maxLength,
                              forward, maxPaths);

  std::list<Sequence> sequences;

  for (const auto &trace : traces) {
    Sequence sequence;
    sequence.reserve(trace.size() * EXPECTED_BASES_PER_NODE);
    if (forward) {
      for (auto it = trace.begin(); it != trace.end(); it++) {
        sequence += it->treeSequence();
      }
    } else {
      for (auto it = trace.rbegin(); it != trace.rend(); it++) {
        sequence += it->treeSequence();
      }
    }
    sequences.push_back(sequence);
  }

  return sequences;
}
