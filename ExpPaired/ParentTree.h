#ifndef PARENT_TREE
#define PARENT_TREE

#include <map>
#include "CommonDefs.h"
#include "AssemblyAlgorithms.h"
#include "NodeScores.h"

typedef std::pair<PackedSeq, PackedSeq> PSPair;
typedef std::map<PackedSeq, PSeqSet> RelationshipTable;

class ParentTree
{
	public:
		ParentTree(const PackedSeq& rootNode, extDirection dir, ISequenceCollection* pSC, int maxDepth);
		
		void addScore(const PackedSeq& node, double weight);
		TreeScore getScore(const PackedSeq& node);
		void print(int maxDepth);
		void printRootChildren();
		PSeqSet getRootChildren();
		
	private:
	
		void generateUpdateSet(const PackedSeq& startNode, PSeqSet& updateSet);
		
		void addToTable(RelationshipTable& table, const PackedSeq& key, const PackedSeq& value);
		void explore(const PackedSeq& parent, int depth);
		void printRecursive(const PackedSeq& node, int depth, int maxDepth);
	
		PackedSeq m_root;
		ISequenceCollection* m_pSC;
		extDirection m_dir;
		int m_maxDepth;
		
		// TODO: Use a proper tree
		RelationshipTable m_childToParent;
		RelationshipTable m_parentToChild;
		NodeScores m_nodeScores;
};

#endif
