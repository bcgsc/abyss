#ifndef _ASSEMBLY_COUNTERS_H_
#define _ASSEMBLY_COUNTERS_H_

#include "Common/IOUtil.h"
#include <iostream>
#include <cassert>
#include <string>

namespace BloomDBG {

	/**
	 * Counters for tracking assembly statistics and producing
	 * progress messages.
	 */
	struct AssemblyCounters
	{
		/** reads consisting entirely of solid k-mers */
		size_t solidReads;
		size_t allKmersVisited;
		size_t allBranchKmersVisited;
		size_t nonBranchingReads;
		size_t nonBranchingContigs;
		size_t readsProcessed;
		size_t basesAssembled;
		size_t contigID;

		AssemblyCounters() : solidReads(0), allKmersVisited(0),
			allBranchKmersVisited(0), nonBranchingReads(0),
			nonBranchingContigs(0), readsProcessed(0),
			basesAssembled(0), contigID(0) {}

		/** serialize counters as a TSV table */
		friend std::ostream& operator<<(std::ostream& out,
			const AssemblyCounters& o)
		{
			/* write headers */

			out << "solid_reads"
				<< '\t' << "processed_reads"
				<< '\t' << "bases_assembled"
				<< '\t' << "next_contig_id"
				<< '\n';

			/* write data */

			out << o.solidReads
				<< '\t' << o.readsProcessed
				<< '\t' << o.basesAssembled
				<< '\t' << o.contigID
				<< '\n';

			return out;
		}

		/** deserialize counters from a TSV table  */
		friend std::istream& operator>>(std::istream& in,
			AssemblyCounters& o)
		{
			/* ignore header line */

			std::string headers;
			getline(in, headers);

			/* read data */

			in >> o.solidReads
				>> expect("\t") >> o.readsProcessed
				>> expect("\t") >> o.basesAssembled
				>> expect("\t") >> o.contigID
				>> expect("\n");

			return in;
		}
	};

} /* end of BloomDBG namespace */

#endif
