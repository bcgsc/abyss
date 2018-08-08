#ifndef _ASSEMBLY_STREAMS_H_
#define _ASSEMBLY_STREAMS_H_

namespace BloomDBG {

	/** Bundles together input and output streams used during assembly */
	template <typename InputStreamT>
	struct AssemblyStreams
	{
		/** input reads stream */
		InputStreamT&  in;
		/** main FASTA output */
		std::ostream& out;
		/** duplicated FASTA output for checkpointing */
		std::ostream& checkpointOut;
		/** trace file output for debugging */
		std::ostream& traceOut;
		/** outcomes of processing each read */
		std::ostream& readLogOut;

		AssemblyStreams(InputStreamT& in, std::ostream& out,
			std::ostream& checkpointOut, std::ostream& traceOut,
			std::ostream& readLogOut) :
			in(in), out(out), checkpointOut(checkpointOut),
			traceOut(traceOut), readLogOut(readLogOut) {}
	};

}

#endif
