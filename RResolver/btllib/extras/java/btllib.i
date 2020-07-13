%module btllib

%{
#include "btllib/index_queue.hpp"
#include "btllib/seq_reader.hpp"
#include "btllib/nthash.hpp"
#include "btllib/data_stream.hpp"
#include "btllib/graph.hpp"
#include "btllib/util.hpp"
#include "btllib/status.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/seq.hpp"
#include "btllib/bloom_filter.hpp"
%}

%include <java.swg>
%include <various.i>
%include <std_string.i>

%include "extra.i"

%include "btllib/index_queue.hpp"
%include "btllib/seq_reader.hpp"
%include "btllib/nthash.hpp"
%include "btllib/data_stream.hpp"
%include "btllib/graph.hpp"
%include "btllib/util.hpp"
%include "btllib/status.hpp"
%include "btllib/seq_writer.hpp"
%include "btllib/seq.hpp"
%include "btllib/bloom_filter.hpp"
