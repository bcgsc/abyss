%module btllib

%{
#define SWIG_FILE_WITH_INIT

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

%include <pyprimtypes.swg>
%include <pyopers.swg>
%include <std_common.i>
%include <cstring.i>
%include <std_string.i>
%include <exception.i>
%include <std_string.i>
%include <std_iostream.i>
%include <std_vector.i>
%include <carrays.i>

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
