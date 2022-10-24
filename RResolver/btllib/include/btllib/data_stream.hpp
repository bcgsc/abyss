#ifndef BTLLIB_DATA_STREAM_HPP
#define BTLLIB_DATA_STREAM_HPP

#include "btllib/process_pipeline.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <memory>
#include <string>
#include <vector>

namespace btllib {

struct Datatype
{
  std::vector<std::string> prefixes;
  std::vector<std::string> suffixes;
  std::vector<std::string> cmds_check_existence;
  std::vector<std::string> read_cmds;
  std::vector<std::string> write_cmds;
  std::vector<std::string> append_cmds;
};

extern const Datatype DATATYPES[12];

class DataStream
{
public:
  enum Operation
  {
    READ,
    WRITE,
    APPEND,
    CLOSE
  };

  DataStream(const std::string& path, Operation op);
  DataStream(const DataStream&) = delete;
  DataStream(DataStream&&) = delete;

  DataStream& operator=(const DataStream&) = delete;
  DataStream& operator=(DataStream&&) = delete;

  ~DataStream() { close(); }
  void close();

  FILE* operator*() const { return file; }
  FILE* operator->() const { return file; }
  operator FILE*() const { return file; }

protected:
  std::string streampath;
  Operation op;
  FILE* file = nullptr;
  std::atomic<bool> closed{ false };
  std::unique_ptr<ProcessPipeline> pipeline;
};

class DataSource : public DataStream
{

public:
  DataSource(const std::string& path)
    : DataStream(path, READ)
  {
  }
};

class DataSink : public DataStream
{

public:
  DataSink(const std::string& path, bool append = false)
    : DataStream(path, append ? APPEND : WRITE)
  {
  }
};

} // namespace btllib

#endif
