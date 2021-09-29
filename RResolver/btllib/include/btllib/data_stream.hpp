#ifndef BTLLIB_DATA_STREAM_HPP
#define BTLLIB_DATA_STREAM_HPP

#include "process_pipeline.hpp"
#include "status.hpp"
#include "util.hpp"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <csignal>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
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

// clang-format off
static const Datatype DATATYPES[] {
  { { "http://", "https://", "ftp://" }, {}, { "command -v wget" }, { "wget -O-" }, { "" }, { "" } }, { {}, { ".url" }, { "command -v wget" }, { "wget -O- -i" }, { "" }, { "" } }, { {}, { ".ar" }, { "command -v ar" }, { "ar -p" }, { "" }, { "" } }, { {}, { ".tar" }, { "command -v tar" }, { "tar -xOf" }, { "" }, { "" } }, { {}, { ".tgz" }, { "command -v tar" }, { "tar -zxOf" }, { "" }, { "" } }, { {}, { ".gz", ".z" }, { "command -v pigz",
  "command -v gzip" }, { "pigz -dc", "gzip -dc" }, { "pigz >", "gzip >" }, {
  "pigz >>", "gzip >>" } }, { {}, { ".bz2" }, { "command -v bzip2" }, {
  "bunzip2 -dc" }, { "bzip2 >" }, { "bzip2 >>" } }, { {}, { ".xz" }, {
  "command -v xz" }, { "unxz -dc" }, { "xz -T0 >" }, { "xz -T0 >>" } }, { {},
  { ".7z" }, { "command -v 7z" }, { "7z -so e" }, { "7z -si a" }, { "7z -si a" } }, { {}, { ".zip" }, { "command -v zip" }, { "unzip -p" }, { "" }, {
  "" } }, { {}, { ".lrz" }, { "command -v lrzip" }, { "lrzip -q -d -o -" }, {
  "lrzip -q >" }, { "" } }, { {}, { ".bam", ".cram" }, { "command -v samtools" }, { "samtools view -h" }, { "samtools -Sb - >" }, { "samtools -Sb - >>" } },
};
// clang-format on

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
  bool closed = false;
  std::unique_ptr<ProcessPipeline> pipeline;
};

class DataSource : public DataStream
{

public:
  DataSource(const std::string& path)
    : DataStream(path, READ)
  {}
};

class DataSink : public DataStream
{

public:
  DataSink(const std::string& path, bool append = false)
    : DataStream(path, append ? APPEND : WRITE)
  {}
};

static inline std::string
get_pipeline_cmd(const std::string& path, DataStream::Operation op);

inline DataStream::DataStream(const std::string& path, Operation op)
  : streampath(path)
  , op(op)
{
  if (path == "-") {
    if (op == READ) {
      file = stdin;
    } else {
      file = stdout;
    }
  } else {
    pipeline = std::unique_ptr<ProcessPipeline>(
      new ProcessPipeline(get_pipeline_cmd(path, op)));
    if (op == READ) {
      file = pipeline->out;
    } else {
      file = pipeline->in;
    }
  }
}

inline void
DataStream::close()
{
  if (!closed) {
    if (streampath != "-") {
      pipeline->end();
    }
    closed = true;
  }
}

static inline std::string
get_datatype_cmd(const std::string& path,
                 const Datatype& datatype,
                 DataStream::Operation op)
{
  bool found_cmd = false;
  int cmd_idx = 0;
  for (const auto& existence_cmd : datatype.cmds_check_existence) {
    const pid_t pid = fork();
    if (pid == 0) {
      int null_fd = open("/dev/null", O_WRONLY, 0);
      dup2(null_fd, STDOUT_FILENO);
      dup2(null_fd, STDERR_FILENO);
      close(null_fd);

      execlp("sh", "sh", "-c", existence_cmd.c_str(), NULL);
      log_error("exec failed: sh -c \"" + existence_cmd + "\'");
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    } else {
      check_error(pid == -1, "Error on fork.");
      int status;
      check_error(waitpid(pid, &status, 0) != pid,
                  "waitpid error: " + get_strerror());
      if (!(WIFSIGNALED(status)) &&
          ((WIFEXITED(status)) && (WEXITSTATUS(status) == 0))) { // NOLINT
        found_cmd = true;
        break;
      }
    }
    cmd_idx++;
  }

  check_error(!found_cmd,
              "Filetype recognized for '" + path +
                "', but no tool available to work with it.");

  std::string cmd;
  switch (op) {
    case DataStream::Operation::READ:
      cmd = datatype.read_cmds[cmd_idx];
      break;
    case DataStream::Operation::WRITE:
      cmd = datatype.write_cmds[cmd_idx];
      break;
    case DataStream::Operation::APPEND:
      cmd = datatype.append_cmds[cmd_idx];
      break;
    default:
      log_error("Invalid operation");
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
  check_error(cmd.empty(),
              "Filetype recognized for '" + path +
                "', but no tool available to work with it.");
  return cmd;
}

static inline std::vector<std::string>
peel_datatype(const std::string& path, DataStream::Operation op)
{
  std::string default_cmd = "cat";
  if (op == DataStream::Operation::WRITE) {
    default_cmd += " >";
  } else if (op == DataStream::Operation::APPEND) {
    default_cmd += " >>";
  }

  std::string path_trimmed = path;
  std::vector<std::string> cmd_layers;
  for (;;) {
    bool found_datatype = false;
    for (const auto& datatype : DATATYPES) {
      size_t trim_start = 0, trim_end = 0;
      bool this_datatype = false;
      for (const auto& prefix : datatype.prefixes) {
        if (startswith(path_trimmed, prefix)) {
          this_datatype = true;
          trim_start += prefix.size();
          break;
        }
      }
      for (const auto& suffix : datatype.suffixes) {
        if (endswith(path_trimmed, suffix)) {
          this_datatype = true;
          trim_end += suffix.size();
          break;
        }
      }

      if (this_datatype) {
        found_datatype = true;
        cmd_layers.push_back(get_datatype_cmd(path_trimmed, datatype, op));
        path_trimmed.erase(0, trim_start);
        path_trimmed.erase(path_trimmed.size() - trim_end);
      }
    }
    if (!found_datatype) {
      break;
    }
  }
  if (cmd_layers.empty()) {
    cmd_layers.push_back(default_cmd);
  }
  if (op == DataStream::Operation::WRITE ||
      op == DataStream::Operation::APPEND) {
    std::reverse(cmd_layers.begin(), cmd_layers.end());
  }

  return cmd_layers;
}

static inline std::string
form_string_cmd(std::vector<std::string>& cmd_layers,
                DataStream::Operation op,
                const std::string& path)
{
  std::string result_cmd;
  for (size_t i = 0; i < cmd_layers.size(); i++) {
    auto& cmd = cmd_layers[i];
    if (op == DataStream::Operation::WRITE ||
        op == DataStream::Operation::APPEND) {
      if (i == cmd_layers.size() - 1) {
        if (cmd.back() == '>') {
          if (path == "-") {
            while (cmd.back() == '>' || cmd.back() == ' ') {
              cmd.pop_back();
            }
          } else {
            cmd += path;
          }
        } else {
          cmd += " ";
          cmd += path;
        }
      } else {
        if (cmd.back() == '>') {
          while (cmd.back() == '>' || cmd.back() == ' ') {
            cmd.pop_back();
          }
        } else {
          cmd += " -";
        }
      }
    } else {
      if (i == 0) {
        cmd += " ";
        cmd += path;
      } else {
        cmd += " -";
      }
    }
    if (i > 0) {
      result_cmd += " | ";
    }
    result_cmd += cmd;
  }

  check_error(result_cmd.empty(),
              (op == DataStream::Operation::READ ? "Error loading from "
                                                 : "Error saving to ") +
                path);

  check_error(result_cmd == "cat" || result_cmd == "cat -",
              "Attempting to create a pipeline on stdin or stdout which is a "
              "redundant operation.");

  return result_cmd;
}

static inline std::string
get_pipeline_cmd(const std::string& path, DataStream::Operation op)
{
  auto cmd_layers = peel_datatype(path, op);
  return form_string_cmd(cmd_layers, op, path);
}

} // namespace btllib

#endif
