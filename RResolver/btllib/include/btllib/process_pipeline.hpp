#ifndef BTLLIB_PROCESS_PIPELINE_HPP
#define BTLLIB_PROCESS_PIPELINE_HPP

#include <atomic>
#include <cstdio>
#include <string>

namespace btllib {

using PipeId = unsigned long;
using PipelineId = unsigned long;

/**
 * Run a process pipeline and obtain the stdin of the first
 * and stdout of the last process.
 */
class ProcessPipeline
{

public:
  /**
   * Runs a process or a pipeline of processes in the background.
   *
   * @param cmd The command to execute and obtain stdin and stdout for.
   * A number of commands can be chained and piped (with the |
   * operator) in which case stdin of the first command stdout
   * of the last command are available.
   */
  ProcessPipeline(const std::string& cmd);
  ~ProcessPipeline() { end(); }

  void close_in();
  void close_out();

  void end();

  FILE *in = nullptr, *out = nullptr;
  std::atomic<bool> in_closed{ false };
  std::atomic<bool> out_closed{ false };

private:
  PipelineId id = 0;
  std::atomic<bool> ended{ false };
};

} // namespace btllib

#endif
