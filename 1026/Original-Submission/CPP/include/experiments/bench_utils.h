#ifndef CP_CALS_BENCH_UTILS_H
#define CP_CALS_BENCH_UTILS_H

namespace cals::bench
{
  const static int ITERATIONS = 5;

  struct BenchResult
  {
    double time;
    double first_time;
    double second_time;
    u_int64_t flops;
    u_int64_t memops;
  };

  void initialize_csv(std::string &file_name, const std::string &sep = ";")
  {
    auto file = std::ofstream(file_name, std::ios::out);
    file << "MODE" << sep
         << "TENSOR_MODES" << sep
         << "FLOPS" << sep
         << "MEMOPS" << sep
         << "RANK" << sep
         << "THREADS" << sep
         << "TIME" << sep
         << "F_TIME" << sep
         << "S_TIME" << std::endl;
  }

  void append_csv(std::string &file_name, int mode, std::string &modes_string, int rank, int threads,
                  BenchResult &br, const std::string &sep = ";")
  {
    auto file = std::ofstream(file_name, std::ios::app);
    file << mode << sep
         << modes_string << sep
         << br.flops << sep
         << br.memops << sep
         << rank << sep
         << threads << sep
         << br.time << sep
         << br.first_time << sep
         << br.second_time << std::endl;
  }

  enum EXTERNAL_MTTKRP_METHODS
  {
    CTF = cals::mttkrp::MTTKRP_METHOD::AUTO + 1,
    PLANC
  };
}

#endif //CP_CALS_BENCH_UTILS_H
