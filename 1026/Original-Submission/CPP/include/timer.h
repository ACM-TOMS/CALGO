#ifndef CALS_TIMER_H
#define CALS_TIMER_H

#include <chrono>
#include <string>

namespace cals
{
  class Timer
  {
    std::chrono::high_resolution_clock::time_point t0{};
    double t{std::numeric_limits<double>::max()};

  public:
    Timer() = default;

    ~Timer() = default;

    inline void start()
    { t0 = std::chrono::high_resolution_clock::now(); };

    inline void stop()
    {
      auto t1 = std::chrono::high_resolution_clock::now();
      t = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / 1e9;
    };

    inline void reset()
    { t = 0.0; };

    inline double get_time() const
    { return (t == std::numeric_limits<double>::max()) ? 0.0 : t; };
  };

  struct MttkrpTimers
  {
    enum TIMERS
    {
      MT_KRP = 0,
      MT_GEMM,
      TS_GEMM,
      TS_GEMV,
      LENGTH
    };
    std::string names[TIMERS::LENGTH] =
        {
            "MT_KRP",
            "MT_GEMM",
            "TS_GEMM",
            "TS_GEMV"
        };
    Timer timers[TIMERS::LENGTH];

    Timer &operator[](int timer)
    { return timers[timer]; }
  };

  struct ModeTimers
  {
    enum TIMERS
    {
      MTTKRP = 0,
      UPDATE,
      LENGTH
    };
    std::string names[TIMERS::LENGTH] =
        {
            "TOTAL_MTTKRP",
            "UPDATE"
        };
    Timer timers[TIMERS::LENGTH];

    Timer &operator[](int timer)
    { return timers[timer]; }
  };

  struct AlsTimers
  {
    enum TIMERS
    {
      ITERATION = 0,
      DEFRAGMENTATION,
      ERROR,
      G_COPY,
      LENGTH
    };
    std::string names[TIMERS::LENGTH] =
        {
            "ITERATION",
            "DEFRAGMENTATION",
            "ERROR",
            "G_COPY"
        };

    Timer timers[TIMERS::LENGTH];

    Timer &operator[](int timer)
    { return timers[timer]; }
  };

} // namespace cals

#endif // CALS_TIMER_H
