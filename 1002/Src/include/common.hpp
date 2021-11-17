//
// common.hpp
//
// Common types/functions/macros
//

#ifndef COMMON_HPP
#define COMMON_HPP

#include <algorithm>
#include <chrono>
#include <iostream>
#include <sstream>
#include <vector>

#include <omp.h>

// types

# ifdef ENABLE_LONG
typedef signed long int sval_t;
typedef unsigned long int uval_t;
# else
typedef signed int sval_t;
typedef unsigned int uval_t;
# endif

typedef signed long int slval_t;
typedef unsigned long int ulval_t;
typedef unsigned int node_t;

struct Arc {
    node_t tail;
    node_t head;
    uval_t upper;
};

struct ResArc {
    node_t head;
    uval_t rcap;
    sval_t cost;
    ResArc* pair;
};

struct GraphColoring {
    std::vector<node_t> colors;
    node_t num_colors;
};

struct Problem {
    std::vector<Arc> arcs;  // original arcs (not necessarily ordered)
    node_t s, t;            // source/sink nodes
    node_t n;               // number of nodes
    ulval_t m;              // number of arcs
};

// functions

inline std::vector<node_t>
join_vectors(const std::vector<std::vector<node_t>>& vs)
{
    node_t n = 0;

    std::vector<node_t> begs(vs.size());
    for (size_t i = 0; i < vs.size(); ++i) {
        begs[i] = n;
        n += vs[i].size();
    }

    std::vector<node_t> result(n);
    for (size_t i = 0; i < vs.size(); ++i) {
        for (size_t j = 0; j < vs[i].size(); ++j) {
            result[begs[i] + j] = vs[i][j];
        }
    }

    return result;
}

// macros

template <typename T>
inline std::ostream&
operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (const auto& e : v) {
        os << e << ", ";
    }
    os << "]";
    return os;
}

template <typename T>
inline void
_log(T first)
{ std::cerr << first << std::endl; }

template <typename T, typename... Ts>
inline void
_log(T first, Ts... rest)
{
    std::cerr << first;
    _log(rest...);
}

# ifdef ENABLE_LOG
#     define LOG(...)          _log(__FILE__, ":", __LINE__, ": log: ", __VA_ARGS__)
#     define LOGF(format, ...) fprintf(stderr, "%s:%d: log: " format "\n", __FILE__, __LINE__, __VA_ARGS__)
# else
#     define LOG(...)          ((void)0)
#     define LOGF(...)         ((void)0)
# endif

# ifdef ENABLE_INFO
#     define INFO(...)          _log("info: ", __VA_ARGS__)
#     define INFOF(format, ...) fprintf(stderr, "info: " format "\n", __VA_ARGS__)
# else
#     define INFO(...)          ((void)0)
#     define INFOF(...)         ((void)0)
# endif

struct PushRelabelStat {
    ulval_t pushes_1;
    ulval_t pushes_2;
    ulval_t relabels_1;
    ulval_t relabels_2;
    ulval_t discharges_1;
    ulval_t discharges_2;
    ulval_t global_relabels_1;
    ulval_t global_relabels_2;
};

PushRelabelStat pr_stat{0, 0, 0, 0, 0, 0, 0, 0};

inline std::ostream&
operator<<(std::ostream& os, const PushRelabelStat& stat)
{
    os << "c\n"
       << "c # pushes          : " << stat.pushes_1 << " + "
                                   << stat.pushes_2 << " = "
                                   << stat.pushes_1 + stat.pushes_2 << "\n"
       << "c # relabels        : " << stat.relabels_1 << " + "
                                   << stat.relabels_2 << " = "
                                   << stat.relabels_1 + stat.relabels_2 << "\n"
       << "c # discharges      : " << stat.discharges_1 << " + "
                                   << stat.discharges_2 << " = "
                                   << stat.discharges_1 + stat.discharges_2 << "\n"
       << "c # global relabels : " << stat.global_relabels_1 << " + "
                                   << stat.global_relabels_2 << " = "
                                   << stat.global_relabels_1 + stat.global_relabels_2 << std::endl;

    return os;
}

struct ColoredParallelPushRelabelStat {
    ulval_t pushes_1;
    ulval_t pushes_2;
    ulval_t relabels_1;
    ulval_t relabels_2;
    ulval_t discharges_1;
    ulval_t discharges_2;
    ulval_t global_relabels_1;
    ulval_t global_relabels_2;
    ulval_t colors;
    ulval_t color_ticks_1;
    ulval_t color_ticks_2;
};

ColoredParallelPushRelabelStat cppr_stat{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

inline std::ostream&
operator<<(std::ostream& os, const ColoredParallelPushRelabelStat& stat)
{
    os << "c\n"
       << "c # pushes          : " << stat.pushes_1 << " + "
                                   << stat.pushes_2 << " = "
                                   << stat.pushes_1 + stat.pushes_2 << "\n"
       << "c # relabels        : " << stat.relabels_1 << " + "
                                   << stat.relabels_2 << " = "
                                   << stat.relabels_1 + stat.relabels_2 << "\n"
       << "c # discharges      : " << stat.discharges_1 << " + "
                                   << stat.discharges_2 << " = "
                                   << stat.discharges_1 + stat.discharges_2 << "\n"
       << "c # global relabels : " << stat.global_relabels_1 << " + "
                                   << stat.global_relabels_2 << " = "
                                   << stat.global_relabels_1 + stat.global_relabels_2 << "\n"
       << "c # colors          : " << stat.colors << "\n"
       << "c # color ticks     : " << stat.color_ticks_1 << " + "
                                   << stat.color_ticks_2 << " = "
                                   << stat.color_ticks_1 + stat.color_ticks_2 << std::endl;

    return os;
}

# ifdef ENABLE_STAT
#     define STAT_INC(x)    ++(x)
#     define STAT_SET(x, y) (x) = (y)
#     define STAT_OUT(o, c) (o) << (c)
# else
#     define STAT_INC(x)    ((void)0)
#     define STAT_SET(x, y) ((void)0)
#     define STAT_OUT(o, c) ((void)0)
# endif

std::chrono::time_point<std::chrono::system_clock> beg_time, beg_tot_time;

std::chrono::duration<double> read_time;
std::chrono::duration<double> convert_residual_time;
std::chrono::duration<double> check_solution_time;
std::chrono::duration<double> write_time;

struct PushRelabelTime {
    std::chrono::duration<double> discharge_time_1;
    std::chrono::duration<double> discharge_time_2;
    std::chrono::duration<double> global_relabel_time_1;
    std::chrono::duration<double> global_relabel_time_2;
};

PushRelabelTime pr_time;

inline std::ostream&
operator<<(std::ostream& os, const PushRelabelTime& time)
{
    auto init_time = read_time + convert_residual_time;
    auto solve_time_1 = time.discharge_time_1 + time.global_relabel_time_1;
    auto solve_time_2 = time.discharge_time_2 + time.global_relabel_time_2;
    auto report_time = check_solution_time + write_time;
    auto total_time = init_time + solve_time_1 + solve_time_2 + report_time;

    os << "c\n"
       << "c total measured time       : " << total_time.count() << "\n"
       << "c - init time               : " << init_time.count() << "\n"
       << "c   - read time             : " << read_time.count() << "\n"
       << "c   - convert residual time : " << convert_residual_time.count() << "\n"
       << "c - solve time              : " << solve_time_1.count() << " + "
                                           << solve_time_2.count() << " = "
                                           << solve_time_1.count() +
                                              solve_time_2.count() << "\n"
       << "c   - discharge time        : " << time.discharge_time_1.count() << " + "
                                           << time.discharge_time_2.count() << " = "
                                           << time.discharge_time_1.count() +
                                              time.discharge_time_2.count() << "\n"
       << "c   - global relabel time   : " << time.global_relabel_time_1.count() << " + "
                                           << time.global_relabel_time_2.count() << " = "
                                           << time.global_relabel_time_1.count() +
                                              time.global_relabel_time_2.count() << "\n"
       << "c - report time             : " << report_time.count() << "\n"
       << "c   - check solution time   : " << check_solution_time.count() << "\n"
       << "c   - write time            : " << write_time.count() << std::endl;

    return os;
}

struct ColoredParallelPushRelabelTime {
    std::chrono::duration<double> color_graph_time;
    std::chrono::duration<double> discharge_time_1;
    std::chrono::duration<double> discharge_time_2;
    std::chrono::duration<double> global_relabel_time_1;
    std::chrono::duration<double> global_relabel_time_2;
};

ColoredParallelPushRelabelTime cppr_time;

inline std::ostream&
operator<<(std::ostream& os, const ColoredParallelPushRelabelTime& time)
{
    auto init_time = read_time + convert_residual_time;
    auto solve_time_1 = time.color_graph_time + time.discharge_time_1 + time.global_relabel_time_1;
    auto solve_time_2 = time.discharge_time_2 + time.global_relabel_time_2;
    auto report_time = check_solution_time + write_time;
    auto total_time = init_time + solve_time_1 + solve_time_2 + report_time;

    os << "c\n"
       << "c total measured time           : " << total_time.count() << "\n"
       << "c - init time                   : " << init_time.count() << "\n"
       << "c   - read time                 : " << read_time.count() << "\n"
       << "c   - convert residual time     : " << convert_residual_time.count() << "\n"
       << "c - solve time                  : " << solve_time_1.count() << " + "
                                               << solve_time_2.count() << " = "
                                               << solve_time_1.count() +
                                                  solve_time_2.count() << "\n"
       << "c   - color graph time          : " << time.color_graph_time.count() << "\n"
       << "c   - discharge time (par)      : " << time.discharge_time_1.count() << " + "
                                               << time.discharge_time_2.count() << " = "
                                               << time.discharge_time_1.count() +
                                                  time.discharge_time_2.count() << "\n"
       << "c   - global relabel time (par) : " << time.global_relabel_time_1.count() << " + "
                                               << time.global_relabel_time_2.count() << " = "
                                               << time.global_relabel_time_1.count() +
                                                  time.global_relabel_time_2.count() << "\n"
       << "c - report time                 : " << report_time.count() << "\n"
       << "c   - check solution time       : " << check_solution_time.count() << "\n"
       << "c   - write time                : " << write_time.count() << std::endl;

    return os;
}

# ifdef ENABLE_TIME
#     define TIME_BEG            beg_time = std::chrono::system_clock::now()
#     define TIME_END(t)         (t) = std::chrono::system_clock::now() - beg_time
#     define TIME_BEG_TOT        beg_tot_time = std::chrono::system_clock::now()
#     define TIME_END_TOT(t, t0) (t) = (std::chrono::system_clock::now() - beg_tot_time) - (t0)
#     define TIME_ADD(t)         (t) += std::chrono::system_clock::now() - beg_time
#     define TIME_SET(t1, t2)    (t1) = (t2)
#     define TIME_OUT(o, t)      (o) << (t)
# else
#     define TIME_BEG            ((void)0)
#     define TIME_END(t)         ((void)0)
#     define TIME_BEG_TOT        ((void)0)
#     define TIME_END_TOT(t, t0) ((void)0)
#     define TIME_ADD(t)         ((void)0)
#     define TIME_SET(t1, t2)    ((void)0)
#     define TIME_OUT(o, t)      ((void)0)
# endif

template <typename Graph>
inline void
check_solution_max(std::ostream& os, Graph& graph)
{
    TIME_BEG;
    INFO("checking solution..");

    node_t n = graph.num_nodes();

    std::vector<ulval_t> excs(n, 0);
    for (node_t i = 0; i < n; ++i) {
        for (auto& arc : graph[i]) {
            if (arc.cost > 0) {
                node_t j = arc.head;
                uval_t f = arc.pair->rcap;
                excs[i] -= f;
                excs[j] += f;
            }
        }
    }

    // flow conservation
    bool feasible = true;
    for (node_t i = 0; i < n; ++i) {
        if (excs[i] && i != graph.s && i != graph.t) {
            LOG("excess of node ", i, " is ", excs[i]);
            feasible = false;
        }
    }

    TIME_END(check_solution_time);

    if (feasible) {
        os << "c\nc solution is feasible" << std::endl;
    } else {
        os << "c\nc solution is NOT feasible (flow conservation is violated)" << std::endl;
    }
}

#endif // COMMON_HPP
