//
// push_relabel.hpp
//
// Push relabel algorithm for maximum flow problem
//

#ifndef PUSH_RELABEL_HPP
#define PUSH_RELABEL_HPP

#include <limits>
#include <queue>
#include <set>
#include <unordered_set>

#include <omp.h>

#include "common.hpp"
#include "color_queue.hpp"
#include "multi_queue.hpp"
#include "greedy_coloring.hpp"

//
// (reverse) breadth-first search (sequential)
//
template <typename Graph>
inline void
global_relabel(Graph& graph,
        std::vector<uval_t>& dist,
        std::queue<node_t>& q,
        bool forward)
{
    node_t n = graph.num_nodes();

    node_t dmax = forward ? n : (2 * n);

    std::fill(dist.begin(), dist.end(), dmax);

    dist[graph.s] = n;
    dist[graph.t] = 0;

    if (forward) {
        q.push(graph.t);
    } else {
        q.push(graph.s);
    }

    while (!q.empty()) {
        node_t i = q.front();
        q.pop();
        for (auto& arc : graph[i]) {
            if (arc.pair->rcap == 0) { continue; }

            node_t j = arc.head;

            if (dist[j] == dmax) {
                dist[j] = dist[i] + 1;
                q.push(j);
            }
        }
    }
}

//
// (reverse) breadth-first search (parallel)
//
template <typename Graph>
inline void
global_relabel(Graph& graph,
        std::vector<uval_t>& dist,
        MultiQueue& q,
        bool forward)
{
    node_t n = graph.num_nodes();

    node_t dmax = forward ? n : (2 * n);

    # pragma omp parallel for
    for (node_t i = 0; i < n; ++i) {
        dist[i] = dmax;
    }

    dist[graph.s] = n;
    dist[graph.t] = 0;

    if (forward) {
        q.push(graph.t);
    } else {
        q.push(graph.s);
    }

    node_t level = forward ? 0 : n;

    while (!q.empty()) {
        auto elems = q.pop_elems();

        ++level;

        # pragma omp parallel for
        for (size_t i = 0; i < elems.size(); ++i) {
            for (auto& arc : graph[elems[i]]) {
                if (arc.pair->rcap == 0) { continue; }

                node_t j = arc.head;

                if (dist[j] == dmax) {
                    uval_t old;

                    # pragma omp atomic capture
                    { old = dist[j]; dist[j] = level; }

                    if (old == dmax) {
                        q.push(arc.head, omp_get_thread_num());
                    }
                }
            }
        }
    }
}

inline void
inc_excs(std::queue<node_t>& acts,
        std::vector<ulval_t>& excs,
        uval_t f,
        node_t j,
        node_t t)
{
    ulval_t old = excs[j];

    excs[j] += f;

    if (old == 0 && j != t) {
        acts.push(j);
    }
}

inline void
inc_excs(ColorQueue& acts,
        std::vector<ulval_t>& excs,
        uval_t f,
        node_t j,
        node_t t)
{
    ulval_t old;

    # pragma omp atomic capture
    { old = excs[j]; excs[j] += f; }

    if (old == 0 && j != t) {
        acts.push(j, omp_get_thread_num());
    }
}

template <typename Graph, typename Queue>
inline void
discharge(Graph& graph,
        std::vector<ulval_t>& excs,
        std::vector<uval_t>& dist,
        Queue& acts,
        node_t i,
        ulval_t& num_pushes,
        ulval_t& num_relabels,
        bool forward)
{
    (void)num_pushes;

    while (true) {
        if (forward && dist[i] >= graph.num_nodes()) { return; }

        uval_t dmin = std::numeric_limits<uval_t>::max();

        for (auto& arc : graph[i]) {
            if (arc.rcap == 0) { continue; }

            node_t j = arc.head;

            if (dist[i] == dist[j] + 1) {
                uval_t f = std::min<ulval_t>(excs[i], arc.rcap);

                LOG("push ", f, " from ", i, " to ", j);
                STAT_INC(num_pushes);

                arc.rcap -= f;
                arc.pair->rcap += f;

                excs[i] -= f;

                inc_excs(acts, excs, f, j, graph.t);

                if (excs[i] == 0) { return; }
            } else {
                dmin = std::min(dmin, dist[j]);
            }
        }

        LOG("relabel node ", i, " to ", dmin + 1);
        ++num_relabels;

        dist[i] = dmin + 1;
    }
}

template <typename Graph>
inline void
process_queue(Graph& graph,
        std::vector<ulval_t>& excs,
        std::vector<uval_t>& dist,
        std::queue<node_t>& acts,
        ulval_t& num_pushes,
        ulval_t& num_relabels,
        ulval_t& num_discharges,
        ulval_t& num_color_ticks,
        bool forward)
{
    (void)num_color_ticks;

    node_t i = acts.front();
    acts.pop();

    LOG("select ", i, " as the active node");
    STAT_INC(num_discharges);

    discharge(graph, excs, dist, acts, i, num_pushes, num_relabels, forward);
}

template <typename Graph>
inline void
process_queue(Graph& graph,
        std::vector<ulval_t>& excs,
        std::vector<uval_t>& dist,
        ColorQueue& acts,
        ulval_t& num_pushes,
        ulval_t& num_relabels,
        ulval_t& num_discharges,
        ulval_t& num_color_ticks,
        bool forward)
{
    auto elems = acts.pop_elems();

    LOG("select ", elems, " as active nodes");
    STAT_SET(num_discharges, num_discharges + elems.size());
    STAT_INC(num_color_ticks);

    # pragma omp parallel for reduction(+:num_pushes) reduction(+:num_relabels)
    for (size_t i = 0; i < elems.size(); ++i) {
        discharge(graph, excs, dist, acts, elems[i], num_pushes, num_relabels, forward);
    }
}

template <typename Graph>
inline void
init_queue(std::queue<node_t>&, Graph&, ulval_t&)
{ }

template <typename Graph>
inline void
init_queue(ColorQueue& queue, Graph& graph, ulval_t& num_colors)
{
    TIME_BEG;
    INFO("coloring graph..");

    auto coloring = greedy_coloring<Graph>(graph);
    queue.set_colors(coloring);
    STAT_SET(num_colors, coloring.num_colors);

    TIME_END(cppr_time.color_graph_time);
}

template <typename Graph>
inline void
init_gr_queue(std::queue<node_t>&, Graph&)
{ }

template <typename Graph>
inline void
init_gr_queue(MultiQueue& queue, Graph& graph)
{
    node_t n = graph.num_nodes();
    queue.reserve(n);
}

template <typename Graph, typename Queue, typename GRQueue>
inline ulval_t
push_relabel(Graph& graph, bool mincut)
{
    ulval_t num_colors = 0;

    Queue acts;
    init_queue(acts, graph, num_colors);

    STAT_SET(cppr_stat.colors, num_colors);

    TIME_BEG_TOT;
    INFO("running algorithm (phase 1)..");

    node_t n = graph.num_nodes();

    std::vector<uval_t> dist(n, 0);
    std::vector<ulval_t> excs(n, 0);

    dist[graph.s] = n;

    bool forward = true;

    ulval_t num_pushes = 0;
    ulval_t num_relabels = n;
    ulval_t num_discharges = 0;
    ulval_t num_global_relabels = 0;
    ulval_t num_color_ticks = 0;

    for (auto& arc : graph[graph.s]) {
        node_t j = arc.head;
        uval_t f = arc.rcap;

        LOG("initial push ", f, " from ", graph.s, " to ", j);
        STAT_INC(num_pushes);

        arc.rcap -= f;
        arc.pair->rcap += f;

        excs[graph.s] -= f;
        excs[j] += f;

        if (excs[j] == f && j != graph.t) {
            acts.push(j);
        }
    }

    GRQueue grq;
    init_gr_queue(grq, graph);

    while (!acts.empty()) {
        if (num_relabels >= n) {
            TIME_BEG;

            LOG("global relabeling..");
            STAT_INC(num_global_relabels);
            STAT_SET(pr_stat.relabels_1, pr_stat.relabels_1 + num_relabels);

            num_relabels = 0;

            global_relabel(graph, dist, grq, forward);

            TIME_ADD(pr_time.global_relabel_time_1);
        }

        process_queue(graph,
                excs,
                dist,
                acts,
                num_pushes,
                num_relabels,
                num_discharges,
                num_color_ticks,
                forward);
    }

    STAT_SET(pr_stat.pushes_1, num_pushes);
    STAT_SET(pr_stat.relabels_1, pr_stat.relabels_1 + num_relabels - n);
    STAT_SET(pr_stat.discharges_1, num_discharges);
    STAT_SET(pr_stat.global_relabels_1, num_global_relabels);

    STAT_SET(cppr_stat.pushes_1, pr_stat.pushes_1);
    STAT_SET(cppr_stat.relabels_1, pr_stat.relabels_1);
    STAT_SET(cppr_stat.discharges_1, pr_stat.relabels_1);
    STAT_SET(cppr_stat.global_relabels_1, pr_stat.global_relabels_1);

    STAT_SET(cppr_stat.color_ticks_1, num_color_ticks);

    TIME_END_TOT(pr_time.discharge_time_1, pr_time.global_relabel_time_1);

    TIME_SET(cppr_time.discharge_time_1, pr_time.discharge_time_1);
    TIME_SET(cppr_time.global_relabel_time_1, pr_time.global_relabel_time_1);

    if (mincut) { return excs[graph.t]; }

    TIME_BEG_TOT;
    INFO("running algorithm (phase 2)..");

    forward = false;

    num_pushes = 0;
    num_relabels = n;
    num_discharges = 0;
    num_global_relabels = 0;
    num_color_ticks = 0;

    for (node_t i = 0; i < n; ++i) {
        if (i == graph.s || i == graph.t) { continue; }
        if (excs[i] > 0) {
            acts.push(i);
        }
    }

    while (!acts.empty()) {
        if (num_relabels >= n) {
            TIME_BEG;

            LOG("global relabeling..");
            STAT_INC(num_global_relabels);
            STAT_SET(pr_stat.relabels_2, pr_stat.relabels_2 + num_relabels);

            num_relabels = 0;

            global_relabel(graph, dist, grq, forward);

            TIME_ADD(pr_time.global_relabel_time_2);
        }

        process_queue(graph,
                excs,
                dist,
                acts,
                num_pushes,
                num_relabels,
                num_discharges,
                num_color_ticks,
                forward);
    }

    STAT_SET(pr_stat.pushes_2, num_pushes);
    STAT_SET(pr_stat.relabels_2, pr_stat.relabels_2 + num_relabels - n);
    STAT_SET(pr_stat.discharges_2, num_discharges);
    STAT_SET(pr_stat.global_relabels_2, num_global_relabels);

    STAT_SET(cppr_stat.pushes_2, pr_stat.pushes_2);
    STAT_SET(cppr_stat.relabels_2, pr_stat.relabels_2);
    STAT_SET(cppr_stat.discharges_2, pr_stat.relabels_2);
    STAT_SET(cppr_stat.global_relabels_2, pr_stat.global_relabels_2);

    STAT_SET(cppr_stat.color_ticks_2, num_color_ticks);

    TIME_END_TOT(pr_time.discharge_time_2, pr_time.global_relabel_time_2);

    TIME_SET(cppr_time.discharge_time_2, pr_time.discharge_time_2);
    TIME_SET(cppr_time.global_relabel_time_2, pr_time.global_relabel_time_2);

    return excs[graph.t];
}

#endif // PUSH_RELABEL_HPP
