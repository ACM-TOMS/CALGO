//
// list_graph.hpp
//
// Adjacency list graph representation
//

#ifndef LIST_GRAPH_HPP
#define LIST_GRAPH_HPP

#include "common.hpp"

struct ListGraph {
    std::vector<std::vector<ResArc>> nodes;

    ulval_t m;

    node_t s, t;

    ListGraph(const Problem& problem);

    node_t
    num_nodes()
    { return nodes.size(); }

    ulval_t
    num_arcs()
    { return m; }

    std::vector<ResArc>&
    operator[](node_t i)
    { return nodes[i]; }
};

ListGraph::ListGraph(const Problem& problem)
    : nodes(problem.n)
    , m(problem.m)
    , s(problem.s)
    , t(problem.t)
{
    TIME_BEG;
    INFO("converting to residual network..");

    node_t n = nodes.size();

    std::vector<node_t> lens(n);
    for (auto& arc : problem.arcs) {
        ++lens[arc.head];
        ++lens[arc.tail];
    }

    // required to prevent reallocations which would change memory addresses
    for (node_t i = 0; i < n; ++i) {
        nodes[i].reserve(lens[i]);
    }

    for (auto& arc : problem.arcs) {
        node_t i = arc.tail;
        node_t j = arc.head;
        uval_t u = arc.upper;
        nodes[i].push_back(ResArc{j, u, 1, nullptr});
        nodes[j].push_back(ResArc{i, 0, -1, nullptr});
        nodes[i].back().pair = &nodes[j].back();
        nodes[j].back().pair = &nodes[i].back();
    }

    TIME_END(convert_residual_time);
}

#endif // LIST_GRAPH_HPP
