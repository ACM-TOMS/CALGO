//
// star_graph.hpp
//
// Star graph representation
//

#ifndef STAR_GRAPH_HPP
#define STAR_GRAPH_HPP

#include "common.hpp"

struct StarNode {
    std::vector<ResArc>::iterator it_beg;
    std::vector<ResArc>::iterator it_end;

    std::vector<ResArc>::iterator
    begin()
    { return it_beg; }

    std::vector<ResArc>::iterator
    end()
    { return it_end; }
};

struct StarGraph {
    std::vector<StarNode> nodes;
    std::vector<ResArc> arcs;

    ulval_t m;

    node_t s, t;

    StarGraph(const Problem& problem);

    node_t
    num_nodes()
    { return nodes.size(); }

    ulval_t
    num_arcs()
    { return m; }

    StarNode&
    operator[](node_t i)
    { return nodes[i]; }
};

StarGraph::StarGraph(const Problem& problem)
    : nodes(problem.n)
    , m(problem.m)
    , s(problem.s)
    , t(problem.t)
{
    TIME_BEG;
    INFO("converting to residual network..");

    node_t n = problem.n;

    std::vector<node_t> lens(n);
    for (ulval_t i = 0; i < m; ++i) {
        ++lens[problem.arcs[i].head];
        ++lens[problem.arcs[i].tail];
    }

    nodes.resize(n);

    ulval_t size = 0;
    std::vector<node_t> begs(n);
    for (node_t i = 0; i < n; ++i) {
        begs[i] = size;
        size += lens[i];
    }

    arcs.resize(size);

    for (node_t i = 0; i < n; ++i) {
        nodes[i].it_beg = arcs.begin() + begs[i];
        nodes[i].it_end = arcs.begin() + begs[i] + lens[i];
    }

    for (ulval_t curr = 0; curr < m; ++curr) {
        node_t i = problem.arcs[curr].tail;
        node_t j = problem.arcs[curr].head;
        uval_t u = problem.arcs[curr].upper;
        arcs[begs[i]] = ResArc{j, u, 1, nullptr};
        arcs[begs[j]] = ResArc{i, 0, -1, nullptr};
        arcs[begs[i]].pair = &arcs[begs[j]];
        arcs[begs[j]].pair = &arcs[begs[i]];
        ++begs[i];
        ++begs[j];
    }

    TIME_END(convert_residual_time);
}

#endif // STAR_GRAPH_HPP
