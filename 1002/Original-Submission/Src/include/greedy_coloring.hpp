//
// greedy_coloring.hpp
//
// Greedy vertex coloring algorithm
//

#ifndef GREEDY_COLORING_HPP
#define GREEDY_COLORING_HPP

#include <vector>

#include "common.hpp"

template <typename Graph>
inline GraphColoring
greedy_coloring(Graph& graph)
{
    node_t n = graph.num_nodes();

    node_t num_colors = 0;

    node_t blank = std::numeric_limits<node_t>::max();

    std::vector<node_t> colors(n, blank);

    std::vector<char> mark(n, false);

    for (node_t i = 0; i < n; ++i) {
        for (auto& arc : graph[i]) {
            node_t j = arc.head;
            node_t c = colors[j];
            if (c != blank) {
                mark[c] = true;
            }
        }

        bool found = false;
        for (node_t c = 0; c < num_colors; ++c) {
            if (mark[c] == false) {
                colors[i] = c;
                found = true;
                break;
            }
        }

        for (node_t c = 0; c < num_colors; ++c) {
            mark[c] = false;
        }

        if (found) { continue; }

        colors[i] = num_colors++;
    }

    return GraphColoring{colors, num_colors};
}

#endif // GREEDY_COLORING_HPP
