//
// color_queue.hpp
//
// 2d array of arrays as a queue with colors and threads
//

#ifndef COLOR_QUEUE_HPP
#define COLOR_QUEUE_HPP

#include <algorithm>
#include <vector>

#include <omp.h>

#include "common.hpp"

class ColorQueue {
public:
    ColorQueue()
        : curr_color(0)
    { }

    void
    set_colors(const GraphColoring& coloring);

    void
    push(node_t i, int tid = 0)
    { elems[colors[i]][tid].push_back(i); }

    std::vector<node_t>
    pop_elems();

    bool
    empty();

private:
    node_t curr_color;
    std::vector<node_t> colors;
    std::vector<std::vector<std::vector<node_t>>> elems;
};

inline void
ColorQueue::set_colors(const GraphColoring& coloring)
{
    colors = coloring.colors;

    elems = std::vector<std::vector<std::vector<node_t>>>(
            coloring.num_colors,
            std::vector<std::vector<node_t>>(omp_get_max_threads()));

    node_t size = colors.size() / (omp_get_max_threads() * coloring.num_colors);

    for (auto& vv : elems) {
        for (auto& v : vv) {
            v.reserve(size);
        }
    }
}

inline std::vector<node_t>
ColorQueue::pop_elems()
{
    std::vector<node_t> curr;
    while (curr.size() == 0) {
        ++curr_color;
        curr_color %= elems.size();
        curr = join_vectors(elems[curr_color]);
    }

    // INFO("number of active vertices : ", curr.size());

    for (auto& v : elems[curr_color]) {
        v.clear();
    }

    return curr;
}

inline bool
ColorQueue::empty()
{
    for (auto& vv : elems) {
        for (auto& v : vv) {
            if (!v.empty()) { return false; }
        }
    }

    return true;
}

#endif // COLOR_QUEUE_HPP
