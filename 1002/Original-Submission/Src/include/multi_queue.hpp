//
// multi_queue.hpp
//
// An array of arrays as a queue with threads
//

#ifndef MULTI_QUEUE_HPP
#define MULTI_QUEUE_HPP

#include <algorithm>
#include <vector>

#include <omp.h>

#include "common.hpp"

class MultiQueue {
public:
    MultiQueue()
        : elems(omp_get_max_threads())
    { }

    void
    reserve(node_t n);

    void
    push(node_t i, int tid = 0)
    { elems[tid].push_back(i); }

    std::vector<node_t>
    pop_elems();

    bool
    empty();

private:
    std::vector<std::vector<node_t>> elems;
};

inline void
MultiQueue::reserve(node_t n)
{
    node_t size = n / omp_get_max_threads();

    for (auto& v : elems) {
        v.reserve(size);
    }
}

inline std::vector<node_t>
MultiQueue::pop_elems()
{
    auto curr = join_vectors(elems);

    for (auto& v : elems) {
        v.clear();
    }

    return curr;
}

inline bool
MultiQueue::empty()
{
    for (auto& v : elems) {
        if (!v.empty()) { return false; }
    }

    return true;
}

#endif // MULTI_QUEUE_HPP
