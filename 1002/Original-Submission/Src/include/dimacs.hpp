//
// dimacs.hpp
//
// I/O functions for DIMACS format
//

#ifndef DIMACS_HPP
#define DIMACS_HPP

#include <cassert>
#include <iostream>
#include <sstream>

#include "common.hpp"
#include "list_graph.hpp"
#include "star_graph.hpp"

inline Problem
read_dimacs(std::istream& inp)
{
    TIME_BEG;
    INFO("reading problem..");

    Problem problem;

    char c;
    std::string p;

    node_t node;
    char type;

    node_t tail, head;
    uval_t upper;

    std::string line;
    while (getline(inp, line)) {
        std::istringstream ss(line);
        ss >> c;
        switch (c) {
        case 'p': // problem line
            ss >> p >> problem.n >> problem.m;
            break;
        }
        if (inp.peek() == 'n') { break; }
    }

    problem.arcs.reserve(problem.m);

    if (p == "max") {
        while (getline(inp, line)) {
            std::istringstream ss(line);
            ss >> c;
            switch (c) {
            case 'n': // node descriptor
                ss >> node >> type;
                switch (type) {
                case 's':
                    problem.s = node - 1;
                    break;
                case 't':
                    problem.t = node - 1;
                    break;
                }
                break;
            case 'a': // arc descriptor
                ss >> tail >> head >> upper;
                problem.arcs.push_back(Arc{tail - 1, head - 1, upper});
                break;
            }
        }
    }

    TIME_END(read_time);

    return problem;
}

template <typename Graph>
inline void
print_dimacs(Graph& graph)
{
    TIME_BEG;
    INFO("writing solution..");

    node_t n = graph.num_nodes();

    for (node_t i = 0; i < n; ++i) {
        for (auto& arc : graph[i]) {
            if (arc.cost > 0) {
                std::cout << "f "
                          << i + 1 << " "
                          << arc.head + 1 << " "
                          << arc.pair->rcap << "\n";
            }
        }
    }

    std::cout << std::flush;

    TIME_END(write_time);
}

#endif // DIMACS_HPP
