//
// push_relabel_test.cpp
//

#include <cassert>
#include <fstream>
#include <queue>

#include "dimacs.hpp"
#include "push_relabel.hpp"

using std::ifstream;
using std::queue;
using std::string;

template <typename Graph, typename Queue, typename GRQueue>
void
solve_problem(const string& fname, ulval_t flow)
{
    ifstream inp(fname);
    auto graph = Graph(read_dimacs(inp));
    auto result = push_relabel<Graph, Queue, GRQueue>(graph, false);
    assert(result == flow);
}

int
main()
{
    solve_problem<ListGraph, queue<node_t>, queue<node_t>>("test/sample/inp.max", 15);
    solve_problem<StarGraph, queue<node_t>, queue<node_t>>("test/sample/inp.max", 15);

    solve_problem<ListGraph, ColorQueue, MultiQueue>("test/sample/inp.max", 15);
    solve_problem<StarGraph, ColorQueue, MultiQueue>("test/sample/inp.max", 15);

    return 0;
}
