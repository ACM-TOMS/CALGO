//
// color_queue_test.cpp
//

#include <cassert>
#include <vector>

#include "color_queue.hpp"

using std::vector;

int
main()
{
    ColorQueue queue;

    auto colors = vector<node_t>{0, 1, 1, 0, 0, 1};

    auto coloring = GraphColoring{colors, 2};

    queue.set_colors(coloring);

    queue.push(2, 0);
    queue.push(3, 1);
    queue.push(1, 2);
    queue.push(0, 0);
    queue.push(5, 1);
    queue.push(4, 2);

    vector<node_t> elems;

    elems = queue.pop_elems();
    for (auto elem : elems) {
        assert(colors[elem] == colors[elems[0]]);
    }

    elems = queue.pop_elems();
    for (auto elem : elems) {
        assert(colors[elem] == colors[elems[0]]);
    }

    return 0;
}
