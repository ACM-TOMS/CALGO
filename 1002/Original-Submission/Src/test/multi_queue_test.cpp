//
// multi_queue_test.cpp
//

#include <cassert>
#include <vector>

#include "multi_queue.hpp"

using std::vector;

int
main()
{
    MultiQueue queue;

    vector<node_t> elems;

    queue.push(2, 0);

    elems = queue.pop_elems();
    assert(elems[0] == 2);

    queue.push(3, 0);
    queue.push(5, 1);

    elems = queue.pop_elems();
    assert(elems.size() == 2);

    queue.push(3, 0);
    queue.push(1, 1);
    queue.push(2, 2);

    elems = queue.pop_elems();
    assert(elems.size() == 3);

    return 0;
}
