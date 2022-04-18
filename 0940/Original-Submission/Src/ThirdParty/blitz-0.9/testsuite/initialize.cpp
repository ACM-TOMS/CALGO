#include "testsuite.h"
#include <blitz/vector.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Vector<int> a(6);
    a = 1, 4, 6, 3, 7, 8;

    BZTEST(a(0) == 1);
    BZTEST(a(1) == 4);
    BZTEST(a(2) == 6);
    BZTEST(a(3) == 3);
    BZTEST(a(4) == 7);
    BZTEST(a(5) == 8);

    Vector<int> b = a(Range(0,4,2));
    b = -1, -2, -3;

    BZTEST(a(0) == -1);
    BZTEST(a(1) == 4);
    BZTEST(a(2) == -2);
    BZTEST(a(3) == 3);
    BZTEST(a(4) == -3);
    BZTEST(a(5) == 8);

    Vector<int> c = a.reverse();
    c = 0;
    BZTEST(a(0) == 0);
    BZTEST(a(1) == 0);
    BZTEST(a(2) == 0);
    BZTEST(a(3) == 0);
    BZTEST(a(4) == 0);
    BZTEST(a(5) == 0);

    c = 8, 7, 3, 6, 4, 1;
    BZTEST(a(0) == 1);
    BZTEST(a(1) == 4);
    BZTEST(a(2) == 6);
    BZTEST(a(3) == 3);
    BZTEST(a(4) == 7);
    BZTEST(a(5) == 8);

    return 0;
}

