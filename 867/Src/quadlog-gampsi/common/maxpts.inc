*
*     Parameter variables.  A prime value differing from a power of
*     two by at least the cache line size may help to reduce cache
*     collisions in expressions involving different arrays with the
*     same subscript.
*
      INTEGER             maxpts
      PARAMETER           (maxpts = 1049)
