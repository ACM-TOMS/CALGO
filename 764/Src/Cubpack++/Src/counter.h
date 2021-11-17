/////////////////////////////////////////////////////////
//                                                     //
//    Cubpack++                                        //
//                                                     //
//        A Package For Automatic Cubature             //
//                                                     //
//        Authors : Ronald Cools                       //
//                  Dirk Laurie                        //
//                  Luc Pluym                          //
//                                                     //
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// File : counter.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   28 Mar 1996     V0.1h(long instead of int)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS Counter
// ---------------------------
//
// BASECLASSES:
//   None
//
//
// PURPOSE:
//   Counter is a pure virtual base class for implementing
//   stop-watch-like classes.
//
// TEMPLATES:
//   None
// METHODS:
//   CONSTRUCTORS:
//     None
//
//   SELECTORS:
//     1) virtual unsigned long Read()=0
//     --------------------------
//     Read reads out the counter
//   MODIFIERS:
//     1) virtual void Start()=0
//     -------------------------
//     Start() makes the stop-watch running. Successive calls
//     of Start() have no effect.
//     2) virtual void Stop()=0
//     -------------------------
//     Stops the stop-watch. After Stop(), Read() will always
//     read out the same value, until the stop-watch is started
//     again. Stopping a stopped stop-watch has no effect.
//     3) virtual void Reset()=0
//     -------------------------
//     After a call of Reset(), Read() will read out zero,
//     until the stop-watch is started again. A running
//     stop-watch should not be Reset().
//
//   OPERATORS:
//     None
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

#ifndef COUNTER_H
#define COUNTER_H

//////////////////////////////////////////////
class Counter
  {
  public:

  virtual void Start()=0;
  virtual void Stop()=0;
  virtual void Reset()=0;
  virtual unsigned long  Read()=0;
  };
//////////////////////////////////////////////

#endif
