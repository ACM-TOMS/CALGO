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
// File : eval_ctr.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   28 Mar 1996     V0.1h(long instead of int)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS EvaluationCounter
// -------------------------------------
//
// BASECLASSES:
//   Counter
//
//
// PURPOSE:
//     EvaluationCounter acts as a kind of stop-watch
//     counting the number of Integrand evaluations
//     during the time it's running.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) EvaluationCounter()
//     ----------------------
//     Default constructor. The stop-watch isn't running after
//     construction.
//   SELECTORS:
//     1) unsigned long Read()
//     ----------------
//     Read reads out the number of
//     IntegrandEvaluations made since Start() was
//     called. This value remains constant after Stop() is
//     called.
//   MODIFIERS:
//     1) void Start()
//     ---------------
//     Start() makes the stop-watch running. Successive calls
//     of Start() have no effect.
//     2) void Stop()
//     ---------------
//     Stops the stop-watch. After Stop(), Read() will always
//     read out the same value, until the stop-watch is started
//     again. Stopping a stopped stop-watch has no effect.
//     3) void Reset()
//     ---------------
//     After a call of Reset(), Read() will read out zero,
//     until the stop-watch is started again. A running
//     stop-watch should not be Reset().
//
//     4) void Reset(unsigned long value)
//     ----------------------------------
//     allows you to reset the counter to some other
//     value than zero. Counting will start from this value.
//
//   OPERATORS:
//     None
//   SPECIAL:
//     None
///////////////////////////////////////////////////////////


#ifndef EVAL_CTR_H
#define EVAL_CTR_H
/////////////////////////////////
#include <counter.h>
#include <boolean.h>

////////////////////////////////////////////////
class EvaluationCounter : public Counter
  {
  public:

  EvaluationCounter();
  void Start();
  void Stop();
  void Reset();
  void Reset(unsigned long);
  unsigned long Read();

  private:

  long Strt;
  long End;
  unsigned long Bias;
  Boolean Running;
  };

////////////////////////////////////////////////
#endif
