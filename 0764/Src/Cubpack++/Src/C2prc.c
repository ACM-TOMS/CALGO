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
//File C2prc.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(typedef introduced and long lines split)
////////////////////////////////////////////////////////////
#include <C2prc.h>
#include <C2rule13.h>
#include <C2dv2.h>
#include <C2dv4.h>
#include <error.h>
#include <atomic.h>
////////////////////////////////////////////////////////
Pointer < Parallelogram_Processor::RuleParallelogram > 
  Parallelogram_Processor::TheRule = new Parallelogram_Rule13;
Pointer < Parallelogram_Processor::SameShapeDivisorParallelogram > 
  Parallelogram_Processor::TheDivisor2 = new Parallelogram_Divide2;
Pointer < Parallelogram_Processor::SameShapeDivisorParallelogram > 
  Parallelogram_Processor::TheDivisor4 = new Parallelogram_Divide4;
////////////////////////////////////////////////////////
Parallelogram_Processor::Parallelogram_Processor()
  :TimesCalled(0),
   Diffs(2)
   {
   }
/////////////////////////////////////////////////////////////
Parallelogram_Processor*
Parallelogram_Processor::Descendant()
const
  {
  Parallelogram_Processor* r = new Parallelogram_Processor;
  return r;
  }
////////////////////////////////////////////////
void
Parallelogram_Processor::Process( Stack<AtomicRegion>&  Offspring)
  {
  TimesCalled ++;
  if (TimesCalled == 1)
    {
    TheRule->ApplyWithDiffs(LocalIntegrand(),Geometry(),Integral(),
                   AbsoluteError(),Diffs);
    Offspring.MakeEmpty();
    return;
    };
  if(TimesCalled == 2)
    {
    real NewVolume
          = Geometry().Volume()/2;
    Stack<Parallelogram> Parts;
    Vector<unsigned int> DiffOrder(Diffs.Size());
    const real difffac = 1.81818, difftreshold = 0.00004;
    if (max(Diffs[0],Diffs[1]) < difftreshold)
      {
      TheDivisor4->Apply(Geometry(),Parts,DiffOrder);
      NewVolume /=2;
      }
    else if (Diffs[0]>difffac*Diffs[1])
      {
      DiffOrder[0] = 1 ; DiffOrder[1] = 2;
      TheDivisor2->Apply (Geometry(),Parts,DiffOrder);
      }
    else if (Diffs[1]>difffac*Diffs[0])
      {
      DiffOrder[0] = 2; DiffOrder[1] =1;
      TheDivisor2->Apply (Geometry(),Parts,DiffOrder);
      }
    else
      {
      TheDivisor4->Apply(Geometry(),Parts,DiffOrder);
      NewVolume /=2;
      };
    unsigned int N = Parts.Size();
    for (unsigned int i =0;i<N;i++)
      {
      Parallelogram* g = Parts.Pop();
      g->Volume(NewVolume);
      Processor<Parallelogram>* p = Descendant();
      Atomic<Parallelogram>* a = new Atomic<Parallelogram>(g,p);
      a->LocalIntegrand(&LocalIntegrand());
      Offspring.Push(a);
      };
    return;
    };
   Error(TimesCalled > 2,
     "Parallelogram_Processor : more than two calls of Process()");
   }
///////////////////////////////////////////////
Processor<Parallelogram>*
Parallelogram_Processor::NewCopy()
const
  {
  return new Parallelogram_Processor(*this);
  }
///////////////////////////////////////////////
