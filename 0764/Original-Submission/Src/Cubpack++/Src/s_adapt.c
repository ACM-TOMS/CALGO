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
//////////////////////////////////////////
//File s_adapt.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(typedef introduced, code from .c to .h)
//////////////////////////////////////////////

#include <s_adapt.h>
#include <atomic.h>
#include <error.h>
#include <pointer.h>

/////////////////////////////////////////////

//template <class GEOMETRY>
//SimpleAdaptive<GEOMETRY>::SimpleAdaptive(const Pointer<RuleGEOMETRY>& R,
//  const Pointer<SameShapeDivisorGEOMETRY>& D)
//    :Processor<GEOMETRY>(),
//    TimesCalled(0),
//    Diffs(2),
//    TheRule(R),
//    TheDivisor(D)
//  {
//  }
///////////////////////////////////////////////

template <class GEOMETRY>
SimpleAdaptive<GEOMETRY>::SimpleAdaptive(Rule<GEOMETRY>* R,
  SameShapeDivisor<GEOMETRY>* D)
    :Processor<GEOMETRY>(),
    TimesCalled(0),
    Diffs(2),
    TheRule(R),
    TheDivisor(D)
  {
  }
///////////////////////////////////////////////
template <class GEOMETRY>
SimpleAdaptive<GEOMETRY>*
SimpleAdaptive<GEOMETRY>::Descendant()
const
  {
  SimpleAdaptive<GEOMETRY>* r = new
        SimpleAdaptive<GEOMETRY>(TheRule,TheDivisor);
  return r;
  }
////////////////////////////////////////////////
template <class GEOMETRY>
void
SimpleAdaptive<GEOMETRY>::Process( Stack<AtomicRegion>&  Offspring)
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
          = Geometry().Volume()/TheDivisor->NumberOfParts();
    Stack<GEOMETRY> Parts;
    Vector<unsigned int> DiffOrder(Diffs.Size());
    if (Diffs[0]>Diffs[1])
      {
      DiffOrder[0] = 1 ; DiffOrder[1] = 2;
      }
    else
      {
      DiffOrder[0] = 2; DiffOrder[1] =1;
      };
    TheDivisor->Apply (Geometry(),Parts,DiffOrder);
    unsigned int N = Parts.Size();
    for (unsigned int i =0;i<N;i++)
      {
      GEOMETRY* g = Parts.Pop();
      g->Volume(NewVolume);
      Processor<GEOMETRY>* p = Descendant();
      Atomic<GEOMETRY>* a = new Atomic<GEOMETRY>(g,p);
      a->LocalIntegrand(&LocalIntegrand());
      Offspring.Push(a);
      };
    return;
    };
   Error(TimesCalled > 2,
     "SimpleAdaptive : more than two calls of Process()");
   }
///////////////////////////////////////////////
template <class GEOMETRY>
Processor<GEOMETRY>*
SimpleAdaptive<GEOMETRY>::NewCopy()
const
  {
  return new SimpleAdaptive<GEOMETRY>(*this);
  }
/////////////////////////////////////////////

