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
//File S2adapt.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   30 Jan 1996     V0.1g(no longer compare signed and unsigned)
//////////////////////////////////////////
#include <S2adapt.h>
#include <point.h>
#include <real.h>
#include <polC2itf.h>
#include <atomic.h>

////////////////////////////////////////////////
void
CircleAdaptive::Process(Stack<AtomicRegion>& Offspring)
  {
  Circle& C = Geometry();
  static unsigned int mxgen=2;
  //if (mxgen == 0)
    //{
      //cerr<<"How many times may I apply the circle rule? ";
      //cin>>mxgen;
    //};
  TimesCalled++;
  if(TimesCalled ==1)
    {
    if (GenerationNumber < mxgen)
      {
       TheRule->Apply(LocalIntegrand(),C,Integral(),AbsoluteError());
      }
    else
      {
       Point rv (C.Radius(),0);
       Point bp = C.Center()+rv;
       AtomicRegion* A ;
       A= (AtomicRegion*) POLAR_RECTANGLE(C.Center(),bp,bp);
       A->LocalIntegrand(&(LocalIntegrand()));
       Offspring.Push(A);
      };
    }
  else
    {
    const real CircleRatio = 0.44721359549995793928;
    Point m1(C.Radius()*CircleRatio ,0.);
    Point m2(C.Radius() , 0.);
    Point m3(0.,C.Radius());
    Point m4(0.,C.Radius()*CircleRatio);
    Point origin=C.Center();

    AtomicRegion*
        t1=(AtomicRegion*) POLAR_RECTANGLE(origin+m1,origin+m2,origin+m3);
    AtomicRegion*
        t2=(AtomicRegion*) POLAR_RECTANGLE(origin+m4,origin+m3,origin-m2);
    AtomicRegion*
        t3=(AtomicRegion*) POLAR_RECTANGLE(origin-m1,origin-m2,origin-m3);
    AtomicRegion*
        t4=(AtomicRegion*) POLAR_RECTANGLE(origin-m4,origin-m3,origin+m2);

    Offspring.Push(t1);
    Offspring.Push(t2);
    Offspring.Push(t3);
    Offspring.Push(t4);
    Offspring.IteratorReset();
    while (!Offspring.IteratorAtEnd())
      {
      Offspring.IteratorNext()->LocalIntegrand(&LocalIntegrand());
      };
    if (CircleRatio != 0.0)
      {
       Circle* tmp = new Circle(origin,C.Radius()*CircleRatio);
       Atomic<Circle>* a =new Atomic<Circle>(tmp,Descendant());
       a->LocalIntegrand(&LocalIntegrand());
       Offspring.Push(a);
      };
    };
  }
  //////////////////////////////////////////////////
  CircleAdaptive::CircleAdaptive(Rule <Circle>*  R)
    :Processor<Circle>(),TheRule(R),TimesCalled(0),GenerationNumber(0)
    {
    }
/////////////////////////////////////////////////////
CircleAdaptive*
CircleAdaptive::Descendant()
const
  {
  CircleAdaptive * CA = new CircleAdaptive (&(*(TheRule)));
  CA->GenerationNumber = GenerationNumber+1;
  return CA;
  }
//////////////////////////////////////////////////////
Processor<Circle>*
CircleAdaptive::NewCopy()
const
  {
  return new CircleAdaptive(*this);
  }
//////////////////////////////////////////////////
