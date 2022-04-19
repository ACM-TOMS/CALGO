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
//File polC2prc.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <polC2prc.h>
#include <C2toS2.h>
#include <translat.h>
#include <integran.h>
#include <polC2.h>
#include <C2interf.h>
#include <atomic.h>

////////////////////////////////////////////////////////

void
PolarRectangle_Processor::Process( Stack<AtomicRegion>& Offspring)
  {
  PolarRectangle& G = Geometry();
  Point P1(G.InnerRadius(),G.SmallAngle()), P2(G.OuterRadius(),G.SmallAngle()),
        P3(G.OuterRadius(),G.BigAngle());
  AtomicRegion* A ;
  A= (AtomicRegion*) PARALLELOGRAM(P2,P1,P3);
  Integrand I1(LocalIntegrand(),new Translation(G.Center()));
  A->LocalIntegrand(new Integrand(I1, new PolarToRectangular));
  Offspring.Push(A);
  }
/////////////////////////////////////////////////
PolarRectangle_Processor::PolarRectangle_Processor()
  :Processor<PolarRectangle>()
  {
  }
/////////////////////////////////////////////////       
Processor<PolarRectangle>*
PolarRectangle_Processor::NewCopy()
const
  {
  return new PolarRectangle_Processor(*this);
  }
////////////////////////////////////////////////
