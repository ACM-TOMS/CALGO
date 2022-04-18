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
////////////////////////////////////////////////////
//File E2secprc.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <E2secprc.h>
#include <C2toS2.h>
#include <translat.h>
#include <integran.h>
#include <polC2.h>
#include <semstitf.h>
#include <atomic.h>

////////////////////////////////////////////////////////

void
PlaneSector_Processor::Process( Stack<AtomicRegion>& Offspring)
  {
  PlaneSector& G = Geometry();
  Point P1(G.InnerRadius(),G.BigAngle()), P2(G.InnerRadius(),G.SmallAngle());
  AtomicRegion* A ;
  A= (AtomicRegion*) SEMI_INFINITE_STRIP(P1,P2);
  Integrand I1(LocalIntegrand(),new Translation(G.Center()));
  A->LocalIntegrand(new Integrand(I1, new PolarToRectangular));
  Offspring.Push(A);
  }
/////////////////////////////////////////////////
PlaneSector_Processor::PlaneSector_Processor()
  :Processor<PlaneSector>()
  {
  }
/////////////////////////////////////////////////       
Processor<PlaneSector>*
PlaneSector_Processor::NewCopy()
const
  {
  return new PlaneSector_Processor(*this);
  }
////////////////////////////////////////////////
