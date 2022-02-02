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
//File polC2prc.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
#include <gsprc.h>
#include <C2toS2.h>
#include <translat.h>
#include <integran.h>
#include <gs.h>
#include <C2interf.h>
#include <atomic.h>

////////////////////////////////////////////////////////
// LOCAL CLASS
//////////////////////////////////////////
class VariableScaling : public Transformation
  {
  private: real (*B) (real);
  public:

  VariableScaling(real(*b)(real))
    :Transformation(),B(b)
    {};
  void Transform( real& w, Point& p)
    {
    p.R() = B(p.Theta())*p.R();
    w*= B(p.Theta());
    };
  };
///////////////////////////////////////////////////////////

void
GeneralizedSector_Processor::Process( Stack<AtomicRegion>& Offspring)
  {
  GeneralizedSector& G = Geometry();
  Point P1(0,G.Alpha()), P2(0,G.Beta()),
        P3(1,G.Alpha());
  AtomicRegion* A ;
  A= (AtomicRegion*) PARALLELOGRAM(P1,P2,P3);
  Integrand I1(LocalIntegrand(),new Translation(G.Center()));
  Integrand I2(I1, new PolarToRectangular);
  A->LocalIntegrand(new Integrand(I2, new VariableScaling(G.Boundary())));
  Offspring.Push(A);
  }
/////////////////////////////////////////////////
GeneralizedSector_Processor::GeneralizedSector_Processor()
  :Processor<GeneralizedSector>()
  {
  }
/////////////////////////////////////////////////       
Processor<GeneralizedSector>*
GeneralizedSector_Processor::NewCopy()
const
  {
  return new GeneralizedSector_Processor(*this);
  }
/////////////////////////////////////////////////////
