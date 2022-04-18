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
/////////////////////////////////////////////////
//File regcoll.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//////////////////////////////////////////////////

#include <regcoll.h>
#include <error.h>

///////////////////////////////////////////////////
void
REGION_COLLECTION::LocalIntegrand(Integrand* ip)
  {
  SCR_ptr->IteratorReset();
  while(!SCR_ptr->IteratorAtEnd())
    {
    SCR_ptr->IteratorNext()->LocalIntegrand(ip);
    };
  }
///////////////////////////////////////////////////
void
REGION_COLLECTION::LocalIntegrand(Function f)
  {
  LocalIntegrand(new Integrand(f));
  }
///////////////////////////////////////////////////
REGION_COLLECTION&
REGION_COLLECTION::operator+= (const COMPOUND_REGION& C)
  {
  SCR_ptr->Push( C.NewCopy());
  return *this;
  }
//////////////////////////////////////////////////
REGION_COLLECTION
REGION_COLLECTION::operator+ (const COMPOUND_REGION& C)
  {
  REGION_COLLECTION r;
  r+= (*this);
  r+= C;
  return r;
  }
///////////////////////////////////////////////////
REGION_COLLECTION::REGION_COLLECTION()
  :COMPOUND_REGION()
  {
  SCR_ptr = new Stack<COMPOUND_REGION>;
  }
/////////////////////////////////////////////////
REGION_COLLECTION::REGION_COLLECTION(const REGION_COLLECTION& rc)
  :COMPOUND_REGION(rc),SCR_ptr(rc.SCR_ptr)
  {
  }
//////////////////////////////////////////////////////
void
REGION_COLLECTION::Preprocess()
  {
  SCR_ptr->IteratorReset();
  while(!SCR_ptr->IteratorAtEnd())
    {
    COMPOUND_REGION* C =SCR_ptr->IteratorNext();
    real OldIntegral = C->Integral();
    real OldError = C->AbsoluteError();
    C->Preprocess();
    LocalInfo().Integral() += C->Integral()-OldIntegral;
    LocalInfo().AbsoluteError ()+= C->AbsoluteError()-OldError;
    };
  // if all parts are hopeless then the result is
  // hopeless too:
  SCR_ptr->IteratorReset();
  while (!SCR_ptr->IteratorAtEnd())
    {
    COMPOUND_REGION* C = SCR_ptr->IteratorNext();
    if(!C->Hopeless()) return;
    };
  LocalInfo().Hopeless() = True;
  }
 ////////////////////////////////////////////////////
void
REGION_COLLECTION::Improve()
  {
  real MaxCompoundError = 0;
  COMPOUND_REGION* MaxCompoundErrorAtRegion  =0;
  SCR_ptr->IteratorReset();
  while (!SCR_ptr->IteratorAtEnd())
    {
    COMPOUND_REGION* C = SCR_ptr->IteratorNext();
    real e =C->MaxAtomicError();
    if (e > MaxCompoundError && !C->Hopeless())
      {
      MaxCompoundError = e;
      MaxCompoundErrorAtRegion = C;
      };
    };
  LocalInfo().Integral() -=
      MaxCompoundErrorAtRegion->Integral();
  LocalInfo().AbsoluteError() -=
      MaxCompoundErrorAtRegion->AbsoluteError();
  MaxCompoundErrorAtRegion->Process();
  LocalInfo().Integral() +=
      MaxCompoundErrorAtRegion->Integral();
  LocalInfo().AbsoluteError() +=
      MaxCompoundErrorAtRegion->AbsoluteError();
  // if all parts are hopeless then the result is
  // hopeless too:
  SCR_ptr->IteratorReset();
  while (!SCR_ptr->IteratorAtEnd())
    {
    COMPOUND_REGION* C = SCR_ptr->IteratorNext();
    if(!C->Hopeless()) return;
    };
  LocalInfo().Hopeless() = True;
  }
/////////////////////////////////////////////////////
real
REGION_COLLECTION::MaxAtomicError()
const
  {
  real MaxError=0;
  SCR_ptr->IteratorReset();
  while(!SCR_ptr->IteratorAtEnd())
    {
    real m = SCR_ptr->IteratorNext()->MaxAtomicError();
    if(m > MaxError)
      {
      MaxError = m;
      };
    };
  return MaxError;
 }
/////////////////////////////////////////////////////
COMPOUND_REGION*
REGION_COLLECTION::NewCopy()
const
  {
  return new REGION_COLLECTION(*this);
  }
///////////////////////////////////////////////////////
