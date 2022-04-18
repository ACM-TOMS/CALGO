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
//////////////////////////////////////////////////
//File userint.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//    9 Sep 1994     V0.1a(operator added)
//   10 Sep 1994     V0.1b(initialisation of Int_ptr added)
//   25 Jan 1996     V0.1f(typedef introduced, code from .c to .h)
///////////////////////////////////////////////////
#include <userint.h>
#include <error.h>
///////////////////////////////////////////////////
template <class GEOMETRY>
USERINTERFACE<GEOMETRY>::USERINTERFACE()
  :COMPOUND_REGION()
  {
  SAR_ptr = new Stack<AtomicRegion>;
  HAR_ptr = new Heap<AtomicRegion>;
  HopelessAR_ptr = new Stack<AtomicRegion>;
  Int_ptr = (Pointer<Integrand>) (0);
  }
//////////////////////////////////////////////////
template <class GEOMETRY>
USERINTERFACE<GEOMETRY>::USERINTERFACE(const USERINTERFACE<GEOMETRY>& u)
  :COMPOUND_REGION(u)
  {
  SAR_ptr = u.SAR_ptr;
  HAR_ptr = u.HAR_ptr;
  HopelessAR_ptr = u.HopelessAR_ptr;
  Int_ptr = u.Int_ptr;
  }
///////////////////////////////////////////////////
template <class GEOMETRY>
USERINTERFACE<GEOMETRY>::~USERINTERFACE()
  {
  }
///////////////////////////////////////////////////
template <class GEOMETRY>
void
USERINTERFACE<GEOMETRY>::Use(Processor<GEOMETRY>* rp)
  {
  Error(SAR_ptr-> Size() != 1,
  "Attempt to specify region processor for multiple regions"
        );
  Atomic<GEOMETRY>* A =(Atomic<GEOMETRY>*)(SAR_ptr->Top());
  A->Use(rp);
  }
//////////////////////////////////////////////////
template <class GEOMETRY>
void
USERINTERFACE<GEOMETRY> ::StoreAtomic (GEOMETRY* g,
     Processor<GEOMETRY>* p)
  {
  Atomic<GEOMETRY>* A = new Atomic<GEOMETRY>(g,p);
  SAR_ptr->Push(A);
  }
//////////////////////////////////////////////////
//template <class GEOMETRY>
//void
//USERINTERFACE<GEOMETRY> ::StoreAtomic (GEOMETRY* g)
  //{
  //Atomic<GEOMETRY>* A = new Atomic<GEOMETRY>(g);
  //SAR_ptr->Push(A);
  //}
//////////////////////////////////////////////////
template <class GEOMETRY>
void
USERINTERFACE<GEOMETRY>::Preprocess()
  {
  while(!SAR_ptr->Empty())
    {
    AtomicRegion* A;
    A = SAR_ptr->Pop();
    Stack<AtomicRegion> Offspring;
    A->Process(Offspring);
    if( Offspring.Empty())
      {
      if (!A->Hopeless())
        {
        (*HAR_ptr) += A;
        }
      else
        {
        HopelessAR_ptr->Push( A);
        };
      LocalInfo().Integral() += A->Integral();
      LocalInfo().AbsoluteError ()   += A->AbsoluteError();
      }
    else
      {
      delete A;
      SAR_ptr->Merge(Offspring);
      };
    };
  if (HAR_ptr->Empty())
    {
    LocalInfo().Hopeless() = True;
    };
  }
/////////////////////////////////////////////////////
template <class GEOMETRY>
void
USERINTERFACE<GEOMETRY>::Improve()
  {
  // (*this) shouldn't be Hopeless()
    AtomicRegion* A = HAR_ptr->Get();
//HAR_ptr->Print();
//cout << A->Integral()<<" "<<A->AbsoluteError()<<endl;
    LocalInfo().Integral() -= A->Integral();
    LocalInfo().AbsoluteError() -= A->AbsoluteError();
    Stack<AtomicRegion> Offspring;
    A->Process(Offspring);
    if( Offspring.Empty())
      {
      if (!A->Hopeless())
        {
        (*HAR_ptr) += A;
        }
      else
        {
        HopelessAR_ptr->Push( A);
        };
      LocalInfo().Integral() += A->Integral();
      LocalInfo().AbsoluteError()    += A->AbsoluteError();
      }
    else
      {
      delete A;
      SAR_ptr->Merge(Offspring);
      Preprocess();
      };
  if (HAR_ptr->Empty())
    {
    LocalInfo().Hopeless() = True;
    };
  }
//////////////////////////////////////////////////////
template <class GEOMETRY>
void
USERINTERFACE<GEOMETRY>::LocalIntegrand(Integrand* ip)
  {
  if (Int_ptr == (Pointer<Integrand>) (0)) Int_ptr = ip; 
  else
    {
      Error(!(*Int_ptr == *ip),
            "Attempt to modify integrand during integration");
    }
//    Error(!HAR_ptr->Empty(),
//          "Attempt to modify integrand during integration");
      if ( !(SAR_ptr->Empty()) )
        { 
         SAR_ptr->IteratorReset();
         while(!SAR_ptr->IteratorAtEnd())
           {
           SAR_ptr->IteratorNext()->LocalIntegrand(ip);
           };
        }
  }
//////////////////////////////////////////////////////
template <class GEOMETRY>
void
USERINTERFACE<GEOMETRY>::LocalIntegrand(Function f)
  {
  LocalIntegrand(new Integrand(f));
  }
///////////////////////////////////////////////////////
template <class GEOMETRY>
REGION_COLLECTION
USERINTERFACE<GEOMETRY>::operator+(  const COMPOUND_REGION& C)
  {
  REGION_COLLECTION r;
  r += (*this);
  r += C;
  return r;
  }
//////////////////////////////////////////////////////
template <class GEOMETRY>
real
USERINTERFACE<GEOMETRY>::MaxAtomicError()
const
  {
  real MaxError=0;
  if (!HAR_ptr->Empty())
    {
    AtomicRegion* t = HAR_ptr->Look();
    MaxError = t->AbsoluteError();
    };
  return MaxError;
  }
////////////////////////////////////////////////////////
template <class GEOMETRY>
COMPOUND_REGION*
USERINTERFACE<GEOMETRY>::NewCopy()
const
  {
  return new USERINTERFACE<GEOMETRY>(*this);
  }
///////////////////////////////////////////////////////
//template <class GEOMETRY>
//USERINTERFACE<GEOMETRY>::operator AtomicRegion*()
//  {
//  Error (SAR_ptr->Empty(),
//    "converting empty compound region to atomic.");
//  return SAR_ptr->Pop();
//  }
////////////////////////////////////////////////////
