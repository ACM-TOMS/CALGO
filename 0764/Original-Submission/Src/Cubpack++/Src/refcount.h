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
// File : refcount.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// DEFINITION OF CLASS ReferenceCounting
// -------------------------------------------------
//
// BASECLASSES:
//   None
//
// PURPOSE:
//   provides a means to count references to an object.
//   it is to be used preferably together with
//   class Pointer<>, which automatically calls
//   Refer() and Unrefer() wherever necessary.
//
// TEMPLATES:
//   None
//
// METHODS:
//   CONSTRUCTORS:
//     1) ReferenceCounting()
//     ----------------------
//     initialises the reference count to zero.
//
//     2) ReferenceCounting(const ReferenceCounting&)
//     ----------------------------------------------
//     also initialises the count to zero
//
//   SELECTORS:
//     1) unsigned int NumberOfReferences() const
//     ------------------------------------------
//
//   MODIFIERS:
//     1) void Refer()
//     ---------------
//     increments the reference count.
//
//     2) void UnRefer()
//     -----------------
//     decrements the reference count
//
//   OPERATORS:
//     1) ReferenceCounting& operator=
//                          (const ReferenceCounting&)
//     -------------------------------------------------
//     resets the reference count.
//
//   SPECIAL:
//     None
//
/////////////////////////////////////////////////////////
#ifndef  REFCOUNT_H
#define  REFCOUNT_H
///////////////////////////////////////////////////

class ReferenceCounting
  {
  public:

  ReferenceCounting();
  ReferenceCounting(const ReferenceCounting&);
  void Refer();
  void UnRefer();
  ReferenceCounting& operator=(const ReferenceCounting&);
  unsigned int  NumberOfReferences()const;

  private:

  unsigned int numref;
  };

////////////////////////////////////////////////
#endif
