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
// File : templist.h
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////////////
// PURPOSE:
//   to enumerate the compilers that need inclusion of
//   template implementations
/////////////////////////////////////////////////////////
#ifndef TEMPLIST_H
#define TEMPLIST_H

#ifdef __GNUG__
#define TEMPLATEINCLUDE
#endif

#ifdef __TCPLUSPLUS__
#define TEMPLATEINCLUDE
#endif

#ifdef __DECCXX
#define TEMPLATEINCLUDE
#endif

#endif
