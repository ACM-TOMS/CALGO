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
//File C2dv2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(long lines split)
/////////////////////////////////////////////////

#include <C2.h>
#include <point.h>
#include <stack.h>
#include <C2dv2.h>

//////////////////////////////////////////////////
void
Parallelogram_Divide2::Apply (const Parallelogram& t,
  Stack<Parallelogram>& Offspring, 
  const Vector<unsigned int>& DiffOrder)
  {
  Point m12 = (t.Vertex(0)+t.Vertex(DiffOrder[0]))/2;
  Point m34 = t.Vertex(DiffOrder[1])+(t.Vertex(DiffOrder[0])-t.Vertex(0))/2;

  Parallelogram* t1 = new  
    Parallelogram(t.Vertex(0),m12,t.Vertex(DiffOrder[1]));
  Parallelogram* t2 = new  Parallelogram(m12,t.Vertex(DiffOrder[0]),m34);

  Offspring.Push(t1);
  Offspring.Push(t2);
  }
//////////////////////////////////////////////////
Parallelogram_Divide2::Parallelogram_Divide2()
  :SameShapeDivisor<Parallelogram>()
  {
  }
//////////////////////////////////////////////////
