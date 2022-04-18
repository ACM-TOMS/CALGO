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
//File C2dv4.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(long lines split)
/////////////////////////////////////////////////

#include <point.h>
#include <stack.h>
#include <C2dv4.h>

//////////////////////////////////////////////////
void
Parallelogram_Divide4::Apply (const Parallelogram& t,
  Stack<Parallelogram>& Offspring, const Vector<unsigned int>& DiffOrder)
  {
  Point m12 = (t.Vertex(0)+t.Vertex(1))/2;
  Point m13 = (t.Vertex(0)+t.Vertex(2))/2;
  Point m23 = (t.Vertex(2)+t.Vertex(1))/2;
  Point m24 = t.Vertex(1)+(t.Vertex(2)-t.Vertex(0))/2;
  Point m34 = t.Vertex(2)+(t.Vertex(1)-t.Vertex(0))/2;

  Parallelogram* t1 = new  Parallelogram(t.Vertex(0),m12,m13);
  Parallelogram* t2 = new  Parallelogram(m12,t.Vertex(1),m23);
  Parallelogram* t3 = new  Parallelogram(m13,m23,t.Vertex(2));
  Parallelogram* t4 = new  Parallelogram(m23,m24,m34);

  Offspring.Push(t1);
  Offspring.Push(t2);
  Offspring.Push(t3);
  Offspring.Push(t4);
  }
//////////////////////////////////////////////////
Parallelogram_Divide4::Parallelogram_Divide4()
  :SameShapeDivisor<Parallelogram>()
  {
  }
//////////////////////////////////////////////////
