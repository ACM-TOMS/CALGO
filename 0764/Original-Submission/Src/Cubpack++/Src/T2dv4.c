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
//File T2dv4.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
/////////////////////////////////////////////////

#include <point.h>
#include <stack.h>
#include <T2dv4.h>

//////////////////////////////////////////////////
void
Triangle_Divide4::Apply (const Triangle& t,Stack<Triangle>& Offspring,
                         const Vector<unsigned int>& D)
  {

  Point m12 = (t.Vertex(0)+t.Vertex(1))/2;
  Point m13 = (t.Vertex(0)+t.Vertex(2))/2;
  Point m23 = (t.Vertex(2)+t.Vertex(1))/2;


  Triangle* t1 =  new Triangle(t.Vertex(0),m12,m13);
  Triangle* t2 =  new Triangle(t.Vertex(1),m12,m23);
  Triangle* t3 =  new Triangle(t.Vertex(2),m23,m13);
  Triangle* t4 =  new Triangle(m23,m12,m13);
  Offspring.Push(t1);
  Offspring.Push(t2);
  Offspring.Push(t3);
  Offspring.Push(t4);
  }
//////////////////////////////////////////////////
Triangle_Divide4::Triangle_Divide4()
  :SameShapeDivisor<Triangle>()
  {
  }
//////////////////////////////////////////////////
