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
//File T2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////

#include <T2.h>
#include <math.h>

//////////////////////////////////////////////

Triangle::Triangle (const Point& p1,const Point& p2,const Point& p3)
  : Geometry(2),Vertices(3),TheVolumeKnown(False)
  {
  Vertices[0] = p1;
  Vertices[1] = p2;
  Vertices[2] = p3;
  }

//////////////////////////////////////////////
void
Triangle::ComputeVolume()
  {
  TheVolume = fabs( 0.5*
          ((Vertices[0][0]*Vertices[1][1]-Vertices[1][0]*Vertices[0][1])
          +(Vertices[1][0]*Vertices[2][1]-Vertices[2][0]*Vertices[1][1])
          +(Vertices[2][0]*Vertices[0][1]-Vertices[0][0]*Vertices[2][1]))
        );
  TheVolumeKnown = True;
  }

///////////////////////////////////////////////
const Point&
Triangle::Vertex(int i)
const
  {
  return Vertices[i];
  }
//////////////////////////////////////////////////
real
Triangle::Volume()
const
  {
  Triangle* T= (Triangle*) this;
  if (!TheVolumeKnown)
    {
    T->ComputeVolume();
    };
  return  T->TheVolume;
  }
///////////////////////////////////////////////
void
Triangle::Volume(real v)
  {
  TheVolume = v;
  TheVolumeKnown = True;
  }
////////////////////////////////////////////////
//Processor<Triangle>*
//Triangle::DefaultProcessor()
//const
  //{
  //return new  SimpleAdaptive<Triangle>(
               //new Triangle_Rule13,
               //new Triangle_Divide4);
  //}
////////////////////////////////////////////////
