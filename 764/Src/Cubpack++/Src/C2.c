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
//File C2.c
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////

#include <C2.h>
#include <math.h>

//////////////////////////////////////////////

Parallelogram::Parallelogram (const Point& p1,const Point& p2,const Point& p3)
  : Geometry(2),Vertices(3),TheVolumeKnown(False)
  {
  Vertices[0] = p1;
  Vertices[1] = p2;
  Vertices[2] = p3;
  }

//////////////////////////////////////////////
void
Parallelogram::ComputeVolume()
  {
  TheVolume = fabs(
          ((Vertices[0][0]*Vertices[1][1]-Vertices[1][0]*Vertices[0][1])
          +(Vertices[1][0]*Vertices[2][1]-Vertices[2][0]*Vertices[1][1])
          +(Vertices[2][0]*Vertices[0][1]-Vertices[0][0]*Vertices[2][1]))
        );
  TheVolumeKnown = True;
  }

///////////////////////////////////////////////
const Point&
Parallelogram::Vertex(int i)
const
  {
  return Vertices[i];
  }
//////////////////////////////////////////////////
real
Parallelogram::Volume()
const
  {
  Parallelogram* P = (Parallelogram*) this;
  if (!TheVolumeKnown)
    {
    P->ComputeVolume();
    };
  return  P->TheVolume;
  }
///////////////////////////////////////////////
void
Parallelogram::Volume(real v)
  {
  TheVolume = v;
  TheVolumeKnown = True;
  }
////////////////////////////////////////////////
//Processor<Parallelogram>*
//Parallelogram::DefaultProcessor()
//const
  //{
  //return new SimpleAdaptive<Parallelogram>
                //(new Parallelogram_Rule13,
                 //new Parallelogram_Divide4);
  //}
///////////////////////////////////////////////////
