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
//
//File C2rule13.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//   25 Jan 1996     V0.1f(long lines split)
/////////////////////////////////////////////////////////
//
// A cubature formula of degree 13 with 37 points and associated
// null rules.
// The cubature formula, from Rabinowitz & Richter, was recomputed
// to obtain higher accuracy.
//
/////////////////////////////////////////////////////////
#include <real.h>
static int K[]={1,2,3,2};
const int NumberOfPoints=K[0]+4*K[1]+4*K[2]+8*K[3];
const int Orbits=K[0]+K[1]+K[2]+K[3];

static real Type1[]={ 0.9909890363004326469792722978603e+00,
                      0.6283940712305315063814483471116e+00};
static real l1sdl2s = Type1[1]*Type1[1]/Type1[0]/Type1[0];
static real Type2[]={ 0.9194861553393073086142137772149e+00,
                      0.6973201917871173078084506730937e+00,
                      0.3805687186904854497424188074662e+00};
static real Type3[2][2]={{0.9708504361720225062147290554088e+00,
                          0.6390348393207252159077623446225e+00},
                         {0.8623637916722781475018696425693e+00,
                          0.3162277660168700033875075593701e+00}};
static real
Weight[8][8]={{0.2995235559387052215463143056692e+00,
               0.3311006686692356205977471655046e-01,
               0.1802214941550624038355347399683e+00,
               0.3916727896035153300761243260674e-01,
               0.1387748348777288706306435595057e+00,
               0.2268881207335707037147066705814e+00,
               0.3657395765508995601240002438981e-01,
               0.1169047000557533546701746277951e+00},
              {7.610781847149629154716409791983e-2,
               1.486101247399760261471935168346e-1,
              -2.077685631717747007172983323970e-1,
               6.850758313011924198538315395405e-2,
               2.024205813317813585572881715385e-1,
               1.108627473745508429879249169864e-1,
              -1.187411393304862640859204217487e-1,
              -5.208857468077715683772080394959e-2},
              {4.016494861405949747097510013162e-2,
              -1.093132962444079541048635452881e-1,
              -2.270251673633777452624380129694e-1,
               1.231674163356097016086203579325e-2,
              -1.420402526499201540699111172200e-1,
               1.189080551229557928776504129312e-1,
              -4.482039658150474743804189300793e-3,
               1.730383808319875827592824151609e-1},
             {-5.643905795781771973971259866415e-1,
               2.878418073676293225652331648545e-2,
               1.159354231997583294689565314470e-1,
               1.376081498690624477894043101438e-1,
              -7.909780225340130915490382973570e-2,
               1.174335441429478112778176601234e-1,
              -1.107251942334134124782600707843e-1,
               2.094226883312045633400182488252e-2},
             {-2.269001713589584730602581694579e-1,
               2.976190892690301120078774620049e-2,
              -7.440193483272787588251423144751e-2,
              -1.224665989043784131260454301280e-1,
              -4.857910454732976198562745578156e-2,
               2.228157325962656425537280474671e-1,
               1.459764751457503859063666414952e-1,
              -1.211789553452468781539987084682e-1},
             {-3.326760468009974589269134283992e-1,
               1.796655319904795478676993902115e-1,
              -4.389976396805911868560791966472e-2,
              -2.295841771339316497310760908889e-1,
               6.182618387692816082856552878852e-2,
              -1.202703885325137746461829140891e-1,
               5.109536580363550180208564374234e-3,
               1.126062761533095493689566169969e-1},
              {2.290638530086106512999345512401e-1,
               2.702070398116919449911037051753e-1,
              -9.078047988731123605988441792069e-3,
               4.618480310858703283999169489655e-2,
              -2.598231009547631799096616255056e-1,
              -2.518433931146441037986247681820e-2,
              -1.257796993152456033984707367389e-2,
              -2.720818902721190304043617320910e-2},
              {2.746908885094872977794932213372e-1,
              -1.149427039769738298032807785523e-2,
               1.596178537820019535731955591283e-1,
              -2.180626972663360443142752377527e-1,
              -8.711748038292630173597899697063e-3,
               1.902786182960269617633915869710e-1,
              -1.189840649092108827784089292890e-1,
               2.883382565767354162177931122471e-2}};

//
// The error estimator constants
//
static const real crival=0.4 , facmed=8.0;
static const real facopt = facmed/(crival*crival);

#include <C2rule13.h>
#include <error.h>
#include <tools.h>
#include <math.h>
#define sqr(x) ((x)*(x))

void Parallelogram_Rule13::ApplyWithDiffs(Integrand& F,Parallelogram& R,
  real& TheResult,real& TheError, Vector<real>& DiffOrder)
  {
  int i,j,p,Type,nr,number;
  real Null[7];
  real Tres = -1.0;
  real noise,r,r1,r2,r3,Deg5,Deg3,Deg1,Deg7,z1,z2,sumval;
  real D1,D2;
  real F0;
  real FE[2][2][2];
  //     [l1,l2][vertex1,vertex2][plus,minus]
  Point x[8];

  if (Tres <= 0)
  {
  Tres=50*REAL_EPSILON;
  }
  TheResult = 0.0;
  for (i=0;i<7;i++){Null[i]=0;}
  p=0;

  for (Type=0;Type <=3;Type++){
     for (nr=0;nr<K[Type];nr++){
            if (Type == 0) {
               number=1;
               x[0]=(R.Vertex(1)+R.Vertex(2))/2;
            }
            else if (Type == 1) {
               z1=Type1[nr];
               number=4;
               x[0] = (-R.Vertex(0)*z1+R.Vertex(1)*(1+z1)+R.Vertex(2))/2;
               x[1] = (-R.Vertex(0)*z1 +R.Vertex(1)+R.Vertex(2)*(z1+1))/2;
               x[2] = (R.Vertex(0)*z1 +R.Vertex(1)+R.Vertex(2)*(-z1+1))/2;
               x[3] = (R.Vertex(0)*z1+R.Vertex(1)*(1-z1)+R.Vertex(2))/2;
            }
            else if (Type == 2) {
               z1=Type2[nr];
               number=4;
         x[1]=(-2*R.Vertex(0)*z1+R.Vertex(1)*(1+z1)+R.Vertex(2)*(z1+1))/2;
          x[2]=(2*R.Vertex(0)*z1+R.Vertex(1)*(1-z1)+R.Vertex(2)*(1-z1))/2;
               x[3]=(R.Vertex(1)*(1+z1) + R.Vertex(2)*(1-z1))/2;
               x[0]=(R.Vertex(1)*(1-z1) + R.Vertex(2)*(1+z1))/2;
            }
            else{
               z1=Type3[nr][0];
               z2=Type3[nr][1] ;
               number=8;
      x[1]=(R.Vertex(0)*(-z1-z2)+R.Vertex(1)*(1+z2)+R.Vertex(2)*(1+z1))/2;
      x[2]=(R.Vertex(0)*(+z1-z2)+R.Vertex(1)*(1+z2)+R.Vertex(2)*(1-z1))/2;
      x[3]=(R.Vertex(0)*(-z1+z2)+R.Vertex(1)*(1-z2)+R.Vertex(2)*(1+z1))/2;
      x[4]=(R.Vertex(0)*(+z1+z2)+R.Vertex(1)*(1-z2)+R.Vertex(2)*(1-z1))/2;
      x[5]=(R.Vertex(0)*(-z1-z2)+R.Vertex(1)*(1+z1)+R.Vertex(2)*(1+z2))/2;
      x[6]=(R.Vertex(0)*(-z1+z2)+R.Vertex(1)*(1+z1)+R.Vertex(2)*(1-z2))/2;
      x[7]=(R.Vertex(0)*(+z1-z2)+R.Vertex(1)*(1-z1)+R.Vertex(2)*(1+z2))/2;
      x[0]=(R.Vertex(0)*(+z1+z2)+R.Vertex(1)*(1-z1)+R.Vertex(2)*(1-z2))/2;
            }

            if(Type==0)
              {
              sumval = F0 = F(x[0]);
              }
            else if (Type == 1)
              {
              sumval = FE[nr][0][0] = F(x[0]);
              sumval += FE[nr][0][1] = F(x[3]);
              sumval += FE[nr][1][0] = F(x[1]);
              sumval += FE[nr][1][1] = F(x[2]);
              }
            else
              {
            sumval = F(x[0]);
            for (j=1 ; j<number ; j++){ sumval += F(x[j]); }
              }

            TheResult += Weight[0][p]*sumval;
            for (i=1 ; i<=7 ; i++){ Null[i-1] += Weight[i][p]*sumval ;}
            p=p+1;
     }
  }
  //    Compute differences
  D1 = ( FE[1][0][0] +  FE[1][0][1] -2*F0)
      - l1sdl2s * ( FE[0][0][0] +  FE[0][0][1] -2*F0);
  D2 = ( FE[1][1][0] +  FE[1][1][1] -2*F0)
      - l1sdl2s * ( FE[0][1][0] +  FE[0][1][1] -2*F0);
  DiffOrder[0] = fabs(D1);
  DiffOrder[1] = fabs(D2);

  //    Compute errors.
  noise = fabs(TheResult)*Tres;
  Deg7 = sqrt(sqr(Null[0])+sqr(Null[1]));
  if (Deg7 <= noise)
  {
    TheError = noise;
  }
  else
  {
    Deg5 = sqrt(sqr(Null[2])+sqr(Null[3]));
    Deg3 = sqrt(sqr(Null[4])+sqr(Null[5]));
    Deg1 = sqrt(sqr(Null[6])+sqr(Null[5]));
    if(Deg5 != 0)
      {
      r1=Deg7/Deg5;
      }
    else
      {
      r1 = 1;
      }
    if(Deg3 != 0)
      {
      r2=Deg5/Deg3;
      }
    else
      {
      r2 = 1;
      }
    if(Deg1 != 0)
      {
      r3=Deg3/Deg1;
      }
    else
      {
      r3 = 1;
      }
    r=max(r1,max(r2,r3));
    if (r>= 1)
      {
      TheError =facmed*Deg7;
      }
    else
      {
      if (r>=crival)
        {
        TheError = facmed*r*Deg7;
        }
      else
        {
        TheError = facopt *r*r*r*Deg7;
        }
      }
    TheError = max(noise,TheError);
    }
  TheError *= R.Volume()/4;
  TheResult *=R.Volume()/4;
  }
//////////////////////////////////////////////////
Parallelogram_Rule13::Parallelogram_Rule13()
  :Rule<Parallelogram>()
  {
  }
//////////////////////////////////////////////////
