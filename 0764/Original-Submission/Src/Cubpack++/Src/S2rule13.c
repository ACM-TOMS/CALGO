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
//File S2rule13.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
//
// A cubature formula of degree 13 with 36 points and associated
// null rules (Cools & Haegemans).
//
#include <real.h>
static real facmed  = 1.5;
static real crival  = 0.3;
static int K[]={0,3,2,2};
const int NumberOfPoints=K[0]+4*K[1]+4*K[2]+8*K[3];
const int Orbits=K[0]+K[1]+K[2]+K[3];

static real Type1[]=
   {0.283402832348840340260055531925e+00,
    0.649007230578023857142954047869e+00,
    0.920575474036255978410400945265e+00};
static real Type2[]=
   {0.677355106028069632559870394597e+00,
    0.412672216763278882147845015245e+00};
static real Type3[2][2]=
  {{0.335767600829543545504321849367e+00,
    0.938694583835106683089433405007e+00},
   {0.752042776803937943743749637665e+00,
    0.379717011170079696840056017956e+00}};
static real Weight[7][7]=
  {{0.497326951852944153115255529183e-01,
    0.439732090591084987635001867639e-01,
    0.221057969395291114959209696751e-01,
    0.138071954801617944635783199052e-01,
    0.478191639416850908138653549756e-01,
    0.579302409092634509107760741653e-02,
    0.304879456061841994847272004644e-01},
//  Nullrule exact for the first 6 basispolynomials,
     {-0.8073879159066197e-02,
     -0.2498315518214109e-01,
      0.5349170539028515e-01,
      0.3595547716318644e-01,
      0.4320518775944201e-01,
     -0.1967607360410405e-01,
     -0.3012159438174910e-01},
//  Nullrule exact for the first 5 basispolynomials,
     {0.2175100120594341e-01,
     -0.7117820786684134e-01,
      0.1801191272610644e-01,
     -0.1775642346529020e-01,
     -0.2754333027079005e-02,
     -0.1170796597217573e-01,
      0.3767099118575607e-01},
//  Nullrule exact for the first 4 basispolynomials,
     {0.2413156434309756e-01,
      0.2327141313816948e-01,
      0.3438009174891747e-01,
      0.3084063431284037e-01,
     -0.6858211678704531e-01,
     -0.2561913343310581e-01,
      0.3598340055116032e-02},
//  Nullrule exact for the first 3 basispolynomials,
    {-0.4853627844062673e-01,
      0.3645578330411178e-01,
      0.4003943562394422e-01,
     -0.4494496660283691e-01,
      0.1266357696444992e-01,
     -0.2028566592818438e-01,
      0.2244689050366323e-01},
//  Nullrule exact for the first            2 basispolynomials,
    {-0.2979266039117054e-01,
      0.3977474716810517e-02,
     -0.4396076050818155e-01,
      0.5463380585542593e-01,
      0.1777067704103634e-01,
     -0.2876745680547196e-01,
      0.2745318844851159e-01},
//  Nullrule exact for the first            1 basispolynomials,
     {0.6246499865107039e-01,
      0.2624836264288300e-01,
     -0.1903648724133519e-01,
     -0.2649038633167001e-01,
      0.3481279195257287e-01,
     -0.3459307949232926e-01,
     -0.4406560344431272e-02}};
//
// The error estimator constants
//

#include <S2rule13.h>
#include <error.h>
#include <tools.h>
#include <math.h>
//#include <iostream.h>

#define sqr(x) ((x)*(x))

void Circle_Rule13::Apply(Integrand& F,Circle& R,real& TheResult,real& TheError)
  {
  //if (crival<0 || facmed <0)
   // {
   // cerr << " crival,facmed? " << endl;
   // cin >>  crival;
   // cin >> facmed;
   // };
const real facopt = facmed/(crival*crival);
  int i,j,p,Type,nr,number;
  real Null[6];
  real Tres = -1.0;
  real noise,r,r1,r2,Deg5,Deg3,Deg7,sumval;
  Point x[8];

  if (Tres <= 0)
  {
  Tres=50*REAL_EPSILON;
  }
  TheResult = 0.0;
  for (i=0;i<6;i++){Null[i]=0;}
  p=0;

  for (Type=0;Type <=3;Type++){
     for (nr=0;nr<K[Type];nr++){
            if (Type == 0) {
               number=1;
               x[0]=R.Center();
            }
            else if (Type == 1) {
               Point z1(R.Radius()*Type1[nr],0),
                     z2(0,R.Radius()*Type1[nr]);
               number=4;
               x[0] = R.Center() + z1;
               x[1] = R.Center() - z1;
               x[2] = R.Center() + z2;
               x[3] = R.Center() - z2;
            }
            else if (Type == 2) {
               Point z1(R.Radius()*Type2[nr],R.Radius()*Type2[nr]),
                     z2(-R.Radius()*Type2[nr],R.Radius()*Type2[nr]);
               number=4;
               x[0] = R.Center() + z1;
               x[1] = R.Center() - z1;
               x[2] = R.Center() + z2;
               x[3] = R.Center() - z2;
            }
            else{
               Point z1(R.Radius()*Type3[nr][0],R.Radius()*Type3[nr][1]),
                     z2(-R.Radius()*Type3[nr][0],R.Radius()*Type3[nr][1]),
                     z3(R.Radius()*Type3[nr][1],R.Radius()*Type3[nr][0]),
                     z4(-R.Radius()*Type3[nr][1],R.Radius()*Type3[nr][0]);
               number=8;
               x[0] = R.Center() + z1;
               x[1] = R.Center() - z1;
               x[2] = R.Center() + z2;
               x[3] = R.Center() - z2;
               x[4] = R.Center() + z3;
               x[5] = R.Center() - z3;
               x[6] = R.Center() + z4;
               x[7] = R.Center() - z4;
            }

            sumval = F(x[0]);
            for (j=1 ; j<number ; j++){ sumval += F(x[j]); }

            TheResult += Weight[0][p]*sumval;
            for (i=1 ; i<=6 ; i++){ Null[i-1] += Weight[i][p]*sumval ;}
            p=p+1;
     }
  }
//  cout << "Benadering voor cirkel = " <<TheResult*R.Volume()<<endl;
//  cout << "Nulregels:" <<endl;
//  for (i=0 ; i<=5 ; i++){cout<<Null[i]<<endl;}
//  cout <<"---------------------------"<<endl;
  //    Compute errors.
  noise = fabs(TheResult)*Tres;
  Deg7 = fabs(Null[0]);
//  cout <<"Deg7 = "<<Deg7<<endl;
  if (Deg7 <= noise)
  {
    TheError = noise;
  }
  else
  {
    Deg5 = sqrt(sqr(Null[1])+sqr(Null[2]));
    Deg3 = sqrt(sqr(Null[3])+sqr(Null[4]));
//    cout<<"Deg5 = "<<Deg5<<" Deg3 = "<<Deg3<<endl;
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
      r2 = Deg5/Deg3;
      }
    else
      {
      r2 = 1;
      }
//    cout<<"R = "<<r1<<" "<<r2<<endl;
    r = max(r1,r2);
    if (r>= 1)
      {
      TheError =10*Deg7;
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
  TheError *= R.Volume();
  TheResult *= R.Volume();
//  cout <<"=> Error = "<<TheError<<endl;
  }
 /// ///////////////////////////////////////////
Circle_Rule13::Circle_Rule13()
  :Rule<Circle>()
  {
  }
///////////////////////////////////////////
