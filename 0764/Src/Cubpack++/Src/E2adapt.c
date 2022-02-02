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
//File E2.cpp
// History:
//   (date)          (version)
//   19 Aug 1994     V0.1 (first limited distribution)
////////////////////////////////////////////////
#include <math.h>
#include <real.h>
#include <E2adapt.h>
#include <tools.h>
#include <invert.h>
#include <atomic.h>
#include <S2interf.h>
#include <outS2itf.h>
///////////////////////////////////////////////
#define sqr(x) ((x)*(x))
static int K[]={0,3,2,2};
const int NumberOfPoints=K[0]+4*K[1]+4*K[2]+8*K[3];
const int Orbits=K[0]+K[1]+K[2]+K[3];

static real Type1[]=
       {0.5460844152585403666984417e+00 ,
        0.1376486679696350642323050e+01 ,
        0.2358932261980680967828022e+01};
static real Type2[]=
      {0.8462499884480199068814807e+00 ,
       0.1985531339917884990890422e+01};
static real Type3[2][2]=
     {{0.1779246738905319006374908e+01 ,
       0.9318240227761998500475049e+00},
      {0.3030436065368579151664045e+01 ,
       0.9962534261727632104611162e+00}};

static real OrigWeight[7]=
      {0.1418036987406628465242303e+00,
       0.3807850175899473939094384e-01,
       0.1290856883945041937607418e-02,
       0.5772571826901623901546648e-01,
       0.1620590503146543183518239e-03,
       0.5445472363593369147167489e-02,
       0.2411028493987025953260209e-04};

static int d[]={0,3,1,5,2,4,6};

static real Weight[7][7]=
{{0.6002713316445080080142126,
0.7955990316249109629216728,
1.0584888443837658614605899,
0.7595383160846588054989698,
1.3523576929325600784468914,
0.9663531601396869556259704,
1.9895758469520118113052969},
{-0.1714120818975139,
-1.4261615934658296,
20.9491066359021694,
1.7029699666146943,
57.2013795938251143,
-4.3332227154742827,
-175.0624620496226842},
{-0.2883800695213331,
2.5450105556734605,
2.8699527274413806,
0.304671504643251,
138.1913876850488233,
-8.5573830821954439,
-18.7930973279112966},
{-0.250456836419066,
0.3273926337592694,
-83.9551965137636766,
0.4609375440071548,
-365.1974538111883157,
10.4127651350263225,
1537.3605596124529503},
{-0.2595143739988262,
0.1255865390202144,
85.4652862386879934,
-0.2489288730337109,
-712.5513360902265439,
6.3018929219370941,
-645.6559424743654782},
{-0.24820543230912,
-0.2143965818640334,
8.5337065366470588,
-0.2238393816932916,
867.6556565044116888,
7.6084037588845301,
-4712.804729866979117},
{0.3052968227487744,
1.0011209710514362,
-6.3357177267870899,
0.7227067892680122,
-358.0249035190603089,
2.7469614313645364,
-6407.2910382311609777}};

// Radii of rings associated with weights of published formula
static real RingRadius[]=
{0.9151577658134376,
1.2649388389767648,
1.7333214530567136,
2.2616892861784437,
2.6609731353556634,
2.9246248487342039,
6.0979478618219725};  //This is `infinity'. Never used actually

//
// The error estimator constants
//
const real crival=0.5 , facmed=5.0;
const real facopt = facmed/(crival*crival);

//////////////////////////////////////////////////////
void
PlaneAdaptive::Rule(Integrand& F,Plane& R,real& TheResult,real& TheError,real& HalfValueRadius)
  {
  int i,j,p,Type,nr,number;
  real Null[6],OrbitContrib[7],below,above,fraction;
  real Tres = 50*REAL_EPSILON;
  real noise,r,r1,r2,r3,Deg5,Deg3,Deg7,Deg1,sumval;
  real A=1.0/R.ScaleX() , B=1.0/R.ScaleY();
  Point x[8];

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
               Point z1(A*Type1[nr],0),
                     z2(0,B*Type1[nr]);
               number=4;
               x[0] = R.Center() + z1;
               x[1] = R.Center() - z1;
               x[2] = R.Center() + z2;
               x[3] = R.Center() - z2;
            }
            else if (Type == 2) {
               Point z1(A*Type2[nr],B*Type2[nr]),
                     z2(-A*Type2[nr],B*Type2[nr]);
               number=4;
               x[0] = R.Center() + z1;
               x[1] = R.Center() - z1;
               x[2] = R.Center() + z2;
               x[3] = R.Center() - z2;
            }
            else{
               Point z1(A*Type3[nr][0],B*Type3[nr][1]),
                     z2(-A*Type3[nr][0],B*Type3[nr][1]),
                     z3(A*Type3[nr][1],B*Type3[nr][0]),
                     z4(-A*Type3[nr][1],B*Type3[nr][0]);
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

            OrbitContrib[p] = Weight[0][p]*sumval;
            TheResult += Weight[0][p]*sumval;
            for (i=1 ; i<=6 ; i++){ Null[i-1] += Weight[i][p]*sumval ;}
            p=p+1;
     }
  }
  i = -1;   above = 0;
  do{
     i++;
     below = above;
     above +=OrbitContrib[d[i]];
    } while (above < TheResult*0.5);
    // 50% circle is lying between orbits d[i-1] and d[i]
    fraction = (TheResult*0.5-below)/OrbitContrib[d[i]];
    if (d[i] < K[1]+K[2])
       {
        fraction *= 4*OrigWeight[d[i]];
       }
    else
       {
        fraction *= 8*OrigWeight[d[i]];
       }
    if (i == 0)
      {
       HalfValueRadius = sqrt(-log(1-fraction));
      }
    else
      {
       HalfValueRadius =
          sqrt(-log(exp(-RingRadius[i-1]*RingRadius[i-1]) - fraction));
      }

//  cout<<"Bijdrage per orbiet (toenemende afstand):"<<endl;
//  for(i=0; i<=6 ; i++)
//    {cout<<d[i]<<"  "<<OrbitContrib[d[i]]<<" "<<Radius[d[i]]<<endl;}
//  cout<<"HalfValueRadius = "<<HalfValueRadius<<"  past orbit "<<i<<endl;
//  cout << "Benadering voor vlak= " <<TheResult<<endl;
//  cout << "Nulregels:" <<endl;
//  for (i=0 ; i<=5 ; i++){cout<<Null[i]<<endl;}
//  cout <<"---------------------------"<<endl;
  //    Compute errors.
  noise = fabs(TheResult)*Tres;
  Deg7 = fabs(Null[0]);
  if (Deg7 <= noise)
  {
    TheError = noise;
  }
  else
  {
    Deg5 = sqrt(sqr(Null[1])+sqr(Null[2]));
    Deg3 = sqrt(sqr(Null[3])+sqr(Null[4]));
    Deg1 = fabs(Null[5]);
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
      r2=sqrt(Deg7/Deg3);
      }
    else
      {
      r2 = 1;
      }
    if(Deg1 != 0)
      {
      r3=sqrt(sqrt(Deg7/Deg1));
      }
    else
      {
      r3 = 1;
      }
//    cout<<"R = "<<r1<<" "<<r2<<" "<<r3<<endl;
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
  TheResult *= (A*B);
  TheError *= (A*B);
//  cout <<"=> Error = "<<TheError<<endl;
  }

//////////////////////////////////////////////
void
PlaneAdaptive::Process(Stack<AtomicRegion>& Offspring)
  {
  TimesCalled++;
  if (TimesCalled==1)
    {
    Rule(LocalIntegrand(),Geometry(),Integral(),AbsoluteError(),HalfValueRadius);
    }
  else
    {
    const Point& C=Geometry().Center();
    Offspring.Push(( AtomicRegion*) CIRCLE(C,HalfValueRadius));
    Offspring.Push(( AtomicRegion*) OUT_CIRCLE(C,HalfValueRadius));
    Offspring.IteratorReset();
    while (!Offspring.IteratorAtEnd())
      {
      Offspring.IteratorNext()->LocalIntegrand(&LocalIntegrand());
      };
    }
  }
///////////////////////////////////////////////////
PlaneAdaptive::PlaneAdaptive()
  :Processor<Plane>(),TimesCalled(0),HalfValueRadius(0)
  {
  }
////////////////////////////////////////////////////

Processor<Plane>*
PlaneAdaptive::NewCopy()
const
  {
  return new PlaneAdaptive(*this);
  }
//////////////////////////////////////////////////////
