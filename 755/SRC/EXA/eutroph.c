/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include "adouble.h"
#include "adutils.h"
void eutroph(unsigned short tag, double* px,  double* pxp)
{
double IK=0.11 ;
double FRZ=0.3 ;
double EFFUZ=0.6;
double PRITZ=1.0e-3;
double RESP=5.0e-3;
double sinK=5.0e-3;
double PRITA=0.1;
double RZ=1.0e-2;
double K2=4.0e-2;
double K3=5.0e-1;
double KSP=2.0e2;
double KSF=1.0;
double BETA=100.0/1.25;
double ALPHA=0.002;
double TRZ=2;
double EPSP = 0.4;
double FI1 = 230.4;
double FI3=282.8;
double FI4=127.5;
double FI5=141.9;
double p = 40.0;
double DEPTH = 45;
/************* fix controls ********/
double PRFOS = 0.5*p;
double M = 0.1;
double ZMIX = (45+RZ)/2;
double QIV=0.297E-02/3;
/******initialize adoubles *************/
adouble x[5],xp[5];
int i;
trace_on(tag);
for(i=0;i<4;i++)
  x[i]<<= px[i];
adouble T;
T <<= px[4];
xp[4] = 1;
double tdum=0.0;
adouble TEMP=9.5+7.9*sin(T+FI1);
adouble FOTOP = 12.0+4.19*sin(T+280.0);
adouble I=229.0+215.0*sin(T+FI3)+15.3*sin(2.0*T+FI4)+ 21.7*sin(3.0*T+FI5);
adouble PIDI=.8+.25*cos(T)-.12*cos(2.*T);
double MORITZ = 0.075;
double Q = 0.786E6;
double VND = 0.265E9;
double V = VND;
if (T<72) I *= 0.603;
adouble EPS = ALPHA * x[0] + x[3] + EPSP;
adouble temp = I * exp(-EPS*ZMIX);
adouble temp2 = 2*IK*FOTOP;
adouble GROW;
GROW = 1.2*FOTOP/EPS/ZMIX * (1.333 * atan ( I / temp2 )
 -IK*FOTOP / I * log( 1 + pow( (I /temp2 ),2) ) 
 -1.333 * atan ( temp / temp2)
 +IK*FOTOP/temp* log( 1+pow(temp/temp2, 2) )) * x[2] /(KSF+x[2]) 
       * 0.366 * pow(K2,0.52) * exp(0.09*TEMP) * pow(x[0],(1-0.52));
xp[0] = GROW - RESP * TEMP * x[0] - FRZ * x[0] * x[1] - sinK * PIDI * x[0] 
          + (PRITA - x[0]) * Q/VND;
xp[1] = FRZ * x[0] / K2 * x[1] / 1000 * EFFUZ*KSP / KSP+x[0] 
             - RZ * x[1] - MORITZ * x[1] + (PRITZ - x[1] ) * Q/V;
xp[2] = K3 * (-GROW + RESP * TEMP * x[0] + FRZ * x[0] * x[1] *
              (1 - EFFUZ*KSP /(KSP+x[0]) ) + RZ * K2 * 1000 *
              x[1] + MORITZ * K2 * 1000 * x[1] ) + (PRFOS - x[2])* Q/V;
xp[3] = (- x[3] * Q  + BETA * M / TRZ)/VND;
for (i=0;i<4;i++)
  xp[i] >>= pxp[i];
xp[4] >>= tdum;
trace_off();
}

void eutroph(double* px,  double* pxp)
{
double IK=0.11 ;
double FRZ=0.3 ;
double EFFUZ=0.6;
double PRITZ=1.0e-3;
double RESP=5.0e-3;
double sinK=5.0e-3;
double PRITA=0.1;
double RZ=1.0e-2;
double K2=4.0e-2;
double K3=5.0e-1;
double KSP=2.0e2;
double KSF=1.0;
double BETA=100.0/1.25;
double ALPHA=0.002;
double TRZ=2;
double EPSP = 0.4;
double FI1 = 230.4;
double FI3=282.8;
double FI4=127.5;
double FI5=141.9;
double p = 40.0;
double DEPTH = 45;
/************* fix controls ********/
double PRFOS = 0.5*p;
double M = 0.1;
double ZMIX = (45+RZ)/2;
double QIV=0.297E-02/3;
/******initialize doubles *************/
double x[5],xp[5];
int i;
for(i=0;i<4;i++)
  x[i]= px[i];
double T;
T = px[4];
xp[4] = 1;
double TEMP=9.5+7.9*sin(T+FI1);
double FOTOP = 12.0+4.19*sin(T+280.0);
double I=229.0+215.0*sin(T+FI3)+15.3*sin(2.0*T+FI4)+ 21.7*sin(3.0*T+FI5);
double PIDI=.8+.25*cos(T)-.12*cos(2.*T);
double MORITZ = 0.075;
double Q = 0.786E6;
double VND = 0.265E9;
double V = VND;
if (T<72) I *= 0.603;
double EPS = ALPHA * x[0] + x[3] + EPSP;
double temp = I * exp(-EPS*ZMIX);
double temp2 = 2*IK*FOTOP;
double GROW;
GROW = 1.2*FOTOP/EPS/ZMIX * (1.333 * atan ( I / temp2 )
 -IK*FOTOP / I * log( 1 + pow( (I /temp2 ),2) ) 
 -1.333 * atan ( temp / temp2)
 +IK*FOTOP/temp* log( 1+pow(temp/temp2, 2) )) * x[2] /(KSF+x[2]) 
       * 0.366 * pow(K2,0.52) * exp(0.09*TEMP) * pow(x[0],(1-0.52));
xp[0] = GROW - RESP * TEMP * x[0] - FRZ * x[0] * x[1] - sinK * PIDI * x[0] 
          + (PRITA - x[0]) * Q/VND;
xp[1] = FRZ * x[0] / K2 * x[1] / 1000 * EFFUZ*KSP / KSP+x[0] 
             - RZ * x[1] - MORITZ * x[1] + (PRITZ - x[1] ) * Q/V;
xp[2] = K3 * (-GROW + RESP * TEMP * x[0] + FRZ * x[0] * x[1] *
              (1 - EFFUZ*KSP /(KSP+x[0]) ) + RZ * K2 * 1000 *
              x[1] + MORITZ * K2 * 1000 * x[1] ) + (PRFOS - x[2])* Q/V;
xp[3] = (- x[3] * Q  + BETA * M / TRZ)/VND;
for (i=0;i<5;i++)
  pxp[i] = xp[i];
}
