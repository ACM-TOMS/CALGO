#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

const size_t kirby2_N = 151;
const size_t kirby2_P = 5;

/* double kirby2_x0[5] = { 2, -0.1, 0.003, -0.001, 0.00001 }; */

MpIeee kirby2_x0[5] =  { MpIeee( "1.5" ), -MpIeee( "0.15" ), MpIeee( "0.0025" ), -MpIeee( "0.0015" ), MpIeee( "0.00002" ) }; 

MpIeee kirby2_x[5] =  {
  MpIeee( "1.6745063063E+00" ),
  -MpIeee( "1.3927397867E-01" ),
  MpIeee( "2.5961181191E-03" ),
  -MpIeee( "1.7241811870E-03" ),
  MpIeee( "2.1664802578E-05" )
};

MpIeee kirby2_sumsq=  MpIeee( "3.9050739624E+00" );

MpIeee kirby2_sigma[5] =  {
  MpIeee( "8.7989634338E-02" ),
  MpIeee( "4.1182041386E-03" ),
  MpIeee( "4.1856520458E-05" ),
  MpIeee( "5.8931897355E-05" ),
  MpIeee( "2.0129761919E-07" )
};

MpIeee kirby2_F1[151] =  {
       MpIeee( "0.0082E0" ),
       MpIeee( "0.0112E0" ),
       MpIeee( "0.0149E0" ),
       MpIeee( "0.0198E0" ),
       MpIeee( "0.0248E0" ),
       MpIeee( "0.0324E0" ),
       MpIeee( "0.0420E0" ),
       MpIeee( "0.0549E0" ),
       MpIeee( "0.0719E0" ),
       MpIeee( "0.0963E0" ),
       MpIeee( "0.1291E0" ),
       MpIeee( "0.1710E0" ),
       MpIeee( "0.2314E0" ),
       MpIeee( "0.3227E0" ),
       MpIeee( "0.4809E0" ),
       MpIeee( "0.7084E0" ),
       MpIeee( "1.0220E0" ),
       MpIeee( "1.4580E0" ),
       MpIeee( "1.9520E0" ),
       MpIeee( "2.5410E0" ),
       MpIeee( "3.2230E0" ),
       MpIeee( "3.9990E0" ),
       MpIeee( "4.8520E0" ),
       MpIeee( "5.7320E0" ),
       MpIeee( "6.7270E0" ),
       MpIeee( "7.8350E0" ),
       MpIeee( "9.0250E0" ),
      MpIeee( "10.2670E0" ),
      MpIeee( "11.5780E0" ),
      MpIeee( "12.9440E0" ),
      MpIeee( "14.3770E0" ),
      MpIeee( "15.8560E0" ),
      MpIeee( "17.3310E0" ),
      MpIeee( "18.8850E0" ),
      MpIeee( "20.5750E0" ),
      MpIeee( "22.3200E0" ),
      MpIeee( "22.3030E0" ),
      MpIeee( "23.4600E0" ),
      MpIeee( "24.0600E0" ),
      MpIeee( "25.2720E0" ),
      MpIeee( "25.8530E0" ),
      MpIeee( "27.1100E0" ),
      MpIeee( "27.6580E0" ),
      MpIeee( "28.9240E0" ),
      MpIeee( "29.5110E0" ),
      MpIeee( "30.7100E0" ),
      MpIeee( "31.3500E0" ),
      MpIeee( "32.5200E0" ),
      MpIeee( "33.2300E0" ),
      MpIeee( "34.3300E0" ),
      MpIeee( "35.0600E0" ),
      MpIeee( "36.1700E0" ),
      MpIeee( "36.8400E0" ),
      MpIeee( "38.0100E0" ),
      MpIeee( "38.6700E0" ),
      MpIeee( "39.8700E0" ),
      MpIeee( "40.0300E0" ),
      MpIeee( "40.5000E0" ),
      MpIeee( "41.3700E0" ),
      MpIeee( "41.6700E0" ),
      MpIeee( "42.3100E0" ),
      MpIeee( "42.7300E0" ),
      MpIeee( "43.4600E0" ),
      MpIeee( "44.1400E0" ),
      MpIeee( "44.5500E0" ),
      MpIeee( "45.2200E0" ),
      MpIeee( "45.9200E0" ),
      MpIeee( "46.3000E0" ),
      MpIeee( "47.0000E0" ),
      MpIeee( "47.6800E0" ),
      MpIeee( "48.0600E0" ),
      MpIeee( "48.7400E0" ),
      MpIeee( "49.4100E0" ),
      MpIeee( "49.7600E0" ),
      MpIeee( "50.4300E0" ),
      MpIeee( "51.1100E0" ),
      MpIeee( "51.5000E0" ),
      MpIeee( "52.1200E0" ),
      MpIeee( "52.7600E0" ),
      MpIeee( "53.1800E0" ),
      MpIeee( "53.7800E0" ),
      MpIeee( "54.4600E0" ),
      MpIeee( "54.8300E0" ),
      MpIeee( "55.4000E0" ),
      MpIeee( "56.4300E0" ),
      MpIeee( "57.0300E0" ),
      MpIeee( "58.0000E0" ),
      MpIeee( "58.6100E0" ),
      MpIeee( "59.5800E0" ),
      MpIeee( "60.1100E0" ),
      MpIeee( "61.1000E0" ),
      MpIeee( "61.6500E0" ),
      MpIeee( "62.5900E0" ),
      MpIeee( "63.1200E0" ),
      MpIeee( "64.0300E0" ),
      MpIeee( "64.6200E0" ),
      MpIeee( "65.4900E0" ),
      MpIeee( "66.0300E0" ),
      MpIeee( "66.8900E0" ),
      MpIeee( "67.4200E0" ),
      MpIeee( "68.2300E0" ),
      MpIeee( "68.7700E0" ),
      MpIeee( "69.5900E0" ),
      MpIeee( "70.1100E0" ),
      MpIeee( "70.8600E0" ),
      MpIeee( "71.4300E0" ),
      MpIeee( "72.1600E0" ),
      MpIeee( "72.7000E0" ),
      MpIeee( "73.4000E0" ),
      MpIeee( "73.9300E0" ),
      MpIeee( "74.6000E0" ),
      MpIeee( "75.1600E0" ),
      MpIeee( "75.8200E0" ),
      MpIeee( "76.3400E0" ),
      MpIeee( "76.9800E0" ),
      MpIeee( "77.4800E0" ),
      MpIeee( "78.0800E0" ),
      MpIeee( "78.6000E0" ),
      MpIeee( "79.1700E0" ),
      MpIeee( "79.6200E0" ),
      MpIeee( "79.8800E0" ),
      MpIeee( "80.1900E0" ),
      MpIeee( "80.6600E0" ),
      MpIeee( "81.2200E0" ),
      MpIeee( "81.6600E0" ),
      MpIeee( "82.1600E0" ),
      MpIeee( "82.5900E0" ),
      MpIeee( "83.1400E0" ),
      MpIeee( "83.5000E0" ),
      MpIeee( "84.0000E0" ),
      MpIeee( "84.4000E0" ),
      MpIeee( "84.8900E0" ),
      MpIeee( "85.2600E0" ),
      MpIeee( "85.7400E0" ),
      MpIeee( "86.0700E0" ),
      MpIeee( "86.5400E0" ),
      MpIeee( "86.8900E0" ),
      MpIeee( "87.3200E0" ),
      MpIeee( "87.6500E0" ),
      MpIeee( "88.1000E0" ),
      MpIeee( "88.4300E0" ),
      MpIeee( "88.8300E0" ),
      MpIeee( "89.1200E0" ),
      MpIeee( "89.5400E0" ),
      MpIeee( "89.8500E0" ),
      MpIeee( "90.2500E0" ),
      MpIeee( "90.5500E0" ),
      MpIeee( "90.9300E0" ),
      MpIeee( "91.2000E0" ),
      MpIeee( "91.5500E0" ),
      MpIeee( "92.2000E0" )
};


MpIeee kirby2_F0[151] =  {
    MpIeee( "9.65E0" ),
   MpIeee( "10.74E0" ),
   MpIeee( "11.81E0" ),
   MpIeee( "12.88E0" ),
   MpIeee( "14.06E0" ),
   MpIeee( "15.28E0" ),
   MpIeee( "16.63E0" ),
   MpIeee( "18.19E0" ),
   MpIeee( "19.88E0" ),
   MpIeee( "21.84E0" ),
   MpIeee( "24.00E0" ),
   MpIeee( "26.25E0" ),
   MpIeee( "28.86E0" ),
   MpIeee( "31.85E0" ),
   MpIeee( "35.79E0" ),
   MpIeee( "40.18E0" ),
   MpIeee( "44.74E0" ),
   MpIeee( "49.53E0" ),
   MpIeee( "53.94E0" ),
   MpIeee( "58.29E0" ),
   MpIeee( "62.63E0" ),
   MpIeee( "67.03E0" ),
   MpIeee( "71.25E0" ),
   MpIeee( "75.22E0" ),
   MpIeee( "79.33E0" ),
   MpIeee( "83.56E0" ),
   MpIeee( "87.75E0" ),
   MpIeee( "91.93E0" ),
   MpIeee( "96.10E0" ),
  MpIeee( "100.28E0" ),
  MpIeee( "104.46E0" ),
  MpIeee( "108.66E0" ),
  MpIeee( "112.71E0" ),
  MpIeee( "116.88E0" ),
  MpIeee( "121.33E0" ),
  MpIeee( "125.79E0" ),
  MpIeee( "125.79E0" ),
  MpIeee( "128.74E0" ),
  MpIeee( "130.27E0" ),
  MpIeee( "133.33E0" ),
  MpIeee( "134.79E0" ),
  MpIeee( "137.93E0" ),
  MpIeee( "139.33E0" ),
  MpIeee( "142.46E0" ),
  MpIeee( "143.90E0" ),
  MpIeee( "146.91E0" ),
  MpIeee( "148.51E0" ),
  MpIeee( "151.41E0" ),
  MpIeee( "153.17E0" ),
  MpIeee( "155.97E0" ),
  MpIeee( "157.76E0" ),
  MpIeee( "160.56E0" ),
  MpIeee( "162.30E0" ),
  MpIeee( "165.21E0" ),
  MpIeee( "166.90E0" ),
  MpIeee( "169.92E0" ),
  MpIeee( "170.32E0" ),
  MpIeee( "171.54E0" ),
  MpIeee( "173.79E0" ),
  MpIeee( "174.57E0" ),
  MpIeee( "176.25E0" ),
  MpIeee( "177.34E0" ),
  MpIeee( "179.19E0" ),
  MpIeee( "181.02E0" ),
  MpIeee( "182.08E0" ),
  MpIeee( "183.88E0" ),
  MpIeee( "185.75E0" ),
  MpIeee( "186.80E0" ),
  MpIeee( "188.63E0" ),
  MpIeee( "190.45E0" ),
  MpIeee( "191.48E0" ),
  MpIeee( "193.35E0" ),
  MpIeee( "195.22E0" ),
  MpIeee( "196.23E0" ),
  MpIeee( "198.05E0" ),
  MpIeee( "199.97E0" ),
  MpIeee( "201.06E0" ),
  MpIeee( "202.83E0" ),
  MpIeee( "204.69E0" ),
  MpIeee( "205.86E0" ),
  MpIeee( "207.58E0" ),
  MpIeee( "209.50E0" ),
  MpIeee( "210.65E0" ),
  MpIeee( "212.33E0" ),
  MpIeee( "215.43E0" ),
  MpIeee( "217.16E0" ),
  MpIeee( "220.21E0" ),
  MpIeee( "221.98E0" ),
  MpIeee( "225.06E0" ),
  MpIeee( "226.79E0" ),
  MpIeee( "229.92E0" ),
  MpIeee( "231.69E0" ),
  MpIeee( "234.77E0" ),
  MpIeee( "236.60E0" ),
  MpIeee( "239.63E0" ),
  MpIeee( "241.50E0" ),
  MpIeee( "244.48E0" ),
  MpIeee( "246.40E0" ),
  MpIeee( "249.35E0" ),
  MpIeee( "251.32E0" ),
  MpIeee( "254.22E0" ),
  MpIeee( "256.24E0" ),
  MpIeee( "259.11E0" ),
  MpIeee( "261.18E0" ),
  MpIeee( "264.02E0" ),
  MpIeee( "266.13E0" ),
  MpIeee( "268.94E0" ),
  MpIeee( "271.09E0" ),
  MpIeee( "273.87E0" ),
  MpIeee( "276.08E0" ),
  MpIeee( "278.83E0" ),
  MpIeee( "281.08E0" ),
  MpIeee( "283.81E0" ),
  MpIeee( "286.11E0" ),
  MpIeee( "288.81E0" ),
  MpIeee( "291.08E0" ),
  MpIeee( "293.75E0" ),
  MpIeee( "295.99E0" ),
  MpIeee( "298.64E0" ),
  MpIeee( "300.84E0" ),
  MpIeee( "302.02E0" ),
  MpIeee( "303.48E0" ),
  MpIeee( "305.65E0" ),
  MpIeee( "308.27E0" ),
  MpIeee( "310.41E0" ),
  MpIeee( "313.01E0" ),
  MpIeee( "315.12E0" ),
  MpIeee( "317.71E0" ),
  MpIeee( "319.79E0" ),
  MpIeee( "322.36E0" ),
  MpIeee( "324.42E0" ),
  MpIeee( "326.98E0" ),
  MpIeee( "329.01E0" ),
  MpIeee( "331.56E0" ),
  MpIeee( "333.56E0" ),
  MpIeee( "336.10E0" ),
  MpIeee( "338.08E0" ),
  MpIeee( "340.60E0" ),
  MpIeee( "342.57E0" ),
  MpIeee( "345.08E0" ),
  MpIeee( "347.02E0" ),
  MpIeee( "349.52E0" ),
  MpIeee( "351.44E0" ),
  MpIeee( "353.93E0" ),
  MpIeee( "355.83E0" ),
  MpIeee( "358.32E0" ),
  MpIeee( "360.20E0" ),
  MpIeee( "362.67E0" ),
  MpIeee( "364.53E0" ),
  MpIeee( "367.00E0" ),
  MpIeee( "371.30E0" )
};


int
 kirby2_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee b[5];
  size_t i;

  for (i = 0; i < 5; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < 151; i++)
    {
      MpIeee x=  kirby2_F0[i];
      MpIeee y=  ((b[0] + x* (b[1]  + x * b[2]))
                  / (MpIeee( "1" ) + x*(b[3]  + x *b[4])));
      gsl_vector_set (f, i, kirby2_F1[i] - y);
    }

  return GSL_SUCCESS;
}

int
 kirby2_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee b[5];
  size_t i;

  for (i = 0; i < 5; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < 151; i++)
    {
      MpIeee x=  kirby2_F0[i];
      MpIeee u=  (b[0] + x*(b[1] + x*b[2]));
      MpIeee v=  (MpIeee( "1" ) + x*(b[3] + x*b[4]));
      gsl_matrix_set (df, i, 0, -1/v);
      gsl_matrix_set (df, i, 1, -x/v);
      gsl_matrix_set (df, i, 2, -x*x/v);
      gsl_matrix_set (df, i, 3, x*u/(v*v));
      gsl_matrix_set (df, i, 4, x*x*u/(v*v));
    }

  return GSL_SUCCESS;
}

int
 kirby2_fdf(const gsl_vector * x, void *params,
           gsl_vector * f, gsl_matrix * df)
{
  kirby2_f (x, params, f);
  kirby2_df (x, params, df);

  return GSL_SUCCESS;
}
