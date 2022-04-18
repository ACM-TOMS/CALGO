#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/fermi_dirac.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_fermi_dirac.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#define locEPS  (1000.0*GSL_DBL_EPSILON)


/* Chebyshev fit for F_{1}(t);  -1 < t < 1, -1 < x < 1
 */
static MpIeee fd_1_a_data[22] =  {
  MpIeee( "1.8949340668482264365" ),
  MpIeee( "0.7237719066890052793" ),
  MpIeee( "0.1250000000000000000" ),
  MpIeee( "0.0101065196435973942" ),
  MpIeee( "0.0" ),
 -MpIeee( "0.0000600615242174119" ),
  MpIeee( "0.0" ),
  MpIeee( "6.816528764623e-7" ),
  MpIeee( "0.0" ),
 -MpIeee( "9.5895779195e-9" ),
  MpIeee( "0.0" ),
  MpIeee( "1.515104135e-10" ),
  MpIeee( "0.0" ),
 -MpIeee( "2.5785616e-12" ),
  MpIeee( "0.0" ),
  MpIeee( "4.62270e-14" ),
  MpIeee( "0.0" ),
 -MpIeee( "8.612e-16" ),
  MpIeee( "0.0" ),
  MpIeee( "1.65e-17" ),
  MpIeee( "0.0" ),
 -MpIeee( "3.e-19" )
};
static cheb_series fd_1_a_cs = {
  fd_1_a_data,
  21,
  -1, 1,
  12
};


/* Chebyshev fit for F_{1}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
 */
static MpIeee fd_1_b_data[22] =  {
  MpIeee( "10.409136795234611872" ),
  MpIeee( "3.899445098225161947" ),
  MpIeee( "0.513510935510521222" ),
  MpIeee( "0.010618736770218426" ),
 -MpIeee( "0.001584468020659694" ),
  MpIeee( "0.000146139297161640" ),
 -MpIeee( "1.408095734499e-6" ),
 -MpIeee( "2.177993899484e-6" ),
  MpIeee( "3.91423660640e-7" ),
 -MpIeee( "2.3860262660e-8" ),
 -MpIeee( "4.138309573e-9" ),
  MpIeee( "1.283965236e-9" ),
 -MpIeee( "1.39695990e-10" ),
 -MpIeee( "4.907743e-12" ),
  MpIeee( "4.399878e-12" ),
 -MpIeee( "7.17291e-13" ),
  MpIeee( "2.4320e-14" ),
  MpIeee( "1.4230e-14" ),
 -MpIeee( "3.446e-15" ),
  MpIeee( "2.93e-16" ),
  MpIeee( "3.7e-17" ),
 -MpIeee( "1.6e-17" )
};
static cheb_series fd_1_b_cs = {
  fd_1_b_data,
  21,
  -1, 1,
  11
};


/* Chebyshev fit for F_{1}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
 */
static MpIeee fd_1_c_data[23] =  {
  MpIeee( "56.78099449124299762" ),
  MpIeee( "21.00718468237668011" ),
  MpIeee( "2.24592457063193457" ),
  MpIeee( "0.00173793640425994" ),
 -MpIeee( "0.00058716468739423" ),
  MpIeee( "0.00016306958492437" ),
 -MpIeee( "0.00003817425583020" ),
  MpIeee( "7.64527252009e-6" ),
 -MpIeee( "1.31348500162e-6" ),
  MpIeee( "1.9000646056e-7" ),
 -MpIeee( "2.141328223e-8" ),
  MpIeee( "1.23906372e-9" ),
  MpIeee( "2.1848049e-10" ),
 -MpIeee( "1.0134282e-10" ),
  MpIeee( "2.484728e-11" ),
 -MpIeee( "4.73067e-12" ),
  MpIeee( "7.3555e-13" ),
 -MpIeee( "8.740e-14" ),
  MpIeee( "4.85e-15" ),
  MpIeee( "1.23e-15" ),
 -MpIeee( "5.6e-16" ),
  MpIeee( "1.4e-16" ),
 -MpIeee( "3.e-17" )
};
static cheb_series fd_1_c_cs = {
  fd_1_c_data,
  22,
  -1, 1,
  13
};


/* Chebyshev fit for F_{1}(x) / x^2
 * 10 < x < 30 
 * -1 < t < 1
 * t = 1/10 (x-10) - 1 = x/10 - 2
 * x = 10(t+2)
 */
static MpIeee fd_1_d_data[30] =  {
  MpIeee( "1.0126626021151374442" ),
 -MpIeee( "0.0063312525536433793" ),
  MpIeee( "0.0024837319237084326" ),
 -MpIeee( "0.0008764333697726109" ),
  MpIeee( "0.0002913344438921266" ),
 -MpIeee( "0.0000931877907705692" ),
  MpIeee( "0.0000290151342040275" ),
 -MpIeee( "8.8548707259955e-6" ),
  MpIeee( "2.6603474114517e-6" ),
 -MpIeee( "7.891415690452e-7" ),
  MpIeee( "2.315730237195e-7" ),
 -MpIeee( "6.73179452963e-8" ),
  MpIeee( "1.94048035606e-8" ),
 -MpIeee( "5.5507129189e-9" ),
  MpIeee( "1.5766090896e-9" ),
 -MpIeee( "4.449310875e-10" ),
  MpIeee( "1.248292745e-10" ),
 -MpIeee( "3.48392894e-11" ),
  MpIeee( "9.6791550e-12" ),
 -MpIeee( "2.6786240e-12" ),
  MpIeee( "7.388852e-13" ),
 -MpIeee( "2.032828e-13" ),
  MpIeee( "5.58115e-14" ),
 -MpIeee( "1.52987e-14" ),
  MpIeee( "4.1886e-15" ),
 -MpIeee( "1.1458e-15" ),
  MpIeee( "3.132e-16" ),
 -MpIeee( "8.56e-17" ),
  MpIeee( "2.33e-17" ),
 -MpIeee( "5.9e-18" )
};
static cheb_series fd_1_d_cs = {
  fd_1_d_data,
  29,
  -1, 1,
  14
};


/* Chebyshev fit for F_{1}(x) / x^2
 * 30 < x < Inf
 * -1 < t < 1
 * t = 60/x - 1
 * x = 60/(t+1)
 */
static MpIeee fd_1_e_data[10] =  {
  MpIeee( "1.0013707783890401683" ),
  MpIeee( "0.0009138522593601060" ),
  MpIeee( "0.0002284630648400133" ),
 -MpIeee( "1.57e-17" ),
 -MpIeee( "1.27e-17" ),
 -MpIeee( "9.7e-18" ),
 -MpIeee( "6.9e-18" ),
 -MpIeee( "4.6e-18" ),
 -MpIeee( "2.9e-18" ),
 -MpIeee( "1.7e-18" )
};
static cheb_series fd_1_e_cs = {
  fd_1_e_data,
  9,
  -1, 1,
  4
};


/* Chebyshev fit for F_{2}(t);  -1 < t < 1, -1 < x < 1
 */
static MpIeee fd_2_a_data[21] =  {
  MpIeee( "2.1573661917148458336" ),
  MpIeee( "0.8849670334241132182" ),
  MpIeee( "0.1784163467613519713" ),
  MpIeee( "0.0208333333333333333" ),
  MpIeee( "0.0012708226459768508" ),
  MpIeee( "0.0" ),
 -MpIeee( "5.0619314244895e-6" ),
  MpIeee( "0.0" ),
  MpIeee( "4.32026533989e-8" ),
  MpIeee( "0.0" ),
 -MpIeee( "4.870544166e-10" ),
  MpIeee( "0.0" ),
  MpIeee( "6.4203740e-12" ),
  MpIeee( "0.0" ),
 -MpIeee( "9.37424e-14" ),
  MpIeee( "0.0" ),
  MpIeee( "1.4715e-15" ),
  MpIeee( "0.0" ),
 -MpIeee( "2.44e-17" ),
  MpIeee( "0.0" ),
  MpIeee( "4.e-19" )
};
static cheb_series fd_2_a_cs = {
  fd_2_a_data,
  20,
  -1, 1,
  12
};


/* Chebyshev fit for F_{2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
 */
static MpIeee fd_2_b_data[22] =  {
  MpIeee( "16.508258811798623599" ),
  MpIeee( "7.421719394793067988" ),
  MpIeee( "1.458309885545603821" ),
  MpIeee( "0.128773850882795229" ),
  MpIeee( "0.001963612026198147" ),
 -MpIeee( "0.000237458988738779" ),
  MpIeee( "0.000018539661382641" ),
 -MpIeee( "1.92805649479e-7" ),
 -MpIeee( "2.01950028452e-7" ),
  MpIeee( "3.2963497518e-8" ),
 -MpIeee( "1.885817092e-9" ),
 -MpIeee( "2.72632744e-10" ),
  MpIeee( "8.0554561e-11" ),
 -MpIeee( "8.313223e-12" ),
 -MpIeee( "2.24489e-13" ),
  MpIeee( "2.18778e-13" ),
 -MpIeee( "3.4290e-14" ),
  MpIeee( "1.225e-15" ),
  MpIeee( "5.81e-16" ),
 -MpIeee( "1.37e-16" ),
  MpIeee( "1.2e-17" ),
  MpIeee( "1.e-18" )
};
static cheb_series fd_2_b_cs = {
  fd_2_b_data,
  21,
  -1, 1,
  12
};


/* Chebyshev fit for F_{1}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
 */
static MpIeee fd_2_c_data[20] =  {
  MpIeee( "168.87129776686440711" ),
  MpIeee( "81.80260488091659458" ),
  MpIeee( "15.75408505947931513" ),
  MpIeee( "1.12325586765966440" ),
  MpIeee( "0.00059057505725084" ),
 -MpIeee( "0.00016469712946921" ),
  MpIeee( "0.00003885607810107" ),
 -MpIeee( "7.89873660613e-6" ),
  MpIeee( "1.39786238616e-6" ),
 -MpIeee( "2.1534528656e-7" ),
  MpIeee( "2.831510953e-8" ),
 -MpIeee( "2.94978583e-9" ),
  MpIeee( "1.6755082e-10" ),
  MpIeee( "2.234229e-11" ),
 -MpIeee( "1.035130e-11" ),
  MpIeee( "2.41117e-12" ),
 -MpIeee( "4.3531e-13" ),
  MpIeee( "6.447e-14" ),
 -MpIeee( "7.39e-15" ),
  MpIeee( "4.3e-16" )
};
static cheb_series fd_2_c_cs = {
  fd_2_c_data,
  19,
  -1, 1,
  12
};


/* Chebyshev fit for F_{1}(x) / x^3
 * 10 < x < 30 
 * -1 < t < 1
 * t = 1/10 (x-10) - 1 = x/10 - 2
 * x = 10(t+2)
 */
static MpIeee fd_2_d_data[30] =  {
  MpIeee( "0.3459960518965277589" ),
 -MpIeee( "0.00633136397691958024" ),
  MpIeee( "0.00248382959047594408" ),
 -MpIeee( "0.00087651191884005114" ),
  MpIeee( "0.00029139255351719932" ),
 -MpIeee( "0.00009322746111846199" ),
  MpIeee( "0.00002904021914564786" ),
 -MpIeee( "8.86962264810663e-6" ),
  MpIeee( "2.66844972574613e-6" ),
 -MpIeee( "7.9331564996004e-7" ),
  MpIeee( "2.3359868615516e-7" ),
 -MpIeee( "6.824790880436e-8" ),
  MpIeee( "1.981036528154e-8" ),
 -MpIeee( "5.71940426300e-9" ),
  MpIeee( "1.64379426579e-9" ),
 -MpIeee( "4.7064937566e-10" ),
  MpIeee( "1.3432614122e-10" ),
 -MpIeee( "3.823400534e-11" ),
  MpIeee( "1.085771994e-11" ),
 -MpIeee( "3.07727465e-12" ),
  MpIeee( "8.7064848e-13" ),
 -MpIeee( "2.4595431e-13" ),
  MpIeee( "6.938531e-14" ),
 -MpIeee( "1.954939e-14" ),
  MpIeee( "5.50162e-15" ),
 -MpIeee( "1.54657e-15" ),
  MpIeee( "4.3429e-16" ),
 -MpIeee( "1.2178e-16" ),
  MpIeee( "3.394e-17" ),
 -MpIeee( "8.81e-18" )
};
static cheb_series fd_2_d_cs = {
  fd_2_d_data,
  29,
  -1, 1,
  14
};


/* Chebyshev fit for F_{2}(x) / x^3
 * 30 < x < Inf
 * -1 < t < 1
 * t = 60/x - 1
 * x = 60/(t+1)
 */
static MpIeee fd_2_e_data[4] =  {
  MpIeee( "0.3347041117223735227" ),
  MpIeee( "0.00091385225936012645" ),
  MpIeee( "0.00022846306484003205" ),
  MpIeee( "5.2e-19" )
};
static cheb_series fd_2_e_cs = {
  fd_2_e_data,
  3,
  -1, 1,
  3
};


/* Chebyshev fit for F_{-1/2}(t);  -1 < t < 1, -1 < x < 1
 */
static MpIeee fd_mhalf_a_data[20] =  {
  MpIeee( "1.2663290042859741974" ),
  MpIeee( "0.3697876251911153071" ),
  MpIeee( "0.0278131011214405055" ),
 -MpIeee( "0.0033332848565672007" ),
 -MpIeee( "0.0004438108265412038" ),
  MpIeee( "0.0000616495177243839" ),
  MpIeee( "8.7589611449897e-6" ),
 -MpIeee( "1.2622936986172e-6" ),
 -MpIeee( "1.837464037221e-7" ),
  MpIeee( "2.69495091400e-8" ),
  MpIeee( "3.9760866257e-9" ),
 -MpIeee( "5.894468795e-10" ),
 -MpIeee( "8.77321638e-11" ),
  MpIeee( "1.31016571e-11" ),
  MpIeee( "1.9621619e-12" ),
 -MpIeee( "2.945887e-13" ),
 -MpIeee( "4.43234e-14" ),
  MpIeee( "6.6816e-15" ),
  MpIeee( "1.0084e-15" ),
 -MpIeee( "1.561e-16" )
};
static cheb_series fd_mhalf_a_cs = {
  fd_mhalf_a_data,
  19,
  -1, 1,
  12
};


/* Chebyshev fit for F_{-1/2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
 */
static MpIeee fd_mhalf_b_data[20] =  {
  MpIeee( "3.270796131942071484" ),
  MpIeee( "0.5809004935853417887" ),
 -MpIeee( "0.0299313438794694987" ),
 -MpIeee( "0.0013287935412612198" ),
  MpIeee( "0.0009910221228704198" ),
 -MpIeee( "0.0001690954939688554" ),
  MpIeee( "6.5955849946915e-6" ),
  MpIeee( "3.5953966033618e-6" ),
 -MpIeee( "9.430672023181e-7" ),
  MpIeee( "8.75773958291e-8" ),
  MpIeee( "1.06247652607e-8" ),
 -MpIeee( "4.9587006215e-9" ),
  MpIeee( "7.160432795e-10" ),
  MpIeee( "4.5072219e-12" ),
 -MpIeee( "2.3695425e-11" ),
  MpIeee( "4.9122208e-12" ),
 -MpIeee( "2.905277e-13" ),
 -MpIeee( "9.59291e-14" ),
  MpIeee( "3.00028e-14" ),
 -MpIeee( "3.4970e-15" )
};
static cheb_series fd_mhalf_b_cs = {
  fd_mhalf_b_data,
  19,
  -1, 1,
  12
};


/* Chebyshev fit for F_{-1/2}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
 */
static MpIeee fd_mhalf_c_data[25] =  {
  MpIeee( "5.828283273430595507" ),
  MpIeee( "0.677521118293264655" ),
 -MpIeee( "0.043946248736481554" ),
  MpIeee( "0.005825595781828244" ),
 -MpIeee( "0.000864858907380668" ),
  MpIeee( "0.000110017890076539" ),
 -MpIeee( "6.973305225404e-6" ),
 -MpIeee( "1.716267414672e-6" ),
  MpIeee( "8.59811582041e-7" ),
 -MpIeee( "2.33066786976e-7" ),
  MpIeee( "4.8503191159e-8" ),
 -MpIeee( "8.130620247e-9" ),
  MpIeee( "1.021068250e-9" ),
 -MpIeee( "5.3188423e-11" ),
 -MpIeee( "1.9430559e-11" ),
  MpIeee( "8.750506e-12" ),
 -MpIeee( "2.324897e-12" ),
  MpIeee( "4.83102e-13" ),
 -MpIeee( "8.1207e-14" ),
  MpIeee( "1.0132e-14" ),
 -MpIeee( "4.64e-16" ),
 -MpIeee( "2.24e-16" ),
  MpIeee( "9.7e-17" ),
 -MpIeee( "2.6e-17" ),
  MpIeee( "5.e-18" )
};
static cheb_series fd_mhalf_c_cs = {
  fd_mhalf_c_data,
  24,
  -1, 1,
  13
};


/* Chebyshev fit for F_{-1/2}(x) / x^(1/2)
 * 10 < x < 30 
 * -1 < t < 1
 * t = 1/10 (x-10) - 1 = x/10 - 2
 */
static MpIeee fd_mhalf_d_data[30] =  {
  MpIeee( "2.2530744202862438709" ),
  MpIeee( "0.0018745152720114692" ),
 -MpIeee( "0.0007550198497498903" ),
  MpIeee( "0.0002759818676644382" ),
 -MpIeee( "0.0000959406283465913" ),
  MpIeee( "0.0000324056855537065" ),
 -MpIeee( "0.0000107462396145761" ),
  MpIeee( "3.5126865219224e-6" ),
 -MpIeee( "1.1313072730092e-6" ),
  MpIeee( "3.577454162766e-7" ),
 -MpIeee( "1.104926666238e-7" ),
  MpIeee( "3.31304165692e-8" ),
 -MpIeee( "9.5837381008e-9" ),
  MpIeee( "2.6575790141e-9" ),
 -MpIeee( "7.015201447e-10" ),
  MpIeee( "1.747111336e-10" ),
 -MpIeee( "4.04909605e-11" ),
  MpIeee( "8.5104999e-12" ),
 -MpIeee( "1.5261885e-12" ),
  MpIeee( "1.876851e-13" ),
  MpIeee( "1.00574e-14" ),
 -MpIeee( "1.82002e-14" ),
  MpIeee( "8.6634e-15" ),
 -MpIeee( "3.2058e-15" ),
  MpIeee( "1.0572e-15" ),
 -MpIeee( "3.259e-16" ),
  MpIeee( "9.60e-17" ),
 -MpIeee( "2.74e-17" ),
  MpIeee( "7.6e-18" ),
 -MpIeee( "1.9e-18" )
};
static cheb_series fd_mhalf_d_cs = {
  fd_mhalf_d_data,
  29,
  -1, 1,
  15
};


/* Chebyshev fit for F_{1/2}(t);  -1 < t < 1, -1 < x < 1
 */
static MpIeee fd_half_a_data[23] =  {
  MpIeee( "1.7177138871306189157" ),
  MpIeee( "0.6192579515822668460" ),
  MpIeee( "0.0932802275119206269" ),
  MpIeee( "0.0047094853246636182" ),
 -MpIeee( "0.0004243667967864481" ),
 -MpIeee( "0.0000452569787686193" ),
  MpIeee( "5.2426509519168e-6" ),
  MpIeee( "6.387648249080e-7" ),
 -MpIeee( "8.05777004848e-8" ),
 -MpIeee( "1.04290272415e-8" ),
  MpIeee( "1.3769478010e-9" ),
  MpIeee( "1.847190359e-10" ),
 -MpIeee( "2.51061890e-11" ),
 -MpIeee( "3.4497818e-12" ),
  MpIeee( "4.784373e-13" ),
  MpIeee( "6.68828e-14" ),
 -MpIeee( "9.4147e-15" ),
 -MpIeee( "1.3333e-15" ),
  MpIeee( "1.898e-16" ),
  MpIeee( "2.72e-17" ),
 -MpIeee( "3.9e-18" ),
 -MpIeee( "6.e-19" ),
  MpIeee( "1.e-19" )
};
static cheb_series fd_half_a_cs = {
  fd_half_a_data,
  22,
  -1, 1,
  11
};


/* Chebyshev fit for F_{1/2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
 */
static MpIeee fd_half_b_data[20] =  {
  MpIeee( "7.651013792074984027" ),
  MpIeee( "2.475545606866155737" ),
  MpIeee( "0.218335982672476128" ),
 -MpIeee( "0.007730591500584980" ),
 -MpIeee( "0.000217443383867318" ),
  MpIeee( "0.000147663980681359" ),
 -MpIeee( "0.000021586361321527" ),
  MpIeee( "8.07712735394e-7" ),
  MpIeee( "3.28858050706e-7" ),
 -MpIeee( "7.9474330632e-8" ),
  MpIeee( "6.940207234e-9" ),
  MpIeee( "6.75594681e-10" ),
 -MpIeee( "3.10200490e-10" ),
  MpIeee( "4.2677233e-11" ),
 -MpIeee( "2.1696e-14" ),
 -MpIeee( "1.170245e-12" ),
  MpIeee( "2.34757e-13" ),
 -MpIeee( "1.4139e-14" ),
 -MpIeee( "3.864e-15" ),
  MpIeee( "1.202e-15" )
};
static cheb_series fd_half_b_cs = {
  fd_half_b_data,
  19,
  -1, 1,
  12
};


/* Chebyshev fit for F_{1/2}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
 */
static MpIeee fd_half_c_data[23] =  {
  MpIeee( "29.584339348839816528" ),
  MpIeee( "8.808344283250615592" ),
  MpIeee( "0.503771641883577308" ),
 -MpIeee( "0.021540694914550443" ),
  MpIeee( "0.002143341709406890" ),
 -MpIeee( "0.000257365680646579" ),
  MpIeee( "0.000027933539372803" ),
 -MpIeee( "1.678525030167e-6" ),
 -MpIeee( "2.78100117693e-7" ),
  MpIeee( "1.35218065147e-7" ),
 -MpIeee( "3.3740425009e-8" ),
  MpIeee( "6.474834942e-9" ),
 -MpIeee( "1.009678978e-9" ),
  MpIeee( "1.20057555e-10" ),
 -MpIeee( "6.636314e-12" ),
 -MpIeee( "1.710566e-12" ),
  MpIeee( "7.75069e-13" ),
 -MpIeee( "1.97973e-13" ),
  MpIeee( "3.9414e-14" ),
 -MpIeee( "6.374e-15" ),
  MpIeee( "7.77e-16" ),
 -MpIeee( "4.0e-17" ),
 -MpIeee( "1.4e-17" )
};
static cheb_series fd_half_c_cs = {
  fd_half_c_data,
  22,
  -1, 1,
  13
};


/* Chebyshev fit for F_{1/2}(x) / x^(3/2)
 * 10 < x < 30 
 * -1 < t < 1
 * t = 1/10 (x-10) - 1 = x/10 - 2
 */
static MpIeee fd_half_d_data[30] =  {
  MpIeee( "1.5116909434145508537" ),
 -MpIeee( "0.0036043405371630468" ),
  MpIeee( "0.0014207743256393359" ),
 -MpIeee( "0.0005045399052400260" ),
  MpIeee( "0.0001690758006957347" ),
 -MpIeee( "0.0000546305872688307" ),
  MpIeee( "0.0000172223228484571" ),
 -MpIeee( "5.3352603788706e-6" ),
  MpIeee( "1.6315287543662e-6" ),
 -MpIeee( "4.939021084898e-7" ),
  MpIeee( "1.482515450316e-7" ),
 -MpIeee( "4.41552276226e-8" ),
  MpIeee( "1.30503160961e-8" ),
 -MpIeee( "3.8262599802e-9" ),
  MpIeee( "1.1123226976e-9" ),
 -MpIeee( "3.204765534e-10" ),
  MpIeee( "9.14870489e-11" ),
 -MpIeee( "2.58778946e-11" ),
  MpIeee( "7.2550731e-12" ),
 -MpIeee( "2.0172226e-12" ),
  MpIeee( "5.566891e-13" ),
 -MpIeee( "1.526247e-13" ),
  MpIeee( "4.16121e-14" ),
 -MpIeee( "1.12933e-14" ),
  MpIeee( "3.0537e-15" ),
 -MpIeee( "8.234e-16" ),
  MpIeee( "2.215e-16" ),
 -MpIeee( "5.95e-17" ),
  MpIeee( "1.59e-17" ),
 -MpIeee( "4.0e-18" )
};
static cheb_series fd_half_d_cs = {
  fd_half_d_data,
  29,
  -1, 1,
  15
};



/* Chebyshev fit for F_{3/2}(t);  -1 < t < 1, -1 < x < 1
 */
static MpIeee fd_3half_a_data[20] =  {
  MpIeee( "2.0404775940601704976" ),
  MpIeee( "0.8122168298093491444" ),
  MpIeee( "0.1536371165644008069" ),
  MpIeee( "0.0156174323847845125" ),
  MpIeee( "0.0005943427879290297" ),
 -MpIeee( "0.0000429609447738365" ),
 -MpIeee( "3.8246452994606e-6" ),
  MpIeee( "3.802306180287e-7" ),
  MpIeee( "4.05746157593e-8" ),
 -MpIeee( "4.5530360159e-9" ),
 -MpIeee( "5.306873139e-10" ),
  MpIeee( "6.37297268e-11" ),
  MpIeee( "7.8403674e-12" ),
 -MpIeee( "9.840241e-13" ),
 -MpIeee( "1.255952e-13" ),
  MpIeee( "1.62617e-14" ),
  MpIeee( "2.1318e-15" ),
 -MpIeee( "2.825e-16" ),
 -MpIeee( "3.78e-17" ),
  MpIeee( "5.1e-18" )
};
static cheb_series fd_3half_a_cs = {
  fd_3half_a_data,
  19,
  -1, 1,
  11
};


/* Chebyshev fit for F_{3/2}(3/2(t+1) + 1);  -1 < t < 1, 1 < x < 4
 */
static MpIeee fd_3half_b_data[22] =  {
  MpIeee( "13.403206654624176674" ),
  MpIeee( "5.574508357051880924" ),
  MpIeee( "0.931228574387527769" ),
  MpIeee( "0.054638356514085862" ),
 -MpIeee( "0.001477172902737439" ),
 -MpIeee( "0.000029378553381869" ),
  MpIeee( "0.000018357033493246" ),
 -MpIeee( "2.348059218454e-6" ),
  MpIeee( "8.3173787440e-8" ),
  MpIeee( "2.6826486956e-8" ),
 -MpIeee( "6.011244398e-9" ),
  MpIeee( "4.94345981e-10" ),
  MpIeee( "3.9557340e-11" ),
 -MpIeee( "1.7894930e-11" ),
  MpIeee( "2.348972e-12" ),
 -MpIeee( "1.2823e-14" ),
 -MpIeee( "5.4192e-14" ),
  MpIeee( "1.0527e-14" ),
 -MpIeee( "6.39e-16" ),
 -MpIeee( "1.47e-16" ),
  MpIeee( "4.5e-17" ),
 -MpIeee( "5.e-18" )
};
static cheb_series fd_3half_b_cs = {
  fd_3half_b_data,
  21,
  -1, 1,
  12
};


/* Chebyshev fit for F_{3/2}(3(t+1) + 4);  -1 < t < 1, 4 < x < 10
 */
static MpIeee fd_3half_c_data[21] =  {
  MpIeee( "101.03685253378877642" ),
  MpIeee( "43.62085156043435883" ),
  MpIeee( "6.62241373362387453" ),
  MpIeee( "0.25081415008708521" ),
 -MpIeee( "0.00798124846271395" ),
  MpIeee( "0.00063462245101023" ),
 -MpIeee( "0.00006392178890410" ),
  MpIeee( "6.04535131939e-6" ),
 -MpIeee( "3.4007683037e-7" ),
 -MpIeee( "4.072661545e-8" ),
  MpIeee( "1.931148453e-8" ),
 -MpIeee( "4.46328355e-9" ),
  MpIeee( "7.9434717e-10" ),
 -MpIeee( "1.1573569e-10" ),
  MpIeee( "1.304658e-11" ),
 -MpIeee( "7.4114e-13" ),
 -MpIeee( "1.4181e-13" ),
  MpIeee( "6.491e-14" ),
 -MpIeee( "1.597e-14" ),
  MpIeee( "3.05e-15" ),
 -MpIeee( "4.8e-16" )
};
static cheb_series fd_3half_c_cs = {
  fd_3half_c_data,
  20,
  -1, 1,
  12
};


/* Chebyshev fit for F_{3/2}(x) / x^(5/2)
 * 10 < x < 30 
 * -1 < t < 1
 * t = 1/10 (x-10) - 1 = x/10 - 2
 */
static MpIeee fd_3half_d_data[25] =  {
  MpIeee( "0.6160645215171852381" ),
 -MpIeee( "0.0071239478492671463" ),
  MpIeee( "0.0027906866139659846" ),
 -MpIeee( "0.0009829521424317718" ),
  MpIeee( "0.0003260229808519545" ),
 -MpIeee( "0.0001040160912910890" ),
  MpIeee( "0.0000322931223232439" ),
 -MpIeee( "9.8243506588102e-6" ),
  MpIeee( "2.9420132351277e-6" ),
 -MpIeee( "8.699154670418e-7" ),
  MpIeee( "2.545460071999e-7" ),
 -MpIeee( "7.38305056331e-8" ),
  MpIeee( "2.12545670310e-8" ),
 -MpIeee( "6.0796532462e-9" ),
  MpIeee( "1.7294556741e-9" ),
 -MpIeee( "4.896540687e-10" ),
  MpIeee( "1.380786037e-10" ),
 -MpIeee( "3.88057305e-11" ),
  MpIeee( "1.08753212e-11" ),
 -MpIeee( "3.0407308e-12" ),
  MpIeee( "8.485626e-13" ),
 -MpIeee( "2.364275e-13" ),
  MpIeee( "6.57636e-14" ),
 -MpIeee( "1.81807e-14" ),
  MpIeee( "4.6884e-15" )
};
static cheb_series fd_3half_d_cs = {
  fd_3half_d_data,
  24,
  -1, 1,
  16
};


/* Goano's modification of the Levin-u implementation.
 * This is a simplification of the original WHIZ algorithm.
 * See [Fessler et al., ACM Toms 9, 346 (1983)].
 */
static
int
 fd_whiz(const MpIeee term, const int iterm,
        MpIeee * qnum, MpIeee * qden,
        MpIeee * result, MpIeee * s)
{
  if(iterm == 0) *s = MpIeee( "0.0" );

  *s += term;

  qden[iterm] = MpIeee( "1.0" )/(term*(iterm+MpIeee( "1.0" ))*(iterm+MpIeee( "1.0" )));
  qnum[iterm] = *s * qden[iterm];

  if(iterm > 0) {
    MpIeee factor=  MpIeee( "1.0" );
    MpIeee ratio=  iterm/(iterm+MpIeee( "1.0" ));
    int  j;
    for(j=iterm-1; j>=0; j--) {
      MpIeee c=  factor * (j+MpIeee( "1.0" )) / (iterm+MpIeee( "1.0" ));
      factor *= ratio;
      qden[j] = qden[j+1] - c * qden[j];
      qnum[j] = qnum[j+1] - c * qnum[j];
    }
  }

  *result = qnum[0] / qden[0];
  return GSL_SUCCESS;
}


/* Handle case of integer j <= -2.
 */
static
int
 fd_nint(const int j, const MpIeee x, gsl_sf_result * result)
{
/*    const int nsize = 100 + 1; */
  enum {
    nsize = 100+1
  };
  MpIeee qcoeff[nsize];

  if(j >= -1) {
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_ESANITY);
  }
  else if(j < -(nsize)) {
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EUNIMPL);
  }
  else {
    MpIeee a;MpIeee  p;MpIeee  f;
    int  i;int   k;
    int  n=  -(j+1);

    qcoeff[1] = MpIeee( "1.0" );

    for(k=2; k<=n; k++) {
      qcoeff[k] = -qcoeff[k-1];
      for(i=k-1; i>=2; i--) {
        qcoeff[i] = i*qcoeff[i] - (k-(i-MpIeee( "1" )))*qcoeff[i-1];
      }
    }

    if(x >= 0.0) {
      a = exp(-x);
      f = qcoeff[1];
      for(i=2; i<=n; i++) {
        f = f*a + qcoeff[i];
      }
    }
    else {
      a = exp(x);
      f = qcoeff[n];
      for(i=n-1; i>=1; i--) {
        f = f*a + qcoeff[i];
      }
    }

    p = gsl_sf_pow_int(MpIeee( "1.0" )+a, j);
    result->val = f*a*p;
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(f*a*p);
    return GSL_SUCCESS;
  }
}


/* x < 0
 */
static
int
 fd_neg(const MpIeee j, const MpIeee x, gsl_sf_result * result)
{
  enum {
    itmax = 100,
    qsize = 100+1
  };
/*    const int itmax = 100; */
/*    const int qsize = 100 + 1; */
  MpIeee qnum[qsize];MpIeee  qden[qsize];

  if(x < GSL_LOG_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < -1.0 && x < -fabs(j+1.0)) {
    /* Simple series implementation. Avoid the
     * complexity and extra work of the series
     * acceleration method below.
     */
    MpIeee ex=  exp(x);
    MpIeee term=  ex;
    MpIeee sum=  term;
    int  n;
    for(n=2; n<100; n++) {
      MpIeee rat=  (n-MpIeee( "1.0" ))/n;
      MpIeee p=  pow(rat, j+MpIeee( "1.0" ));
      term *= -ex * p;
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }
    result->val = sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(sum);
    return GSL_SUCCESS;
  }
  else {
    MpIeee s;
    MpIeee xn=  x;
    MpIeee ex=  -exp(x);
    MpIeee enx=  -ex;
    MpIeee f=  MpIeee( "0.0" );
    MpIeee f_previous;
    int  jterm;
    for(jterm=0; jterm<=itmax; jterm++) {
      MpIeee p=  pow(jterm+MpIeee( "1.0" ), j+MpIeee( "1.0" ));
      MpIeee term=  enx/p;
      f_previous = f;
      fd_whiz(term, jterm, qnum, qden, &f, &s);
      xn += x;
      if(fabs(f-f_previous) < fabs(f)*2.0*GSL_DBL_EPSILON || xn < GSL_LOG_DBL_MIN) break;
      enx *= ex;
    }

    result->val  = f;
    result->err  = fabs(f-f_previous);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(f);

    if(jterm == itmax)
      GSL_ERROR ("error", GSL_EMAXITER);
    else
      return GSL_SUCCESS;
  }
}


/* asymptotic expansion
 * j + 2.0 > 0.0
 */
static
int
 fd_asymp(const MpIeee j, const MpIeee x, gsl_sf_result * result)
{
  const int j_integer = ( fabs(j - floor(j+0.5)) < 100.0*GSL_DBL_EPSILON );
  const int itmax = 200;
  gsl_sf_result lg;
  int  stat_lg=  gsl_sf_lngamma_e(j + 2.0, &lg);
  MpIeee seqn_val=  MpIeee( "0.5" );
  MpIeee seqn_err=  MpIeee( "0.0" );
  MpIeee xm2=  (MpIeee( "1.0" )/x)/x;
  MpIeee xgam=  MpIeee( "1.0" );
  MpIeee add=  GSL_DBL_MAX;
  MpIeee cos_term;
  MpIeee ln_x;
  MpIeee ex_term_1;
  MpIeee ex_term_2;
  gsl_sf_result fneg;
  gsl_sf_result ex_arg;
  gsl_sf_result ex;
  int  stat_fneg;
  int  stat_e;
  int  n;
  for(n=1; n<=itmax; n++) {
    MpIeee add_previous=  add;
    gsl_sf_result eta;
    gsl_sf_eta_int_e(2*n, &eta);
    xgam = xgam * xm2 * (j + MpIeee( "1.0" ) - (MpIeee( "2" )*n-MpIeee( "2" ))) * (j + MpIeee( "1.0" ) - (MpIeee( "2" )*n-MpIeee( "1" )));
    add  = eta.val * xgam;
    if(!j_integer && fabs(add) > fabs(add_previous)) break;
    if(fabs(add/seqn_val) < GSL_DBL_EPSILON) break;
    seqn_val += add;
    seqn_err += MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs(add);
  }
  seqn_err += fabs(add);

  stat_fneg = fd_neg(j, -x, &fneg);
  ln_x = log(x);
  ex_term_1 = (j+MpIeee( "1.0" ))*ln_x;
  ex_term_2 = lg.val;
  ex_arg.val = ex_term_1 - ex_term_2; /*(j+1.0)*ln_x - lg.val; */
  ex_arg.err = GSL_DBL_EPSILON*(fabs(ex_term_1) + fabs(ex_term_2)) + lg.err;
  stat_e    = gsl_sf_exp_err_e(ex_arg.val, ex_arg.err, &ex);
  cos_term  = cos(j*M_PI);
  result->val  = cos_term * fneg.val + 2.0 * seqn_val * ex.val;
  result->err  = fabs(2.0 * ex.err * seqn_val);
  result->err += fabs(2.0 * ex.val * seqn_err);
  result->err += fabs(cos_term) * fneg.err;
  result->err += 4.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_ERROR_SELECT_3(stat_e, stat_fneg, stat_lg);
}


/* Series evaluation for small x, generic j.
 * [Goano (8)]
 */
#if 0
static
int
 fd_series(const MpIeee j, const MpIeee x, MpIeee * result)
{
  const int nmax = 1000;
  int  n;
  MpIeee sum=  MpIeee( "0.0" );
  MpIeee prev;
  MpIeee pow_factor=  MpIeee( "1.0" );
  MpIeee eta_factor;
  gsl_sf_eta_e(j + 1.0, &eta_factor);
  prev = pow_factor * eta_factor;
  sum += prev;
  for(n=1; n<nmax; n++) {
    MpIeee term;
    gsl_sf_eta_e(j+1.0-n, &eta_factor);
    pow_factor *= x/n;
    term = pow_factor * eta_factor;
    sum += term;
    if(fabs(term/sum) < GSL_DBL_EPSILON && fabs(prev/sum) < GSL_DBL_EPSILON) break;
    prev = term;
  }

  *result = sum;
  return GSL_SUCCESS;
}
#endif /* 0 */


/* Series evaluation for small x > 0, integer j > 0; x < Pi.
 * [Goano (8)]
 */
static
int
 fd_series_int(const int j, const MpIeee x, gsl_sf_result * result)
{
  int  n;
  MpIeee sum=  MpIeee( "0.0" );
  MpIeee del;
  MpIeee pow_factor=  MpIeee( "1.0" );
  gsl_sf_result eta_factor;
  gsl_sf_eta_int_e(j + 1, &eta_factor);
  del  = pow_factor * eta_factor.val;
  sum += del;

  /* Sum terms where the argument
   * of eta() is positive.
   */
  for(n=1; n<=j+2; n++) {
    gsl_sf_eta_int_e(j+1-n, &eta_factor);
    pow_factor *= x/n;
    del  = pow_factor * eta_factor.val;
    sum += del;
    if(fabs(del/sum) < GSL_DBL_EPSILON) break;
  }

  /* Now sum the terms where eta() is negative.
   * The argument of eta() must be odd as well,
   * so it is convenient to transform the series
   * as follows:
   *
   * Sum[ eta(j+1-n) x^n / n!, {n,j+4,Infinity}]
   * = x^j / j! Sum[ eta(1-2m) x^(2m) j! / (2m+j)! , {m,2,Infinity}]
   *
   * We do not need to do this sum if j is large enough.
   */
  if(j < 32) {
    int  m;
    gsl_sf_result jfact;
    MpIeee sum2;
    MpIeee pre2;

    gsl_sf_fact_e((unsigned int)j, &jfact);
    pre2 = gsl_sf_pow_int(x, j) / jfact.val;

    gsl_sf_eta_int_e(-3, &eta_factor);
    pow_factor = x*x*x*x / ((j+MpIeee( "4" ))*(j+MpIeee( "3" ))*(j+MpIeee( "2" ))*(j+MpIeee( "1" )));
    sum2 = eta_factor.val * pow_factor;

    for(m=3; m<24; m++) {
      gsl_sf_eta_int_e(1-2*m, &eta_factor);
      pow_factor *= x*x / ((j+MpIeee( "2" )*m)*(j+MpIeee( "2" )*m-MpIeee( "1" )));
      sum2 += eta_factor.val * pow_factor;
    }

    sum += pre2 * sum2;
  }

  result->val = sum;
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(sum);

  return GSL_SUCCESS;
}


/* series of hypergeometric functions for integer j > 0, x > 0
 * [Goano (7)]
 */
static
int
 fd_UMseries_int(const int j, const MpIeee x, gsl_sf_result * result)
{
  const int nmax = 2000;
  MpIeee pre;
  MpIeee lnpre_val;
  MpIeee lnpre_err;
  MpIeee sum_even_val=  MpIeee( "1.0" );
  MpIeee sum_even_err=  MpIeee( "0.0" );
  MpIeee sum_odd_val=  MpIeee( "0.0" );
  MpIeee sum_odd_err=  MpIeee( "0.0" );
  int  stat_sum;
  int  stat_e;
  int  stat_h=  GSL_SUCCESS;
  int  n;

  if(x < 500.0 && j < 80) {
    MpIeee p=  gsl_sf_pow_int(x, j+MpIeee( "1" ));
    gsl_sf_result g;
    gsl_sf_fact_e(j+1, &g); /* Gamma(j+2) */
    lnpre_val = MpIeee( "0.0" );
    lnpre_err = MpIeee( "0.0" );
    pre   = p/g.val;
  }
  else {
    MpIeee lnx=  log(x);
    gsl_sf_result lg;
    gsl_sf_lngamma_e(j + 2.0, &lg);
    lnpre_val = (j+MpIeee( "1.0" ))*lnx - lg.val;
    lnpre_err = MpIeee( "2.0" ) * GSL_DBL_EPSILON * fabs((j+MpIeee( "1.0" ))*lnx) + lg.err;
    pre = MpIeee( "1.0" );
  }

  /* Add up the odd terms of the sum.
   */
  for(n=1; n<nmax; n+=2) {
    MpIeee del_val;
    MpIeee del_err;
    gsl_sf_result U;
    gsl_sf_result M;
    int  stat_h_U=  gsl_sf_hyperg_U_int_e(1, j+2, n*x, &U);
    int  stat_h_F=  gsl_sf_hyperg_1F1_int_e(1, j+2, -n*x, &M);
    stat_h = GSL_ERROR_SELECT_3(stat_h, stat_h_U, stat_h_F);
    del_val = ((j+MpIeee( "1.0" ))*U.val - M.val);
    del_err = (fabs(j+MpIeee( "1.0" ))*U.err + M.err);
    sum_odd_val += del_val;
    sum_odd_err += del_err;
    if(fabs(del_val/sum_odd_val) < GSL_DBL_EPSILON) break;
  }

  /* Add up the even terms of the sum.
   */
  for(n=2; n<nmax; n+=2) {
    MpIeee del_val;
    MpIeee del_err;
    gsl_sf_result U;
    gsl_sf_result M;
    int  stat_h_U=  gsl_sf_hyperg_U_int_e(1, j+2, n*x, &U);
    int  stat_h_F=  gsl_sf_hyperg_1F1_int_e(1, j+2, -n*x, &M);
    stat_h = GSL_ERROR_SELECT_3(stat_h, stat_h_U, stat_h_F);
    del_val = ((j+MpIeee( "1.0" ))*U.val - M.val);
    del_err = (fabs(j+MpIeee( "1.0" ))*U.err + M.err);
    sum_even_val -= del_val;
    sum_even_err += del_err;
    if(fabs(del_val/sum_even_val) < GSL_DBL_EPSILON) break;
  }

  stat_sum = ( n >= nmax ? GSL_EMAXITER : GSL_SUCCESS );
  stat_e   = gsl_sf_exp_mult_err_e(lnpre_val, lnpre_err,
                                      pre*(sum_even_val + sum_odd_val),
                                      pre*(sum_even_err + sum_odd_err),
                                      result);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_3(stat_e, stat_h, stat_sum);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* [Goano (4)] */
int  gsl_sf_fermi_dirac_m1_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < 0.0) {
    const MpIeee ex=  exp(x);
    result->val = ex/(1.0+ex);
    result->err = 2.0 * (fabs(x) + 1.0) * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee ex=  exp(-x);
    result->val = 1.0/(1.0 + ex);
    result->err = 2.0 * GSL_DBL_EPSILON * (x + 1.0) * ex;
    return GSL_SUCCESS;
  }
}


/* [Goano (3)] */
int  gsl_sf_fermi_dirac_0_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < -5.0) {
    MpIeee ex=  exp(x);
    MpIeee ser=  MpIeee( "1.0" ) - ex*(MpIeee( "0.5" ) - ex*(MpIeee( "1.0" )/MpIeee( "3.0" ) - ex*(MpIeee( "1.0" )/MpIeee( "4.0" ) - ex*(MpIeee( "1.0" )/MpIeee( "5.0" ) - ex/MpIeee( "6.0" )))));
    result->val = ex * ser;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 10.0) {
    result->val = log(1.0 + exp(x));
    result->err = fabs(x * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
  else {
    MpIeee ex=  exp(-x);
    result->val = x + ex * (1.0 - 0.5*ex + ex*ex/3.0 - ex*ex*ex/4.0);
    result->err = (x + ex) * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}


int  gsl_sf_fermi_dirac_1_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < -1.0) {
    /* series [Goano (6)]
     */
    MpIeee ex=  exp(x);
    MpIeee term=  ex;
    MpIeee sum=  term;
    int  n;
    for(n=2; n<100 ; n++) {
      MpIeee rat=  (n-MpIeee( "1.0" ))/n;
      term *= -ex * rat * rat;
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }
    result->val = sum;
    result->err = 2.0 * fabs(sum) * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 1.0) {
    return cheb_eval_e(&fd_1_a_cs, x, result);
  }
  else if(x < 4.0) {
    MpIeee t=  MpIeee( "2.0" )/MpIeee( "3.0" )*(x-MpIeee( "1.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_1_b_cs, t, result);
  }
  else if(x < 10.0) {
    MpIeee t=  MpIeee( "1.0" )/MpIeee( "3.0" )*(x-MpIeee( "4.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_1_c_cs, t, result);
  }
  else if(x < 30.0) {
    MpIeee t=  MpIeee( "0.1" )*x - MpIeee( "2.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_1_d_cs, t, &c);
    result->val  = c.val * x*x;
    result->err  = c.err * x*x + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 1.0/GSL_SQRT_DBL_EPSILON) {
    MpIeee t=  MpIeee( "60.0" )/x - MpIeee( "1.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_1_e_cs, t, &c);
    result->val  = c.val * x*x;
    result->err  = c.err * x*x + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_SQRT_DBL_MAX) {
    result->val = 0.5 * x*x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int  gsl_sf_fermi_dirac_2_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < -1.0) {
    /* series [Goano (6)]
     */
    MpIeee ex=  exp(x);
    MpIeee term=  ex;
    MpIeee sum=  term;
    int  n;
    for(n=2; n<100 ; n++) {
      MpIeee rat=  (n-MpIeee( "1.0" ))/n;
      term *= -ex * rat * rat * rat;
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }
    result->val = sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(sum);
    return GSL_SUCCESS;
  }
  else if(x < 1.0) {
    return cheb_eval_e(&fd_2_a_cs, x, result);
  }
  else if(x < 4.0) {
    MpIeee t=  MpIeee( "2.0" )/MpIeee( "3.0" )*(x-MpIeee( "1.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_2_b_cs, t, result);
  }
  else if(x < 10.0) {
    MpIeee t=  MpIeee( "1.0" )/MpIeee( "3.0" )*(x-MpIeee( "4.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_2_c_cs, t, result);
  }
  else if(x < 30.0) {
    MpIeee t=  MpIeee( "0.1" )*x - MpIeee( "2.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_2_d_cs, t, &c);
    result->val  = c.val * x*x*x;
    result->err  = c.err * x*x*x + 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 1.0/GSL_ROOT3_DBL_EPSILON) {
    MpIeee t=  MpIeee( "60.0" )/x - MpIeee( "1.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_2_e_cs, t, &c);
    result->val = c.val * x*x*x;
    result->err = c.err * x*x*x + 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_ROOT3_DBL_MAX) {
    result->val = 1.0/6.0 * x*x*x;
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int  gsl_sf_fermi_dirac_int_e(const int j, const MpIeee x, gsl_sf_result * result)
{
  if(j < -1) {
    return fd_nint(j, x, result);
  }
  else if (j == -1) {
    return gsl_sf_fermi_dirac_m1_e(x, result);
  }
  else if(j == 0) {
    return gsl_sf_fermi_dirac_0_e(x, result);
  }
  else if(j == 1) {
    return gsl_sf_fermi_dirac_1_e(x, result);
  }
  else if(j == 2) {
    return gsl_sf_fermi_dirac_2_e(x, result);
  }
  else if(x < 0.0) {
    return fd_neg(j, x, result);
  }
  else if(x == 0.0) {
    return gsl_sf_eta_int_e(j+1, result);
  }
  else if(x < 1.5) {
    return fd_series_int(j, x, result);
  }
  else {
    gsl_sf_result fasymp;
    int  stat_asymp=  fd_asymp(j, x, &fasymp);

    if(stat_asymp == GSL_SUCCESS) {
      result->val  = fasymp.val;
      result->err  = fasymp.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_asymp;
    }
    else {
      return fd_UMseries_int(j, x, result);
    }
  }
}


int  gsl_sf_fermi_dirac_mhalf_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < -1.0) {
    /* series [Goano (6)]
     */
    MpIeee ex=  exp(x);
    MpIeee term=  ex;
    MpIeee sum=  term;
    int  n;
    for(n=2; n<200 ; n++) {
      MpIeee rat=  (n-MpIeee( "1.0" ))/n;
      term *= -ex * sqrt(rat);
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }
    result->val = sum;
    result->err = 2.0 * fabs(sum) * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 1.0) {
    return cheb_eval_e(&fd_mhalf_a_cs, x, result);
  }
  else if(x < 4.0) {
    MpIeee t=  MpIeee( "2.0" )/MpIeee( "3.0" )*(x-MpIeee( "1.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_mhalf_b_cs, t, result);
  }
  else if(x < 10.0) {
    MpIeee t=  MpIeee( "1.0" )/MpIeee( "3.0" )*(x-MpIeee( "4.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_mhalf_c_cs, t, result);
  }
  else if(x < 30.0) {
    MpIeee rtx=  sqrt(x);
    MpIeee t=  MpIeee( "0.1" )*x - MpIeee( "2.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_mhalf_d_cs, t, &c);
    result->val  = c.val * rtx;
    result->err  = c.err * rtx + 0.5 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return fd_asymp(-0.5, x, result);
  }
}


int  gsl_sf_fermi_dirac_half_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < -1.0) {
    /* series [Goano (6)]
     */
    MpIeee ex=  exp(x);
    MpIeee term=  ex;
    MpIeee sum=  term;
    int  n;
    for(n=2; n<100 ; n++) {
      MpIeee rat=  (n-MpIeee( "1.0" ))/n;
      term *= -ex * rat * sqrt(rat);
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }
    result->val = sum;
    result->err = 2.0 * fabs(sum) * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 1.0) {
    return cheb_eval_e(&fd_half_a_cs, x, result);
  }
  else if(x < 4.0) {
    MpIeee t=  MpIeee( "2.0" )/MpIeee( "3.0" )*(x-MpIeee( "1.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_half_b_cs, t, result);
  }
  else if(x < 10.0) {
    MpIeee t=  MpIeee( "1.0" )/MpIeee( "3.0" )*(x-MpIeee( "4.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_half_c_cs, t, result);
  }
  else if(x < 30.0) {
    MpIeee x32=  x*sqrt(x);
    MpIeee t=  MpIeee( "0.1" )*x - MpIeee( "2.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_half_d_cs, t, &c);
    result->val = c.val * x32;
    result->err = c.err * x32 + 1.5 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return fd_asymp(0.5, x, result);
  }
}


int  gsl_sf_fermi_dirac_3half_e(const MpIeee x, gsl_sf_result * result)
{
  if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else if(x < -1.0) {
    /* series [Goano (6)]
     */
    MpIeee ex=  exp(x);
    MpIeee term=  ex;
    MpIeee sum=  term;
    int  n;
    for(n=2; n<100 ; n++) {
      MpIeee rat=  (n-MpIeee( "1.0" ))/n;
      term *= -ex * rat * rat * sqrt(rat);
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }
    result->val = sum;
    result->err = 2.0 * fabs(sum) * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 1.0) {
    return cheb_eval_e(&fd_3half_a_cs, x, result);
  }
  else if(x < 4.0) {
    MpIeee t=  MpIeee( "2.0" )/MpIeee( "3.0" )*(x-MpIeee( "1.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_3half_b_cs, t, result);
  }
  else if(x < 10.0) {
    MpIeee t=  MpIeee( "1.0" )/MpIeee( "3.0" )*(x-MpIeee( "4.0" )) - MpIeee( "1.0" );
    return cheb_eval_e(&fd_3half_c_cs, t, result);
  }
  else if(x < 30.0) {
    MpIeee x52=  x*x*sqrt(x);
    MpIeee t=  MpIeee( "0.1" )*x - MpIeee( "2.0" );
    gsl_sf_result c;
    cheb_eval_e(&fd_3half_d_cs, t, &c);
    result->val = c.val * x52;
    result->err = c.err * x52 + 2.5 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return fd_asymp(1.5, x, result);
  }
}

/* [Goano p. 222] */
int  gsl_sf_fermi_dirac_inc_0_e(const MpIeee x, const MpIeee b, gsl_sf_result * result)
{
  if(b < 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    MpIeee arg=  b - x;
    gsl_sf_result f0;
    int  status=  gsl_sf_fermi_dirac_0_e(arg, &f0);
    result->val = f0.val - arg;
    result->err = f0.err + GSL_DBL_EPSILON * (fabs(x) + fabs(b));
    return status;
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_fermi_dirac_m1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_m1_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_0(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_0_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_1(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_1_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_2(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_2_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_int(const int j, const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_int_e(j, x, &result));
}

MpIeee gsl_sf_fermi_dirac_mhalf(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_mhalf_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_half(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_half_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_3half(const MpIeee x)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_3half_e(x, &result));
}

MpIeee gsl_sf_fermi_dirac_inc_0(const MpIeee x, const MpIeee b)
{
  EVAL_RESULT(gsl_sf_fermi_dirac_inc_0_e(x, b, &result));
}
