#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/airy.c
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
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_airy.h>

#include "error.h"
#include "check.h"

#include "chebyshev.h"
#include "cheb_eval_mode.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* chebyshev expansions for Airy modulus and phase
   based on SLATEC r9aimp()

 Series for AM21       on the interval -1.25000D-01 to  0.
                                        with weighted error   2.89E-17
                                         log weighted error  16.54
                               significant figures required  14.15
                                    decimal places required  17.34

 Series for ATH1       on the interval -1.25000D-01 to  0.
                                        with weighted error   2.53E-17
                                         log weighted error  16.60
                               significant figures required  15.15
                                    decimal places required  17.38

 Series for AM22       on the interval -1.00000D+00 to -1.25000D-01
                                        with weighted error   2.99E-17
                                         log weighted error  16.52
                               significant figures required  14.57
                                    decimal places required  17.28

 Series for ATH2       on the interval -1.00000D+00 to -1.25000D-01
                                        with weighted error   2.57E-17
                                         log weighted error  16.59
                               significant figures required  15.07
                                    decimal places required  17.34
*/

static MpIeee am21_data[37] =  {
  MpIeee( "0.0065809191761485" ),
  MpIeee( "0.0023675984685722" ),
  MpIeee( "0.0001324741670371" ),
  MpIeee( "0.0000157600904043" ),
  MpIeee( "0.0000027529702663" ),
  MpIeee( "0.0000006102679017" ),
  MpIeee( "0.0000001595088468" ),
  MpIeee( "0.0000000471033947" ),
  MpIeee( "0.0000000152933871" ),
  MpIeee( "0.0000000053590722" ),
  MpIeee( "0.0000000020000910" ),
  MpIeee( "0.0000000007872292" ),
  MpIeee( "0.0000000003243103" ),
  MpIeee( "0.0000000001390106" ),
  MpIeee( "0.0000000000617011" ),
  MpIeee( "0.0000000000282491" ),
  MpIeee( "0.0000000000132979" ),
  MpIeee( "0.0000000000064188" ),
  MpIeee( "0.0000000000031697" ),
  MpIeee( "0.0000000000015981" ),
  MpIeee( "0.0000000000008213" ),
  MpIeee( "0.0000000000004296" ),
  MpIeee( "0.0000000000002284" ),
  MpIeee( "0.0000000000001232" ),
  MpIeee( "0.0000000000000675" ),
  MpIeee( "0.0000000000000374" ),
  MpIeee( "0.0000000000000210" ),
  MpIeee( "0.0000000000000119" ),
  MpIeee( "0.0000000000000068" ),
  MpIeee( "0.0000000000000039" ),
  MpIeee( "0.0000000000000023" ),
  MpIeee( "0.0000000000000013" ),
  MpIeee( "0.0000000000000008" ),
  MpIeee( "0.0000000000000005" ),
  MpIeee( "0.0000000000000003" ),
  MpIeee( "0.0000000000000001" ),
  MpIeee( "0.0000000000000001" )
};
static cheb_series am21_cs = {
  am21_data,
  36,
  -1, 1,
  20
};

static MpIeee ath1_data[36] =  {
  -MpIeee( "0.07125837815669365" ),
  -MpIeee( "0.00590471979831451" ),
  -MpIeee( "0.00012114544069499" ),
  -MpIeee( "0.00000988608542270" ),
  -MpIeee( "0.00000138084097352" ),
  -MpIeee( "0.00000026142640172" ),
  -MpIeee( "0.00000006050432589" ),
  -MpIeee( "0.00000001618436223" ),
  -MpIeee( "0.00000000483464911" ),
  -MpIeee( "0.00000000157655272" ),
  -MpIeee( "0.00000000055231518" ),
  -MpIeee( "0.00000000020545441" ),
  -MpIeee( "0.00000000008043412" ),
  -MpIeee( "0.00000000003291252" ),
  -MpIeee( "0.00000000001399875" ),
  -MpIeee( "0.00000000000616151" ),
  -MpIeee( "0.00000000000279614" ),
  -MpIeee( "0.00000000000130428" ),
  -MpIeee( "0.00000000000062373" ),
  -MpIeee( "0.00000000000030512" ),
  -MpIeee( "0.00000000000015239" ),
  -MpIeee( "0.00000000000007758" ),
  -MpIeee( "0.00000000000004020" ),
  -MpIeee( "0.00000000000002117" ),
  -MpIeee( "0.00000000000001132" ),
  -MpIeee( "0.00000000000000614" ),
  -MpIeee( "0.00000000000000337" ),
  -MpIeee( "0.00000000000000188" ),
  -MpIeee( "0.00000000000000105" ),
  -MpIeee( "0.00000000000000060" ),
  -MpIeee( "0.00000000000000034" ),
  -MpIeee( "0.00000000000000020" ),
  -MpIeee( "0.00000000000000011" ),
  -MpIeee( "0.00000000000000007" ),
  -MpIeee( "0.00000000000000004" ),
  -MpIeee( "0.00000000000000002" )
};
static cheb_series ath1_cs = {
  ath1_data,
  35,
  -1, 1,
  15
};

static MpIeee am22_data[33] =  {
 -MpIeee( "0.01562844480625341" ),
  MpIeee( "0.00778336445239681" ),
  MpIeee( "0.00086705777047718" ),
  MpIeee( "0.00015696627315611" ),
  MpIeee( "0.00003563962571432" ),
  MpIeee( "0.00000924598335425" ),
  MpIeee( "0.00000262110161850" ),
  MpIeee( "0.00000079188221651" ),
  MpIeee( "0.00000025104152792" ),
  MpIeee( "0.00000008265223206" ),
  MpIeee( "0.00000002805711662" ),
  MpIeee( "0.00000000976821090" ),
  MpIeee( "0.00000000347407923" ),
  MpIeee( "0.00000000125828132" ),
  MpIeee( "0.00000000046298826" ),
  MpIeee( "0.00000000017272825" ),
  MpIeee( "0.00000000006523192" ),
  MpIeee( "0.00000000002490471" ),
  MpIeee( "0.00000000000960156" ),
  MpIeee( "0.00000000000373448" ),
  MpIeee( "0.00000000000146417" ),
  MpIeee( "0.00000000000057826" ),
  MpIeee( "0.00000000000022991" ),
  MpIeee( "0.00000000000009197" ),
  MpIeee( "0.00000000000003700" ),
  MpIeee( "0.00000000000001496" ),
  MpIeee( "0.00000000000000608" ),
  MpIeee( "0.00000000000000248" ),
  MpIeee( "0.00000000000000101" ),
  MpIeee( "0.00000000000000041" ),
  MpIeee( "0.00000000000000017" ),
  MpIeee( "0.00000000000000007" ),
  MpIeee( "0.00000000000000002" )
};
static cheb_series am22_cs = {
  am22_data,
  32,
  -1, 1,
  15
};

static MpIeee ath2_data[32] =  {
   MpIeee( "0.00440527345871877" ),
  -MpIeee( "0.03042919452318455" ),
  -MpIeee( "0.00138565328377179" ),
  -MpIeee( "0.00018044439089549" ),
  -MpIeee( "0.00003380847108327" ),
  -MpIeee( "0.00000767818353522" ),
  -MpIeee( "0.00000196783944371" ),
  -MpIeee( "0.00000054837271158" ),
  -MpIeee( "0.00000016254615505" ),
  -MpIeee( "0.00000005053049981" ),
  -MpIeee( "0.00000001631580701" ),
  -MpIeee( "0.00000000543420411" ),
  -MpIeee( "0.00000000185739855" ),
  -MpIeee( "0.00000000064895120" ),
  -MpIeee( "0.00000000023105948" ),
  -MpIeee( "0.00000000008363282" ),
  -MpIeee( "0.00000000003071196" ),
  -MpIeee( "0.00000000001142367" ),
  -MpIeee( "0.00000000000429811" ),
  -MpIeee( "0.00000000000163389" ),
  -MpIeee( "0.00000000000062693" ),
  -MpIeee( "0.00000000000024260" ),
  -MpIeee( "0.00000000000009461" ),
  -MpIeee( "0.00000000000003716" ),
  -MpIeee( "0.00000000000001469" ),
  -MpIeee( "0.00000000000000584" ),
  -MpIeee( "0.00000000000000233" ),
  -MpIeee( "0.00000000000000093" ),
  -MpIeee( "0.00000000000000037" ),
  -MpIeee( "0.00000000000000015" ),
  -MpIeee( "0.00000000000000006" ),
  -MpIeee( "0.00000000000000002" )
};
static cheb_series ath2_cs = {
  ath2_data,
  31,
  -1, 1,
  16
};


/* Airy modulus and phase for x < -1 */
static
int
 airy_mod_phase(const MpIeee x, gsl_mode_t mode, gsl_sf_result * mod, gsl_sf_result * phase)
{
  gsl_sf_result result_m;
  gsl_sf_result result_p;
  MpIeee m;MpIeee  p;
  MpIeee sqx;

  if(x < -2.0) {
    MpIeee z=  MpIeee( "16.0" )/(x*x*x) + MpIeee( "1.0" );
    cheb_eval_mode_e(&am21_cs, z, mode, &result_m);
    cheb_eval_mode_e(&ath1_cs, z, mode, &result_p);
  }
  else if(x <= -1.0) {
    MpIeee z=  (MpIeee( "16.0" )/(x*x*x) + MpIeee( "9.0" ))/MpIeee( "7.0" );
    cheb_eval_mode_e(&am22_cs, z, mode, &result_m);
    cheb_eval_mode_e(&ath2_cs, z, mode, &result_p);
  }
  else {
    mod->val = 0.0;
    mod->err = 0.0;
    phase->val = 0.0;
    phase->err = 0.0;
    GSL_ERROR ("x is greater than 1.0", GSL_EDOM);
  }

  m =  MpIeee( "0.3125" ) + result_m.val;
  p = -MpIeee( "0.625" )  + result_p.val;

  sqx = sqrt(-x);

  mod->val   = sqrt(m/sqx);
  mod->err  = fabs(mod->val) * (GSL_DBL_EPSILON + fabs(result_m.err/result_m.val));
  phase->val = M_PI_4 - x*sqx * p;
  phase->err = fabs(phase->val) * (GSL_DBL_EPSILON + fabs(result_p.err/result_p.val));

  return GSL_SUCCESS;
}



/* Chebyshev series for Ai(x) with x in [-1,1]
   based on SLATEC ai(x)

 series for aif        on the interval -1.00000d+00 to  1.00000d+00
                                   with weighted error   1.09e-19
                                    log weighted error  18.96
                          significant figures required  17.76
                               decimal places required  19.44

 series for aig        on the interval -1.00000d+00 to  1.00000d+00
                                   with weighted error   1.51e-17
                                    log weighted error  16.82
                          significant figures required  15.19
                               decimal places required  17.27
 */
static MpIeee ai_data_f[9] =  {
  -MpIeee( "0.03797135849666999750" ),
   MpIeee( "0.05919188853726363857" ),
   MpIeee( "0.00098629280577279975" ),
   MpIeee( "0.00000684884381907656" ),
   MpIeee( "0.00000002594202596219" ),
   MpIeee( "0.00000000006176612774" ),
   MpIeee( "0.00000000000010092454" ),
   MpIeee( "0.00000000000000012014" ),
   MpIeee( "0.00000000000000000010" )
};
static cheb_series aif_cs = {
  ai_data_f,
  8,
  -1, 1,
  8
};

static MpIeee ai_data_g[8] =  {
   MpIeee( "0.01815236558116127" ),
   MpIeee( "0.02157256316601076" ),
   MpIeee( "0.00025678356987483" ),
   MpIeee( "0.00000142652141197" ),
   MpIeee( "0.00000000457211492" ),
   MpIeee( "0.00000000000952517" ),
   MpIeee( "0.00000000000001392" ),
   MpIeee( "0.00000000000000001" )
};
static cheb_series aig_cs = {
  ai_data_g,
  7,
  -1, 1,
  7
};


/* Chebvyshev series for Bi(x) with x in [-1,1]
   based on SLATEC bi(x)

 series for bif        on the interval -1.00000d+00 to  1.00000d+00
                                        with weighted error   1.88e-19
                                         log weighted error  18.72
                               significant figures required  17.74
                                    decimal places required  19.20

 series for big        on the interval -1.00000d+00 to  1.00000d+00
                                        with weighted error   2.61e-17
                                         log weighted error  16.58
                               significant figures required  15.17
                                    decimal places required  17.03
 */
static MpIeee data_bif[9] =  {
  -MpIeee( "0.01673021647198664948" ),
   MpIeee( "0.10252335834249445610" ),
   MpIeee( "0.00170830925073815165" ),
   MpIeee( "0.00001186254546774468" ),
   MpIeee( "0.00000004493290701779" ),
   MpIeee( "0.00000000010698207143" ),
   MpIeee( "0.00000000000017480643" ),
   MpIeee( "0.00000000000000020810" ),
   MpIeee( "0.00000000000000000018" )
};
static cheb_series bif_cs = {
  data_bif,
  8,
  -1, 1,
  8
};

static MpIeee data_big[8] =  {
   MpIeee( "0.02246622324857452" ),
   MpIeee( "0.03736477545301955" ),
   MpIeee( "0.00044476218957212" ),
   MpIeee( "0.00000247080756363" ),
   MpIeee( "0.00000000791913533" ),
   MpIeee( "0.00000000001649807" ),
   MpIeee( "0.00000000000002411" ),
   MpIeee( "0.00000000000000002" )
};
static cheb_series big_cs = {
  data_big,
  7,
  -1, 1,
  7
};


/* Chebyshev series for Bi(x) with x in [1,8]
   based on SLATEC bi(x)
 */
static MpIeee data_bif2[10] =  {
  MpIeee( "0.0998457269381604100" ),
  MpIeee( "0.4786249778630055380" ),
  MpIeee( "0.0251552119604330118" ),
  MpIeee( "0.0005820693885232645" ),
  MpIeee( "0.0000074997659644377" ),
  MpIeee( "0.0000000613460287034" ),
  MpIeee( "0.0000000003462753885" ),
  MpIeee( "0.0000000000014288910" ),
  MpIeee( "0.0000000000000044962" ),
  MpIeee( "0.0000000000000000111" )
};
static cheb_series bif2_cs = {
  data_bif2,
  9,
  -1, 1,
  9
};

static MpIeee data_big2[10] =  {
  MpIeee( "0.033305662145514340" ),
  MpIeee( "0.161309215123197068" ),
  MpIeee( "0.0063190073096134286" ),
  MpIeee( "0.0001187904568162517" ),
  MpIeee( "0.0000013045345886200" ),
  MpIeee( "0.0000000093741259955" ),
  MpIeee( "0.0000000000474580188" ),
  MpIeee( "0.0000000000001783107" ),
  MpIeee( "0.0000000000000005167" ),
  MpIeee( "0.0000000000000000011" )
};
static cheb_series big2_cs = {
  data_big2,
  9,
  -1, 1,
  9
};


/* chebyshev for Ai(x) asymptotic factor 
   based on SLATEC aie()

 Series for AIP        on the interval  0.          to  1.00000D+00
                   with weighted error   5.10E-17
                    log weighted error  16.29
          significant figures required  14.41
               decimal places required  17.06
               
 [GJ] Sun Apr 19 18:14:31 EDT 1998
 There was something wrong with these coefficients. I was getting
 errors after 3 or 4 digits. So I recomputed this table. Now I get
 double precision agreement with Mathematica. But it does not seem
 possible that the small differences here would account for the
 original discrepancy. There must have been something wrong with my
 original usage...
*/
static MpIeee data_aip[36] =  {
 -MpIeee( "0.0187519297793867540198" ),
 -MpIeee( "0.0091443848250055004725" ),
  MpIeee( "0.0009010457337825074652" ),
 -MpIeee( "0.0001394184127221491507" ),
  MpIeee( "0.0000273815815785209370" ),
 -MpIeee( "0.0000062750421119959424" ),
  MpIeee( "0.0000016064844184831521" ),
 -MpIeee( "0.0000004476392158510354" ),
  MpIeee( "0.0000001334635874651668" ),
 -MpIeee( "0.0000000420735334263215" ),
  MpIeee( "0.0000000139021990246364" ),
 -MpIeee( "0.0000000047831848068048" ),
  MpIeee( "0.0000000017047897907465" ),
 -MpIeee( "0.0000000006268389576018" ),
  MpIeee( "0.0000000002369824276612" ),
 -MpIeee( "0.0000000000918641139267" ),
  MpIeee( "0.0000000000364278543037" ),
 -MpIeee( "0.0000000000147475551725" ),
  MpIeee( "0.0000000000060851006556" ),
 -MpIeee( "0.0000000000025552772234" ),
  MpIeee( "0.0000000000010906187250" ),
 -MpIeee( "0.0000000000004725870319" ),
  MpIeee( "0.0000000000002076969064" ),
 -MpIeee( "0.0000000000000924976214" ),
  MpIeee( "0.0000000000000417096723" ),
 -MpIeee( "0.0000000000000190299093" ),
  MpIeee( "0.0000000000000087790676" ),
 -MpIeee( "0.0000000000000040927557" ),
  MpIeee( "0.0000000000000019271068" ),
 -MpIeee( "0.0000000000000009160199" ),
  MpIeee( "0.0000000000000004393567" ),
 -MpIeee( "0.0000000000000002125503" ),
  MpIeee( "0.0000000000000001036735" ),
 -MpIeee( "0.0000000000000000509642" ),
  MpIeee( "0.0000000000000000252377" ),
 -MpIeee( "0.0000000000000000125793" )
/*
  -.0187519297793868
  -.0091443848250055,
   .0009010457337825,
  -.0001394184127221,
   .0000273815815785,
  -.0000062750421119,
   .0000016064844184,
  -.0000004476392158,
   .0000001334635874,
  -.0000000420735334,
   .0000000139021990,
  -.0000000047831848,
   .0000000017047897,
  -.0000000006268389,
   .0000000002369824,
  -.0000000000918641,
   .0000000000364278,
  -.0000000000147475,
   .0000000000060851,
  -.0000000000025552,
   .0000000000010906,
  -.0000000000004725,
   .0000000000002076,
  -.0000000000000924,
   .0000000000000417,
  -.0000000000000190,
   .0000000000000087,
  -.0000000000000040,
   .0000000000000019,
  -.0000000000000009,
   .0000000000000004,
  -.0000000000000002,
   .0000000000000001,
  -.0000000000000000
*/
};
static cheb_series aip_cs = {
  data_aip,
  35,
  -1, 1,
  17
};


/* chebyshev for Bi(x) asymptotic factor 
   based on SLATEC bie()

 Series for BIP        on the interval  1.25000D-01 to  3.53553D-01
                   with weighted error   1.91E-17
                    log weighted error  16.72
          significant figures required  15.35
               decimal places required  17.41

 Series for BIP2       on the interval  0.          to  1.25000D-01
                   with weighted error   1.05E-18
                    log weighted error  17.98
          significant figures required  16.74
               decimal places required  18.71
*/
static MpIeee data_bip[24] =  {
  -MpIeee( "0.08322047477943447" ),
   MpIeee( "0.01146118927371174" ),
   MpIeee( "0.00042896440718911" ),
  -MpIeee( "0.00014906639379950" ),
  -MpIeee( "0.00001307659726787" ),
   MpIeee( "0.00000632759839610" ),
  -MpIeee( "0.00000042226696982" ),
  -MpIeee( "0.00000019147186298" ),
   MpIeee( "0.00000006453106284" ),
  -MpIeee( "0.00000000784485467" ),
  -MpIeee( "0.00000000096077216" ),
   MpIeee( "0.00000000070004713" ),
  -MpIeee( "0.00000000017731789" ),
   MpIeee( "0.00000000002272089" ),
   MpIeee( "0.00000000000165404" ),
  -MpIeee( "0.00000000000185171" ),
   MpIeee( "0.00000000000059576" ),
  -MpIeee( "0.00000000000012194" ),
   MpIeee( "0.00000000000001334" ),
   MpIeee( "0.00000000000000172" ),
  -MpIeee( "0.00000000000000145" ),
   MpIeee( "0.00000000000000049" ),
  -MpIeee( "0.00000000000000011" ),
   MpIeee( "0.00000000000000001" )
};
static cheb_series bip_cs = {
  data_bip,
  23,
  -1, 1,
  14
};

static MpIeee data_bip2[29] =  {    
  -MpIeee( "0.113596737585988679" ),
   MpIeee( "0.0041381473947881595" ),
   MpIeee( "0.0001353470622119332" ),
   MpIeee( "0.0000104273166530153" ),
   MpIeee( "0.0000013474954767849" ),
   MpIeee( "0.0000001696537405438" ),
  -MpIeee( "0.0000000100965008656" ),
  -MpIeee( "0.0000000167291194937" ),
  -MpIeee( "0.0000000045815364485" ),
   MpIeee( "0.0000000003736681366" ),
   MpIeee( "0.0000000005766930320" ),
   MpIeee( "0.0000000000621812650" ),
  -MpIeee( "0.0000000000632941202" ),
  -MpIeee( "0.0000000000149150479" ),
   MpIeee( "0.0000000000078896213" ),
   MpIeee( "0.0000000000024960513" ),
  -MpIeee( "0.0000000000012130075" ),
  -MpIeee( "0.0000000000003740493" ),
   MpIeee( "0.0000000000002237727" ),
   MpIeee( "0.0000000000000474902" ),
  -MpIeee( "0.0000000000000452616" ),
  -MpIeee( "0.0000000000000030172" ),
   MpIeee( "0.0000000000000091058" ),
  -MpIeee( "0.0000000000000009814" ),
  -MpIeee( "0.0000000000000016429" ),
   MpIeee( "0.0000000000000005533" ),
   MpIeee( "0.0000000000000002175" ),
  -MpIeee( "0.0000000000000001737" ),
  -MpIeee( "0.0000000000000000010" )
};
static cheb_series bip2_cs = {
  data_bip2,
  28,
  -1, 1,
  10
};


/* assumes x >= 1.0 */
inline static int 
 airy_aie(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  MpIeee sqx=  sqrt(x);
  MpIeee z=  MpIeee( "2.0" )/(x*sqx) - MpIeee( "1.0" );
  MpIeee y=  sqrt(sqx);
  gsl_sf_result result_c;
  cheb_eval_mode_e(&aip_cs, z, mode, &result_c);
  result->val = (0.28125 + result_c.val)/y;
  result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

/* assumes x >= 2.0 */
static int  airy_bie(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee ATR=   8.7506905708484345;
  const MpIeee BTR=  -2.0938363213560543;

  if(x < 4.0) {
    MpIeee sqx=  sqrt(x);
    MpIeee z=  ATR/(x*sqx) + BTR;
    MpIeee y=  sqrt(sqx);
    gsl_sf_result result_c;
    cheb_eval_mode_e(&bip_cs, z, mode, &result_c);
    result->val = (0.625 + result_c.val)/y;
    result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    MpIeee sqx=  sqrt(x);
    MpIeee z=  MpIeee( "16.0" )/(x*sqx) - MpIeee( "1.0" );
    MpIeee y=  sqrt(sqx);
    gsl_sf_result result_c;
    cheb_eval_mode_e(&bip2_cs, z, mode, &result_c);
    result->val = (0.625 + result_c.val)/y;
    result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_airy_Ai_e(const MpIeee x, const gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result cos_result;
    int  stat_mp=  airy_mod_phase(x, mode, &mod, &theta);
    int  stat_cos=  gsl_sf_cos_err_e(theta.val, theta.err, &cos_result);
    result->val  = mod.val * cos_result.val;
    result->err  = fabs(mod.val * cos_result.err) + fabs(cos_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_cos);
  }
  else if(x <= 1.0) {
    const MpIeee z=  x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs, z, mode, &result_c1);
    result->val  = 0.375 + (result_c0.val - x*(0.25 + result_c1.val));
    result->err  = result_c0.err + fabs(x*result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee x32=  x * sqrt(x);
    MpIeee s=  exp(-MpIeee( "2.0" )*x32/MpIeee( "3.0" ));
    gsl_sf_result result_aie;
    int  stat_aie=  airy_aie(x, mode, &result_aie);
    result->val  = result_aie.val * s;
    result->err  = result_aie.err * s + result->val * x32 * GSL_DBL_EPSILON;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return stat_aie;
  }
}


int
 gsl_sf_airy_Ai_scaled_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result cos_result;
    int  stat_mp=  airy_mod_phase(x, mode, &mod, &theta);
    int  stat_cos=  gsl_sf_cos_err_e(theta.val, theta.err, &cos_result);
    result->val  = mod.val * cos_result.val;
    result->err  = fabs(mod.val * cos_result.err) + fabs(cos_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_cos);
  }
  else if(x <= 1.0) {
    const MpIeee z=  x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs, z, mode, &result_c1);
    result->val  = 0.375 + (result_c0.val - x*(0.25 + result_c1.val));
    result->err  = result_c0.err + fabs(x*result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    if(x > 0.0) {
      const MpIeee scale=  exp(2.0/3.0 * sqrt(z));
      result->val *= scale;
      result->err *= scale;
    }

    return GSL_SUCCESS;
  }
  else {
    return airy_aie(x, mode, result);
  }
}


int  gsl_sf_airy_Bi_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result sin_result;
    int  stat_mp=  airy_mod_phase(x, mode, &mod, &theta);
    int  stat_sin=  gsl_sf_sin_err_e(theta.val, theta.err, &sin_result);
    result->val  = mod.val * sin_result.val;
    result->err  = fabs(mod.val * sin_result.err) + fabs(sin_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_sin);
  }
  else if(x < 1.0) {
    const MpIeee z=  x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big_cs, z, mode, &result_c1);
    result->val  = 0.625 + result_c0.val + x*(0.4375 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    const MpIeee z=  (2.0*x*x*x - 9.0)/7.0;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c1);
    result->val  = 1.125 + result_c0.val + x*(0.625 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee y=  2.0*x*sqrt(x)/3.0;
    const MpIeee s=  exp(y);

    if(y > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR(result);
    }
    else {
      gsl_sf_result result_bie;
      int  stat_bie=  airy_bie(x, mode, &result_bie);
      result->val  = result_bie.val * s;
      result->err  = result_bie.err * s + fabs(1.5*y * (GSL_DBL_EPSILON * result->val));
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return stat_bie;
    }
  }
}


int
 gsl_sf_airy_Bi_scaled_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result sin_result;
    int  stat_mp=  airy_mod_phase(x, mode, &mod, &theta);
    int  stat_sin=  gsl_sf_sin_err_e(theta.val, theta.err, &sin_result);
    result->val  = mod.val * sin_result.val;
    result->err  = fabs(mod.val * sin_result.err) + fabs(sin_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_sin);
  }
  else if(x < 1.0) {
    const MpIeee z=  x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big_cs, z, mode, &result_c1);
    result->val  = 0.625 + result_c0.val + x*(0.4375 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    if(x > 0.0) {
      const MpIeee scale=  exp(-2.0/3.0 * sqrt(z));
      result->val *= scale;
      result->err *= scale;
    }
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    const MpIeee x3=  x*x*x;
    const MpIeee z=  (2.0*x3 - 9.0)/7.0;
    const MpIeee s=  exp(-2.0/3.0 * sqrt(x3));
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c1);
    result->val  = s * (1.125 + result_c0.val + x*(0.625 + result_c1.val));
    result->err  = s * (result_c0.err + fabs(x * result_c1.err));
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return airy_bie(x, mode, result);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_airy_Ai(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_e(x, mode, &result));
}

MpIeee gsl_sf_airy_Ai_scaled(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_scaled_e(x, mode, &result));
}

MpIeee gsl_sf_airy_Bi(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_e(x, mode, &result));
}

MpIeee gsl_sf_airy_Bi_scaled(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_scaled_e(x, mode, &result));
}


