#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/airy_der.c
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
#include <gsl/gsl_sf_airy.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval_mode.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on SLATEC aide.f, bide.f, aid.f, bid.f, r9admp.f */
 
/* 
 series for aif on the interval -1.00000e+00 to  1.00000e+00
                                        with weighted error   5.22e-18
                                         log weighted error  17.28
                               significant figures required  16.01
                                    decimal places required  17.73
*/
static MpIeee aif_data[8] =  {
   MpIeee( "0.10527461226531408809" ),
   MpIeee( "0.01183613628152997844" ),
   MpIeee( "0.00012328104173225664" ),
   MpIeee( "0.00000062261225638140" ),
   MpIeee( "0.00000000185298887844" ),
   MpIeee( "0.00000000000363328873" ),
   MpIeee( "0.00000000000000504622" ),
   MpIeee( "0.00000000000000000522" )
};
static cheb_series aif_cs = {
  aif_data,
  7,
  -1, 1,
  7
};

/*
 series for aig on the interval -1.00000e+00 to  1.00000e+00
                                        with weighted error   3.14e-19
                                         log weighted error  18.50
                               significant figures required  17.44
                                    decimal places required  18.98
*/
static MpIeee aig_data[9] =  {
   MpIeee( "0.021233878150918666852" ),
   MpIeee( "0.086315930335214406752" ),
   MpIeee( "0.001797594720383231358" ),
   MpIeee( "0.000014265499875550693" ),
   MpIeee( "0.000000059437995283683" ),
   MpIeee( "0.000000000152403366479" ),
   MpIeee( "0.000000000000264587660" ),
   MpIeee( "0.000000000000000331562" ),
   MpIeee( "0.000000000000000000314" )
};
static cheb_series aig_cs = {
  aig_data,
  8,
  -1, 1,
  8
};

/*
 series for aip2 on the interval  0.00000e+00 to  1.25000e-01
                                        with weighted error   2.15e-17
                                         log weighted error  16.67
                               significant figures required  14.27
                                    decimal places required  17.26
*/
static MpIeee aip2_data[15] =  {
    MpIeee( "0.0065457691989713757" ),
    MpIeee( "0.0023833724120774592" ),
   -MpIeee( "0.0000430700770220586" ),
    MpIeee( "0.0000015629125858629" ),
   -MpIeee( "0.0000000815417186163" ),
    MpIeee( "0.0000000054103738057" ),
   -MpIeee( "0.0000000004284130883" ),
    MpIeee( "0.0000000000389497963" ),
   -MpIeee( "0.0000000000039623161" ),
    MpIeee( "0.0000000000004428184" ),
   -MpIeee( "0.0000000000000536297" ),
    MpIeee( "0.0000000000000069650" ),
   -MpIeee( "0.0000000000000009620" ),
    MpIeee( "0.0000000000000001403" ),
   -MpIeee( "0.0000000000000000215" )
};
static cheb_series aip2_cs = {
  aip2_data,
  14,
  -1, 1,
  9
};

/*
 series for aip1 on the interval  1.25000e-01 to  1.00000e+00
                                        with weighted error   2.60e-17
                                         log weighted error  16.58
                               significant figures required  14.91
                                    decimal places required  17.28
*/
static MpIeee aip1_data[25] =  {
    MpIeee( "0.0358865097808301538" ),
    MpIeee( "0.0114668575627764899" ),
   -MpIeee( "0.0007592073583861400" ),
    MpIeee( "0.0000869517610893841" ),
   -MpIeee( "0.0000128237294298592" ),
    MpIeee( "0.0000022062695681038" ),
   -MpIeee( "0.0000004222295185921" ),
    MpIeee( "0.0000000874686415726" ),
   -MpIeee( "0.0000000192773588418" ),
    MpIeee( "0.0000000044668460054" ),
   -MpIeee( "0.0000000010790108052" ),
    MpIeee( "0.0000000002700029447" ),
   -MpIeee( "0.0000000000696480108" ),
    MpIeee( "0.0000000000184489907" ),
   -MpIeee( "0.0000000000050027817" ),
    MpIeee( "0.0000000000013852243" ),
   -MpIeee( "0.0000000000003908218" ),
    MpIeee( "0.0000000000001121536" ),
   -MpIeee( "0.0000000000000326862" ),
    MpIeee( "0.0000000000000096619" ),
   -MpIeee( "0.0000000000000028935" ),
    MpIeee( "0.0000000000000008770" ),
   -MpIeee( "0.0000000000000002688" ),
    MpIeee( "0.0000000000000000832" ),
   -MpIeee( "0.0000000000000000260" )
};
static cheb_series aip1_cs = {
  aip1_data,
  24,
  -1, 1,
  14
};


/*
 series for bif on the interval -1.00000e+00 to  1.00000e+00
                                        with weighted error   9.05e-18
                                         log weighted error  17.04
                               significant figures required  15.83
                                    decimal places required  17.49
*/
static MpIeee bif_data[8] =  {
   MpIeee( "0.1153536790828570243" ),
   MpIeee( "0.0205007894049192875" ),
   MpIeee( "0.0002135290278902876" ),
   MpIeee( "0.0000010783960614677" ),
   MpIeee( "0.0000000032094708833" ),
   MpIeee( "0.0000000000062930407" ),
   MpIeee( "0.0000000000000087403" ),
   MpIeee( "0.0000000000000000090" )
};
static cheb_series bif_cs = {
  bif_data,
  7,
  -1, 1,
  7
};

/*
 series for big on the interval -1.00000e+00 to  1.00000e+00
                                        with weighted error   5.44e-19
                                         log weighted error  18.26
                               significant figures required  17.46
                                    decimal places required  18.74
*/
static MpIeee big_data[9] =  {
   -MpIeee( "0.097196440416443537390" ),
    MpIeee( "0.149503576843167066571" ),
    MpIeee( "0.003113525387121326042" ),
    MpIeee( "0.000024708570579821297" ),
    MpIeee( "0.000000102949627731379" ),
    MpIeee( "0.000000000263970373987" ),
    MpIeee( "0.000000000000458279271" ),
    MpIeee( "0.000000000000000574283" ),
    MpIeee( "0.000000000000000000544" )
};
static cheb_series big_cs = {
  big_data,
  8,
  -1, 1,
  8
};

/*
 series for bif2 on the interval  1.00000e+00 to  8.00000e+00
                                        with weighted error   3.82e-19
                                         log weighted error  18.42
                               significant figures required  17.68
                                    decimal places required  18.92
*/
static MpIeee bif2_data[10] =  {
   MpIeee( "0.323493987603522033521" ),
   MpIeee( "0.086297871535563559139" ),
   MpIeee( "0.002994025552655397426" ),
   MpIeee( "0.000051430528364661637" ),
   MpIeee( "0.000000525840250036811" ),
   MpIeee( "0.000000003561751373958" ),
   MpIeee( "0.000000000017146864007" ),
   MpIeee( "0.000000000000061663520" ),
   MpIeee( "0.000000000000000171911" ),
   MpIeee( "0.000000000000000000382" )
};
static cheb_series bif2_cs = {
  bif2_data,
  9,
  -1, 1,
  9
};

/*
 series for big2 on the interval  1.00000e+00 to  8.00000e+00
                                        with weighted error   3.35e-17
                                         log weighted error  16.48
                               significant figures required  16.52
                                    decimal places required  16.98
*/
static MpIeee big2_data[10] =  {
   MpIeee( "1.6062999463621294578" ),
   MpIeee( "0.7449088819876088652" ),
   MpIeee( "0.0470138738610277380" ),
   MpIeee( "0.0012284422062548239" ),
   MpIeee( "0.0000173222412256624" ),
   MpIeee( "0.0000001521901652368" ),
   MpIeee( "0.0000000009113560249" ),
   MpIeee( "0.0000000000039547918" ),
   MpIeee( "0.0000000000000130017" ),
   MpIeee( "0.0000000000000000335" )
};
static cheb_series big2_cs = {
  big2_data,
  9,
  -1, 1,
  9
};

/*
 series for bip2 on the interval  0.00000e+00 to  1.25000e-01
                                        with weighted error   2.07e-18
                                         log weighted error  17.69
                               significant figures required  16.51
                                    decimal places required  18.42
*/
static MpIeee bip2_data[29] =  {
    -MpIeee( "0.13269705443526630495" ),
    -MpIeee( "0.00568443626045977481" ),
    -MpIeee( "0.00015643601119611610" ),
    -MpIeee( "0.00001136737203679562" ),
    -MpIeee( "0.00000143464350991284" ),
    -MpIeee( "0.00000018098531185164" ),
     MpIeee( "0.00000000926177343611" ),
     MpIeee( "0.00000001710005490721" ),
     MpIeee( "0.00000000476698163504" ),
    -MpIeee( "0.00000000035195022023" ),
    -MpIeee( "0.00000000058890614316" ),
    -MpIeee( "0.00000000006678499608" ),
     MpIeee( "0.00000000006395565102" ),
     MpIeee( "0.00000000001554529427" ),
    -MpIeee( "0.00000000000792397000" ),
    -MpIeee( "0.00000000000258326243" ),
     MpIeee( "0.00000000000121655048" ),
     MpIeee( "0.00000000000038707207" ),
    -MpIeee( "0.00000000000022487045" ),
    -MpIeee( "0.00000000000004953477" ),
     MpIeee( "0.00000000000004563782" ),
     MpIeee( "0.00000000000000332998" ),
    -MpIeee( "0.00000000000000921750" ),
     MpIeee( "0.00000000000000094157" ),
     MpIeee( "0.00000000000000167154" ),
    -MpIeee( "0.00000000000000055134" ),
    -MpIeee( "0.00000000000000022369" ),
     MpIeee( "0.00000000000000017487" ),
     MpIeee( "0.00000000000000000207" )
};
static cheb_series bip2_cs = {
  bip2_data,
  28,
  -1, 1,
  14
};

/*
 series for bip1 on the interval  1.25000e-01 to  3.53553e-01
                                        with weighted error   1.86e-17
                                         log weighted error  16.73
                               significant figures required  15.67
                                    decimal places required  17.42
*/
static MpIeee bip1_data[24] =  {
   -MpIeee( "0.1729187351079553719" ),
   -MpIeee( "0.0149358492984694364" ),
   -MpIeee( "0.0005471104951678566" ),
    MpIeee( "0.0001537966292958408" ),
    MpIeee( "0.0000154353476192179" ),
   -MpIeee( "0.0000065434113851906" ),
    MpIeee( "0.0000003728082407879" ),
    MpIeee( "0.0000002072078388189" ),
   -MpIeee( "0.0000000658173336470" ),
    MpIeee( "0.0000000074926746354" ),
    MpIeee( "0.0000000011101336884" ),
   -MpIeee( "0.0000000007265140553" ),
    MpIeee( "0.0000000001782723560" ),
   -MpIeee( "0.0000000000217346352" ),
   -MpIeee( "0.0000000000020302035" ),
    MpIeee( "0.0000000000019311827" ),
   -MpIeee( "0.0000000000006044953" ),
    MpIeee( "0.0000000000001209450" ),
   -MpIeee( "0.0000000000000125109" ),
   -MpIeee( "0.0000000000000019917" ),
    MpIeee( "0.0000000000000015154" ),
   -MpIeee( "0.0000000000000004977" ),
    MpIeee( "0.0000000000000001155" ),
   -MpIeee( "0.0000000000000000186" )
};
static cheb_series bip1_cs = {
  bip1_data,
  23,
  -1, 1,
  13
};

/*
 series for an22 on the interval -1.00000e+00 to -1.25000e-01
                                        with weighted error   3.30e-17
                                         log weighted error  16.48
                               significant figures required  14.95
                                    decimal places required  17.24
*/
static MpIeee an22_data[33] =  {
    MpIeee( "0.0537418629629794329" ),
   -MpIeee( "0.0126661435859883193" ),
   -MpIeee( "0.0011924334106593007" ),
   -MpIeee( "0.0002032327627275655" ),
   -MpIeee( "0.0000446468963075164" ),
   -MpIeee( "0.0000113359036053123" ),
   -MpIeee( "0.0000031641352378546" ),
   -MpIeee( "0.0000009446708886149" ),
   -MpIeee( "0.0000002966562236472" ),
   -MpIeee( "0.0000000969118892024" ),
   -MpIeee( "0.0000000326822538653" ),
   -MpIeee( "0.0000000113144618964" ),
   -MpIeee( "0.0000000040042691002" ),
   -MpIeee( "0.0000000014440333684" ),
   -MpIeee( "0.0000000005292853746" ),
   -MpIeee( "0.0000000001967763374" ),
   -MpIeee( "0.0000000000740800096" ),
   -MpIeee( "0.0000000000282016314" ),
   -MpIeee( "0.0000000000108440066" ),
   -MpIeee( "0.0000000000042074801" ),
   -MpIeee( "0.0000000000016459150" ),
   -MpIeee( "0.0000000000006486827" ),
   -MpIeee( "0.0000000000002574095" ),
   -MpIeee( "0.0000000000001027889" ),
   -MpIeee( "0.0000000000000412846" ),
   -MpIeee( "0.0000000000000166711" ),
   -MpIeee( "0.0000000000000067657" ),
   -MpIeee( "0.0000000000000027585" ),
   -MpIeee( "0.0000000000000011296" ),
   -MpIeee( "0.0000000000000004645" ),
   -MpIeee( "0.0000000000000001917" ),
   -MpIeee( "0.0000000000000000794" ),
   -MpIeee( "0.0000000000000000330" )
};
static cheb_series an22_cs = {
  an22_data,
  32,
  -1, 1,
  18
};

/*
 series for an21 on the interval -1.25000e-01 to -1.56250e-02
                                        with weighted error   3.43e-17
                                         log weighted error  16.47
                               significant figures required  14.48
                                    decimal places required  17.16
*/
static MpIeee an21_data[24] =  {
    MpIeee( "0.0198313155263169394" ),
   -MpIeee( "0.0029376249067087533" ),
   -MpIeee( "0.0001136260695958196" ),
   -MpIeee( "0.0000100554451087156" ),
   -MpIeee( "0.0000013048787116563" ),
   -MpIeee( "0.0000002123881993151" ),
   -MpIeee( "0.0000000402270833384" ),
   -MpIeee( "0.0000000084996745953" ),
   -MpIeee( "0.0000000019514839426" ),
   -MpIeee( "0.0000000004783865344" ),
   -MpIeee( "0.0000000001236733992" ),
   -MpIeee( "0.0000000000334137486" ),
   -MpIeee( "0.0000000000093702824" ),
   -MpIeee( "0.0000000000027130128" ),
   -MpIeee( "0.0000000000008075954" ),
   -MpIeee( "0.0000000000002463214" ),
   -MpIeee( "0.0000000000000767656" ),
   -MpIeee( "0.0000000000000243883" ),
   -MpIeee( "0.0000000000000078831" ),
   -MpIeee( "0.0000000000000025882" ),
   -MpIeee( "0.0000000000000008619" ),
   -MpIeee( "0.0000000000000002908" ),
   -MpIeee( "0.0000000000000000993" ),
   -MpIeee( "0.0000000000000000343" )
};
static cheb_series an21_cs = {
  an21_data,
  23,
  -1, 1,
  12
};

/*
 series for an20 on the interval -1.56250e-02 to  0.00000e+00
                                        with weighted error   4.41e-17
                                         log weighted error  16.36
                               significant figures required  14.16
                                    decimal places required  16.96
*/
static MpIeee an20_data[16] =  {
    MpIeee( "0.0126732217145738027" ),
   -MpIeee( "0.0005212847072615621" ),
   -MpIeee( "0.0000052672111140370" ),
   -MpIeee( "0.0000001628202185026" ),
   -MpIeee( "0.0000000090991442687" ),
   -MpIeee( "0.0000000007438647126" ),
   -MpIeee( "0.0000000000795494752" ),
   -MpIeee( "0.0000000000104050944" ),
   -MpIeee( "0.0000000000015932426" ),
   -MpIeee( "0.0000000000002770648" ),
   -MpIeee( "0.0000000000000535343" ),
   -MpIeee( "0.0000000000000113062" ),
   -MpIeee( "0.0000000000000025772" ),
   -MpIeee( "0.0000000000000006278" ),
   -MpIeee( "0.0000000000000001621" ),
   -MpIeee( "0.0000000000000000441" )
};
static cheb_series an20_cs = {
  an20_data,
  15,
  -1, 1,
  8
};

/*
 series for aph2 on the interval -1.00000e+00 to -1.25000e-01
                                        with weighted error   2.94e-17
                                         log weighted error  16.53
                               significant figures required  15.58
                                    decimal places required  17.28
*/
static MpIeee aph2_data[32] =  {
   -MpIeee( "0.2057088719781465107" ),
    MpIeee( "0.0422196961357771922" ),
    MpIeee( "0.0020482560511207275" ),
    MpIeee( "0.0002607800735165006" ),
    MpIeee( "0.0000474824268004729" ),
    MpIeee( "0.0000105102756431612" ),
    MpIeee( "0.0000026353534014668" ),
    MpIeee( "0.0000007208824863499" ),
    MpIeee( "0.0000002103236664473" ),
    MpIeee( "0.0000000644975634555" ),
    MpIeee( "0.0000000205802377264" ),
    MpIeee( "0.0000000067836273921" ),
    MpIeee( "0.0000000022974015284" ),
    MpIeee( "0.0000000007961306765" ),
    MpIeee( "0.0000000002813860610" ),
    MpIeee( "0.0000000001011749057" ),
    MpIeee( "0.0000000000369306738" ),
    MpIeee( "0.0000000000136615066" ),
    MpIeee( "0.0000000000051142751" ),
    MpIeee( "0.0000000000019351689" ),
    MpIeee( "0.0000000000007393607" ),
    MpIeee( "0.0000000000002849792" ),
    MpIeee( "0.0000000000001107281" ),
    MpIeee( "0.0000000000000433412" ),
    MpIeee( "0.0000000000000170801" ),
    MpIeee( "0.0000000000000067733" ),
    MpIeee( "0.0000000000000027017" ),
    MpIeee( "0.0000000000000010835" ),
    MpIeee( "0.0000000000000004367" ),
    MpIeee( "0.0000000000000001769" ),
    MpIeee( "0.0000000000000000719" ),
    MpIeee( "0.0000000000000000294" )
};
static cheb_series aph2_cs = {
  aph2_data,
  31,
  -1, 1,
  16
};

/*
 series for aph1 on the interval -1.25000e-01 to -1.56250e-02
                                        with weighted error   6.38e-17
                                         log weighted error  16.20
                               significant figures required  14.91
                                    decimal places required  16.87
*/
static MpIeee aph1_data[22] =  {
  -MpIeee( "0.1024172908077571694" ),
   MpIeee( "0.0071697275146591248" ),
   MpIeee( "0.0001209959363122329" ),
   MpIeee( "0.0000073361512841220" ),
   MpIeee( "0.0000007535382954272" ),
   MpIeee( "0.0000001041478171741" ),
   MpIeee( "0.0000000174358728519" ),
   MpIeee( "0.0000000033399795033" ),
   MpIeee( "0.0000000007073075174" ),
   MpIeee( "0.0000000001619187515" ),
   MpIeee( "0.0000000000394539982" ),
   MpIeee( "0.0000000000101192282" ),
   MpIeee( "0.0000000000027092778" ),
   MpIeee( "0.0000000000007523806" ),
   MpIeee( "0.0000000000002156369" ),
   MpIeee( "0.0000000000000635283" ),
   MpIeee( "0.0000000000000191757" ),
   MpIeee( "0.0000000000000059143" ),
   MpIeee( "0.0000000000000018597" ),
   MpIeee( "0.0000000000000005950" ),
   MpIeee( "0.0000000000000001934" ),
   MpIeee( "0.0000000000000000638" )
};
static cheb_series aph1_cs = {
  aph1_data,
  21,
  -1, 1,
  10
};

/*
 series for aph0 on the interval -1.56250e-02 to  0.00000e+00
                                        with weighted error   2.29e-17
                                         log weighted error  16.64
                               significant figures required  15.27
                                    decimal places required  17.23
*/
static MpIeee aph0_data[15] =  {
 -MpIeee( "0.0855849241130933257" ),
  MpIeee( "0.0011214378867065261" ),
  MpIeee( "0.0000042721029353664" ),
  MpIeee( "0.0000000817607381483" ),
  MpIeee( "0.0000000033907645000" ),
  MpIeee( "0.0000000002253264423" ),
  MpIeee( "0.0000000000206284209" ),
  MpIeee( "0.0000000000023858763" ),
  MpIeee( "0.0000000000003301618" ),
  MpIeee( "0.0000000000000527010" ),
  MpIeee( "0.0000000000000094555" ),
  MpIeee( "0.0000000000000018709" ),
  MpIeee( "0.0000000000000004024" ),
  MpIeee( "0.0000000000000000930" ),
  MpIeee( "0.0000000000000000229" )
};
static cheb_series aph0_cs = {
  aph0_data,
  14,
  -1, 1,
  7
};


static
int
 airy_deriv_mod_phase(const MpIeee x, gsl_mode_t mode,
                     gsl_sf_result * ampl, gsl_sf_result * phi)
{
  const MpIeee pi34=  2.356194490192344928847;
  gsl_sf_result result_a;
  gsl_sf_result result_p;
  MpIeee a;MpIeee  p;
  MpIeee sqx;
  MpIeee x32;

  if(x <= -4.0) {
    MpIeee z=  MpIeee( "128.0" )/(x*x*x) + MpIeee( "1.0" );
    cheb_eval_mode_e(&an20_cs, z, mode, &result_a);
    cheb_eval_mode_e(&aph0_cs, z, mode, &result_p);
  }
  else if(x <= -2.0) {
    MpIeee z=  (MpIeee( "128.0" )/(x*x*x) + MpIeee( "9.0" )) / MpIeee( "7.0" );
    cheb_eval_mode_e(&an21_cs, z, mode, &result_a);
    cheb_eval_mode_e(&aph1_cs, z, mode, &result_p);
  }
  else if(x <= -1.0) {
    MpIeee z=  (MpIeee( "16.0" )/(x*x*x) + MpIeee( "9.0" )) / MpIeee( "7.0" );
    cheb_eval_mode_e(&an22_cs, z, mode, &result_a);
    cheb_eval_mode_e(&aph2_cs, z, mode, &result_p);
  }
  else {
    ampl->val = 0.0;
    ampl->err = 0.0;
    phi->val  = 0.0;
    phi->err  = 0.0;
    GSL_ERROR ("x is greater than 1.0", GSL_EDOM);
  }

  a =  MpIeee( "0.3125" ) + result_a.val;
  p = -MpIeee( "0.625" )  + result_p.val;
 
  sqx = sqrt(-x);
  x32   = x*sqx;

  ampl->val = sqrt(a * sqx);
  ampl->err = fabs(ampl->val) * (GSL_DBL_EPSILON + fabs(result_a.err/result_a.val));
  phi->val  = pi34 - x * sqx * p;
  phi->err = fabs(phi->val) * (GSL_DBL_EPSILON + fabs(result_p.err/result_p.val));

  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
 gsl_sf_airy_Ai_deriv_scaled_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int  status_ap=  airy_deriv_mod_phase(x, mode, &a, &p);
    MpIeee c=  cos(p.val);
    result->val  = a.val * c;
    result->err  = fabs(result->val * p.err) + fabs(c * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x <= 1.0) {
    const MpIeee x3=  x*x*x;
    const MpIeee x2=  x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs, x3, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs, x3, mode, &result_c1);

    result->val  = (x2*(0.125 + result_c0.val) - result_c1.val) - 0.25;
    result->err  = fabs(x2*result_c0.val) + result_c1.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    if(x > GSL_ROOT3_DBL_EPSILON*GSL_ROOT3_DBL_EPSILON) {
      /* scale only if x is positive */
      MpIeee s=  exp(MpIeee( "2.0" )*x*sqrt(x)/MpIeee( "3.0" ));
      result->val *= s;
      result->err *= s;
    }

    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const MpIeee sqrtx=  sqrt(x);
    const MpIeee z=  (16.0/(x*sqrtx) - 9.0)/7.0;
    const MpIeee s=  sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&aip1_cs, z, mode, &result_c0);
    result->val  = -(0.28125 + result_c0.val) * s;
    result->err  = result_c0.err * s;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee sqrtx=  sqrt(x);
    const MpIeee z=  16.0/(x*sqrtx) - 1.0;
    const MpIeee s=  sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&aip2_cs, z, mode, &result_c0);
    result->val  = -(0.28125 + result_c0.val) * s;
    result->err  = result_c0.err * s;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_airy_Ai_deriv_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int  status_ap=  airy_deriv_mod_phase(x, mode, &a, &p);
    MpIeee c=  cos(p.val);
    result->val  = a.val * c;
    result->err  = fabs(result->val * p.err) + fabs(c * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x < 1.0) {
    const MpIeee x3=  x*x*x;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&aif_cs, x3, mode, &result_c1);
    cheb_eval_mode_e(&aig_cs, x3, mode, &result_c2);
    result->val  = (x*x*(0.125 + result_c1.val) - result_c2.val) - 0.25;
    result->err  = fabs(x*x*result_c1.err) + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x*x*x < 9.0/4.0 * GSL_LOG_DBL_MIN*GSL_LOG_DBL_MIN) {
    gsl_sf_result result_aps;
    const MpIeee arg=  -2.0*x*sqrt(x)/3.0;
    const int stat_a = gsl_sf_airy_Ai_deriv_scaled_e(x, mode, &result_aps);
    const int stat_e = gsl_sf_exp_mult_err_e(arg, 1.5*fabs(arg*GSL_DBL_EPSILON),
                                                result_aps.val, result_aps.err,
                                                result);
    return GSL_ERROR_SELECT_2(stat_e, stat_a);
  }
  else {
    UNDERFLOW_ERROR(result);
  }
}


int
 gsl_sf_airy_Bi_deriv_scaled_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  const MpIeee atr=   8.7506905708484345;   /* 16./(sqrt(8)-1) */
  const MpIeee btr=  -2.0938363213560543;   /* -(sqrt(8)+1)/(sqrt(8)-1) */

  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int  status_ap=  airy_deriv_mod_phase(x, mode, &a, &p);
    MpIeee s=  sin(p.val);
    result->val  = a.val * s;
    result->err  = fabs(result->val * p.err) + fabs(s * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x < 1.0) {
    const MpIeee x3=  x*x*x;
    const MpIeee x2=  x*x;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&bif_cs, x3, mode, &result_c1);
    cheb_eval_mode_e(&big_cs, x3, mode, &result_c2);
    result->val  = x2 * (result_c1.val + 0.25) + result_c2.val + 0.5;
    result->err  = x2 * result_c1.err + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    if(x > GSL_ROOT3_DBL_EPSILON*GSL_ROOT3_DBL_EPSILON) {
      /* scale only if x is positive */
      const MpIeee s=  exp(-2.0*x*sqrt(x)/3.0);
      result->val *= s;
      result->err *= s;
    }

    return GSL_SUCCESS;
  }
  else if(x < 2.0) {
    const MpIeee z=  (2.0*x*x*x - 9.0) / 7.0;
    const MpIeee s=  exp(-2.0*x*sqrt(x)/3.0);
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c1);
    result->val  = s * (x*x * (0.25 + result_c0.val) + 0.5 + result_c1.val);
    result->err  = s * (x*x * result_c0.err + result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 4.0) {
    const MpIeee sqrtx=  sqrt(x);
    const MpIeee z=  atr/(x*sqrtx) + btr;
    const MpIeee s=  sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&bip1_cs, z, mode, &result_c0);
    result->val  = s * (0.625 + result_c0.val);
    result->err  = s * result_c0.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const MpIeee sqrtx=  sqrt(x);
    const MpIeee z=  16.0/(x*sqrtx) - 1.0;
    const MpIeee s=  sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&bip2_cs, z, mode, &result_c0);
    result->val  = s * (0.625 + result_c0.val);
    result->err  = s * result_c0.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
 gsl_sf_airy_Bi_deriv_e(const MpIeee x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int  status_ap=  airy_deriv_mod_phase(x, mode, &a, &p);
    MpIeee s=  sin(p.val);
    result->val  = a.val * s;
    result->err  = fabs(result->val * p.err) + fabs(s * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x < 1.0) {
    const MpIeee x3=  x*x*x;
    const MpIeee x2=  x*x;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&bif_cs, x3, mode, &result_c1);
    cheb_eval_mode_e(&big_cs, x3, mode, &result_c2);
    result->val  = x2 * (result_c1.val + 0.25) + result_c2.val + 0.5;
    result->err  = x2 * result_c1.err + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 2.0) {
    const MpIeee z=  (2.0*x*x*x - 9.0) / 7.0;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c1);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c2);
    result->val  = x*x * (result_c1.val + 0.25) + result_c2.val + 0.5;
    result->err  = x*x * result_c1.err + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_ROOT3_DBL_MAX*GSL_ROOT3_DBL_MAX) {
    gsl_sf_result result_bps;
    const MpIeee arg=  2.0*(x*sqrt(x)/3.0);
    int  stat_b=  gsl_sf_airy_Bi_deriv_scaled_e(x, mode, &result_bps);
    int  stat_e=  gsl_sf_exp_mult_err_e(arg, 1.5*fabs(arg*GSL_DBL_EPSILON),
                                          result_bps.val, result_bps.err,
                                          result);
    return GSL_ERROR_SELECT_2(stat_e, stat_b);
  }
  else {
    OVERFLOW_ERROR(result);
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_airy_Ai_deriv_scaled(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_deriv_scaled_e(x, mode, &result));
}

MpIeee gsl_sf_airy_Ai_deriv(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_deriv_e(x, mode, &result));
}

MpIeee gsl_sf_airy_Bi_deriv_scaled(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_deriv_scaled_e(x, mode, &result));
}

MpIeee gsl_sf_airy_Bi_deriv(const MpIeee x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_deriv_e(x, mode, &result));
}
