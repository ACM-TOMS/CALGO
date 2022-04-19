#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/zeta.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_zeta.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#define LogTwoPi_  1.8378770664093454835606594728111235279723


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* chebyshev fit for (s(t)-1)Zeta[s(t)]
 * s(t)= (t+1)/2
 * -1 <= t <= 1
 */
static MpIeee zeta_xlt1_data[14] =  {
  MpIeee( "1.48018677156931561235192914649" ),
  MpIeee( "0.25012062539889426471999938167" ),
  MpIeee( "0.00991137502135360774243761467" ),
 -MpIeee( "0.00012084759656676410329833091" ),
 -MpIeee( "4.7585866367662556504652535281e-06" ),
  MpIeee( "2.2229946694466391855561441361e-07" ),
 -MpIeee( "2.2237496498030257121309056582e-09" ),
 -MpIeee( "1.0173226513229028319420799028e-10" ),
  MpIeee( "4.3756643450424558284466248449e-12" ),
 -MpIeee( "6.2229632593100551465504090814e-14" ),
 -MpIeee( "6.6116201003272207115277520305e-16" ),
  MpIeee( "4.9477279533373912324518463830e-17" ),
 -MpIeee( "1.0429819093456189719660003522e-18" ),
  MpIeee( "6.9925216166580021051464412040e-21" ),
};
static cheb_series zeta_xlt1_cs = {
  zeta_xlt1_data,
  13,
  -1, 1,
  8
};

/* chebyshev fit for (s(t)-1)Zeta[s(t)]
 * s(t)= (19t+21)/2
 * -1 <= t <= 1
 */
static MpIeee zeta_xgt1_data[30] =  {
  MpIeee( "19.3918515726724119415911269006" ),
   MpIeee( "9.1525329692510756181581271500" ),
   MpIeee( "0.2427897658867379985365270155" ),
  -MpIeee( "0.1339000688262027338316641329" ),
   MpIeee( "0.0577827064065028595578410202" ),
  -MpIeee( "0.0187625983754002298566409700" ),
   MpIeee( "0.0039403014258320354840823803" ),
  -MpIeee( "0.0000581508273158127963598882" ),
  -MpIeee( "0.0003756148907214820704594549" ),
   MpIeee( "0.0001892530548109214349092999" ),
  -MpIeee( "0.0000549032199695513496115090" ),
   MpIeee( "8.7086484008939038610413331863e-6" ),
   MpIeee( "6.4609477924811889068410083425e-7" ),
  -MpIeee( "9.6749773915059089205835337136e-7" ),
   MpIeee( "3.6585400766767257736982342461e-7" ),
  -MpIeee( "8.4592516427275164351876072573e-8" ),
   MpIeee( "9.9956786144497936572288988883e-9" ),
   MpIeee( "1.4260036420951118112457144842e-9" ),
  -MpIeee( "1.1761968823382879195380320948e-9" ),
   MpIeee( "3.7114575899785204664648987295e-10" ),
  -MpIeee( "7.4756855194210961661210215325e-11" ),
   MpIeee( "7.8536934209183700456512982968e-12" ),
   MpIeee( "9.9827182259685539619810406271e-13" ),
  -MpIeee( "7.5276687030192221587850302453e-13" ),
   MpIeee( "2.1955026393964279988917878654e-13" ),
  -MpIeee( "4.1934859852834647427576319246e-14" ),
   MpIeee( "4.6341149635933550715779074274e-15" ),
   MpIeee( "2.3742488509048340106830309402e-16" ),
  -MpIeee( "2.7276516388124786119323824391e-16" ),
   MpIeee( "7.8473570134636044722154797225e-17" )
};
static cheb_series zeta_xgt1_cs = {
  zeta_xgt1_data,
  29,
  -1, 1,
  17
};


/* chebyshev fit for Ln[Zeta[s(t)] - 1 - 2^(-s(t))]
 * s(t)= 10 + 5t
 * -1 <= t <= 1; 5 <= s <= 15
 */
static MpIeee zetam1_inter_data[24] =  {
  -MpIeee( "21.7509435653088483422022339374" ),
  -MpIeee( "5.63036877698121782876372020472" ),
   MpIeee( "0.0528041358684229425504861579635" ),
  -MpIeee( "0.0156381809179670789342700883562" ),
   MpIeee( "0.00408218474372355881195080781927" ),
  -MpIeee( "0.0010264867349474874045036628282" ),
   MpIeee( "0.000260469880409886900143834962387" ),
  -MpIeee( "0.0000676175847209968878098566819447" ),
   MpIeee( "0.0000179284472587833525426660171124" ),
  -MpIeee( "4.83238651318556188834107605116e-6" ),
   MpIeee( "1.31913788964999288471371329447e-6" ),
  -MpIeee( "3.63760500656329972578222188542e-7" ),
   MpIeee( "1.01146847513194744989748396574e-7" ),
  -MpIeee( "2.83215225141806501619105289509e-8" ),
   MpIeee( "7.97733710252021423361012829496e-9" ),
  -MpIeee( "2.25850168553956886676250696891e-9" ),
   MpIeee( "6.42269392950164306086395744145e-10" ),
  -MpIeee( "1.83363861846127284505060843614e-10" ),
   MpIeee( "5.25309763895283179960368072104e-11" ),
  -MpIeee( "1.50958687042589821074710575446e-11" ),
   MpIeee( "4.34997545516049244697776942981e-12" ),
  -MpIeee( "1.25597782748190416118082322061e-12" ),
   MpIeee( "3.61280740072222650030134104162e-13" ),
  -MpIeee( "9.66437239205745207188920348801e-14" )
}; 
static cheb_series zetam1_inter_cs = {
  zetam1_inter_data,
  22,
  -1, 1,
  12
};



/* assumes s >= 0 and s != 1.0 */
inline
static int
 riemann_zeta_sgt0(MpIeee s, gsl_sf_result * result)
{
  if(s < MpIeee( "1.0" )) {
    gsl_sf_result c;
    cheb_eval_e(&zeta_xlt1_cs, 2.0*s - 1.0, &c);
    result->val = c.val / (s - 1.0);
    result->err = c.err / fabs(s-1.0) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(s <= MpIeee( "20.0" )) {
    MpIeee x=  (MpIeee( "2.0" )*s - MpIeee( "21.0" ))/MpIeee( "19.0" );
    gsl_sf_result c;
    cheb_eval_e(&zeta_xgt1_cs, x, &c);
    result->val = c.val / (s - 1.0);
    result->err = c.err / (s - 1.0) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee f2=  MpIeee( "1.0" ) - pow(MpIeee( "2.0" ),-s);
    MpIeee f3=  MpIeee( "1.0" ) - pow(MpIeee( "3.0" ),-s);
    MpIeee f5=  MpIeee( "1.0" ) - pow(MpIeee( "5.0" ),-s);
    MpIeee f7=  MpIeee( "1.0" ) - pow(MpIeee( "7.0" ),-s);
    result->val = 1.0/(f2*f3*f5*f7);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


inline
static int
 riemann_zeta1ms_slt0(MpIeee s, gsl_sf_result * result)
{
  if(s > -MpIeee( "19.0" )) {
    MpIeee x=  (-MpIeee( "19" ) - MpIeee( "2.0" )*s)/MpIeee( "19.0" );
    gsl_sf_result c;
    cheb_eval_e(&zeta_xgt1_cs, x, &c);
    result->val = c.val / (-s);
    result->err = c.err / (-s) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    MpIeee f2=  MpIeee( "1.0" ) - pow(MpIeee( "2.0" ),-(MpIeee( "1.0" )-s));
    MpIeee f3=  MpIeee( "1.0" ) - pow(MpIeee( "3.0" ),-(MpIeee( "1.0" )-s));
    MpIeee f5=  MpIeee( "1.0" ) - pow(MpIeee( "5.0" ),-(MpIeee( "1.0" )-s));
    MpIeee f7=  MpIeee( "1.0" ) - pow(MpIeee( "7.0" ),-(MpIeee( "1.0" )-s));
    result->val = 1.0/(f2*f3*f5*f7);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/* works for 5 < s < 15*/
static int
 riemann_zeta_minus_1_intermediate_s(MpIeee s, gsl_sf_result * result)
{
  MpIeee t=  (s - MpIeee( "10.0" ))/MpIeee( "5.0" );
  gsl_sf_result c;
  cheb_eval_e(&zetam1_inter_cs, t, &c);
  result->val = exp(c.val) + pow(2.0, -s);
  result->err = (c.err + 2.0*GSL_DBL_EPSILON)*result->val;
  return GSL_SUCCESS;
}


/* assumes s is large and positive
 * write: zeta(s) - 1 = zeta(s) * (1 - 1/zeta(s))
 * and expand a few terms of the product formula to evaluate 1 - 1/zeta(s)
 *
 * works well for s > 15
 */
static int
 riemann_zeta_minus1_large_s(MpIeee s, gsl_sf_result * result)
{
  MpIeee a=  pow( MpIeee( "2.0" ),-s);
  MpIeee b=  pow( MpIeee( "3.0" ),-s);
  MpIeee c=  pow( MpIeee( "5.0" ),-s);
  MpIeee d=  pow( MpIeee( "7.0" ),-s);
  MpIeee e=  pow(MpIeee( "11.0" ),-s);
  MpIeee f=  pow(MpIeee( "13.0" ),-s);
  MpIeee t1=  a + b + c + d + e + f;
  MpIeee t2=  a*(b+c+d+e+f) + b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f;
  /*
  double t3 = a*(b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f) +
              b*(c*(d+e+f) + d*(e+f) + e*f) +
              c*(d*(e+f) + e*f) +
              d*e*f;
  double t4 = a*(b*(c*(d + e + f) + d*(e + f) + e*f) + c*(d*(e+f) + e*f) + d*e*f) +
              b*(c*(d*(e+f) + e*f) + d*e*f) +
              c*d*e*f;
  double t5 = b*c*d*e*f + a*c*d*e*f+ a*b*d*e*f+ a*b*c*e*f+ a*b*c*d*f+ a*b*c*d*e;
  double t6 = a*b*c*d*e*f;
  */
  MpIeee numt=  t1 - t2 /* + t3 - t4 + t5 - t6 */;
  MpIeee zeta=  MpIeee( "1.0" )/((MpIeee( "1.0" )-a)*(MpIeee( "1.0" )-b)*(MpIeee( "1.0" )-c)*(MpIeee( "1.0" )-d)*(MpIeee( "1.0" )-e)*(MpIeee( "1.0" )-f));
  result->val = numt*zeta;
  result->err = (15.0/s + 1.0) * 6.0*GSL_DBL_EPSILON*result->val;
  return GSL_SUCCESS;
}


#if 0
/* zeta(n) */
#define ZETA_POS_TABLE_NMAX   100
static MpIeee zeta_pos_int_table_OLD[ZETA_POS_TABLE_NMAX+1] =  {
 -MpIeee( "0.50000000000000000000000000000" ),       /* zeta(0) */
  MpIeee( "0.0" ) /* FIXME: DirectedInfinity() */,   /* zeta(1) */
  MpIeee( "1.64493406684822643647241516665" ),       /* ...     */
  MpIeee( "1.20205690315959428539973816151" ),
  MpIeee( "1.08232323371113819151600369654" ),
  MpIeee( "1.03692775514336992633136548646" ),
  MpIeee( "1.01734306198444913971451792979" ),
  MpIeee( "1.00834927738192282683979754985" ),
  MpIeee( "1.00407735619794433937868523851" ),
  MpIeee( "1.00200839282608221441785276923" ),
  MpIeee( "1.00099457512781808533714595890" ),
  MpIeee( "1.00049418860411946455870228253" ),
  MpIeee( "1.00024608655330804829863799805" ),
  MpIeee( "1.00012271334757848914675183653" ),
  MpIeee( "1.00006124813505870482925854511" ),
  MpIeee( "1.00003058823630702049355172851" ),
  MpIeee( "1.00001528225940865187173257149" ),
  MpIeee( "1.00000763719763789976227360029" ),
  MpIeee( "1.00000381729326499983985646164" ),
  MpIeee( "1.00000190821271655393892565696" ),
  MpIeee( "1.00000095396203387279611315204" ),
  MpIeee( "1.00000047693298678780646311672" ),
  MpIeee( "1.00000023845050272773299000365" ),
  MpIeee( "1.00000011921992596531107306779" ),
  MpIeee( "1.00000005960818905125947961244" ),
  MpIeee( "1.00000002980350351465228018606" ),
  MpIeee( "1.00000001490155482836504123466" ),
  MpIeee( "1.00000000745071178983542949198" ),
  MpIeee( "1.00000000372533402478845705482" ),
  MpIeee( "1.00000000186265972351304900640" ),
  MpIeee( "1.00000000093132743241966818287" ),
  MpIeee( "1.00000000046566290650337840730" ),
  MpIeee( "1.00000000023283118336765054920" ),
  MpIeee( "1.00000000011641550172700519776" ),
  MpIeee( "1.00000000005820772087902700889" ),
  MpIeee( "1.00000000002910385044497099687" ),
  MpIeee( "1.00000000001455192189104198424" ),
  MpIeee( "1.00000000000727595983505748101" ),
  MpIeee( "1.00000000000363797954737865119" ),
  MpIeee( "1.00000000000181898965030706595" ),
  MpIeee( "1.00000000000090949478402638893" ),
  MpIeee( "1.00000000000045474737830421540" ),
  MpIeee( "1.00000000000022737368458246525" ),
  MpIeee( "1.00000000000011368684076802278" ),
  MpIeee( "1.00000000000005684341987627586" ),
  MpIeee( "1.00000000000002842170976889302" ),
  MpIeee( "1.00000000000001421085482803161" ),
  MpIeee( "1.00000000000000710542739521085" ),
  MpIeee( "1.00000000000000355271369133711" ),
  MpIeee( "1.00000000000000177635684357912" ),
  MpIeee( "1.00000000000000088817842109308" ),
  MpIeee( "1.00000000000000044408921031438" ),
  MpIeee( "1.00000000000000022204460507980" ),
  MpIeee( "1.00000000000000011102230251411" ),
  MpIeee( "1.00000000000000005551115124845" ),
  MpIeee( "1.00000000000000002775557562136" ),
  MpIeee( "1.00000000000000001387778780973" ),
  MpIeee( "1.00000000000000000693889390454" ),
  MpIeee( "1.00000000000000000346944695217" ),
  MpIeee( "1.00000000000000000173472347605" ),
  MpIeee( "1.00000000000000000086736173801" ),
  MpIeee( "1.00000000000000000043368086900" ),
  MpIeee( "1.00000000000000000021684043450" ),
  MpIeee( "1.00000000000000000010842021725" ),
  MpIeee( "1.00000000000000000005421010862" ),
  MpIeee( "1.00000000000000000002710505431" ),
  MpIeee( "1.00000000000000000001355252716" ),
  MpIeee( "1.00000000000000000000677626358" ),
  MpIeee( "1.00000000000000000000338813179" ),
  MpIeee( "1.00000000000000000000169406589" ),
  MpIeee( "1.00000000000000000000084703295" ),
  MpIeee( "1.00000000000000000000042351647" ),
  MpIeee( "1.00000000000000000000021175824" ),
  MpIeee( "1.00000000000000000000010587912" ),
  MpIeee( "1.00000000000000000000005293956" ),
  MpIeee( "1.00000000000000000000002646978" ),
  MpIeee( "1.00000000000000000000001323489" ),
  MpIeee( "1.00000000000000000000000661744" ),
  MpIeee( "1.00000000000000000000000330872" ),
  MpIeee( "1.00000000000000000000000165436" ),
  MpIeee( "1.00000000000000000000000082718" ),
  MpIeee( "1.00000000000000000000000041359" ),
  MpIeee( "1.00000000000000000000000020680" ),
  MpIeee( "1.00000000000000000000000010340" ),
  MpIeee( "1.00000000000000000000000005170" ),
  MpIeee( "1.00000000000000000000000002585" ),
  MpIeee( "1.00000000000000000000000001292" ),
  MpIeee( "1.00000000000000000000000000646" ),
  MpIeee( "1.00000000000000000000000000323" ),
  MpIeee( "1.00000000000000000000000000162" ),
  MpIeee( "1.00000000000000000000000000081" ),
  MpIeee( "1.00000000000000000000000000040" ),
  MpIeee( "1.00000000000000000000000000020" ),
  MpIeee( "1.00000000000000000000000000010" ),
  MpIeee( "1.00000000000000000000000000005" ),
  MpIeee( "1.00000000000000000000000000003" ),
  MpIeee( "1.00000000000000000000000000001" ),
  MpIeee( "1.00000000000000000000000000001" ),
  MpIeee( "1.00000000000000000000000000000" ),
  MpIeee( "1.00000000000000000000000000000" ),
  MpIeee( "1.00000000000000000000000000000" )
};
#endif /* 0 */


/* zeta(n) - 1 */
#define ZETA_POS_TABLE_NMAX   100
static MpIeee zetam1_pos_int_table[ZETA_POS_TABLE_NMAX+1] =  {
 -MpIeee( "1.5" ),                               /* zeta(0) */
  MpIeee( "0.0" ),       /* FIXME: Infinity */   /* zeta(1) - 1 */
  MpIeee( "0.644934066848226436472415166646" ),  /* zeta(2) - 1 */
  MpIeee( "0.202056903159594285399738161511" ),
  MpIeee( "0.082323233711138191516003696541" ),
  MpIeee( "0.036927755143369926331365486457" ),
  MpIeee( "0.017343061984449139714517929790" ),
  MpIeee( "0.008349277381922826839797549849" ),
  MpIeee( "0.004077356197944339378685238508" ),
  MpIeee( "0.002008392826082214417852769232" ),
  MpIeee( "0.000994575127818085337145958900" ),
  MpIeee( "0.000494188604119464558702282526" ),
  MpIeee( "0.000246086553308048298637998047" ),
  MpIeee( "0.000122713347578489146751836526" ),
  MpIeee( "0.000061248135058704829258545105" ),
  MpIeee( "0.000030588236307020493551728510" ),
  MpIeee( "0.000015282259408651871732571487" ),
  MpIeee( "7.6371976378997622736002935630e-6" ),
  MpIeee( "3.8172932649998398564616446219e-6" ),
  MpIeee( "1.9082127165539389256569577951e-6" ),
  MpIeee( "9.5396203387279611315203868344e-7" ),
  MpIeee( "4.7693298678780646311671960437e-7" ),
  MpIeee( "2.3845050272773299000364818675e-7" ),
  MpIeee( "1.1921992596531107306778871888e-7" ),
  MpIeee( "5.9608189051259479612440207935e-8" ),
  MpIeee( "2.9803503514652280186063705069e-8" ),
  MpIeee( "1.4901554828365041234658506630e-8" ),
  MpIeee( "7.4507117898354294919810041706e-9" ),
  MpIeee( "3.7253340247884570548192040184e-9" ),
  MpIeee( "1.8626597235130490064039099454e-9" ),
  MpIeee( "9.3132743241966818287176473502e-10" ),
  MpIeee( "4.6566290650337840729892332512e-10" ),
  MpIeee( "2.3283118336765054920014559759e-10" ),
  MpIeee( "1.1641550172700519775929738354e-10" ),
  MpIeee( "5.8207720879027008892436859891e-11" ),
  MpIeee( "2.9103850444970996869294252278e-11" ),
  MpIeee( "1.4551921891041984235929632245e-11" ),
  MpIeee( "7.2759598350574810145208690123e-12" ),
  MpIeee( "3.6379795473786511902372363558e-12" ),
  MpIeee( "1.8189896503070659475848321007e-12" ),
  MpIeee( "9.0949478402638892825331183869e-13" ),
  MpIeee( "4.5474737830421540267991120294e-13" ),
  MpIeee( "2.2737368458246525152268215779e-13" ),
  MpIeee( "1.1368684076802278493491048380e-13" ),
  MpIeee( "5.6843419876275856092771829675e-14" ),
  MpIeee( "2.8421709768893018554550737049e-14" ),
  MpIeee( "1.4210854828031606769834307141e-14" ),
  MpIeee( "7.1054273952108527128773544799e-15" ),
  MpIeee( "3.5527136913371136732984695340e-15" ),
  MpIeee( "1.7763568435791203274733490144e-15" ),
  MpIeee( "8.8817842109308159030960913863e-16" ),
  MpIeee( "4.4408921031438133641977709402e-16" ),
  MpIeee( "2.2204460507980419839993200942e-16" ),
  MpIeee( "1.1102230251410661337205445699e-16" ),
  MpIeee( "5.5511151248454812437237365905e-17" ),
  MpIeee( "2.7755575621361241725816324538e-17" ),
  MpIeee( "1.3877787809725232762839094906e-17" ),
  MpIeee( "6.9388939045441536974460853262e-18" ),
  MpIeee( "3.4694469521659226247442714961e-18" ),
  MpIeee( "1.7347234760475765720489729699e-18" ),
  MpIeee( "8.6736173801199337283420550673e-19" ),
  MpIeee( "4.3368086900206504874970235659e-19" ),
  MpIeee( "2.1684043449972197850139101683e-19" ),
  MpIeee( "1.0842021724942414063012711165e-19" ),
  MpIeee( "5.4210108624566454109187004043e-20" ),
  MpIeee( "2.7105054312234688319546213119e-20" ),
  MpIeee( "1.3552527156101164581485233996e-20" ),
  MpIeee( "6.7762635780451890979952987415e-21" ),
  MpIeee( "3.3881317890207968180857031004e-21" ),
  MpIeee( "1.6940658945097991654064927471e-21" ),
  MpIeee( "8.4703294725469983482469926091e-22" ),
  MpIeee( "4.2351647362728333478622704833e-22" ),
  MpIeee( "2.1175823681361947318442094398e-22" ),
  MpIeee( "1.0587911840680233852265001539e-22" ),
  MpIeee( "5.2939559203398703238139123029e-23" ),
  MpIeee( "2.6469779601698529611341166842e-23" ),
  MpIeee( "1.3234889800848990803094510250e-23" ),
  MpIeee( "6.6174449004244040673552453323e-24" ),
  MpIeee( "3.3087224502121715889469563843e-24" ),
  MpIeee( "1.6543612251060756462299236771e-24" ),
  MpIeee( "8.2718061255303444036711056167e-25" ),
  MpIeee( "4.1359030627651609260093824555e-25" ),
  MpIeee( "2.0679515313825767043959679193e-25" ),
  MpIeee( "1.0339757656912870993284095591e-25" ),
  MpIeee( "5.1698788284564313204101332166e-26" ),
  MpIeee( "2.5849394142282142681277617708e-26" ),
  MpIeee( "1.2924697071141066700381126118e-26" ),
  MpIeee( "6.4623485355705318034380021611e-27" ),
  MpIeee( "3.2311742677852653861348141180e-27" ),
  MpIeee( "1.6155871338926325212060114057e-27" ),
  MpIeee( "8.0779356694631620331587381863e-28" ),
  MpIeee( "4.0389678347315808256222628129e-28" ),
  MpIeee( "2.0194839173657903491587626465e-28" ),
  MpIeee( "1.0097419586828951533619250700e-28" ),
  MpIeee( "5.0487097934144756960847711725e-29" ),
  MpIeee( "2.5243548967072378244674341938e-29" ),
  MpIeee( "1.2621774483536189043753999660e-29" ),
  MpIeee( "6.3108872417680944956826093943e-30" ),
  MpIeee( "3.1554436208840472391098412184e-30" ),
  MpIeee( "1.5777218104420236166444327830e-30" ),
  MpIeee( "7.8886090522101180735205378276e-31" )
};


#define ZETA_NEG_TABLE_NMAX  99
#define ZETA_NEG_TABLE_SIZE  50
static MpIeee zeta_neg_int_table[ZETA_NEG_TABLE_SIZE] =  {
 -MpIeee( "0.083333333333333333333333333333" ),     /* zeta(-1) */
  MpIeee( "0.008333333333333333333333333333" ),     /* zeta(-3) */
 -MpIeee( "0.003968253968253968253968253968" ),     /* ...      */
  MpIeee( "0.004166666666666666666666666667" ),
 -MpIeee( "0.007575757575757575757575757576" ),
  MpIeee( "0.021092796092796092796092796093" ),
 -MpIeee( "0.083333333333333333333333333333" ),
  MpIeee( "0.44325980392156862745098039216" ),
 -MpIeee( "3.05395433027011974380395433027" ),
  MpIeee( "26.4562121212121212121212121212" ),
 -MpIeee( "281.460144927536231884057971014" ),
  MpIeee( "3607.5105463980463980463980464" ),
 -MpIeee( "54827.583333333333333333333333" ),
  MpIeee( "974936.82385057471264367816092" ),
 -MpIeee( "2.0052695796688078946143462272e+07" ),
  MpIeee( "4.7238486772162990196078431373e+08" ),
 -MpIeee( "1.2635724795916666666666666667e+10" ),
  MpIeee( "3.8087931125245368811553022079e+11" ),
 -MpIeee( "1.2850850499305083333333333333e+13" ),
  MpIeee( "4.8241448354850170371581670362e+14" ),
 -MpIeee( "2.0040310656516252738108421663e+16" ),
  MpIeee( "9.1677436031953307756992753623e+17" ),
 -MpIeee( "4.5979888343656503490437943262e+19" ),
  MpIeee( "2.5180471921451095697089023320e+21" ),
 -MpIeee( "1.5001733492153928733711440151e+23" ),
  MpIeee( "9.6899578874635940656497942895e+24" ),
 -MpIeee( "6.7645882379292820990945242302e+26" ),
  MpIeee( "5.0890659468662289689766332916e+28" ),
 -MpIeee( "4.1147288792557978697665486068e+30" ),
  MpIeee( "3.5666582095375556109684574609e+32" ),
 -MpIeee( "3.3066089876577576725680214670e+34" ),
  MpIeee( "3.2715634236478716264211227016e+36" ),
 -MpIeee( "3.4473782558278053878256455080e+38" ),
  MpIeee( "3.8614279832705258893092720200e+40" ),
 -MpIeee( "4.5892974432454332168863989006e+42" ),
  MpIeee( "5.7775386342770431824884825688e+44" ),
 -MpIeee( "7.6919858759507135167410075972e+46" ),
  MpIeee( "1.0813635449971654696354033351e+49" ),
 -MpIeee( "1.6029364522008965406067102346e+51" ),
  MpIeee( "2.5019479041560462843656661499e+53" ),
 -MpIeee( "4.1067052335810212479752045004e+55" ),
  MpIeee( "7.0798774408494580617452972433e+57" ),
 -MpIeee( "1.2804546887939508790190849756e+60" ),
  MpIeee( "2.4267340392333524078020892067e+62" ),
 -MpIeee( "4.8143218874045769355129570066e+64" ),
  MpIeee( "9.9875574175727530680652777408e+66" ),
 -MpIeee( "2.1645634868435185631335136160e+69" ),
  MpIeee( "4.8962327039620553206849224516e+71" ),    /* ...        */
 -MpIeee( "1.1549023923963519663954271692e+74" ),    /* zeta(-97)  */
  MpIeee( "2.8382249570693706959264156336e+76" )     /* zeta(-99)  */
};


/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static MpIeee hzeta_c[15] =  {
  MpIeee( "1.00000000000000000000000000000" ),
  MpIeee( "0.083333333333333333333333333333" ),
 -MpIeee( "0.00138888888888888888888888888889" ),
  MpIeee( "0.000033068783068783068783068783069" ),
 -MpIeee( "8.2671957671957671957671957672e-07" ),
  MpIeee( "2.0876756987868098979210090321e-08" ),
 -MpIeee( "5.2841901386874931848476822022e-10" ),
  MpIeee( "1.3382536530684678832826980975e-11" ),
 -MpIeee( "3.3896802963225828668301953912e-13" ),
  MpIeee( "8.5860620562778445641359054504e-15" ),
 -MpIeee( "2.1748686985580618730415164239e-16" ),
  MpIeee( "5.5090028283602295152026526089e-18" ),
 -MpIeee( "1.3954464685812523340707686264e-19" ),
  MpIeee( "3.5347070396294674716932299778e-21" ),
 -MpIeee( "8.9535174270375468504026113181e-23" )
};

#define ETA_POS_TABLE_NMAX  100
static MpIeee eta_pos_int_table[ETA_POS_TABLE_NMAX+1] =  {
MpIeee( "0.50000000000000000000000000000" ),  /* eta(0) */
M_LN2,                            /* eta(1) */
MpIeee( "0.82246703342411321823620758332" ),  /* ...    */
MpIeee( "0.90154267736969571404980362113" ),
MpIeee( "0.94703282949724591757650323447" ),
MpIeee( "0.97211977044690930593565514355" ),
MpIeee( "0.98555109129743510409843924448" ),
MpIeee( "0.99259381992283028267042571313" ),
MpIeee( "0.99623300185264789922728926008" ),
MpIeee( "0.99809429754160533076778303185" ),
MpIeee( "0.99903950759827156563922184570" ),
MpIeee( "0.99951714349806075414409417483" ),
MpIeee( "0.99975768514385819085317967871" ),
MpIeee( "0.99987854276326511549217499282" ),
MpIeee( "0.99993917034597971817095419226" ),
MpIeee( "0.99996955121309923808263293263" ),
MpIeee( "0.99998476421490610644168277496" ),
MpIeee( "0.99999237829204101197693787224" ),
MpIeee( "0.99999618786961011347968922641" ),
MpIeee( "0.99999809350817167510685649297" ),
MpIeee( "0.99999904661158152211505084256" ),
MpIeee( "0.99999952325821554281631666433" ),
MpIeee( "0.99999976161323082254789720494" ),
MpIeee( "0.99999988080131843950322382485" ),
MpIeee( "0.99999994039889239462836140314" ),
MpIeee( "0.99999997019885696283441513311" ),
MpIeee( "0.99999998509923199656878766181" ),
MpIeee( "0.99999999254955048496351585274" ),
MpIeee( "0.99999999627475340010872752767" ),
MpIeee( "0.99999999813736941811218674656" ),
MpIeee( "0.99999999906868228145397862728" ),
MpIeee( "0.99999999953434033145421751469" ),
MpIeee( "0.99999999976716989595149082282" ),
MpIeee( "0.99999999988358485804603047265" ),
MpIeee( "0.99999999994179239904531592388" ),
MpIeee( "0.99999999997089618952980952258" ),
MpIeee( "0.99999999998544809143388476396" ),
MpIeee( "0.99999999999272404460658475006" ),
MpIeee( "0.99999999999636202193316875550" ),
MpIeee( "0.99999999999818101084320873555" ),
MpIeee( "0.99999999999909050538047887809" ),
MpIeee( "0.99999999999954525267653087357" ),
MpIeee( "0.99999999999977262633369589773" ),
MpIeee( "0.99999999999988631316532476488" ),
MpIeee( "0.99999999999994315658215465336" ),
MpIeee( "0.99999999999997157829090808339" ),
MpIeee( "0.99999999999998578914539762720" ),
MpIeee( "0.99999999999999289457268000875" ),
MpIeee( "0.99999999999999644728633373609" ),
MpIeee( "0.99999999999999822364316477861" ),
MpIeee( "0.99999999999999911182158169283" ),
MpIeee( "0.99999999999999955591079061426" ),
MpIeee( "0.99999999999999977795539522974" ),
MpIeee( "0.99999999999999988897769758908" ),
MpIeee( "0.99999999999999994448884878594" ),
MpIeee( "0.99999999999999997224442439010" ),
MpIeee( "0.99999999999999998612221219410" ),
MpIeee( "0.99999999999999999306110609673" ),
MpIeee( "0.99999999999999999653055304826" ),
MpIeee( "0.99999999999999999826527652409" ),
MpIeee( "0.99999999999999999913263826204" ),
MpIeee( "0.99999999999999999956631913101" ),
MpIeee( "0.99999999999999999978315956551" ),
MpIeee( "0.99999999999999999989157978275" ),
MpIeee( "0.99999999999999999994578989138" ),
MpIeee( "0.99999999999999999997289494569" ),
MpIeee( "0.99999999999999999998644747284" ),
MpIeee( "0.99999999999999999999322373642" ),
MpIeee( "0.99999999999999999999661186821" ),
MpIeee( "0.99999999999999999999830593411" ),
MpIeee( "0.99999999999999999999915296705" ),
MpIeee( "0.99999999999999999999957648353" ),
MpIeee( "0.99999999999999999999978824176" ),
MpIeee( "0.99999999999999999999989412088" ),
MpIeee( "0.99999999999999999999994706044" ),
MpIeee( "0.99999999999999999999997353022" ),
MpIeee( "0.99999999999999999999998676511" ),
MpIeee( "0.99999999999999999999999338256" ),
MpIeee( "0.99999999999999999999999669128" ),
MpIeee( "0.99999999999999999999999834564" ),
MpIeee( "0.99999999999999999999999917282" ),
MpIeee( "0.99999999999999999999999958641" ),
MpIeee( "0.99999999999999999999999979320" ),
MpIeee( "0.99999999999999999999999989660" ),
MpIeee( "0.99999999999999999999999994830" ),
MpIeee( "0.99999999999999999999999997415" ),
MpIeee( "0.99999999999999999999999998708" ),
MpIeee( "0.99999999999999999999999999354" ),
MpIeee( "0.99999999999999999999999999677" ),
MpIeee( "0.99999999999999999999999999838" ),
MpIeee( "0.99999999999999999999999999919" ),
MpIeee( "0.99999999999999999999999999960" ),
MpIeee( "0.99999999999999999999999999980" ),
MpIeee( "0.99999999999999999999999999990" ),
MpIeee( "0.99999999999999999999999999995" ),
MpIeee( "0.99999999999999999999999999997" ),
MpIeee( "0.99999999999999999999999999999" ),
MpIeee( "0.99999999999999999999999999999" ),
MpIeee( "1.00000000000000000000000000000" ),
MpIeee( "1.00000000000000000000000000000" ),
MpIeee( "1.00000000000000000000000000000" ),
};


#define ETA_NEG_TABLE_NMAX  99
#define ETA_NEG_TABLE_SIZE  50
static MpIeee eta_neg_int_table[ETA_NEG_TABLE_SIZE] =  {
 MpIeee( "0.25000000000000000000000000000" ),   /* eta(-1) */
-MpIeee( "0.12500000000000000000000000000" ),   /* eta(-3) */
 MpIeee( "0.25000000000000000000000000000" ),   /* ...      */
-MpIeee( "1.06250000000000000000000000000" ),
 MpIeee( "7.75000000000000000000000000000" ),
-MpIeee( "86.3750000000000000000000000000" ),
 MpIeee( "1365.25000000000000000000000000" ),
-MpIeee( "29049.0312500000000000000000000" ),
 MpIeee( "800572.750000000000000000000000" ),
-MpIeee( "2.7741322625000000000000000000e+7" ),
 MpIeee( "1.1805291302500000000000000000e+9" ),
-MpIeee( "6.0523980051687500000000000000e+10" ),
 MpIeee( "3.6794167785377500000000000000e+12" ),
-MpIeee( "2.6170760990658387500000000000e+14" ),
 MpIeee( "2.1531418140800295250000000000e+16" ),
-MpIeee( "2.0288775575173015930156250000e+18" ),
 MpIeee( "2.1708009902623770590275000000e+20" ),
-MpIeee( "2.6173826968455814932120125000e+22" ),
 MpIeee( "3.5324148876863877826668602500e+24" ),
-MpIeee( "5.3042033406864906641493838981e+26" ),
 MpIeee( "8.8138218364311576767253114668e+28" ),
-MpIeee( "1.6128065107490778547354654864e+31" ),
 MpIeee( "3.2355470001722734208527794569e+33" ),
-MpIeee( "7.0876727476537493198506645215e+35" ),
 MpIeee( "1.6890450341293965779175629389e+38" ),
-MpIeee( "4.3639690731216831157655651358e+40" ),
 MpIeee( "1.2185998827061261322605065672e+43" ),
-MpIeee( "3.6670584803153006180101262324e+45" ),
 MpIeee( "1.1859898526302099104271449748e+48" ),
-MpIeee( "4.1120769493584015047981746438e+50" ),
 MpIeee( "1.5249042436787620309090168687e+53" ),
-MpIeee( "6.0349693196941307074572991901e+55" ),
 MpIeee( "2.5437161764210695823197691519e+58" ),
-MpIeee( "1.1396923802632287851130360170e+61" ),
 MpIeee( "5.4180861064753979196802726455e+63" ),
-MpIeee( "2.7283654799994373847287197104e+66" ),
 MpIeee( "1.4529750514918543238511171663e+69" ),
-MpIeee( "8.1705519371067450079777183386e+71" ),
 MpIeee( "4.8445781606678367790247757259e+74" ),
-MpIeee( "3.0246694206649519336179448018e+77" ),
 MpIeee( "1.9858807961690493054169047970e+80" ),
-MpIeee( "1.3694474620720086994386818232e+83" ),
 MpIeee( "9.9070382984295807826303785989e+85" ),
-MpIeee( "7.5103780796592645925968460677e+88" ),
 MpIeee( "5.9598418264260880840077992227e+91" ),
-MpIeee( "4.9455988887500020399263196307e+94" ),
 MpIeee( "4.2873596927020241277675775935e+97" ),
-MpIeee( "3.8791952037716162900707994047e+100" ),
 MpIeee( "3.6600317773156342245401829308e+103" ),
-MpIeee( "3.5978775704117283875784869570e+106" )    /* eta(-99)  */
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int  gsl_sf_hzeta_e(const MpIeee s, const MpIeee q, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s <= 1.0 || q <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    const MpIeee max_bits=  54.0;
    const MpIeee ln_term0=  -s * log(q);  

    if(ln_term0 < GSL_LOG_DBL_MIN + 1.0) {
      UNDERFLOW_ERROR(result);
    }
    else if(ln_term0 > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR (result);
    }
    else if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25)) {
      result->val = pow(q, -s);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(s > 0.5*max_bits && q < 1.0) {
      const MpIeee p1=  pow(q, -s);
      const MpIeee p2=  pow(q/(1.0+q), s);
      const MpIeee p3=  pow(q/(2.0+q), s);
      result->val = p1 * (1.0 + p2 + p3);
      result->err = GSL_DBL_EPSILON * (0.5*s + 2.0) * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Euler-Maclaurin summation formula 
       * [Moshier, p. 400, with several typo corrections]
       */
      const int jmax = 12;
      const int kmax = 10;
      int  j;int   k;
      const MpIeee pmax=  pow(kmax + q, -s);
      MpIeee scp=  s;
      MpIeee pcp=  pmax / (kmax + q);
      MpIeee ans=  pmax*((kmax+q)/(s-MpIeee( "1.0" )) + MpIeee( "0.5" ));

      for(k=0; k<kmax; k++) {
        ans += pow(k + q, -s);
      }

      for(j=0; j<=jmax; j++) {
        MpIeee delta=  hzeta_c[j+1] * scp * pcp;
        ans += delta;
        if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
        scp *= (s+MpIeee( "2" )*j+MpIeee( "1" ))*(s+MpIeee( "2" )*j+MpIeee( "2" ));
        pcp /= (kmax + q)*(kmax + q);
      }

      result->val = ans;
      result->err = 2.0 * (jmax + 1.0) * GSL_DBL_EPSILON * fabs(ans);
      return GSL_SUCCESS;
    }
  }
}


int  gsl_sf_zeta_e(const MpIeee s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s == 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(s >= 0.0) {
    return riemann_zeta_sgt0(s, result);
  }
  else {
    /* reflection formula, [Abramowitz+Stegun, 23.2.5] */

    gsl_sf_result zeta_one_minus_s;
    const int stat_zoms = riemann_zeta1ms_slt0(s, &zeta_one_minus_s);
    const MpIeee sin_term=  sin(0.5*M_PI*s)/M_PI;

    if(sin_term == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(s > -170) {
      /* We have to be careful about losing digits
       * in calculating pow(2 Pi, s). The gamma
       * function is fine because we were careful
       * with that implementation.
       * We keep an array of (2 Pi)^(10 n).
       */
      const MpIeee twopi_pow[18] =  { 1.0,
                                     9.589560061550901348e+007,
                                     9.195966217409212684e+015,
                                     8.818527036583869903e+023,
                                     8.456579467173150313e+031,
                                     8.109487671573504384e+039,
                                     7.776641909496069036e+047,
                                     7.457457466828644277e+055,
                                     7.151373628461452286e+063,
                                     6.857852693272229709e+071,
                                     6.576379029540265771e+079,
                                     6.306458169130020789e+087,
                                     6.047615938853066678e+095,
                                     5.799397627482402614e+103,
                                     5.561367186955830005e+111,
                                     5.333106466365131227e+119,
                                     5.114214477385391780e+127,
                                     4.904306689854036836e+135
                                    };
      const int n = floor((-s)/10.0);
      const MpIeee fs=  s + 10.0*n;
      const MpIeee p=  pow(2.0*M_PI, fs) / twopi_pow[n];

      gsl_sf_result g;
      const int stat_g = gsl_sf_gamma_e(1.0-s, &g);
      result->val  = p * g.val * sin_term * zeta_one_minus_s.val;
      result->err  = fabs(p * g.val * sin_term) * zeta_one_minus_s.err;
      result->err += fabs(p * sin_term * zeta_one_minus_s.val) * g.err;
      result->err += GSL_DBL_EPSILON * (fabs(s)+2.0) * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_g, stat_zoms);
    }
    else {
      /* The actual zeta function may or may not
       * overflow here. But we have no easy way
       * to calculate it when the prefactor(s)
       * overflow. Trying to use log's and exp
       * is no good because we loose a couple
       * digits to the exp error amplification.
       * When we gather a little more patience,
       * we can implement something here. Until
       * then just give up.
       */
      OVERFLOW_ERROR(result);
    }
  }
}


int  gsl_sf_zeta_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 0) {
    if(!GSL_IS_ODD(n)) {
      result->val = 0.0; /* exactly zero at even negative integers */
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ZETA_NEG_TABLE_NMAX) {
      result->val = zeta_neg_int_table[-(n+1)/2];
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      return gsl_sf_zeta_e((double)n, result);
    }
  }
  else if(n == 1){
    DOMAIN_ERROR(result);
  }
  else if(n <= ZETA_POS_TABLE_NMAX){
    result->val = 1.0 + zetam1_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}


int  gsl_sf_zetam1_e(const MpIeee s, gsl_sf_result * result)
{
  if(s <= 5.0)
  {
    int  stat=  gsl_sf_zeta_e(s, result);
    result->val = result->val - 1.0;
    return stat;
  }
  else if(s < 15.0)
  {
    return riemann_zeta_minus_1_intermediate_s(s, result);
  }
  else
  {
    return riemann_zeta_minus1_large_s(s, result);
  }
}


int  gsl_sf_zetam1_int_e(const int n, gsl_sf_result * result)
{
  if(n < 0) {
    if(!GSL_IS_ODD(n)) {
      result->val = 0.0; /* exactly zero at even negative integers */
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ZETA_NEG_TABLE_NMAX) {
      result->val = zeta_neg_int_table[-(n+1)/2] - 1.0;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      return gsl_sf_zeta_e((double)n, result);
    }
  }
  else if(n == 1){
    DOMAIN_ERROR(result);
  }
  else if(n <= ZETA_POS_TABLE_NMAX){
    result->val = zetam1_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return gsl_sf_zetam1_e(n, result);
  }
}


int  gsl_sf_eta_int_e(int  n, gsl_sf_result * result)
{
  if(n > ETA_POS_TABLE_NMAX) {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(n >= 0) {
    result->val = eta_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* n < 0 */

    if(!GSL_IS_ODD(n)) {
      /* exactly zero at even negative integers */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ETA_NEG_TABLE_NMAX) {
      result->val = eta_neg_int_table[-(n+1)/2];
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result z;
      gsl_sf_result p;
      int  stat_z=  gsl_sf_zeta_int_e(n, &z);
      int  stat_p=  gsl_sf_exp_e((1.0-n)*M_LN2, &p);
      int  stat_m=  gsl_sf_multiply_e(-p.val, z.val, result);
      result->err  = fabs(p.err * (M_LN2*(1.0-n)) * z.val) + z.err * fabs(p.val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_3(stat_m, stat_p, stat_z);
    }
  }
}


int  gsl_sf_eta_e(const MpIeee s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s > 100.0) {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(fabs(s-1.0) < 10.0*GSL_ROOT5_DBL_EPSILON) {
    MpIeee del=  s-MpIeee( "1.0" );
    MpIeee c0=  M_LN2;
    MpIeee c1=  M_LN2 * (M_EULER - MpIeee( "0.5" )*M_LN2);
    MpIeee c2=  -MpIeee( "0.0326862962794492996" );
    MpIeee c3=   MpIeee( "0.0015689917054155150" );
    MpIeee c4=   MpIeee( "0.00074987242112047532" );
    result->val = c0 + del * (c1 + del * (c2 + del * (c3 + del * c4)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result z;
    gsl_sf_result p;
    int  stat_z=  gsl_sf_zeta_e(s, &z);
    int  stat_p=  gsl_sf_exp_e((1.0-s)*M_LN2, &p);
    int  stat_m=  gsl_sf_multiply_e(1.0-p.val, z.val, result);
    result->err  = fabs(p.err * (M_LN2*(1.0-s)) * z.val) + z.err * fabs(p.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_3(stat_m, stat_p, stat_z);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

MpIeee gsl_sf_zeta(const MpIeee s)
{
  EVAL_RESULT(gsl_sf_zeta_e(s, &result));
}

MpIeee gsl_sf_hzeta(const MpIeee s, const MpIeee a)
{
  EVAL_RESULT(gsl_sf_hzeta_e(s, a, &result));
}

MpIeee gsl_sf_zeta_int(const int s)
{
  EVAL_RESULT(gsl_sf_zeta_int_e(s, &result));
}

MpIeee gsl_sf_zetam1(const MpIeee s)
{
  EVAL_RESULT(gsl_sf_zetam1_e(s, &result));
}

MpIeee gsl_sf_zetam1_int(const int s)
{
  EVAL_RESULT(gsl_sf_zetam1_int_e(s, &result));
}

MpIeee gsl_sf_eta_int(const int s)
{
  EVAL_RESULT(gsl_sf_eta_int_e(s, &result));
}

MpIeee gsl_sf_eta(const MpIeee s)
{
  EVAL_RESULT(gsl_sf_eta_e(s, &result));
}
