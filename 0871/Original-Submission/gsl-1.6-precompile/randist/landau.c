#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* randist/landau.c
 *
 * Copyright (C) 2001, 2004 David Morrison
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

/* Adapted from the CERN library routines DENLAN, RANLAN, and DISLAN
 * as described in http://consult.cern.ch/shortwrups/g110/top.html.
 * Original author: K.S. K\"olbig.
 *
 * The distribution is given by the complex path integral,
 *
 *  p(x) = (1/(2 pi i)) \int_{c-i\inf}^{c+i\inf} ds exp(s log(s) + x s) 
 *
 * which can be converted into a real integral over [0,+\inf]
 *
 *  p(x) = (1/pi) \int_0^\inf dt \exp(-t log(t) - x t) sin(pi t)
 *
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

MpIeee gsl_ran_landau_pdf(const MpIeee x)
{
  static MpIeee P1[5] = 
    {
      MpIeee( "0.4259894875E0" ), -MpIeee( "0.1249762550E0" ), MpIeee( "0.3984243700E-1" ),
      -MpIeee( "0.6298287635E-2" ), MpIeee( "0.1511162253E-2" )
    };
  static MpIeee P2[5] = 
    {
      MpIeee( "0.1788541609E0" ), MpIeee( "0.1173957403E0" ), MpIeee( "0.1488850518E-1" ),
      -MpIeee( "0.1394989411E-2" ), MpIeee( "0.1283617211E-3" )
    };
  static MpIeee P3[5] = 
    {
      MpIeee( "0.1788544503E0" ), MpIeee( "0.9359161662E-1" ), MpIeee( "0.6325387654E-2" ),
      MpIeee( "0.6611667319E-4" ), -MpIeee( "0.2031049101E-5" )
    };
  static MpIeee P4[5] = 
    {
      MpIeee( "0.9874054407E0" ), MpIeee( "0.1186723273E3" ), MpIeee( "0.8492794360E3" ),
      -MpIeee( "0.7437792444E3" ), MpIeee( "0.4270262186E3" )
    };
  static MpIeee P5[5] = 
    {
      MpIeee( "0.1003675074E1" ), MpIeee( "0.1675702434E3" ), MpIeee( "0.4789711289E4" ),
      MpIeee( "0.2121786767E5" ), -MpIeee( "0.2232494910E5" )
    };
  static MpIeee P6[5] = 
    {
      MpIeee( "0.1000827619E1" ), MpIeee( "0.6649143136E3" ), MpIeee( "0.6297292665E5" ),
      MpIeee( "0.4755546998E6" ), -MpIeee( "0.5743609109E7" )
    };

  static MpIeee Q1[5] = 
    {
      MpIeee( "1.0" ), -MpIeee( "0.3388260629E0" ), MpIeee( "0.9594393323E-1" ),
      -MpIeee( "0.1608042283E-1" ), MpIeee( "0.3778942063E-2" )
    };
  static MpIeee Q2[5] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.7428795082E0" ), MpIeee( "0.3153932961E0" ),
      MpIeee( "0.6694219548E-1" ), MpIeee( "0.8790609714E-2" )
    };
  static MpIeee Q3[5] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.6097809921E0" ), MpIeee( "0.2560616665E0" ),
      MpIeee( "0.4746722384E-1" ), MpIeee( "0.6957301675E-2" )
    };
  static MpIeee Q4[5] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.1068615961E3" ), MpIeee( "0.3376496214E3" ),
      MpIeee( "0.2016712389E4" ), MpIeee( "0.1597063511E4" )
    };
  static MpIeee Q5[5] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.1569424537E3" ), MpIeee( "0.3745310488E4" ),
      MpIeee( "0.9834698876E4" ), MpIeee( "0.6692428357E5" )
    };
  static MpIeee Q6[5] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.6514101098E3" ), MpIeee( "0.5697473333E5" ),
      MpIeee( "0.1659174725E6" ), -MpIeee( "0.2815759939E7" )
    };

  static MpIeee A1[3] = 
    {
      MpIeee( "0.4166666667E-1" ), -MpIeee( "0.1996527778E-1" ), MpIeee( "0.2709538966E-1" )
    };
  static MpIeee A2[2] = 
    {
      -MpIeee( "0.1845568670E1" ), -MpIeee( "0.4284640743E1" )
    };

  MpIeee U;MpIeee  V;MpIeee  DENLAN;

  V = x;
  if (V < -MpIeee( "5.5" ))
    {
      U = exp(V + MpIeee( "1.0" ));
      DENLAN = MpIeee( "0.3989422803" ) * (exp( -MpIeee( "1" ) / U) / sqrt(U)) *
        (MpIeee( "1" ) + (A1[0] + (A1[1] + A1[2] * U) * U) * U);
    }
  else if (V < -MpIeee( "1" ))
    {
      U = exp( -V - MpIeee( "1" ));
      DENLAN = exp( -U) * sqrt(U) *
        (P1[0] + (P1[1] + (P1[2] + (P1[3] + P1[4] * V) * V) * V) * V) /
        (Q1[0] + (Q1[1] + (Q1[2] + (Q1[3] + Q1[4] * V) * V) * V) * V);
    }
  else if (V < MpIeee( "1" ))
    {
      DENLAN = (P2[0] + (P2[1] + (P2[2] + (P2[3] + P2[4] * V) * V) * V) * V) /
        (Q2[0] + (Q2[1] + (Q2[2] + (Q2[3] + Q2[4] * V) * V) * V) * V);
    }
  else if (V < MpIeee( "5" ))
    {
      DENLAN = (P3[0] + (P3[1] + (P3[2] + (P3[3] + P3[4] * V) * V) * V) * V) /
        (Q3[0] + (Q3[1] + (Q3[2] + (Q3[3] + Q3[4] * V) * V) * V) * V);
    }
  else if (V < MpIeee( "12" ))
    {
      U = MpIeee( "1" ) / V;
      DENLAN = U * U *
        (P4[0] + (P4[1] + (P4[2] + (P4[3] + P4[4] * U) * U) * U) * U) /
        (Q4[0] + (Q4[1] + (Q4[2] + (Q4[3] + Q4[4] * U) * U) * U) * U);
    }
  else if (V < MpIeee( "50" ))
    {
      U = MpIeee( "1" ) / V;
      DENLAN = U * U *
        (P5[0] + (P5[1] + (P5[2] + (P5[3] + P5[4] * U) * U) * U) * U) /
        (Q5[0] + (Q5[1] + (Q5[2] + (Q5[3] + Q5[4] * U) * U) * U) * U);
    }
  else if (V < MpIeee( "300" ))
    {
      U = MpIeee( "1" ) / V;
      DENLAN = U * U *
        (P6[0] + (P6[1] + (P6[2] + (P6[3] + P6[4] * U) * U) * U) * U) /
        (Q6[0] + (Q6[1] + (Q6[2] + (Q6[3] + Q6[4] * U) * U) * U) * U);
    }
  else
    {
      U = MpIeee( "1" ) / (V - V * log(V) / (V + MpIeee( "1" )));
      DENLAN = U * U * (MpIeee( "1" ) + (A2[0] + A2[1] * U) * U);
    }

  return DENLAN;
}

#if 0 /* Not needed yet */
/* This function is a translation from the original Fortran of the
 * CERN library routine DISLAN, the integral from -inf to x of the
 * Landau p.d.f.
 */
static
MpIeee gsl_ran_landau_dislan(const MpIeee x)
{
  static MpIeee P1[5] = 
    {
      MpIeee( "0.2514091491E0" ), -MpIeee( "0.6250580444E-1" ),
      MpIeee( "0.1458381230E-1" ), -MpIeee( "0.2108817737E-2" ),
      MpIeee( "0.7411247290E-3" )
    };

  static MpIeee P2[4] = 
    {
      MpIeee( "0.2868328584E0" ), MpIeee( "0.3564363231E0" ),
      MpIeee( "0.1523518695E0" ), MpIeee( "0.2251304883E-1" )
    };

  static MpIeee P3[4] = 
    {
      MpIeee( "0.2868329066E0" ), MpIeee( "0.3003828436E0" ),
      MpIeee( "0.9950951941E-1" ), MpIeee( "0.8733827185E-2" )
    };

  static MpIeee P4[4] = 
    {
      MpIeee( "0.1000351630E1" ), MpIeee( "0.4503592498E1" ),
      MpIeee( "0.1085883880E2" ), MpIeee( "0.7536052269E1" )
    };

  static MpIeee P5[4] = 
    {
      MpIeee( "0.1000006517E1" ), MpIeee( "0.4909414111E2" ),
      MpIeee( "0.8505544753E2" ), MpIeee( "0.1532153455E3" )
    };

  static MpIeee P6[4] = 
    {
      MpIeee( "0.1000000983E1" ), MpIeee( "0.1329868456E3" ),
      MpIeee( "0.9162149244E3" ), -MpIeee( "0.9605054274E3" )
    };

  static MpIeee Q1[5] = 
    {
      MpIeee( "1.0" ), -MpIeee( "0.5571175625E-2" ),
      MpIeee( "0.6225310236E-1" ), -MpIeee( "0.3137378427E-2" ),
      MpIeee( "0.1931496439E-2" )
    };

  static MpIeee Q2[4] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.6191136137E0" ),
      MpIeee( "0.1720721448E0" ), MpIeee( "0.2278594771E-1" )
    };

  static MpIeee Q3[4] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.4237190502E0" ),
      MpIeee( "0.1095631512E0" ), MpIeee( "0.8693851567E-2" )
    };

  static MpIeee Q4[4] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.5539969678E1" ),
      MpIeee( "0.1933581111E2" ), MpIeee( "0.2721321508E2" )
    };

  static MpIeee Q5[4] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.5009928881E2" ),
      MpIeee( "0.1399819104E3" ), MpIeee( "0.4200002909E3" )
    };

  static MpIeee Q6[4] = 
    {
      MpIeee( "1.0" ), MpIeee( "0.1339887843E3" ),
      MpIeee( "0.1055990413E4" ), MpIeee( "0.5532224619E3" )
    };

  static MpIeee A1[3] = 
    {
      -MpIeee( "0.4583333333E0" ), MpIeee( "0.6675347222E0" ), -MpIeee( "0.1641741416E1" )
    };

  static MpIeee A2[3] = 
    {
      MpIeee( "1.0" ), -MpIeee( "0.4227843351E0" ), -MpIeee( "0.2043403138E1" )
    };

  MpIeee U;MpIeee  V;MpIeee  DISLAN;

  V = x;
  if (V < -MpIeee( "5.5" ))
    {
      U = exp(V + MpIeee( "1" ));
      DISLAN = MpIeee( "0.3989422803" ) * exp( -MpIeee( "1" ) / U) * sqrt(U) *
               (MpIeee( "1" ) + (A1[0] + (A1[1] + A1[2] * U) * U) * U);
    }
  else if (V < -MpIeee( "1" ))
    {
      U = exp( -V - MpIeee( "1" ));
      DISLAN = (exp( -U) / sqrt(U)) *
               (P1[0] + (P1[1] + (P1[2] + (P1[3] + P1[4] * V) * V) * V) * V) /
               (Q1[0] + (Q1[1] + (Q1[2] + (Q1[3] + Q1[4] * V) * V) * V) * V);
    }
  else if (V < MpIeee( "1" ))
    {
      DISLAN = (P2[0] + (P2[1] + (P2[2] + P2[3] * V) * V) * V) /
               (Q2[0] + (Q2[1] + (Q2[2] + Q2[3] * V) * V) * V);
    }
  else if (V < MpIeee( "4" ))
    {
      DISLAN = (P3[0] + (P3[1] + (P3[2] + P3[3] * V) * V) * V) /
               (Q3[0] + (Q3[1] + (Q3[2] + Q3[3] * V) * V) * V);
    }
  else if (V < MpIeee( "12" ))
    {
      U = MpIeee( "1" ) / V;
      DISLAN = (P4[0] + (P4[1] + (P4[2] + P4[3] * U) * U) * U) /
               (Q4[0] + (Q4[1] + (Q4[2] + Q4[3] * U) * U) * U);
    }
  else if (V < MpIeee( "50" ))
    {
      U = MpIeee( "1" ) / V;
      DISLAN = (P5[0] + (P5[1] + (P5[2] + P5[3] * U) * U) * U) /
               (Q5[0] + (Q5[1] + (Q5[2] + Q5[3] * U) * U) * U);
    }
  else if (V < MpIeee( "300" ))
    {
      U = MpIeee( "1" ) / V;
      DISLAN = (P6[0] + (P6[1] + (P6[2] + P6[3] * U) * U) * U) /
               (Q6[0] + (Q6[1] + (Q6[2] + Q6[3] * U) * U) * U);
    }
  else
    {
      U = MpIeee( "1" ) / (V - V * log(V) / (V + MpIeee( "1" )));
      DISLAN = MpIeee( "1" ) - (A2[0] + (A2[1] + A2[2] * U) * U) * U;
    }

  return DISLAN;
}
#endif

MpIeee gsl_ran_landau(const gsl_rng * r)
{
  static MpIeee F[983] = 
    {
      MpIeee( "0.0000000" ),   /* Add empty element [0] to account for difference 
                      between C and Fortran convention for lower bound. */
      MpIeee( "00.000000" ), MpIeee( "00.000000" ), MpIeee( "00.000000" ), MpIeee( "00.000000" ), MpIeee( "00.000000" ),
      -MpIeee( "2.244733" ), -MpIeee( "2.204365" ), -MpIeee( "2.168163" ), -MpIeee( "2.135219" ), -MpIeee( "2.104898" ),
      -MpIeee( "2.076740" ), -MpIeee( "2.050397" ), -MpIeee( "2.025605" ), -MpIeee( "2.002150" ), -MpIeee( "1.979866" ),
      -MpIeee( "1.958612" ), -MpIeee( "1.938275" ), -MpIeee( "1.918760" ), -MpIeee( "1.899984" ), -MpIeee( "1.881879" ),
      -MpIeee( "1.864385" ), -MpIeee( "1.847451" ), -MpIeee( "1.831030" ), -MpIeee( "1.815083" ), -MpIeee( "1.799574" ),
      -MpIeee( "1.784473" ), -MpIeee( "1.769751" ), -MpIeee( "1.755383" ), -MpIeee( "1.741346" ), -MpIeee( "1.727620" ),
      -MpIeee( "1.714187" ), -MpIeee( "1.701029" ), -MpIeee( "1.688130" ), -MpIeee( "1.675477" ), -MpIeee( "1.663057" ),
      -MpIeee( "1.650858" ), -MpIeee( "1.638868" ), -MpIeee( "1.627078" ), -MpIeee( "1.615477" ), -MpIeee( "1.604058" ),
      -MpIeee( "1.592811" ), -MpIeee( "1.581729" ), -MpIeee( "1.570806" ), -MpIeee( "1.560034" ), -MpIeee( "1.549407" ),
      -MpIeee( "1.538919" ), -MpIeee( "1.528565" ), -MpIeee( "1.518339" ), -MpIeee( "1.508237" ), -MpIeee( "1.498254" ),
      -MpIeee( "1.488386" ), -MpIeee( "1.478628" ), -MpIeee( "1.468976" ), -MpIeee( "1.459428" ), -MpIeee( "1.449979" ),
      -MpIeee( "1.440626" ), -MpIeee( "1.431365" ), -MpIeee( "1.422195" ), -MpIeee( "1.413111" ), -MpIeee( "1.404112" ),
      -MpIeee( "1.395194" ), -MpIeee( "1.386356" ), -MpIeee( "1.377594" ), -MpIeee( "1.368906" ), -MpIeee( "1.360291" ),
      -MpIeee( "1.351746" ), -MpIeee( "1.343269" ), -MpIeee( "1.334859" ), -MpIeee( "1.326512" ), -MpIeee( "1.318229" ),
      -MpIeee( "1.310006" ), -MpIeee( "1.301843" ), -MpIeee( "1.293737" ), -MpIeee( "1.285688" ), -MpIeee( "1.277693" ),
      -MpIeee( "1.269752" ), -MpIeee( "1.261863" ), -MpIeee( "1.254024" ), -MpIeee( "1.246235" ), -MpIeee( "1.238494" ),
      -MpIeee( "1.230800" ), -MpIeee( "1.223153" ), -MpIeee( "1.215550" ), -MpIeee( "1.207990" ), -MpIeee( "1.200474" ),
      -MpIeee( "1.192999" ), -MpIeee( "1.185566" ), -MpIeee( "1.178172" ), -MpIeee( "1.170817" ), -MpIeee( "1.163500" ),
      -MpIeee( "1.156220" ), -MpIeee( "1.148977" ), -MpIeee( "1.141770" ), -MpIeee( "1.134598" ), -MpIeee( "1.127459" ),
      -MpIeee( "1.120354" ), -MpIeee( "1.113282" ), -MpIeee( "1.106242" ), -MpIeee( "1.099233" ), -MpIeee( "1.092255" ),
      -MpIeee( "1.085306" ), -MpIeee( "1.078388" ), -MpIeee( "1.071498" ), -MpIeee( "1.064636" ), -MpIeee( "1.057802" ),
      -MpIeee( "1.050996" ), -MpIeee( "1.044215" ), -MpIeee( "1.037461" ), -MpIeee( "1.030733" ), -MpIeee( "1.024029" ),
      -MpIeee( "1.017350" ), -MpIeee( "1.010695" ), -MpIeee( "1.004064" ), -MpIeee( "0.997456" ), -MpIeee( "0.990871" ),
      -MpIeee( "0.984308" ), -MpIeee( "0.977767" ), -MpIeee( "0.971247" ), -MpIeee( "0.964749" ), -MpIeee( "0.958271" ),
      -MpIeee( "0.951813" ), -MpIeee( "0.945375" ), -MpIeee( "0.938957" ), -MpIeee( "0.932558" ), -MpIeee( "0.926178" ),
      -MpIeee( "0.919816" ), -MpIeee( "0.913472" ), -MpIeee( "0.907146" ), -MpIeee( "0.900838" ), -MpIeee( "0.894547" ),
      -MpIeee( "0.888272" ), -MpIeee( "0.882014" ), -MpIeee( "0.875773" ), -MpIeee( "0.869547" ), -MpIeee( "0.863337" ),
      -MpIeee( "0.857142" ), -MpIeee( "0.850963" ), -MpIeee( "0.844798" ), -MpIeee( "0.838648" ), -MpIeee( "0.832512" ),
      -MpIeee( "0.826390" ), -MpIeee( "0.820282" ), -MpIeee( "0.814187" ), -MpIeee( "0.808106" ), -MpIeee( "0.802038" ),
      -MpIeee( "0.795982" ), -MpIeee( "0.789940" ), -MpIeee( "0.783909" ), -MpIeee( "0.777891" ), -MpIeee( "0.771884" ),
      -MpIeee( "0.765889" ), -MpIeee( "0.759906" ), -MpIeee( "0.753934" ), -MpIeee( "0.747973" ), -MpIeee( "0.742023" ),
      -MpIeee( "0.736084" ), -MpIeee( "0.730155" ), -MpIeee( "0.724237" ), -MpIeee( "0.718328" ), -MpIeee( "0.712429" ),
      -MpIeee( "0.706541" ), -MpIeee( "0.700661" ), -MpIeee( "0.694791" ), -MpIeee( "0.688931" ), -MpIeee( "0.683079" ),
      -MpIeee( "0.677236" ), -MpIeee( "0.671402" ), -MpIeee( "0.665576" ), -MpIeee( "0.659759" ), -MpIeee( "0.653950" ),
      -MpIeee( "0.648149" ), -MpIeee( "0.642356" ), -MpIeee( "0.636570" ), -MpIeee( "0.630793" ), -MpIeee( "0.625022" ),
      -MpIeee( "0.619259" ), -MpIeee( "0.613503" ), -MpIeee( "0.607754" ), -MpIeee( "0.602012" ), -MpIeee( "0.596276" ),
      -MpIeee( "0.590548" ), -MpIeee( "0.584825" ), -MpIeee( "0.579109" ), -MpIeee( "0.573399" ), -MpIeee( "0.567695" ),
      -MpIeee( "0.561997" ), -MpIeee( "0.556305" ), -MpIeee( "0.550618" ), -MpIeee( "0.544937" ), -MpIeee( "0.539262" ),
      -MpIeee( "0.533592" ), -MpIeee( "0.527926" ), -MpIeee( "0.522266" ), -MpIeee( "0.516611" ), -MpIeee( "0.510961" ),
      -MpIeee( "0.505315" ), -MpIeee( "0.499674" ), -MpIeee( "0.494037" ), -MpIeee( "0.488405" ), -MpIeee( "0.482777" ),
      -MpIeee( "0.477153" ), -MpIeee( "0.471533" ), -MpIeee( "0.465917" ), -MpIeee( "0.460305" ), -MpIeee( "0.454697" ),
      -MpIeee( "0.449092" ), -MpIeee( "0.443491" ), -MpIeee( "0.437893" ), -MpIeee( "0.432299" ), -MpIeee( "0.426707" ),
      -MpIeee( "0.421119" ), -MpIeee( "0.415534" ), -MpIeee( "0.409951" ), -MpIeee( "0.404372" ), -MpIeee( "0.398795" ),
      -MpIeee( "0.393221" ), -MpIeee( "0.387649" ), -MpIeee( "0.382080" ), -MpIeee( "0.376513" ), -MpIeee( "0.370949" ),
      -MpIeee( "0.365387" ), -MpIeee( "0.359826" ), -MpIeee( "0.354268" ), -MpIeee( "0.348712" ), -MpIeee( "0.343157" ),
      -MpIeee( "0.337604" ), -MpIeee( "0.332053" ), -MpIeee( "0.326503" ), -MpIeee( "0.320955" ), -MpIeee( "0.315408" ),
      -MpIeee( "0.309863" ), -MpIeee( "0.304318" ), -MpIeee( "0.298775" ), -MpIeee( "0.293233" ), -MpIeee( "0.287692" ),
      -MpIeee( "0.282152" ), -MpIeee( "0.276613" ), -MpIeee( "0.271074" ), -MpIeee( "0.265536" ), -MpIeee( "0.259999" ),
      -MpIeee( "0.254462" ), -MpIeee( "0.248926" ), -MpIeee( "0.243389" ), -MpIeee( "0.237854" ), -MpIeee( "0.232318" ),
      -MpIeee( "0.226783" ), -MpIeee( "0.221247" ), -MpIeee( "0.215712" ), -MpIeee( "0.210176" ), -MpIeee( "0.204641" ),
      -MpIeee( "0.199105" ), -MpIeee( "0.193568" ), -MpIeee( "0.188032" ), -MpIeee( "0.182495" ), -MpIeee( "0.176957" ),
      -MpIeee( "0.171419" ), -MpIeee( "0.165880" ), -MpIeee( "0.160341" ), -MpIeee( "0.154800" ), -MpIeee( "0.149259" ),
      -MpIeee( "0.143717" ), -MpIeee( "0.138173" ), -MpIeee( "0.132629" ), -MpIeee( "0.127083" ), -MpIeee( "0.121537" ),
      -MpIeee( "0.115989" ), -MpIeee( "0.110439" ), -MpIeee( "0.104889" ), -MpIeee( "0.099336" ), -MpIeee( "0.093782" ),
      -MpIeee( "0.088227" ), -MpIeee( "0.082670" ), -MpIeee( "0.077111" ), -MpIeee( "0.071550" ), -MpIeee( "0.065987" ),
      -MpIeee( "0.060423" ), -MpIeee( "0.054856" ), -MpIeee( "0.049288" ), -MpIeee( "0.043717" ), -MpIeee( "0.038144" ),
      -MpIeee( "0.032569" ), -MpIeee( "0.026991" ), -MpIeee( "0.021411" ), -MpIeee( "0.015828" ), -MpIeee( "0.010243" ),
      -MpIeee( "0.004656" ), MpIeee( "00.000934" ), MpIeee( "00.006527" ), MpIeee( "00.012123" ), MpIeee( "00.017722" ),
      MpIeee( "00.023323" ), MpIeee( "00.028928" ), MpIeee( "00.034535" ), MpIeee( "00.040146" ), MpIeee( "00.045759" ),
      MpIeee( "00.051376" ), MpIeee( "00.056997" ), MpIeee( "00.062620" ), MpIeee( "00.068247" ), MpIeee( "00.073877" ),
      MpIeee( "00.079511" ), MpIeee( "00.085149" ), MpIeee( "00.090790" ), MpIeee( "00.096435" ), MpIeee( "00.102083" ),
      MpIeee( "00.107736" ), MpIeee( "00.113392" ), MpIeee( "00.119052" ), MpIeee( "00.124716" ), MpIeee( "00.130385" ),
      MpIeee( "00.136057" ), MpIeee( "00.141734" ), MpIeee( "00.147414" ), MpIeee( "00.153100" ), MpIeee( "00.158789" ),
      MpIeee( "00.164483" ), MpIeee( "00.170181" ), MpIeee( "00.175884" ), MpIeee( "00.181592" ), MpIeee( "00.187304" ),
      MpIeee( "00.193021" ), MpIeee( "00.198743" ), MpIeee( "00.204469" ), MpIeee( "00.210201" ), MpIeee( "00.215937" ),
      MpIeee( "00.221678" ), MpIeee( "00.227425" ), MpIeee( "00.233177" ), MpIeee( "00.238933" ), MpIeee( "00.244696" ),
      MpIeee( "00.250463" ), MpIeee( "00.256236" ), MpIeee( "00.262014" ), MpIeee( "00.267798" ), MpIeee( "00.273587" ),
      MpIeee( "00.279382" ), MpIeee( "00.285183" ), MpIeee( "00.290989" ), MpIeee( "00.296801" ), MpIeee( "00.302619" ),
      MpIeee( "00.308443" ), MpIeee( "00.314273" ), MpIeee( "00.320109" ), MpIeee( "00.325951" ), MpIeee( "00.331799" ),
      MpIeee( "00.337654" ), MpIeee( "00.343515" ), MpIeee( "00.349382" ), MpIeee( "00.355255" ), MpIeee( "00.361135" ),
      MpIeee( "00.367022" ), MpIeee( "00.372915" ), MpIeee( "00.378815" ), MpIeee( "00.384721" ), MpIeee( "00.390634" ),
      MpIeee( "00.396554" ), MpIeee( "00.402481" ), MpIeee( "00.408415" ), MpIeee( "00.414356" ), MpIeee( "00.420304" ),
      MpIeee( "00.426260" ), MpIeee( "00.432222" ), MpIeee( "00.438192" ), MpIeee( "00.444169" ), MpIeee( "00.450153" ),
      MpIeee( "00.456145" ), MpIeee( "00.462144" ), MpIeee( "00.468151" ), MpIeee( "00.474166" ), MpIeee( "00.480188" ),
      MpIeee( "00.486218" ), MpIeee( "00.492256" ), MpIeee( "00.498302" ), MpIeee( "00.504356" ), MpIeee( "00.510418" ),
      MpIeee( "00.516488" ), MpIeee( "00.522566" ), MpIeee( "00.528653" ), MpIeee( "00.534747" ), MpIeee( "00.540850" ),
      MpIeee( "00.546962" ), MpIeee( "00.553082" ), MpIeee( "00.559210" ), MpIeee( "00.565347" ), MpIeee( "00.571493" ),
      MpIeee( "00.577648" ), MpIeee( "00.583811" ), MpIeee( "00.589983" ), MpIeee( "00.596164" ), MpIeee( "00.602355" ),
      MpIeee( "00.608554" ), MpIeee( "00.614762" ), MpIeee( "00.620980" ), MpIeee( "00.627207" ), MpIeee( "00.633444" ),
      MpIeee( "00.639689" ), MpIeee( "00.645945" ), MpIeee( "00.652210" ), MpIeee( "00.658484" ), MpIeee( "00.664768" ),
      MpIeee( "00.671062" ), MpIeee( "00.677366" ), MpIeee( "00.683680" ), MpIeee( "00.690004" ), MpIeee( "00.696338" ),
      MpIeee( "00.702682" ), MpIeee( "00.709036" ), MpIeee( "00.715400" ), MpIeee( "00.721775" ), MpIeee( "00.728160" ),
      MpIeee( "00.734556" ), MpIeee( "00.740963" ), MpIeee( "00.747379" ), MpIeee( "00.753807" ), MpIeee( "00.760246" ),
      MpIeee( "00.766695" ), MpIeee( "00.773155" ), MpIeee( "00.779627" ), MpIeee( "00.786109" ), MpIeee( "00.792603" ),
      MpIeee( "00.799107" ), MpIeee( "00.805624" ), MpIeee( "00.812151" ), MpIeee( "00.818690" ), MpIeee( "00.825241" ),
      MpIeee( "00.831803" ), MpIeee( "00.838377" ), MpIeee( "00.844962" ), MpIeee( "00.851560" ), MpIeee( "00.858170" ),
      MpIeee( "00.864791" ), MpIeee( "00.871425" ), MpIeee( "00.878071" ), MpIeee( "00.884729" ), MpIeee( "00.891399" ),
      MpIeee( "00.898082" ), MpIeee( "00.904778" ), MpIeee( "00.911486" ), MpIeee( "00.918206" ), MpIeee( "00.924940" ),
      MpIeee( "00.931686" ), MpIeee( "00.938446" ), MpIeee( "00.945218" ), MpIeee( "00.952003" ), MpIeee( "00.958802" ),
      MpIeee( "00.965614" ), MpIeee( "00.972439" ), MpIeee( "00.979278" ), MpIeee( "00.986130" ), MpIeee( "00.992996" ),
      MpIeee( "00.999875" ), MpIeee( "01.006769" ), MpIeee( "01.013676" ), MpIeee( "01.020597" ), MpIeee( "01.027533" ),
      MpIeee( "01.034482" ), MpIeee( "01.041446" ), MpIeee( "01.048424" ), MpIeee( "01.055417" ), MpIeee( "01.062424" ),
      MpIeee( "01.069446" ), MpIeee( "01.076482" ), MpIeee( "01.083534" ), MpIeee( "01.090600" ), MpIeee( "01.097681" ),
      MpIeee( "01.104778" ), MpIeee( "01.111889" ), MpIeee( "01.119016" ), MpIeee( "01.126159" ), MpIeee( "01.133316" ),
      MpIeee( "01.140490" ), MpIeee( "01.147679" ), MpIeee( "01.154884" ), MpIeee( "01.162105" ), MpIeee( "01.169342" ),
      MpIeee( "01.176595" ), MpIeee( "01.183864" ), MpIeee( "01.191149" ), MpIeee( "01.198451" ), MpIeee( "01.205770" ),
      MpIeee( "01.213105" ), MpIeee( "01.220457" ), MpIeee( "01.227826" ), MpIeee( "01.235211" ), MpIeee( "01.242614" ),
      MpIeee( "01.250034" ), MpIeee( "01.257471" ), MpIeee( "01.264926" ), MpIeee( "01.272398" ), MpIeee( "01.279888" ),
      MpIeee( "01.287395" ), MpIeee( "01.294921" ), MpIeee( "01.302464" ), MpIeee( "01.310026" ), MpIeee( "01.317605" ),
      MpIeee( "01.325203" ), MpIeee( "01.332819" ), MpIeee( "01.340454" ), MpIeee( "01.348108" ), MpIeee( "01.355780" ),
      MpIeee( "01.363472" ), MpIeee( "01.371182" ), MpIeee( "01.378912" ), MpIeee( "01.386660" ), MpIeee( "01.394429" ),
      MpIeee( "01.402216" ), MpIeee( "01.410024" ), MpIeee( "01.417851" ), MpIeee( "01.425698" ), MpIeee( "01.433565" ),
      MpIeee( "01.441453" ), MpIeee( "01.449360" ), MpIeee( "01.457288" ), MpIeee( "01.465237" ), MpIeee( "01.473206" ),
      MpIeee( "01.481196" ), MpIeee( "01.489208" ), MpIeee( "01.497240" ), MpIeee( "01.505293" ), MpIeee( "01.513368" ),
      MpIeee( "01.521465" ), MpIeee( "01.529583" ), MpIeee( "01.537723" ), MpIeee( "01.545885" ), MpIeee( "01.554068" ),
      MpIeee( "01.562275" ), MpIeee( "01.570503" ), MpIeee( "01.578754" ), MpIeee( "01.587028" ), MpIeee( "01.595325" ),
      MpIeee( "01.603644" ), MpIeee( "01.611987" ), MpIeee( "01.620353" ), MpIeee( "01.628743" ), MpIeee( "01.637156" ),
      MpIeee( "01.645593" ), MpIeee( "01.654053" ), MpIeee( "01.662538" ), MpIeee( "01.671047" ), MpIeee( "01.679581" ),
      MpIeee( "01.688139" ), MpIeee( "01.696721" ), MpIeee( "01.705329" ), MpIeee( "01.713961" ), MpIeee( "01.722619" ),
      MpIeee( "01.731303" ), MpIeee( "01.740011" ), MpIeee( "01.748746" ), MpIeee( "01.757506" ), MpIeee( "01.766293" ),
      MpIeee( "01.775106" ), MpIeee( "01.783945" ), MpIeee( "01.792810" ), MpIeee( "01.801703" ), MpIeee( "01.810623" ),
      MpIeee( "01.819569" ), MpIeee( "01.828543" ), MpIeee( "01.837545" ), MpIeee( "01.846574" ), MpIeee( "01.855631" ),
      MpIeee( "01.864717" ), MpIeee( "01.873830" ), MpIeee( "01.882972" ), MpIeee( "01.892143" ), MpIeee( "01.901343" ),
      MpIeee( "01.910572" ), MpIeee( "01.919830" ), MpIeee( "01.929117" ), MpIeee( "01.938434" ), MpIeee( "01.947781" ),
      MpIeee( "01.957158" ), MpIeee( "01.966566" ), MpIeee( "01.976004" ), MpIeee( "01.985473" ), MpIeee( "01.994972" ),
      MpIeee( "02.004503" ), MpIeee( "02.014065" ), MpIeee( "02.023659" ), MpIeee( "02.033285" ), MpIeee( "02.042943" ),
      MpIeee( "02.052633" ), MpIeee( "02.062355" ), MpIeee( "02.072110" ), MpIeee( "02.081899" ), MpIeee( "02.091720" ),
      MpIeee( "02.101575" ), MpIeee( "02.111464" ), MpIeee( "02.121386" ), MpIeee( "02.131343" ), MpIeee( "02.141334" ),
      MpIeee( "02.151360" ), MpIeee( "02.161421" ), MpIeee( "02.171517" ), MpIeee( "02.181648" ), MpIeee( "02.191815" ),
      MpIeee( "02.202018" ), MpIeee( "02.212257" ), MpIeee( "02.222533" ), MpIeee( "02.232845" ), MpIeee( "02.243195" ),
      MpIeee( "02.253582" ), MpIeee( "02.264006" ), MpIeee( "02.274468" ), MpIeee( "02.284968" ), MpIeee( "02.295507" ),
      MpIeee( "02.306084" ), MpIeee( "02.316701" ), MpIeee( "02.327356" ), MpIeee( "02.338051" ), MpIeee( "02.348786" ),
      MpIeee( "02.359562" ), MpIeee( "02.370377" ), MpIeee( "02.381234" ), MpIeee( "02.392131" ), MpIeee( "02.403070" ),
      MpIeee( "02.414051" ), MpIeee( "02.425073" ), MpIeee( "02.436138" ), MpIeee( "02.447246" ), MpIeee( "02.458397" ),
      MpIeee( "02.469591" ), MpIeee( "02.480828" ), MpIeee( "02.492110" ), MpIeee( "02.503436" ), MpIeee( "02.514807" ),
      MpIeee( "02.526222" ), MpIeee( "02.537684" ), MpIeee( "02.549190" ), MpIeee( "02.560743" ), MpIeee( "02.572343" ),
      MpIeee( "02.583989" ), MpIeee( "02.595682" ), MpIeee( "02.607423" ), MpIeee( "02.619212" ), MpIeee( "02.631050" ),
      MpIeee( "02.642936" ), MpIeee( "02.654871" ), MpIeee( "02.666855" ), MpIeee( "02.678890" ), MpIeee( "02.690975" ),
      MpIeee( "02.703110" ), MpIeee( "02.715297" ), MpIeee( "02.727535" ), MpIeee( "02.739825" ), MpIeee( "02.752168" ),
      MpIeee( "02.764563" ), MpIeee( "02.777012" ), MpIeee( "02.789514" ), MpIeee( "02.802070" ), MpIeee( "02.814681" ),
      MpIeee( "02.827347" ), MpIeee( "02.840069" ), MpIeee( "02.852846" ), MpIeee( "02.865680" ), MpIeee( "02.878570" ),
      MpIeee( "02.891518" ), MpIeee( "02.904524" ), MpIeee( "02.917588" ), MpIeee( "02.930712" ), MpIeee( "02.943894" ),
      MpIeee( "02.957136" ), MpIeee( "02.970439" ), MpIeee( "02.983802" ), MpIeee( "02.997227" ), MpIeee( "03.010714" ),
      MpIeee( "03.024263" ), MpIeee( "03.037875" ), MpIeee( "03.051551" ), MpIeee( "03.065290" ), MpIeee( "03.079095" ),
      MpIeee( "03.092965" ), MpIeee( "03.106900" ), MpIeee( "03.120902" ), MpIeee( "03.134971" ), MpIeee( "03.149107" ),
      MpIeee( "03.163312" ), MpIeee( "03.177585" ), MpIeee( "03.191928" ), MpIeee( "03.206340" ), MpIeee( "03.220824" ),
      MpIeee( "03.235378" ), MpIeee( "03.250005" ), MpIeee( "03.264704" ), MpIeee( "03.279477" ), MpIeee( "03.294323" ),
      MpIeee( "03.309244" ), MpIeee( "03.324240" ), MpIeee( "03.339312" ), MpIeee( "03.354461" ), MpIeee( "03.369687" ),
      MpIeee( "03.384992" ), MpIeee( "03.400375" ), MpIeee( "03.415838" ), MpIeee( "03.431381" ), MpIeee( "03.447005" ),
      MpIeee( "03.462711" ), MpIeee( "03.478500" ), MpIeee( "03.494372" ), MpIeee( "03.510328" ), MpIeee( "03.526370" ),
      MpIeee( "03.542497" ), MpIeee( "03.558711" ), MpIeee( "03.575012" ), MpIeee( "03.591402" ), MpIeee( "03.607881" ),
      MpIeee( "03.624450" ), MpIeee( "03.641111" ), MpIeee( "03.657863" ), MpIeee( "03.674708" ), MpIeee( "03.691646" ),
      MpIeee( "03.708680" ), MpIeee( "03.725809" ), MpIeee( "03.743034" ), MpIeee( "03.760357" ), MpIeee( "03.777779" ),
      MpIeee( "03.795300" ), MpIeee( "03.812921" ), MpIeee( "03.830645" ), MpIeee( "03.848470" ), MpIeee( "03.866400" ),
      MpIeee( "03.884434" ), MpIeee( "03.902574" ), MpIeee( "03.920821" ), MpIeee( "03.939176" ), MpIeee( "03.957640" ),
      MpIeee( "03.976215" ), MpIeee( "03.994901" ), MpIeee( "04.013699" ), MpIeee( "04.032612" ), MpIeee( "04.051639" ),
      MpIeee( "04.070783" ), MpIeee( "04.090045" ), MpIeee( "04.109425" ), MpIeee( "04.128925" ), MpIeee( "04.148547" ),
      MpIeee( "04.168292" ), MpIeee( "04.188160" ), MpIeee( "04.208154" ), MpIeee( "04.228275" ), MpIeee( "04.248524" ),
      MpIeee( "04.268903" ), MpIeee( "04.289413" ), MpIeee( "04.310056" ), MpIeee( "04.330832" ), MpIeee( "04.351745" ),
      MpIeee( "04.372794" ), MpIeee( "04.393982" ), MpIeee( "04.415310" ), MpIeee( "04.436781" ), MpIeee( "04.458395" ),
      MpIeee( "04.480154" ), MpIeee( "04.502060" ), MpIeee( "04.524114" ), MpIeee( "04.546319" ), MpIeee( "04.568676" ),
      MpIeee( "04.591187" ), MpIeee( "04.613854" ), MpIeee( "04.636678" ), MpIeee( "04.659662" ), MpIeee( "04.682807" ),
      MpIeee( "04.706116" ), MpIeee( "04.729590" ), MpIeee( "04.753231" ), MpIeee( "04.777041" ), MpIeee( "04.801024" ),
      MpIeee( "04.825179" ), MpIeee( "04.849511" ), MpIeee( "04.874020" ), MpIeee( "04.898710" ), MpIeee( "04.923582" ),
      MpIeee( "04.948639" ), MpIeee( "04.973883" ), MpIeee( "04.999316" ), MpIeee( "05.024942" ), MpIeee( "05.050761" ),
      MpIeee( "05.076778" ), MpIeee( "05.102993" ), MpIeee( "05.129411" ), MpIeee( "05.156034" ), MpIeee( "05.182864" ),
      MpIeee( "05.209903" ), MpIeee( "05.237156" ), MpIeee( "05.264625" ), MpIeee( "05.292312" ), MpIeee( "05.320220" ),
      MpIeee( "05.348354" ), MpIeee( "05.376714" ), MpIeee( "05.405306" ), MpIeee( "05.434131" ), MpIeee( "05.463193" ),
      MpIeee( "05.492496" ), MpIeee( "05.522042" ), MpIeee( "05.551836" ), MpIeee( "05.581880" ), MpIeee( "05.612178" ),
      MpIeee( "05.642734" ), MpIeee( "05.673552" ), MpIeee( "05.704634" ), MpIeee( "05.735986" ), MpIeee( "05.767610" ),
      MpIeee( "05.799512" ), MpIeee( "05.831694" ), MpIeee( "05.864161" ), MpIeee( "05.896918" ), MpIeee( "05.929968" ),
      MpIeee( "05.963316" ), MpIeee( "05.996967" ), MpIeee( "06.030925" ), MpIeee( "06.065194" ), MpIeee( "06.099780" ),
      MpIeee( "06.134687" ), MpIeee( "06.169921" ), MpIeee( "06.205486" ), MpIeee( "06.241387" ), MpIeee( "06.277630" ),
      MpIeee( "06.314220" ), MpIeee( "06.351163" ), MpIeee( "06.388465" ), MpIeee( "06.426130" ), MpIeee( "06.464166" ),
      MpIeee( "06.502578" ), MpIeee( "06.541371" ), MpIeee( "06.580553" ), MpIeee( "06.620130" ), MpIeee( "06.660109" ),
      MpIeee( "06.700495" ), MpIeee( "06.741297" ), MpIeee( "06.782520" ), MpIeee( "06.824173" ), MpIeee( "06.866262" ),
      MpIeee( "06.908795" ), MpIeee( "06.951780" ), MpIeee( "06.995225" ), MpIeee( "07.039137" ), MpIeee( "07.083525" ),
      MpIeee( "07.128398" ), MpIeee( "07.173764" ), MpIeee( "07.219632" ), MpIeee( "07.266011" ), MpIeee( "07.312910" ),
      MpIeee( "07.360339" ), MpIeee( "07.408308" ), MpIeee( "07.456827" ), MpIeee( "07.505905" ), MpIeee( "07.555554" ),
      MpIeee( "07.605785" ), MpIeee( "07.656608" ), MpIeee( "07.708035" ), MpIeee( "07.760077" ), MpIeee( "07.812747" ),
      MpIeee( "07.866057" ), MpIeee( "07.920019" ), MpIeee( "07.974647" ), MpIeee( "08.029953" ), MpIeee( "08.085952" ),
      MpIeee( "08.142657" ), MpIeee( "08.200083" ), MpIeee( "08.258245" ), MpIeee( "08.317158" ), MpIeee( "08.376837" ),
      MpIeee( "08.437300" ), MpIeee( "08.498562" ), MpIeee( "08.560641" ), MpIeee( "08.623554" ), MpIeee( "08.687319" ),
      MpIeee( "08.751955" ), MpIeee( "08.817481" ), MpIeee( "08.883916" ), MpIeee( "08.951282" ), MpIeee( "09.019600" ),
      MpIeee( "09.088889" ), MpIeee( "09.159174" ), MpIeee( "09.230477" ), MpIeee( "09.302822" ), MpIeee( "09.376233" ),
      MpIeee( "09.450735" ), MpIeee( "09.526355" ), MpIeee( "09.603118" ), MpIeee( "09.681054" ), MpIeee( "09.760191" ),
      MpIeee( "09.840558" ), MpIeee( "09.922186" ), MpIeee( "10.005107" ), MpIeee( "10.089353" ), MpIeee( "10.174959" ),
      MpIeee( "10.261958" ), MpIeee( "10.350389" ), MpIeee( "10.440287" ), MpIeee( "10.531693" ), MpIeee( "10.624646" ),
      MpIeee( "10.719188" ), MpIeee( "10.815362" ), MpIeee( "10.913214" ), MpIeee( "11.012789" ), MpIeee( "11.114137" ),
      MpIeee( "11.217307" ), MpIeee( "11.322352" ), MpIeee( "11.429325" ), MpIeee( "11.538283" ), MpIeee( "11.649285" ),
      MpIeee( "11.762390" ), MpIeee( "11.877664" ), MpIeee( "11.995170" ), MpIeee( "12.114979" ), MpIeee( "12.237161" ),
      MpIeee( "12.361791" ), MpIeee( "12.488946" ), MpIeee( "12.618708" ), MpIeee( "12.751161" ), MpIeee( "12.886394" ),
      MpIeee( "13.024498" ), MpIeee( "13.165570" ), MpIeee( "13.309711" ), MpIeee( "13.457026" ), MpIeee( "13.607625" ),
      MpIeee( "13.761625" ), MpIeee( "13.919145" ), MpIeee( "14.080314" ), MpIeee( "14.245263" ), MpIeee( "14.414134" ),
      MpIeee( "14.587072" ), MpIeee( "14.764233" ), MpIeee( "14.945778" ), MpIeee( "15.131877" ), MpIeee( "15.322712" ),
      MpIeee( "15.518470" ), MpIeee( "15.719353" ), MpIeee( "15.925570" ), MpIeee( "16.137345" ), MpIeee( "16.354912" ),
      MpIeee( "16.578520" ), MpIeee( "16.808433" ), MpIeee( "17.044929" ), MpIeee( "17.288305" ), MpIeee( "17.538873" ),
      MpIeee( "17.796967" ), MpIeee( "18.062943" ), MpIeee( "18.337176" ), MpIeee( "18.620068" ), MpIeee( "18.912049" ),
      MpIeee( "19.213574" ), MpIeee( "19.525133" ), MpIeee( "19.847249" ), MpIeee( "20.180480" ), MpIeee( "20.525429" ),
      MpIeee( "20.882738" ), MpIeee( "21.253102" ), MpIeee( "21.637266" ), MpIeee( "22.036036" ), MpIeee( "22.450278" ),
      MpIeee( "22.880933" ), MpIeee( "23.329017" ), MpIeee( "23.795634" ), MpIeee( "24.281981" ), MpIeee( "24.789364" ),
      MpIeee( "25.319207" ), MpIeee( "25.873062" ), MpIeee( "26.452634" ), MpIeee( "27.059789" ), MpIeee( "27.696581" ),
      MpIeee( "28.365274" ), MpIeee( "29.068370" ), MpIeee( "29.808638" ), MpIeee( "30.589157" ), MpIeee( "31.413354" ),
      MpIeee( "32.285060" ), MpIeee( "33.208568" ), MpIeee( "34.188705" ), MpIeee( "35.230920" ), MpIeee( "36.341388" ),
      MpIeee( "37.527131" ), MpIeee( "38.796172" ), MpIeee( "40.157721" ), MpIeee( "41.622399" ), MpIeee( "43.202525" ),
      MpIeee( "44.912465" ), MpIeee( "46.769077" ), MpIeee( "48.792279" ), MpIeee( "51.005773" ), MpIeee( "53.437996" ),
      MpIeee( "56.123356" ), MpIeee( "59.103894" )
    };
  MpIeee X;MpIeee  U;MpIeee  V;MpIeee  RANLAN;
  int  I;

  X = gsl_rng_uniform_pos(r);
  U = MpIeee( "1000.0" ) * X;
  I = U.toInt();
  U = U - I;

  if (I >= 70 && I <= 800)
    {
      RANLAN = F[I] + U * (F[I + 1] - F[I]);
    }
  else if (I >= 7 && I <= 980)
    {
      RANLAN = F[I] 
        + U * (F[I + 1] - F[I] 
               - MpIeee( "0.25" ) * (MpIeee( "1" ) - U) * (F[I + 2] - F[I + 1] - F[I] + F[I - 1]));
    }
  else if (I < 7)
    {
      V = log(X);
      U = MpIeee( "1" ) / V;
      RANLAN = ((MpIeee( "0.99858950" ) + (MpIeee( "3.45213058E1" ) + MpIeee( "1.70854528E1" ) * U) * U) /
                (MpIeee( "1" ) + (MpIeee( "3.41760202E1" ) + MpIeee( "4.01244582" ) * U) * U)) *
               ( -log( -MpIeee( "0.91893853" ) - V) - MpIeee( "1" ));
    }
  else
    {
      U = MpIeee( "1" ) - X;
      V = U * U;
      if (X <= MpIeee( "0.999" ))
        {
          RANLAN = (MpIeee( "1.00060006" ) + MpIeee( "2.63991156E2" ) * U + MpIeee( "4.37320068E3" ) * V) /
                   ((MpIeee( "1" ) + MpIeee( "2.57368075E2" ) * U + MpIeee( "3.41448018E3" ) * V) * U);
        }
      else
        {
          RANLAN = (MpIeee( "1.00001538" ) + MpIeee( "6.07514119E3" ) * U + MpIeee( "7.34266409E5" ) * V) /
                   ((MpIeee( "1" ) + MpIeee( "6.06511919E3" ) * U + MpIeee( "6.94021044E5" ) * V) * U);
        }
    }

  return RANLAN;
}

