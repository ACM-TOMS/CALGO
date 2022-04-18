//-*- C++ -*-
#ifndef PHYSCONST_H_
#define PHYSCONST_H_

#include <math.h>
#include "mathconst.h"

const double  ANG=1.0/0.52917724924,
  CENTIMETER = 1.0e+8/0.52917724924;
const double ANG2 = ANG*ANG, ANG4 = ANG2*ANG2, ANG8 = ANG4*ANG4,
  ANG10 = ANG2*ANG8, ANG12 = ANG4*ANG8;

const double NAVOGADRO = 6.022137e+23;

// Masakatsu ITO always uses atomic unit.

const double SEC = 1.0e+0/2.4188843341e-17,
  PICOSEC = SEC*1e-12,
  FEMTOSEC = SEC*1e-15;

const double JMOL = 3.808799e-7,
  EV = 3.674931e-2,
  KCALMOL = 1.593601e-3;

// Boltzmann constant
const double KB = 3.166830e-6;

// Atomic Masses  ref. Rikagaku Jiten 4'th edition
const double AMU = 1.6605655e-27/9.10938970e-31,
  AWH1  = 1.007825*AMU,
  AWC12 = 12.00000*AMU,
  AWN14 = 14.0*AMU, // to be CORRECTED!!
  AWO16 = 15.994915*AMU,
  AWAR40= 39.96238*AMU;

const double COULOMB = 1.0 / 1.602177e-19;

// Lennard Jones Potential Parameters
//
// ref. W.D.Cornel, P.Cieplak, C.I.Bayly, I.R.Gould,
//      K.M.Merz,Jr., D.M.Ferguson, D.C.Spellmeyer, T.Fox,
//      J.W.Caldwell and P.A.Kollman, 
//      J.Am.Chem.Soc, 117, 5179 (1995)
//      Table 14. Van Der Waals Parameters

const double VDW_R_TO_SIGAM = 1.0/pow(2.0,1.0/6.0),
  AMBER_VDW_SIGMA_C = 1.9080 * VDW_R_TO_SIGAM * ANG, // Carbony Carbon
  AMBER_VDW_EPS_C = 0.0860 * KCALMOL,
  AMBER_VDW_SIGMA_CM = 1.9080 * VDW_R_TO_SIGAM * ANG, // sp2 Carbon
  AMBER_VDW_EPS_CM = 0.0860 * KCALMOL,
  AMBER_VDW_SIGMA_CT = 1.9080 * VDW_R_TO_SIGAM * ANG, // sp3 Carbon
  AMBER_VDW_EPS_CT = 0.1094 * KCALMOL,
  AMBER_VDW_SIGMA_HA = 1.4590 * VDW_R_TO_SIGAM * ANG, // Aromatic Hydrogen
  AMBER_VDW_EPS_HA = 0.0150 * KCALMOL,
  AMBER_VDW_SIGMA_HC = 1.4870 * VDW_R_TO_SIGAM * ANG, // Aliphatic Hydrogen
  AMBER_VDW_EPS_HC = 0.0157 * KCALMOL,
  AMBER_VDW_SIGMA_O = 1.6612 * VDW_R_TO_SIGAM * ANG, // Amide Oxygen
  AMBER_VDW_EPS_O = 0.2100 * KCALMOL;

const static double 
AMBER_VDW_SIGMA_CM_HA = ( AMBER_VDW_SIGMA_CM + AMBER_VDW_SIGMA_HA )/2.0,
  AMBER_VDW_EPS_CM_HA = sqrt( AMBER_VDW_EPS_CM * AMBER_VDW_EPS_HA),
  AMBER_VDW_SIGMA_HC_O = ( AMBER_VDW_SIGMA_HC + AMBER_VDW_SIGMA_O )/2.0,
  AMBER_VDW_EPS_HC_O = sqrt( AMBER_VDW_EPS_HC * AMBER_VDW_EPS_O );

const static double 
AMBER_VDW_A_CM = 4.0*AMBER_VDW_EPS_CM * pow(AMBER_VDW_SIGMA_CM,12.0),
  AMBER_VDW_B_CM = 4.0*AMBER_VDW_EPS_CM * pow(AMBER_VDW_SIGMA_CM,6.0),
  AMBER_VDW_A_CM_HA = 4.0*AMBER_VDW_EPS_CM_HA * pow(AMBER_VDW_SIGMA_CM_HA,12.0),
  AMBER_VDW_B_CM_HA = 4.0*AMBER_VDW_EPS_CM_HA * pow(AMBER_VDW_SIGMA_CM_HA,6.0),
  AMBER_VDW_A_HA = 4.0*AMBER_VDW_EPS_HA * pow(AMBER_VDW_SIGMA_HA,12.0),
  AMBER_VDW_B_HA = 4.0*AMBER_VDW_EPS_HA * pow(AMBER_VDW_SIGMA_HA,6.0),
  AMBER_VDW_A_HC_O = 4.0*AMBER_VDW_EPS_HC_O * pow(AMBER_VDW_SIGMA_HC_O,12.0),
  AMBER_VDW_B_HC_O = 4.0*AMBER_VDW_EPS_HC_O * pow(AMBER_VDW_SIGMA_HC_O,6.0);

// Lennard Jones Potential Parameters for Argon
//
// ref. G.C.Maitland, M.Rigby, E.B.Smith and W.A.Wakeham (1981)
//      Intermolecular forces : their origin and determination.
//      Clarendon Press, Oxford.
const double VDW_SIGMA_AR = 3.41*ANG,
  VDW_EPS_AR = 119.8*KB;

const static double VDW_A_AR = 4.0*VDW_EPS_AR * pow(VDW_SIGMA_AR,12.0),
  VDW_B_AR = 4.0*VDW_EPS_AR * pow(VDW_SIGMA_AR,6.0);

#endif // PHYSCONST_H_
