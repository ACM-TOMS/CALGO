C
C------------------------------------------------------
C
C   Problem:    RE_ENTRY
C
C               Re-entry maneuver, i.e. optimal control
C               problem
C
C------------------------------------------------------
C
*     REAL CONSTANT
      RE = 6378.0
      G0 = 9.8065E-3
      RHO0 = 1.25E+9
      BETA1 = 6.9
      C1 = -0.505
      C2 = 0.88
      C3 = 0.52
      MS = 3.0E+8
      PI = 3.141592654
C
*     MACRO ALPHA_T
      ALPHA_T = CA1 + CA2*1.0E-2*T + CA3*1.0E-4*T**2
     /        + CA4*1.0E-6*T**3
     /        + CA5*1.0E-8*T**4
     /        + CA6*1.0E-10*T**5
     /        + CA7*1.0E-12*T**6
C
*     VARIABLE
      CA1, CA2, CA3, CA4, CA5, CA6, CA7, PSI, ZETA, V
      GAMMA, Q, KAPPA, T
C
*     FUNCTION PSIP
      H = ZETA*RE
      R = RE + H
      PSIP = V/R*COS(GAMMA)
C
*     FUNCTION ZETAP
      ZETAP = V/RE*SIN(GAMMA)
C
*     FUNCTION VP
      ALPHA = ALPHA_T*PI/180.0
      CD = C2 + C3*COS(ALPHA)
      RHO = RHO0*EXP(-H/BETA1)
      VP = -0.5*RHO*V**2/MS*CD - G0/(1+ZETA)**2*SIN(GAMMA)
C
*     FUNCTION GAMMAP
      CL = C1*SIN(ALPHA)
      GAMMAP = -0.5*RHO*V/MS*CL + (V/(RE*(1+ZETA))
     /           - G0/(V*(1+ZETA)**2))*COS(GAMMA)
C
*     FUNCTION QP
      QP = SQRT(RHO)*V**3
C
*     FUNCTION KAPPAP
      IF (GAMMA.GT.0) THEN
         KAPPAP = -GAMMA
      ELSE
         KAPPAP = 0.0
      ENDIF
C
*     END
C
