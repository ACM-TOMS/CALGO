@echo off
if "%1"=="" goto 1
copy examples\vb%1.c code\main.c
chdir code
del main.obj
del tst.*
make tst.out
rename tst.out vb%1.out
chdir ..
goto 2
:1
echo ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
echo ³    Available examples for Cubpack++             ³
echo ³                                                 ³
echo ³   vb1      CIRCLE                               ³
echo ³   vb2      PLANE_SECTOR (quadrant)              ³
echo ³   vb3      PLANE_SECTOR (wedge)                 ³
echo ³   vb4      GENERALIZED_RECTANGLE (parabola)     ³
echo ³   vb5      GENERALIZED_RECTANGLE (lemniscate)   ³
echo ³   vb6      SEMI_INFINITE_STRIP                  ³
echo ³   vb7      SEMI_INFINITE_STRIP                  ³
echo ³   vb8      INFINITE_STRIP                       ³
echo ³   vb9      RECTANGLE                            ³
echo ³   vb10     TRIANGLE                             ³
echo ³   vb10     TRIANGLE                             ³
echo ³   vb11     PLANE_SECTOR (wedge)                 ³
echo ³                                                 ³
echo ³  Type 'testcase 1'   etc. at the DOS            ³
echo ³  prompt to activate one of these options.       ³
echo ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
:2

