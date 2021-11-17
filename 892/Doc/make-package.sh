#!/bin/bash
set -x
cd ..
tar czf dispmodule.tgz \
  dispmodule/dispmodule.f90 \
  dispmodule/disp_i1mod.f90 \
  dispmodule/disp_i2mod.f90 \
  dispmodule/disp_i4mod.f90 \
  dispmodule/disp_i8mod.f90 \
  dispmodule/disp_l1mod.f90 \
  dispmodule/disp_r16mod.f90 \
  dispmodule/dispmodule_absoft.f90 \
  dispmodule/disp_i1mod_absoft.f90 \
  dispmodule/dispmodule-userman.doc \
  dispmodule/dispmodule-userman.pdf \
  dispmodule/dispmodule-userman.html \
  dispmodule/dispmodule-userman.txt \
  dispmodule/dispmodule.doc \
  dispmodule/changelog.txt \
  dispmodule/README \
  dispmodule/dispdemo.f90 \
  dispmodule/expected-demo-output.txt \
  dispmodule/test_dispmodule.f90 \
  dispmodule/test_dispmodule_fpp.F90 \
  dispmodule/test_naninf.f90 \
  dispmodule/test_naninf_ieee.f90 \
  dispmodule/test_naninf_mod.f90 \
  dispmodule/mexdispdemo.f90 \
  dispmodule/test-all.sh \
  dispmodule/make-package.sh \
  dispmodule/Makefile \
  dispmodule/portability-evidence.txt
version=`awk '$0~"! Version number"{print $4}' dispmodule/dispmodule.f90`
cp dispmodule.tgz dispmodule-v$version.tgz
cd dispmodule
