FC = epcf90
FFLAGS = -sloppy -C -d5 -g -temp=/tmp -u
TIMER = /nfs/anacreon/d1/hometrh/CVSROOT/Libs/EpcF90/timer.o

admfor: admfor.o src.o
	$(FC) $(FFLAGS) admfor.o src.o -o admfor

adm: adm.o adm_gen.o src.o
	$(FC) $(FFLAGS) adm.o adm_gen.o src.o -o adm

all: chem_rea epe epe_n expl1 fish helical helm_00 helm_05 helm_10 helm_20 \\
     helm_30 helm_40 helm_50 helm_60 helm_70 helm_80 helm_90 hf inertia prtcos \\
     quart_02 quart_10 quart_30 quart_50 re_entry test1 tp1 tp268 tp268_2 \\
     tp295_00 tp295_10 tp295_20 tp295_30 tp295_40 tp295_50 tp295_60 tp295_70 \\
     tp295_80 tp295_90 tp32 tp370 tp45 tp55

chem_rea: Examples/chem_rea.fun Examples/chem_rea.inp admfor adm.o src.o
	cp Examples/chem_rea.inp ADM.INP
	cp Examples/chem_rea.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/chem_rea.res

epe: Examples/epe.fun Examples/epe.inp admfor adm.o src.o
	cp Examples/epe.inp ADM.INP
	cp Examples/epe.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/epe.res

epe_n: Examples/epe_n.fun Examples/epe_n.inp admfor adm.o src.o
	cp Examples/epe_n.inp ADM.INP
	cp Examples/epe_n.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/epe_n.res

expl1: Examples/expl1.fun Examples/expl1.inp admfor adm.o src.o
	cp Examples/expl1.inp ADM.INP
	cp Examples/expl1.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/expl1.res

fish: Examples/fish.fun Examples/fish.inp admfor adm.o src.o
	cp Examples/fish.inp ADM.INP
	cp Examples/fish.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/fish.res

helical: Examples/helical.fun Examples/helical.inp admfor adm.o src.o
	cp Examples/helical.inp ADM.INP
	cp Examples/helical.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helical.res

helm_00: Examples/helm_00.fun Examples/helm_00.inp admfor adm.o src.o
	cp Examples/helm_00.inp ADM.INP
	cp Examples/helm_00.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_00.res

helm_05: Examples/helm_05.fun Examples/helm_05.inp admfor adm.o src.o
	cp Examples/helm_05.inp ADM.INP
	cp Examples/helm_05.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_05.res

helm_10: Examples/helm_10.fun Examples/helm_10.inp admfor adm.o src.o
	cp Examples/helm_10.inp ADM.INP
	cp Examples/helm_10.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_10.res

helm_20: Examples/helm_20.fun Examples/helm_20.inp admfor adm.o src.o
	cp Examples/helm_20.inp ADM.INP
	cp Examples/helm_20.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_20.res

helm_30: Examples/helm_30.fun Examples/helm_30.inp admfor adm.o src.o
	cp Examples/helm_30.inp ADM.INP
	cp Examples/helm_30.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_30.res

helm_40: Examples/helm_40.fun Examples/helm_40.inp admfor adm.o src.o
	cp Examples/helm_40.inp ADM.INP
	cp Examples/helm_40.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_40.res

helm_50: Examples/helm_50.fun Examples/helm_50.inp admfor adm.o src.o
	cp Examples/helm_50.inp ADM.INP
	cp Examples/helm_50.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_50.res

helm_60: Examples/helm_60.fun Examples/helm_60.inp admfor adm.o src.o
	cp Examples/helm_60.inp ADM.INP
	cp Examples/helm_60.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_60.res

helm_70: Examples/helm_70.fun Examples/helm_70.inp admfor adm.o src.o
	cp Examples/helm_70.inp ADM.INP
	cp Examples/helm_70.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_70.res

helm_80: Examples/helm_80.fun Examples/helm_80.inp admfor adm.o src.o
	cp Examples/helm_80.inp ADM.INP
	cp Examples/helm_80.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_80.res

helm_90: Examples/helm_90.fun Examples/helm_90.inp admfor adm.o src.o
	cp Examples/helm_90.inp ADM.INP
	cp Examples/helm_90.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/helm_90.res

hf: Examples/hf.fun Examples/hf.inp admfor adm.o src.o
	cp Examples/hf.inp ADM.INP
	cp Examples/hf.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/hf.res

inertia: Examples/inertia.fun Examples/inertia.inp admfor adm.o src.o
	cp Examples/inertia.inp ADM.INP
	cp Examples/inertia.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/inertia.res

prtcos: Examples/prtcos.fun Examples/prtcos.inp admfor adm.o src.o
	cp Examples/prtcos.inp ADM.INP
	cp Examples/prtcos.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/prtcos.res

quart_02: Examples/quart_02.fun Examples/quart_02.inp admfor adm.o src.o
	cp Examples/quart_02.inp ADM.INP
	cp Examples/quart_02.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/quart_02.res

quart_10: Examples/quart_10.fun Examples/quart_10.inp admfor adm.o src.o
	cp Examples/quart_10.inp ADM.INP
	cp Examples/quart_10.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/quart_10.res

quart_30: Examples/quart_30.fun Examples/quart_30.inp admfor adm.o src.o
	cp Examples/quart_30.inp ADM.INP
	cp Examples/quart_30.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/quart_30.res

quart_50: Examples/quart_50.fun Examples/quart_50.inp admfor adm.o src.o
	cp Examples/quart_50.inp ADM.INP
	cp Examples/quart_50.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/quart_50.res

re_entry: Examples/re_entry.fun Examples/re_entry.inp admfor adm.o src.o
	cp Examples/re_entry.inp ADM.INP
	cp Examples/re_entry.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/re_entry.res

test1: Examples/test1.fun Examples/test1.inp admfor adm.o src.o
	cp Examples/test1.inp ADM.INP
	cp Examples/test1.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/test1.res

tp1: Examples/tp1.fun Examples/tp1.inp admfor adm.o src.o
	cp Examples/tp1.inp ADM.INP
	cp Examples/tp1.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp1.res

tp268: Examples/tp268.fun Examples/tp268.inp admfor adm.o src.o
	cp Examples/tp268.inp ADM.INP
	cp Examples/tp268.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp268.res

tp268_2: Examples/tp268_2.fun Examples/tp268_2.inp admfor adm.o src.o
	cp Examples/tp268_2.inp ADM.INP
	cp Examples/tp268_2.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp268_2.res

tp295_00: Examples/tp295_00.fun Examples/tp295_00.inp admfor adm.o src.o
	cp Examples/tp295_00.inp ADM.INP
	cp Examples/tp295_00.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_00.res

tp295_10: Examples/tp295_10.fun Examples/tp295_10.inp admfor adm.o src.o
	cp Examples/tp295_10.inp ADM.INP
	cp Examples/tp295_10.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_10.res

tp295_20: Examples/tp295_20.fun Examples/tp295_20.inp admfor adm.o src.o
	cp Examples/tp295_20.inp ADM.INP
	cp Examples/tp295_20.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_20.res

tp295_30: Examples/tp295_30.fun Examples/tp295_30.inp admfor adm.o src.o
	cp Examples/tp295_30.inp ADM.INP
	cp Examples/tp295_30.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_30.res

tp295_40: Examples/tp295_40.fun Examples/tp295_40.inp admfor adm.o src.o
	cp Examples/tp295_40.inp ADM.INP
	cp Examples/tp295_40.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_40.res

tp295_50: Examples/tp295_50.fun Examples/tp295_50.inp admfor adm.o src.o
	cp Examples/tp295_50.inp ADM.INP
	cp Examples/tp295_50.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_50.res

tp295_60: Examples/tp295_60.fun Examples/tp295_60.inp admfor adm.o src.o
	cp Examples/tp295_60.inp ADM.INP
	cp Examples/tp295_60.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_60.res

tp295_70: Examples/tp295_70.fun Examples/tp295_70.inp admfor adm.o src.o
	cp Examples/tp295_70.inp ADM.INP
	cp Examples/tp295_70.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_70.res

tp295_80: Examples/tp295_80.fun Examples/tp295_80.inp admfor adm.o src.o
	cp Examples/tp295_80.inp ADM.INP
	cp Examples/tp295_80.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_80.res

tp295_90: Examples/tp295_90.fun Examples/tp295_90.inp admfor adm.o src.o
	cp Examples/tp295_90.inp ADM.INP
	cp Examples/tp295_90.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp295_90.res

tp32: Examples/tp32.fun Examples/tp32.inp admfor adm.o src.o
	cp Examples/tp32.inp ADM.INP
	cp Examples/tp32.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp32.res

tp370: Examples/tp370.fun Examples/tp370.inp admfor adm.o src.o
	cp Examples/tp370.inp ADM.INP
	cp Examples/tp370.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp370.res

tp45: Examples/tp45.fun Examples/tp45.inp admfor adm.o src.o
	cp Examples/tp45.inp ADM.INP
	cp Examples/tp45.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp45.res

tp55: Examples/tp55.fun Examples/tp55.inp admfor adm.o src.o
	cp Examples/tp55.inp ADM.INP
	cp Examples/tp55.fun ADM.FUN
	admfor
	$(FC) $(FFLAGS) adm.o adm_gen.f src.o $(TIMER) -o adm
	adm
	cp ADM.RES Results/tp55.res
