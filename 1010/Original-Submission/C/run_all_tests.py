#!/usr/bin/env python3
import sys, time, os, math 
def del_unuseful():
    os.system('rm -f timing_???_sample?.dat')
    os.system('rm -f P_of_eps*.dat')
    os.system('rm -f F_of_eps_rel-HQRL.dat')
    os.system('rm -f statanalysis.out')
    os.system('rm -f timingtests.out')
AT='accuracytest.out'
SA='statanalysis.out'
TT='timingtests.out'
#number of quartics to generate for statistical analysis
SAT=100000
#executable for accuracy tests
EXE_AT="./bin/accuracytest"
#executable for statistical tests
EXE_ST='./bin/statanalysis'
#executable for timing test of sample B 
EXE_TT_SB='./bin/timingtest'
#executable for timing test of sample F
EXE_TT_SF='./bin/timingtest_sample_F'
#number of quartic to generate for timing tests
NTRIALS=500000
#number of runs to average over in the timing tests
NRUNS=10
del_unuseful()
if (not os.path.exists(EXE_AT)) or (not os.path.exists(EXE_ST)) or (not os.path.exists(EXE_TT_SB)) or (not os.path.exists(EXE_TT_SF)):
    print("Before launching this python script please compile the executables!")
    quit()
def print_error():
    print('run_all_tests [-nruns/-nr <number of runs for timing tests>|-nstat/-ns <number of quartics for statistical tests>|')
    print('-ntimings/-nt <number of quartics for timing tests>]')
    quit()
#parse command line arguments
args=sys.argv
del(args[0])
itargs=iter(args)
for a in itargs:
    if a == '-nstat' or a == '-ns':
        SAT=int(next(itargs))
    elif a == '-nruns' or a == '-nr':
        NRUNS=int(next(itargs))
    elif a == '-ntimings' or a == '-nt':
        SAT=int(next(itargs))
    else:
        print_error()
        quit()
##make clean
##make
############################
##### ACCURACY TESTS #######
############################
print('Performing all accuracy tests...', end='')
sys.stdout.flush()
for i in range(1,25):
    os.system(EXE_AT+' '+str(i)+' >> ' + AT) 
print('done!')
############################
### STATISTICAL ANALYSIS ###
############################
print('Performing all statistical tests:')
sys.stdout.flush()
mapst=[ 'A', 'B', 'C', 'D', 'E', 'F' ]
for i in range(0,6):
    print('Testing sample ' + mapst[i] + ' for all algorithms...',end='')
    sys.stdout.flush()
    os.system(EXE_ST+' '+str(SAT)+' 10 '+str(i)+ ' -1 >> ' + SA)
    print('done!')
print('All statistal tests done!')
os.system('rm -f PE-???-.dat')
############################
###### TIMING TESTS  #######
############################
print('Performing all timing tests:')
def timingtest(execname,ssample,nruns,ntrials):
    mapn=[ 'DRY', 'ODM', 'FLO', 'STR', 'FER', 'FQS', 'HQR', 'SHM' ] 
    for ity in range(0, 7):
        fn='timing_'+mapn[ity]+'_sample'+ssample+'.dat'
        print('Testing '+ mapn[ity]+'\t', end='') 
        sys.stdout.flush()
        with open(fn,'w') as ff:
            for cc in range(1,nruns+1):
                initime=time.time()         
                os.system('nice -n 0 '+execname+' '+str(ity)+' '+str(ntrials) + ' >> ' + TT) 
                endtime=time.time() 
                dtime=endtime-initime
                ff.write(str(cc)+' '+str(dtime)+'\n')
                print('.',end='')
                sys.stdout.flush()
        print('done')
    #calculate averages and error (stddev)
    #estimate dry run first
    ity=0
    fn='timing_'+mapn[ity]+'_sample'+ssample+'.dat'
    sumt=0
    sumsq=0
    cc=0
    with open(fn,'r') as ff:
        lines=ff.readlines()
    for ll in lines:
        lst=ll.strip('\n').split(' ')
        tt=float(lst[1])
        sumt += tt
        sumsq+= tt*tt
        cc += 1
    fcc=float(cc)
    avgt_dry=sumt/fcc
    stdd_dry=sumsq/fcc - (sumt/fcc)*(sumt/fcc)
    oft='timings_sample'+ssample+'.txt'
    with open(oft,'w') as ff:
        ff.write('ALGO'.rjust(4)+' '+'AVG'.rjust(12)+' '+'STDDEV'.rjust(12) + '\n')
        for ity in range(1,7):
            fn='timing_'+mapn[ity]+'_sample'+ssample+'.dat'
            with open(fn,'r') as fft:
                lines=fft.readlines()
            sumt=0  
            sumsq=0
            cc=0
            for ll in lines:
                lst=ll.strip('\n').split(' ')
                tt=float(lst[1])
                sumt += tt
                sumsq+= tt*tt
                cc += 1
            fcc=float(cc)
            avgt=sumt/fcc-avgt_dry
            stdd=sumsq/fcc - (sumt/fcc)*(sumt/fcc)
            stdd += stdd_dry
            ff.write(mapn[ity].rjust(4)+' '+'{:12.6f}'.format(avgt)+' '+'{:12.6f}'.format(math.sqrt(stdd))+'\n')
#timing test of sample B
print('Timing tests of Sample B')
timingtest(EXE_TT_SB,'B',NRUNS,NTRIALS)
#timing test of sample F
print('Timing tests of Sample F')
timingtest(EXE_TT_SF,'F',NRUNS,NTRIALS)
del_unuseful()
print('All timing tests done!')
