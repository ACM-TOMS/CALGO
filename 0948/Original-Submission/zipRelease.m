clc; clear; close all

CURRENT_DIR = pwd;

mkdir tmp
cd tmp

!svn export https://websvn.mcmaster.ca/daesa/DAESA-1.0

cd DAESA-1.0
!rm zipRelease.m
!rm -rf doc

cd src/implementation/

subdir = {'@qla' '@sigma' 'common' 'dataoutput' 'lapdm' 'visualization'};

for i=1:length(subdir)
    
    cd(subdir{i});
    pcode('./*','-inplace');
    !rm -rf *.m
    cd ..
    
end


cd(CURRENT_DIR)
cd tmp

zip('../../DAESA-1.0.zip',{'./DAESA-1.0'})

cd ..
!rm -rf tmp

