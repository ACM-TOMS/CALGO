%% Structure of the Toolbox IPscatt
%
% This file shows the structure of the toolbox |IPscatt|, i.e. dependencies of functions. 
% The structure of files in folders is shown in <../index.html>.
%
% Please note that this structure represents the most important 
% dependencies of functions in case of typical work flows and therefore is 
% not complete.
%
% We distinguish between internal and external functions.
% The external ones are marked in this overview by [external] to make clear
% that they are more or less easy to use outside the typical work flow of
% |IPscatt|.
%
%
%% Main structure
%
% *<start.html>* calls
%
% * <init.html>
% * <setInput.html>
% * <setData.html>
% * <setRecon.html>
% * <recon.html>
%
% Details are below.
%
%% init
% <init.html>
%
%% setInput
% <setInput.html>
%
% * <dirMake.html>
%
%% setData
% <setData.html>
%
% |If seti.expData = 'fresnel'|
%
% * <checkConsisExpData.html>
% * <loadData.html>
% * ---- <readRAWData.html> [external]
% * ---- <setGeomSim.html> (details of substructure are below)
% * ---- <matchIncField.html> [external]
%
% |else|
%
% * <setGeomSim.html> (details of substructure are below)
% * <mimo.html> [external]
% * <addNoise.html> [external]
%
% |end|
%
%% setRecon
% <setRecon.html>
%
% * <setInvType.html>
% * <checkConsisRec.html>
% * <shrinkFuncs.html> [external]
% * <setFuncsPda.html>
%
%% recon
% <recon.html>
%
% * <subplots.html>
% * <plotAndSaveFigures.html>
% * <pdaPlot.html>
% * <minPda.html>
% * ---- <minTolIn.html>
% * ---- <minTolOut.html>
% * ---- <pda.html>
% * ----|---- <mimo.html> [external]
% * ----|---- <opNormNum.html>
% * ----|---- <pdaChoosingStepsizes.html>
% * ---- <mimo.html> [external]
%
%% Substructure of setGeomSim
% <setGeomSim.html>
%
% * <setGrid.html> [external]
% * <setReshapeVecMat.html>
% * <setGridScale.html>
% * ---- <gridUp.html> [external]
% * ---- <gridDown.html> [external]
% * <setIdImagReal.html> [external]
% * <setKernel.html> [external]
% * <expSetup.html>
% * ---- <plotExpSetup.html>
% * ---- <setIncPnts.html>
% * ----|---- <pntsGeometry.html> [external]
% * ----|----|---- <dS2D.html>
% * ----|---- <pntsGeometry3D.html> [external]
% * ----|----|---- <dS3D.html>
% * ---- <setMeasPnts.html>
% * ----|---- <pntsGeometry.html> [external]
% * ----|----|---- <dS2D.html>
% * ----|---- <pntsGeometry3D.html> [external]
% * ----|----|---- <dS3D.html>
% * ---- <setIncField.html>
% * ---- <setMeasKer.html>
% * <setContrast.html> [external]
% * ---- <plotPredefinedContrast.html>
%
%% Structure above of start
%
% * <varalpha.html>
% * <varbeta.html>
% * <vardelta.html>
% * <vartol.html>
%
%% Convenience Functions
%
% * <demo.html>
%
% * <adjOfDer.html> [external]
% * <derivative.html> [external]
% * <forward.html> [external]
%
%% Predefined Contrasts
% All predefined contrasts in the folder |incontrasts| are qualified for
% external use.
%
%% Other supporting files
%
% * <docCreate.html> [external]
%

