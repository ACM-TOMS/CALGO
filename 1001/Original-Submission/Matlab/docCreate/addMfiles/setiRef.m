%% Reference of the structure array seti
%
% Descriptions of fields of structure array seti.
%

%% More About this file
%
% * This file summarizes the fields of the main structure array |seti| for each function in one file. 
% * Therefore we copied the descriptions of the fields seti, that are available in the sections "Input Arguments" and "Output arguments" of the source code documentation of each file (the order bases on <programStructure.html>). Therefore duplicates will occur. Please note that we only summarized the fields of |seti|. For other input or output arguments see the source code documentation of the corresponding file.
% * In the headlines we do not describe the path but the dependence of the function, e.g. |setData:setGeomSim| does mean that the routine |setGeomSim| depends on the function |setData|.
% * The dependence starts at the function |start|.

%%%%%

%% setInput: Input Arguments
%
% All fields are set automatically, if they does not exist:
%
% * seti.dirOutput  :   name of folder in which directories for output 
%                       files and figures are created.
%                       (default: 'output')
% * seti.dirDatetime    :   date and time separated by the character "T",
%                           e.g. 20161006T105735.
% * seti.dirSuffix      :   suffix for dirname (default: '')
% * seti.dirSuffixAdd   :   additional suffix of dirname 
%                           (e.g. in case of varalpha, varbeta, vardelta, gscale; or expData is 'fresnel')
% * seti.dirname        :   dirname of files and figures for current
%                           computation, default:
%                           |seti.dirOutput + / + seti.datetime + _ + seti.inseti + seti.dirSuffix + seti.dirSuffixAdd|
%                           If seti.inseti is empty, 'noinseti' is used in
%                           dirname, e.g.
%                           |seti.dirname =
%                           'output/20161006T102740_noinseti'|.
% Not set automatically:
%
% * seti.expData : see <setData.html>
%
%% setInput: Output Arguments
%
% * seti.dirDatetime    : date and time (format: YYYYMMDDThhmmss)
% * seti.fileSuffix     : suffix name for filename
%
% * seti.inseti     :   name of function including input parameters in folder inseti.
% * seti.closed     :   is always 0 in public version 
%                       (value 1 would mean that closed code is available)
%                        (see also <init.html> for variable |closed|)

%%%%%

%% setInput:dirMake: Input Arguments
%
% seti : structural array
%
% _The following fields can be set by user (otherwise set default values):_
%
% * seti.gscale  : grid scaling on (1) or off (0), see <setGridScale.html>
% * seti.expData : e.g. 'fresnel' if data from Institute Fresnel are used, see <setData.html>
%
% See <setInput.html> for the following fields (they can be defined by user
% too):
%
% * seti.dirOutput 
% * seti.dirDatetime
% * seti.dirSuffix      
% * seti.dirSuffixAdd
% * seti.dirname
%
%% setInput:dirMake: Output Arguments
%
% Several fields in seti, if they was not set as input argument; see above.
%

%%%%%

%% setData: Input Arguments
%
% We only mention explictly the input arguments to deal with experimentally
% measured data from Institute Fresnel.
%
% * seti.expData    :   Set 'fresnel' to load real-world data
%                       (experimentally measured) from Institute Fresnel.
%
% *The following parameters can be set by user only in case of 
%  seti.expData = 'fresnel'*
% 
% For further details see <loadData.html> and <matchIncField.html>.
%
% * seti.fresnelFreq    :   frequency of Fresnel data set in Hz
%                           (default: 5e+09)
% * seti.fresnelFile    :   path to real-world data
%                           (default:
%                           'inexpdata/fresnel_opus_1/twodielTM_8f.exp')
%                           (Make sure to use data with |TM| in the filemname,
%                           which stands for transverse magnetic (TM) 
%                           polarization, because the Helmholtz equation
%                           only models electromagnetic waves in the TM
%                           case, see Section 2 in [1].)
% * seti.nuMax          :   match incident field parameter (default: 7)
% * seti.ampCalc        :   method to compute the coefficients c (1, 2 or 3),
%                           we recommend to use method 1.
%
%% setData: Output Arguments
%
% * seti    :   structure array
%
% For the following output arguments see <setGeomSim.html>:
%
% * seti.tol
%
%
% For the following output arguments see Subfunction |setFigureSettings| in
% <setGeomSim.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.plotFreq       :   default: 1
% * seti.usecbarlim     :   default: 1
% * seti.cbarlim        :   default: [-0.2 1.4]
% * seti.plotPublish    :   default: 0
% * seti.pubFontSize    :   default: 20
% * seti.plotVisible    :   default: 'off'
% * seti.savepng        :   default: 1
% * seti.saveepsc       :   default: 0
% * seti.savefig        :   default: 0
% * seti.plotFreqiPda   :   default: 0
%
% For the following output arguments see <setGrid.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.dim            :   default: 2
% * seti.rCD            :   default: 0.2
% * seti.nCD            :   default: 256
%
% _Are set automatically_:
%
% * seti.h              :   2*rCD/nCD
% * seti.dV
% * seti.grid           :   size: seti.dim x seti.nCD^seti.dim
% * seti.ballMask       :   size: seti.nCD x seti.nCD (logical)
% * seti.nROI
% * seti.gridROI        :   size: seti.dim x seti.nROI^seti.dim
% * seti.ROImask        :   size: nCD x nCD (logical)
% * seti.ballMaskROI    :   size: nROI x nROI (logical)
%
% For the following output arguments see <setReshapeVecMat.html>:
%
% _Are set automatically:_
%
% * seti.GROI
% * seti.GCD
% * seti.G
% * seti.iG
%
% For the following output arguments see <setGridScale.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.gscale     : default: 0
% * seti.nCDinv     :
%
% _Are set automatically:_
%
% * seti.nInv       :
% * seti.hInv       :
% * seti.dVinv      :
% * seti.GInv       :
% * seti.GU         :
% * seti.GD         :
%
% For the following output arguments see <setIdImagReal.html>:
% 
% _Are set automatically:_
%
% * seti.S  :   
% * seti.R  :   
% * seti.I  :   
% * seti.T  :   
%
% For the following output arguments see <setKernel.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.k      : default: 250
% * seti.model  : default: 'helmholtz2D'
%
% _Are set automatically:_
%
% * seti.kHat   : size nCD x nCD
%
% For the following output arguments see <expSetup.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.incType        :   default: 'pointSource'
% * seti.measType       :   default: 'nearField'
% * seti.incNb          :   default: 35
% * seti.measNb         :   default: 35
% * seti.radSrc         :   default: 5
% * seti.radMeas        :   default: 5
% * seti.incPntsType    :   default: 'circle'
% * seti.measPntsType   :   default: 'circle'
%
% _Are set automatically:_
%
% * seti.incPnts        :   size seti.dim x seti.incNb
% * seti.dSInc          :   
% * seti.measPnts       :   size seti.dim x seti.measNb
% * seti.dSMeas         :   
% 
% * seti.incField       :   complex matrix of size 
%                           seti.nROI^seti.dim x seti.incNb
% * seti.measKer        :   complex matrix of size 
%                           seti.measNb x seti.nROI^seti.dim
%
% For the following output arguments see <setContrast.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.contrast   :   default: 'cornerBallSparse2D'
%
% _Are set automatically:_
%
% * seti.qCDexact   :   size: seti.nCD^seti.dim x 1
% * seti.qROIexact  :   size: seti.nROI^seti.dim x 1
%
%
% For the following output argument see <setGeomSim.html>:
%
% Note that seti.mCD unequal 0 is not supported in public version.
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.mCD    :   default: 0
% 
% For the following output arguments see <addNoise.html>:
% 
% _Can be set by user (otherwise default is set)_:
%
% * seti.delta      :   default: 0.01
% * seti.whichNoise :   default: 'normal'
% * seti.seed       :   default: 0
%
% The following output arguments are set in setData (this function):
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.loadFmeas      :   path to mat-file containing variables 
%                           FmeasExact and FmeasDelta.
%                           (default: empty, i.e. '').
% * seti.useFmeasDelta  :   Set it to 1 to use FmeasDelta in loadFmeas
%                           (default: 1).
%
% _Are set automatically:_
%
% * seti.FmeasExact     :   scattered field evaluated on receivers positions 
%                           (complex matrix of size seti.measNb x seti.incNb),
%                           see also <mimo.html>.
% * seti.FmeasDelta     :   FmeasExact with noise (synthetic data), see <addNoise.html>
%                           (complex matrix of size seti.measNb x seti.incNb).
%

%%%%%

%% setData: checkConsisExpData: Input Arguments
%
% * seti    : structure array
%
% * seti.expData    : This field is not required. If it is not set, the
% structure |seti| is not changed. Possible inputs: |'fresnel'|. (Note that
% |'simonetti'| is not supported in the public version of this ode.)
%
% The following input arguments can be set by user, otherwise default values are
% set (in case of |seti.expData = 'fresnel'|):
%
% * |seti.rCD|  : Size of computational domain [-rCD,rCD)^dim
%                 (default: 0.2 m), see <setGrid.html>.
% * |seti.nCD|  : Number of discretization points for each dimension of CD
%                 (default: 256), see <setGrid.html>.
% * |seti.fresnelFreq|  : Frequency of the Fresnel data set (default: 5 GHz
% = |5*1E9|)
% * |seti.fresnelFile|  : path to the file with data from Institute Fresnel
%                         (default: 'inexpdata/fresnel_opus_1/twodielTM_8f.exp').
%
% * |seti.nuMax|    : match incident field using Hankel functions of first kind and orders
%                     $\nu = -\texttt{nuMax}, ..., -1, 0, 1, ...,
%                     \texttt{nuMax}$ (default: 7)
% * |seti.ampCalc|  : method to compute the coefficients c (1, 2 or 3),
%                     (default: 1), see <matchIncField.html>.
%
%
%% setData:checkConsisExpData: Output Arguments
%
% * seti    : structure array
%
% In case of |seti.expData = 'fresnel'| the output arguments of |seti| are
%
% * |seti.dim = 2|  : The dimension of the problem is 2.
% * |seti.incType =  'pointSource'| : The type of incident field is set to point sources.
% * |seti.measType = 'nearField'|   : The measurement type is set to near field data.
%
% Additional fields of |seti| may set, if they was not set by the user (see
% Input Arguments).

%%%%%

%% setData:loadData: Input Arguments
%
% * |seti|      : structure array
%
% The input arguments are the input and output arguments of
% |checkConsisExpData| in case of |seti.expData = 'fresnel'|, i.e.:
%
% * |seti.rCD|  : Size of computational domain [-rCD,rCD)^dim
%                 (default: 0.2 m), see <setGrid.html>.
% * |seti.nCD|  : Number of discretization points for each dimension of CD
%                 (default: 256), see <setGrid.html>.
% * |seti.fresnelFreq|  : Frequency of the Fresnel data set
%                         (default: 5 GHz = |5*1E9|)
% * |seti.fresnelFile|  : path to the file with data from Institute Fresnel
%                         (default: 'inexpdata/fresnel_opus_1/twodielTM_8f.exp').
%
% * |seti.nuMax|    : match incident field using Hankel functions of first kind and orders
%                     $\nu = -\texttt{nuMax}, ..., -1, 0, 1, ...,
%                     \texttt{nuMax}$ (default: 7)
% * |seti.ampCalc|  : method to compute the coefficients c (1, 2 or 3),
%                     (default: 1), see <matchIncField.html>.
%
% * |seti.dim = 2|  : The dimension of the problem is 2.
% * |seti.incType =  'pointSource'| : The type of incident field is set to point sources.
% * |seti.measType = 'nearField'|   : The measurement type is set to near field data.
%
% *Optional Input Arguments*
%
% * dispDepth   : depth of displayed messages (number between 0 and 5).
% * out         : depth of output: no figure (0), plot figure (1), plot and save figure (2).
%
%% setData:loadData: Output Arguments
%
% Most important output arguments, i.e. fields in structure array |seti|.
%
% * |seti.incNb|      :   number of transmitters (number of incident fields)
% * |seti.measNb|     :   number of receivers (number of measurements)
% * |seti.radSrc|     :   radius of circle transmitters are arranged on
% * |seti.radMeas|    :   radius of circle receivers are arranged on
%
% * |seti.k|          :   wave number ($k = 2 \pi f/c$ with frequency f and light velocity c)
% * |seti.FmeasDelta| :   scattered field at receivers' positions (i.e. 'the data') for each transmitter
%                         (uScaRX = uTotRX - uIncRX, i.e. total field minus incident field)
%                         (complex matrix of size seti.measNb x seti.incNb)
% * |seti.incField|   :   Incident field in the region of interest (ROI)
%                         for each transmitter
%                         (complex matrix of size seti.nROI^seti.dim x seti.incNb)

%%%%%

%% setData:loadData:matchIncField: Input Arguments
% 
% * seti.dim        :   dimension of the problem: 2 (3 not available)
% * seti.incNb      :   number of transmitters
% * seti.k          :   wave number
% * seti.incPnts    :   Positions of the transmitters, 
%                       (matrix of size seti.dim x seti.incNb), 
%                       seti.incPnts = [5 -2 3; 0 4 2] 
%                       describes coordinates (5,0), (-2,4), and (3,2).
% * seti.measPnts   :   Positions of the receivers
%                       (matrix of size seti.dim x seti.measNb),
%                       coordinates analogical to seti.incPnts.
% * seti.gridROI    :   grid of region of interest (ROI) 
%                       (in case of region = 'ROI')
%                       (matrix of size seti.dim x seti.nROI^seti.dim)
% * seti.grid       :   grid of computational domain (CD)
%                       (in case of region = 'CD')
%                       (seti.dim x seti.nCD^seti.dim)
%
% * seti.ampCalc    :   method to compute the coefficients c (1, 2 or 3),
%                       we recommend to use method 1 (default)
%                       because it is the most accurate and fast, as shown in [1].
%                       (1: best-approximation to V c = z via MATLAB
%                       function,
%                       2: $c = (\gamma I + V^\ast\ V)^{-1} (V^\ast\
%                       \texttt{uIncRX})$
%                           via linear Tikhonov regularization,
%                       3: Landweber iteration)
% * seti.nuMax      :   match incident field using Hankel functions of first kind and orders
%                       $\nu = -\texttt{nuMax}, ..., -1, 0, 1, ...,
%                       \texttt{nuMax}$ (default: 7)
%
% * uIncRX          :   incident field at receivers' positions for each
%                       transmitter
%                       (complex matrix of size seti.measNb x seti.incNb)
% * region          :   'ROI' or 'CD' (region of interest or computational domain)
%                       (usually we are interested in ROI)
%
%
%% setData:loadData:matchIncField: Output Arguments
%
% No fields of |seti|.
%

%%%%%

%% setData:setGeomSim: Input Arguments
%
% * seti    :   structure array, see below in "Output Arguments" and also <setData.html>.
%
%% setData:setGeomSim: Output Arguments
%
% * seti        :   structure array
%
% Following parameters can be set by user (otherwise default is set):
%
% * seti.tol    :   tolerance for GMRES in <solveLippmannSchwinger.html>
%                   (default: 1E-6)
% * seti.mCD    :   coarse grid size (default: 0). 
%                   If mCD = 0, then two-grid is disabled.
%                   Note that coarse gris is not supported in public version.
%
% * Fields of seti, that are directly processed in setGeomSim, are
% described in this file.
% * Other fields are explained in the corresponding functions
% (<setGrid.html>, <setReshapeVecMat.html>,
% <setGridScale.html>, <setIdImagReal.html>, <setKernel.html>,
% <expSetup.html>, <setContrast.html>).
% An assignment of the fields to the functions is given in <setData.html>.
% 
% *Subfunction: setFigureSettings*
%
% Following parameters can be set by user (otherwise default is set):
%
% * seti.plotFreq       :   (default: 1) 
%                           0: do not plot; 
%                           n > 0: plot each n-th outer iteration.
% * seti.usecbarlim     :   (default: 1) 
%                           0: colorbar limits are set automatically; 
%                           1: set the limits manually (for real and imag colorbars in 2D)
% * seti.cbarlim        :   (default: [-0.2 1.4]) vector [cbarmin,cbarmax].
% * seti.plotPublish    :   (default: 0) 
%                           0 or 1, save figures designed for publications,
%                           i.e. no title, no axis, bigger font size for colorbar.,
%                           see <plot2DstylePublish.html> and <plot3DstylePublish.html>.
% * seti.pubFontSize    :   (default: 20) font size for publishing figures.
% * seti.plotVisible    :   (default: 'off') 
%                           figure visibility by MATLAB figure property
%                           "Visible", see [1].
% * seti.savepng        :   (default: 1) 0 or 1, save figures as *.png.
% * seti.saveepsc       :   (default: 0) 0 or 1, save figures as colored *.eps.
% * seti.savefig        :   (default: 0) 0 or 1, save figures as *.fig.
% * seti.plotFreqiPda   :   (default: 0) Frequency of plots inside pda 
%                           A folder is made for every outer iteration.
%                           (This option is not available in public version.)
%

%%%%%

%% setData:setGeomSim: setGrid: Input Arguments
%
% * seti    :   structure array
%
% If input arguments are not set, default values will be used.
%
% * seti.dim     :  dimension of the problem (2 or 3) (default: 2)
% * seti.nCD     :  number of discretization points for each dimension
%                   of computational domain (CD) (in samples)
%                   (default: 256).
% * seti.rCD     :  Size of computational domain [-rCD,rCD)^dim
%                   (default: 0.2) (in meters)
%
%
%% setData:setGeomSim: setGrid: Output Arguments
%
% * seti    :   structure array
%
% * seti.h          :   length of the infinitesimal element of CD.
% * seti.dV         :   area/volume of the infinitesimal element (pixel/voxel) of CD
% * seti.grid       :   grid of computational domain (CD)
%                       (seti.dim x seti.nCD^seti.dim)
% * seti.ballMask   :   mask (logical matrix of size seti.nCD x seti.nCD)
%                       to restrict (later) the contrast in CD 
%                       to the mathematical sensible region.
% * seti.nROI       :   discretization points for each dimension
%                       of region of interest (ROI) (in samples)
% * seti.gridROI    :   grid of region of interest (ROI) 
%                       (matrix of size seti.dim x seti.nROI^seti.dim)
% * seti.ROImask    :   mask (logical matrix of size seti.nCD x seti.nCD)
%                       to restrict (later) the contrast in CD
%                       to the region of interest (ROI).
% * seti.ballMaskROI :   mask (logical matrix of size seti.nROI x seti.nROI)
%                       to restrict (later) the contrast in ROI
%                       to the mathematical sensible region.
%                       (This mask is currently useless, see "More About".)
%

%%%%%

%% setData:setGeomSim:setReshapeVecMat: Input Arguments
%
% * seti.dim    : see <setGrid.html>.
% * seti.nCD    : see <setGrid.html>.
% * seti.nROI   : see <setGrid.html>.
%
%% setData:setGeomSim:setReshapeVecMat: Output Arguments
%
% * seti.GROI   : reshapes a vector of size seti.nROI^seti.dim 
%                 into a matrix of size seti.nROI x seti.nROI in 2D
%                 and seti.nROI x seti.nROI x seti.nROI in 3D.
% * seti.GCD    : reshapes a vector of size seti.nCD^seti.dim 
%                 into a matrix of size seti.nCD x seti.nCD in 2D
%                 and seti.nCD x seti.nCD x seti.nCD in 3D.
% * seti.G      : same as seti.GROI.
% * seti.iG     : reshapes a matrix into a vector.

%%%%%

%% setData:setGeomSim:setGridScale: Input Arguments
%
% If the fields was not defined, default values are set.
%
% * seti.gscale     :   1: use a down scaled grid for reconstruction;
%                       0: do not use grid scaling (default)
% * seti.nCDinv     :   number of discretization points
%                       _of the downscaled grid_
%                       for each dimension of CD.
%                       (for nCD see <setGrid.html>).
%                       (default: |floor(seti.nCD/2)|).
%
%
%% setData:setGeomSim:setGridScale: Output Arguments
%
% * seti.nInv   :   The analogue to seti.nROI, see <setGrid.html>, for the downscaled grid.
% * seti.hInv   :   The analogue to seti.h, see <setGrid.html>, for the downscaled grid.
% * seti.dVinv  :   The analogoue to seti.dV, see <setGrid.html>, for the downscaled grid.
% * seti.GInv   :   The analogue to seti.G, see <setReshapeVecMat.html>, for the downscaled grid.
%
% * seti.GU     :   function to scale up the grid
%                   (input: down scaled contrast as a vector of size seti.nInv^seti.dim;
%                   output: up scaled contrast as a vector of size seti.nROI^seti.dim)
% * seti.GD     :   function to scale down the grid
%                   (input: up scaled contrast as a vector of size seti.nROI^seti.dim;
%                   output: down scaled contrast as a vector of size seti.nInv^seti.dim)

%%%%%

%% setData:setGeomSim:setIdImagReal: Input Arguments
%
% * seti    :   structure array (no fields are needed)
%
%% setData:setGeomSim:setIdImagReal: Output Arguments
%
% * seti.S  :   function to identify complex with real vectors/matrices
%               (input: complex vector/matrix; output: real vector/matrix)
% * seti.R  :   extracts the originally real (R) part of the output of seti.S
% * seti.I  :   extracts the originally imag (I) part of the output of seti.S
% * seti.T  :   function to identify real with complex vectors/matrices
%               (input: real vector/matrix (R x R); output: complex vector/matrix)

%%%%%

%% setData:setGeomSim:setKernel: Input Arguments
%
% _This arguments was set in rebis before setKernel is called_:
%
% * seti.dim    : see <setGrid.html>.
% * seti.nCD    : see <setGrid.html>.
% * seti.rCD    : see <setGrid.html>.
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.k      : Wave number, default: 250
% * seti.model  : Model of the problem, default: 'helmholtz2D'
%
%% setData:setGeomSim:setKernel: Output Arguments
%
% * seti.kHat   : kernel, complex matrix of size nCD x nCD.
%

%%%%%

%% setData:setGeomSim:expSetup: Input Arguments
%
%% Input Arguments
%
% * |seti|  :   structure array
%
% Several fields in |seti| are required. Because this is an internal
% function we do not list them.
%
% *Optional Input Arguments to differ from default values*
%
% For details of the the following fields look inside subfunction 
% |expSetupCons| below and the functions  <setIncPnts.html>, 
% <setMeasPnts.html>, <pntsGeometry.html> and <pntsGeometry3D.html>:
%
% * seti.incType        :   type of incident field, default: 'pointSource'
%                           ('planeWave' or 'pointSource')
% * seti.measType       :   type of measurement, default: 'nearField'
%                           ('nearField' or 'farField')
% * seti.incNb          :   number of transmitters, default: 35
% * seti.measNb         :   number of receivers, default: 35
% * seti.radSrc         :   radius of circle for transmitters, default: 5
%                           (not necessary in case of incType='planeWave').
% * seti.radMeas        :   radius of circle for receivers, default: 5
%                           (not necessary in case of measType='farField')
% * seti.incPntsType    :   string with type of geometry for transmitters, default: 'circle'
% * seti.measPntsType   :   string with type of geometry for receivers, default: 'circle'
%
%% setData:setGeomSim:expSetup: Output Arguments
%
% * |seti|  :   structure array
%
% The following fields in |seti| are defined. Look in documentation of 
% corresponding functions and subfunctions (subfunctions are in the Section
% Code).
%
% For details of the the following fields see the functions <setIncPnts.html>, 
% <setMeasPnts.html>, <pntsGeometry.html>:
%
% * seti.incPnts        :   coordinates of transmitters, size seti.dim x seti.incNb
% * seti.dSInc          :   Approximation of the infinitesimal element of closed contour with control points.
% * seti.measPnts       :   coordinates of receivers, size seti.dim x seti.measNb
% * seti.dSMeas         :   Approximation of the infinitesimal element of closed contour with control points.
% 
% Further details of the following parameter are in <setIncField.html>:
%
% * seti.incField       :   incident fields, complex matrix of size 
%                           seti.nROI^seti.dim x seti.incNb
%
% Further details of the following parameter are in <seMeasKer.html>:
%
% * seti.measKer        :   measurement kernel on region of interest (ROI) 
%                           for each receiver 
%                           (complex matrix of size seti.measNb x seti.nROI^seti.dim)
%

%%%%%

%% setData:setGeomSim:expSetup:setIncPnts: Input Arguments
%
% * seti.dim            :   dimension of the problem (2 or 3)
% * seti.incType        :   Type of incident field: 'pointSource' (default) or 'planeWave'
%
% *Input Arguments in case of 2D (seti.dim = 2) and point sources (seti.incType = 'pointSource')*
%
% * seti.incPntsType    :   string with type of geometry:
%                           'manually', 'circle' (default), 'square', 'line', 'borehole'
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry.html> (mostly with the prefix inc and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.incPntsType = 'manually'_
%
% * seti.incPnts        :   coordinates of points (real matrix of size 2 x Nb).
%
% _varargin in case of seti.incPntsType = 'circle'_
%
% * seti.incNb             :   Number of points (transmitters)
% * _seti.radSrc_          :   radius of circle
%
% _varargin in case of seti.incPntsType = 'square'_
%
% * seti.incNbEdge      :   incNbEdge+1 points are on one edge
%                           (+1 because a point in a corner belongs to two corners).
%                           (Total number of points in square is |4*incNbEdge|.)
% * seti.incEdgeLength  :   Length of one edge.
%
% _varargin in case of seti.incPntsType = 'line'_
%
% * seti.incNbLine      :   Number of points on line.
% * seti.incLineLength  :   Length of line.
% * seti.incLinePos     :   Distance of line right to origin.
%
% _varargin in case of seti.incPntsType = 'borehole'_
%
% * _seti.boreNbLine_      :   Number of points on line.
% * _seti.boreLineLength_  :   Length of line.
% * _seti.boreLinePos_     :   Distance of line right to origin.
%
% *Input Arguments in case of 3D (seti.dim = 3)*
%
% * seti.incPntsType    :   string with type of geometry:
%                           'manually', 'sphereLatLon', 'sphereFibo' (default).
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry3D.html> (mostly with the prefix inc and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry3D.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.incPntsType = 'manually'_
%
% * seti.incPnts        :   coordinates of points (real matrix of size 3 x Nb).
%
% _varargin in case of seti.incPntsType = 'sphereLatLon' or 'sphereFibo'_
%
% * seti.incNb             :   Number of points (transmitters)
% * _seti.radSrc_          :   radius of sphere
%
%% setData:setGeomSim:expSetup:setIncPnts: Output Arguments
%
% * seti.incPnts    :   Coordinates of points (real matrix of size seti.dim x Nb).
% * seti.incNb      :   Number of points.
% * seti.dSInc      :   Approximation of the infinitesimal element of closed 
%                       contour with control points.
%                       A closed contour does not make sense in the cases
%                       'manually', 'line', and 'borehole', then it is set to 1.
%
% For further information see the corresponding parameters 
% Pnts, Nb, and dS in <pntsGeometry.html> and <pntsGeometry3D.html>.
%

%%%%%

%% setData:setGeomSim:setMeasPnts: Input Arguments
%
% * seti.dim            :   dimension of the problem (2 or 3)
% * seti.measType       :   Type of measurement field: 'nearField' (default) or 'farField'
%
% *Input Arguments in case of 2D (seti.dim = 2) and near field (seti.incType = 'nearField')*
%
% * seti.measPntsType   :   string with type of geometry:
%                           'manually', 'circle' (default), 'square', 'line', 'borehole'
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry.html> (mostly with the prefix meas and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.measPntsType = 'manually'_
%
% * seti.measPnts        :   coordinates of points (real matrix of size 2 x Nb).
%
% _varargin in case of seti.measPntsType = 'circle'_
%
% * seti.measNb             :   Number of points (transmitters)
% * _seti.radMeas_          :   radius of circle
%
% _varargin in case of seti.measPntsType = 'square'_
%
% * seti.measNbEdge      :  measNbEdge+1 points are on one edge
%                           (+1 because a point in a corner belongs to two corners).
%                           (Total number of points in square is |4*measNbEdge|.)
% * seti.measEdgeLength  :  Length of one edge.
%
% _varargin in case of seti.measPntsType = 'line'_
%
% * seti.measNbLine      :   Number of points on line.
% * seti.measLineLength  :   Length of line.
% * seti.measLinePos     :   Distance of line right to origin.
%
% _varargin in case of seti.measPntsType = 'borehole'_
%
% Note: Option 'borehole' is only useful if seti.incPntsType is also 'borhole'.
%
% * _seti.boreNbLine_      :   Number of points on line.
% * _seti.boreLineLength_  :   Length of line.
% * _seti.boreLinePos_     :   Distance of line right to origin.
%
% *Input Arguments in case of 3D (seti.dim = 3)*
%
% * seti.incPntsType    :   string with type of geometry:
%                           'manually', 'sphereLatLon', 'sphereFibo' (default).
%
% It follows a *list of type-depending input parameters* to describe the
% details of the geometry. The names are analog to the parameters in 
% <pntsGeometry3D.html> (mostly with the prefix meas and a big letter,
% exceptions are marked with emphasized).
% Examples are in <pntsGeometry3D.html>.
% If these parameters are not set, default values are set automatically
% (see inside the code).
%
% _varargin in case of seti.measPntsType = 'manually'_
%
% * seti.measPnts        :   coordinates of points (real matrix of size 3 x Nb).
%
% _varargin in case of seti.measPntsType = 'sphereLatLon' or 'sphereFibo'_
%
% * seti.measNb             :   Number of points (receivers)
% * _seti.radMeas_          :   radius of sphere
%               
%% setData:setGeomSim:setMeasPnts: Output Arguments
%
% * seti.measPnts    :   Coordinates of points (real matrix of size seti.dim x Nb).
% * seti.measNb      :   Number of points.
% * seti.dSMeas      :   Approximation of the infinitesimal element of closed 
%                        contour with control points.
%                        A closed contour does not make sense in the cases
%                        'manually', 'line', and 'borehole', then it is set to 1.
%
% For further information see the corresponding parameters 
% Pnts, Nb, and dS in <pntsGeometry.html> and <pntsGeometry3D.html>.
%
%

%%%%%

%% setData:setGeomSim:setIncField: Input Arguments
%
% * seti.incType    :   Type of incident field: 'pointSource' or 'planeWave'
% * seti.incNb      :   number of transmitters
%
% * seti.k      :   wave number, see <setKernel.html>
% * seti.model  :   Model of problem 'helmholtz2D' or 'helmholtz3D',
%                   see <setKernel.html>
%
% * seti.nROI    : discretization points for each dimension
%                  of region of interest (ROI) (in samples),
%                  see <setGrid.html>.
% * seti.gridROI : grid of region of interest (ROI) 
%                  (matrix of size seti.dim x seti.nROI^seti.dim)
%                  see <setGrid.html>.
%
% * seti.incPnts    : coordinates of transmitters
%                     (real matrix of size 2 x seti.incNb), see <setIncPnts.html>.
%
%
%% setData:setGeomSim:setIncField: Output Arguments
%
% * |seti.incField|   :   incident field on region of interest (ROI)
%                         for each transmitter
%                         (complex matrix of size seti.nROI^seti.dim x seti.incNb)
%

%%%%%

%% setData:setGeomSim:setMeasKer: Input Arguments
%
% * seti.measType    :   Type of scattered field: 'nearField' or 'farField'
% * seti.measNb      :   number of receivers
%
% * seti.k      :   wave number, see <setKernel.html>
% * seti.model  :   Model of problem 'helmholtz2D' or 'helmholtz3D',
%                   see <setKernel.html>
%
% * seti.nROI    : discretization points for each dimension
%                  of region of interest (ROI) (in samples),
%                  see <setGrid.html>.
% * seti.gridROI : grid of region of interest (ROI) 
%                  (matrix of size seti.dim x seti.nROI^seti.dim)
%                  see <setGrid.html>.
%
% * seti.measPnts    : coordinates of receivers
%                      (real matrix of size seti.dim x seti.measNb), see <setMeasPnts.html>.
%
%
%% setData:setGeomSim:setMeasKer: Output Arguments
%
% * |seti.measKer|      :   measurement kernel on region of interest (ROI)
%                           for each receiver
%                           (complex matrix of size seti.measNb x seti.nROI^seti.dim)
%

%%%%%

%% setData:setGeomSim:setContrast: Input Arguments
%
% * seti    :   structure array
% 
% * |seti.contrast|   :   name of predefined contrast function
% * |seti.rotation|   :   mathematical positive rotation of the contrast in degrees (only in 2D).
%                         (If the field rotation does not exist, no rotation
%                         is done.)
%
%% setData:setGeomSim:setContrast: Output Arguments
%
% * |seti.qCDexact|     :   predefined contrast evalutated on CD
%                           (computational domain) 
%                           (as vector of size seti.nCD^seti.dim x 1).
% * |seti.qROIexact|    :   predefined contrast evalutated on ROI
%                           (region of interest) 
%                           (as vector of size seti.nROI^seti.dim x 1).

%%%%%

%% setData:mimo: Input Arguments
%
% * |seti|       :    struct, details see below.
%
% *struct |seti|*
%
% As in convenience functions 
% <adjOfDer.html>, <forward.html> and <derivative.html>
% all necessary fields (and more) are set by 
%
%   seti = setGeomSim(seti);
%
% A list of all necessary (and not more) fields is below.
%
% seti in <setGeomSim.html>:
%
% * |seti.mCD|      : internal parameter, set seti.mCD = 0 
%                     (no coarse grid is used)
% * |seti.tol|      : internal parameter, set 1E-6, is used in simo and solveLippmannSchwinger
%
% seti in <setGrid.html>:
%
% * |seti.dim|     : dimension of the problem (2 or 3)
% * |seti.nCD|     : number of discretization points for each dimension
%                    of computational domain (CD)
% * |seti.nROI|    : discretization points for each dimension
%                    of region of interest (ROI)
% * |seti.gridROI| : grid of region of interest (ROI) 
%                    (matrix of size seti.dim x seti.nROI^seti.dim)
% * |seti.ROImask| : biggest square in 2D (cube in 3D) in ball with radius rCD/2
%                    (logical array of size seti.nCD x seti.nCD)
% * |seti.dV|      : pixel/voxel area/volume (number)
% 
% seti in <setKernel.html>:
%
% * |seti.model|   : model (helmholtz2D or helmholtz3D)
% * |seti.k|       : wave number
% * |seti.kHat|    : Fourier coefficients of discretized integral operator
%                    (seti.nCD x seti.nCD)
%
% seti in <expSetup.html> (subfiles <pntsGeometry.html>,
% <setIncField.html>, <setMeasKer.html>):
%
% * |seti.measNb|   : number of receivers (number of measurements)
% * |seti.incNb|    : number of transmitters (number of incident fields)
% * |seti.dSInc|    : approximation of the infinitesimal element of closed contour with control points in case of transmitters
% * |seti.dSMeas|   : approximation of the infinitesimal element of closed contour with control points in case of receivers
% * |seti.incField| : incident field, <setIncField.html>
%                     (matrix of size seti.nROI^seti.dim x seti.incNb)
% * |seti.measKer|  : matrix to compute simulated measurement data from solution, <setMeasKer.html>
%                     (size seti.measNb x seti.nROI^seti.dim)
%
%
%% setData:mimo: Output Arguments
%
% No output fields of |seti|.

%%%%%

%% setData:addNoise: Input Arguments
%
% * seti        : structure array
%
% _The following field of the structure array is required:_
%
% * seti.dSMeas         : Approximation of the infinitesimal element of closed 
%                         contour with control points, see <setMeasPnts.html>.
%
% _The following fields of the structure array are optional (otherwise
% default values are set):_
%
% * seti.delta          : Relative (artificial) noise level of data 
%                         (default: 0.01)
% * seti.whichNoise     : Type of used probability density function to 
%                         noise the data:
%                         'laplace', 'uniform', 'normal' (default).
% * seti.seed           : Number to control the random number generator.
%                         Using a non-negative integer (default: 0).
%
%% setData:addNoise: Output Arguments
%
% * seti    :   structure array
%
% If the fields |delta|, |whichNoise|, and |seed| of |seti| was not set,
% default values are set.
%

%%%%%

%% setRecon: Input Arguments
%
% * seti    :   structure array
%
%% setRecon: Output Arguments
%
% * seti    :   structure array
%
% The added fields in |seti| can be found in the called functions
%
% * <setInvType.html>
% * <checkConsisRec.html>
% * <shrinkFuncs.html>
% * <setFuncsPda.html>

%%%%%

%% setRecon:setInvType: Input Arguments
%
% * |seti|        :   struct seti
%
% *Optional Input Arguments of structure |seti|*
%
% If the fields does not exist, default values are set in this function.
%
% * |seti.invNo|  : Number of reconstruction setting (default: 6)
%                   (invNo = 6 is the only one supported in public version)
% * |seti.alpha|  : Regularization parameter for sparsity penalty 
%                   (default: 500).
% * |seti.beta|   : Regularization parameter for total variation penalty
%                   (default: 1E-5).
% * |seti.tau|    : Discrepancy principle stops at parameter tau * delta
%                   (default: 2.5). (Note that delta is the relative noise
%                   level |seti.delta|.)
%
%% setRecon:setInvType: Output Arguments
%
% * |seti.inv|        : Type of reconstruction method: 'pda' or 'shrinkage'
% * |seti.useWavelet| : Using wavelets (1) or not (0)
%                       (Wavelets are not available in public version)
% * |seti.tF|     : Terms belonging to functional F
%                   (e.g. |seti.tF = {'fd','fg'}|; See table below.)
% * |seti.tG|     : Terms belonging to functional G
%                   (e.g. |seti.tG = {'fs','fp'}|; See table below.)
% * |seti.pNorm|  : Index p for p-Schatten-Norm of measurements to define
%                   discrepancy.
% * _More about pNorm:_   default: pNorm = 2, i.e. essentially Frobenius-Norm
% * _More about pNorm:_   in case of inv = 'pda' with invNo 3 to 6
%                         set pNorm = 2 (automatically) because it is required
% * |seti.qNorm|      : Index q for q-Norm of contrasts to define sparsity penalty
% * _More about qNorm_: Standard: qNorm = 1 -> l^1 minimization.
% * _More about qNorm_: in case of inv = 'pda' set qNorm = 1 (automatically)
%                       because it is required.
% * |seti.p|      : Exponent in discrepancy-term $1/p ||...||^p_2$ with $p = 1$ or $2$
%                   (used e.g. in <setFuncsPda.html>)
%                   (default: |seti.p = 2|).
%
% Note: *p ~= pNorm* in general!
%                   (but: in case of inv = 'shrinkage':
%                   pNorm and qNorm are used as exponents too.)
%

%%%%%

%% setRecon: checkConsisRec: Input Arguments
%
% * |seti|        :   struct seti
%
% *Optional Input Arguments of structure |seti|*
%
% If the fields does not exist, default values are set in this function.
%
% * |seti.useWavelet|   : Use wavelets for reconstruction (default: 0).
%                         Note that wavelets are not supported in the 
%                         public version of this code. Therefore the
%                         wavelet-specific fields |wavWeight|,
%                         |genWhiteNoiseNew|, |adaptWeightsMan|, |wavelet|,
%                         |extension|, |smin| and |samples| are not
%                         documented.
% * |seti.recname|      : What is reconstructed? A |'constrast'| (default) 
%                         or a |'source'|. Influences the title in
%                         <subplots.html>.
% * |seti.physBounds|   : Bounds for real/imaginary part of contrast: 
%                         [reMin, reMax, imMin, imMax] 
%                         (default: [-1,3,0,3])
%                         (consider also the penalty for physical bounds in
%                         <pda.html> Section "More About".)
% * |seti.tolDelta|     : Tolerance for physical bounds' delta-function
%                         (default: 1E-14) in <setFuncsPda.html>.
% * |seti.useDis|       : If it is 1 (default), the discrepancy principle
%                         is used to stop the outer iteration. (Otherwise 0).
% * |seti.nOut|         : Maximal number of outer iterations (default: 30).
%
% If the primal-dual algorithm is used (|seti.inv = 'pda')
% and the following fields does not exist, default values are set in this function.
%
% * |seti.pdaN|         : Number of inner iterations (PDA) (default: 50) (see <pda.html>).
% * |seti.pdaStepsize|  : Method to choose primal and dual stepsizes
%                         (|'fix'| (default) or |'adaptive'|)
%                         see <pdaChoosingStepsizes.html>.
% * |seti.vartheta|     : Parameter in (0,1) used in case of adaptive stepsizes,
%                         see <pdaChoosingStepsizes.html>.
% * |seti.useTolIn|     : If it is set to 1 (default: 0), the inner 
%                         tolerance principle is used to stop the inner 
%                         iterations (primal-dual algorithm), see 
%                         <minPda.html> and <minTolIn.html>.
% * |seti.useTolOut|    : If it is set to 1 (default: 0), the outer 
%                         tolerance principle is used to stop the inner 
%                         iterations, see <minPda.html>, <minTolOut.html>.
%
% In case of |seti.inv = 'pda'| and |seti.useTolIn = 1| and the following
% fields does not exist, default values are set in this function.
%
% * |seti.ThetaStart|   : (default: 0.925), see <minTolIn.html>.
% * |seti.ThetaMax|     : (default: 0.95), see <minTolIn.html>.
% * |seti.TolGamma|     : (default: 0.90), see <minTolIn.html>.
% * |seti.pdaNmax|      : Maximal iteration number, if |seti.useTolIn| or
%                         |seti.useTolOut| are 1 (default : 250).
%
% In case of |seti.inv = 'pda'| and |seti.useTolOut = 1| and the following
% fields does not exist, default values are set in this function.
%
% * |seti.relDisTol|    : Outer tolerance (default: 0.05), see <minTolOut.html>.
% * |seti.pdaNmax|      : Maximal iteration number, if |seti.useTolIn| or
%                         |seti.useTolOut| are 1 (default : 250).
%
% The following fields should have been set in <setInvType.html>:
%
% * |seti.inv|
% * |seti.alpha|
% * |seti.beta|
% * |seti.pNorm|
% * |seti.qNorm|
% * |seti.tau|
%
% See <recon.html> for the following fields of |seti|:
%
% * |seti.saveqROIcomp|
% * |seti.loadqROIcomp|
% * |seti.savedata|
%
% Currently unused fields:
%
% * |seti.lambda|   : Wavelength (default: $\lambda = 2\pi/k$ with 
%                     wave number $k = \texttt{seti.k}$) (currently not used anywhere).
%
%% setRecon:checkConsisRec: Output Argument
%
% * |seti|  : structure array, see "Optional Input Arguments of structure
% |seti|", because default values are set, if fields was not already
% defined.
%

%%%%%

%% setRecon:shrinkFuncs: Input Arguments
%
% * seti    :   structure array
%
% * seti.qNorm  :   in sparsity term $f_\mathrm{spa}$ a norm with index 
%                   |seti.qNorm = 1| is used 
%                   (other values are not supported in public version).
%                   (Index q = 1 does essentially (up to some factor) mean $l^1$ minimization)
% * seti.physBounds     :   contains the physical bounds (min and max)
%                           for real and imaginary part of the contrast as
%                           4-vector with structure 
%                           |[reMin reMax imMin imMax]|.
%
%% setRecon:shrinkFuncs: Output Arguments
%
% * seti.shrkRe = @(x,alpha)    :   extended soft-shrinkage function for real part
% * seti.shrkIm = @(x,alpha)    :   extended soft-shrinkage function for imaginary part

%%%%%

%% setRecon:setFuncsPda: Input Arguments
%
% * seti    :   structure array
%
%% setRecon:setFuncsPda: Output Arguments
%
% * seti    :   structure array
%
% For the fields of the structrual array |seti| see the section "More About".
%
% *--- More About ---*
%
% We only give a overview of the functions. See also <pda.html>.
% Detailed information are in Section 4.5 in [1].
%
% Note that grid scaling is respected in the code, but for 
% simplicitiy is not presented in the formulas.
%
% *Auxiliary quantities*
%
% For $\texttt{seti.S} = T_{\bf{C}\to\bf{R}^2}$ and $\texttt{seti.T} = T_{\bf{R}^2\to\bf{C}}$, see "More About" in <pda.html>.
%
% * $\texttt{seti.Kd} = K_\mathrm{dis} :=
%   T_{\bf{C}\to\bf{R}^2} [\mathcal{F}'(T_{\bf{R}^2\to\bf{C}}\ q)]T_{\bf{R}^2\to\bf{C}}$,
% * $\texttt{seti.vd} = v_\mathrm{dis} := 
%   T_{\bf{C}\to\bf{R}^2}(\mathcal{F}(T_{\bf{R}^2\to\bf{C}}\ q)-F_\mathrm{meas}^\delta)$,
% * $\texttt{seti.Kg} = K_\mathrm{tv}  := 
%   \beta\,T_{\bf{C}\to\bf{R}^2}\nabla T_{\bf{R}^2\to\bf{C}}$
% * $\texttt{seti.vg} = v_\mathrm{tv} :=
%   \beta\, T_{\bf{C}\to\bf{R}^2} \nabla T_{\bf{R}^2\to\bf{C}}\ q$.
%
% Note that the names |Kg| and |vg| were choosen because the appearance of the gradient.
%
% *Parts of functional to minimize*
%
% * $\texttt{seti.fd} = f_\mathrm{dis}(h) := 
%   \frac{1}{2} \|K_\mathrm{dis} h+  v_\mathrm{dis} \|_{\mathrm{dis},\bf{R}}^2$
%   $\quad$ is the discrepancy (linearized problem),
% * $\texttt{seti.fs} = f_\mathrm{spa}(h) := 
%   \alpha \| h + q\|_{\mathrm{spa}, \bf{R}}$
%   $\quad$ is the sparsity penalty,
% * $\texttt{seti.fg} = f_\mathrm{tv}(h) :=
%   \|K_\mathrm{tv} h + v_\mathrm{tv} \|_{\mathrm{tv},\bf{R}}$
%   $\quad$ is the TV-penalty,
% * $\texttt{seti.fp} = f_\mathrm{phy}(h) :=
%   \delta_{[a,b,c,d]}(h + q)$
%   $\quad$ is the penalty for physical bounds.
%
% Note that $\delta_{[a,b,c,d]}(x)$ is a indicator function, i.e. is 0, 
% if all entries of the vector $\mathrm{real}(x)$ are between $a$ and $b$
% and all entries of the vector $\mathrm{imag}(x)$ are between $c$ and $d$;
% and is $\infty$, otherwise.
%
% *Definition of parts of Tikhonov functional*
%
% The Tikhonov functional is |MT = M1 + M2|. 
% Function $F$ is $M_1$, function $G$ is $M_2$ (with $F$ and $G$ from PDA)
% This is done in function |minPda| (<minPda.html>) when function 
% |pda| (<pda.html>) is called.
%
% *Further comments (especially interesting in not public version)*
%
% * In case of inversion by pda, it is always (automatically) set 
%   |seti.pNorm = 2| and |seti.qNorm = 1| (except for invNo = 7), see
%   <setInvType.html>.
% * Exponent p of the discrepancy term, see <setInvType.html>, is usually
%   |seti.p = 2|.
%   In case of shrinkage or pda with wavelets it may differ, see
%   <setInvType.html>.
% * |vs| and |vp| in term $G$: If you choose |vs| and |vp| manually and 
%   |fs| and |fp| belong to term $G$, then |vs = vp| is required.
%

%%%%%

%% recon: Input Arguments
%
% * |seti|    :   structure array
%
% Because this is a internal function we do not explain all fields in seti.
%
% * |seti.loadqROIcomp| : Path to load a mat-file containing the reconstructed 
%                         contrast |qROIcomp| after |iOutStop| outer iterations 
%                         (default: empty, i.\,e. |''| or |[]|, then no data is loaded).
% * |seti.saveqROIcomp| : Optional: Path to save reconstructed contrast
%                         (default: same folder as figures, i.e. 
%                         |sprintf('%s/save_qROIcomp_iOutStop%s.mat',seti.dirname,seti.fileSuffix)|)
%
% * |seti.nOut|         : Maximal number of outer iterations.
%
% * |seti.useDis|       : Use discrepancy principle to stop outer iteration?
%                         (If yes set it to 1.)
% * |seti.tau|          : Tolerance parameter for discrepancy principle
%                         (seti.tau > 1).
%
% * |seti.plotFreq|     : Frequency to plot figures of outer iteration steps
%                         (0: no figures).
% * |seti.savedata|     : If it is set to 1,
%                         relative discrepancies, errors, and differences 
%                         of iterated contrasts are saved
%                         (as |save_dis.mat|, |save_err.mat| and |save_dif.mat|).
%
% * |seti.qROIexact|    : Predefined (exact) contrast
%                         (complex vector of |size seti.nROI|^|seti.dim| x 1).
%
% *Further specific input parameters* for primal-dual algorithm can be found in
% <minPda.html> and <pda.html>.
%
% *Optional Input Argument*
%
% * |dispDepth|     : depth of displayed text.
%                     (from 0 to 5, default sparse: 2; default details: 4; too much information: 5);
% * |out|           : output depth for figures and files: no figure (0), plot figure (1), 
%                     save figures and files (2).
%                     (This argument requires that argument |dispDepth| was
%                     defined.)
%
%% recon: Output Arguments
%
% * |seti|    :   structure array
%
% Because this is a internal function we do not explain all fields in seti.
%
% * |seti.qROIcomp|     : Reconstructed contrast
%                         (complex vector of size |seti.nROI|^|seti.dim| x 1).
% * |seti.iOutStop|     : Stop index of outer iterations.
%
%
% *Discrepancy, error, and difference*
%
% * |seti.dis|        : Relative discrepancy for each outer iteration
%                       (vector of size 1 x seti.nOut).
% * |seti.err|        : Relative error for each outer iteration
%                       (vector of size 1 x seti.nOut).
% * |seti.dif|        : Relative difference of iterated contrast qROI to previously qROI
%                       for each outer iteration
%                       (vector of size 1 x seti.nOut).
%
% Note that |seti.dis(iOut)| is the relative discrepancy 
% _after_ iteration |iOut|. (Analog seti.err(iOut) and seti.dif(iOut).)
%
% *Minimized Tikhonov functional*
%
% * |seti.MTv|        :   Result of minimized Tikhonov functional 
%                         for each outer iteration 
%                         (vector of size 1 x seti.nOut).
% * |seti.M1v|        :   Result of first part of min. functional
%                         for each outer iteration 
%                         (vector of size 1 x seti.nOut).
% * |seti.M2v|        :   Result of second part of min. functional 
%                         for each outer iteration 
%                         (vector of size 1 x seti.nOut).
%
% _More about M1v and M2v_
%
% * In case of primal-dual algorithm (i.e. |seti.inv = 'pda'|)
%   |M1v| is F and |M2| is G, see <pda.html> for F and G.
% * In case of shrinkage (i.e. |seti.inv = 'shrinkage'|)
%   (not available in public version) 
%   |M1| = discrepancy and |M2| = sparsity (setFuncsShrink.m)
%
% *structure array pdas*
%
% * |pdas|        : Structure array
%                   (only in case of pda, i.e. |seti.inv = 'pda'|).
%
% * |pdas.disLin| : Last relative discrepancy of linearized problem in inner iteration of pda
%                   for each outer iteration 
%                   (vector of size 1 x seti.nOut).
% * |pdas.relDis| : Quotient disLin/dis 
%                   for each outer iteration
%                   (rel. discrepancy of linearized problem / rel. discrepancy of the non-linearized problem)
%                   (vector of size 1 x seti.nOut).
% * |pdas.pdaNv|  : Number of inner iterations for each outer iteration
%                   (vector of size 1 x seti.nOut).
% * |pdas.ThetaiOutV| : Inner tolerance for each outer iteration
%                       (in case of inner tolerance principle, see
%                       <minTolIn.html>).
%
% *Further output arguments*
%
% * |seti.disIni|   : Relative discrepancy before 1st iteration.
% * |seti.errIni|   : Relative error before 1st iteration.
%

%%%%%

%% recon: subplots: Input Arguments
%
% * seti      : structural array
% * qROIexact : predefined (i.e. exact) contrast
%               (usually qROIexact = seti.qROIexact, but in other cases 
%               (e.g. 'source') you can put in something else to plot.)
% * qROIcomp  : reconstructed (i.e. computed) contrast
% * iOut      : number of outer iteration
%
% * seti.recname : reconstruction name, e.g. 'contrast' or 'source'.
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
%% recon: subplots: Output Arguments
%
% Output is figure 11 as described in description and in section 
% "Output: Figures" in <start.html>.

%%%%%

%% recon: plotAndSaveFigures: Input Arguments
%
% * seti      : structural array
% * qROIexact : predefined (i.e. exact) contrast
% * qROIcomp  : reconstructed (i.e. computed) contrast
% * iOut      : number of outer iteration
% * out       : Output depth: generate no plots (0),
%               generate plots (1), generate plots and save them (2).
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
% The most important fields for this function in |seti| are:
%
% * seti.dis : relative discrepancies
% * seti.err : relative errors
% * seti.MTv : result of Tikhonov funtional
% * seti.M1v : result of first part of Tikhonov functional
% * seti.M2v : result of second part of Tikhonov functional
%
% Further details of this fields are described in <recon.html>.
%
%% recon: plotAndSaveFigures: Output Arguments
%
% Plots, see <start.html>.
%
% Note that figure 11 is plotted in <subplots.html>, but saved in this
% file.
%

%%%%%

%% recon: pdaPlot: Input Arguments
%
% * iOut    : number of outer iteration
% * seti    : structural array
% * pdas    : structural array with specific results from primal-dual
%             algorithm, see <recon.html> for details.
% * out     : Output depth: generate plots (1), generate plots and save them (2).
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
%% recon: pdaPlot: Output Arguments
%
% Figures are plotted and saved, see <start.html>

%%%%%

%% recon: minPda: Input Arguments
%
% Specific input arguments for primal-dual algorithm are described in
% <pda.html>.
%
% * seti            :   structure array
% * seti.pdaNmax    :   maximal iteration number of inner iteration,
%                       if useTolIn or useTolOut (default : 250)
%
% * pdas        :   structure array for pda specific outputs, 
%                   see below and also <recon.html>.
%
% *structure array pdas*
%
% * pdas.pdaStopInd     :   stop index of inner iteration (i.e. primal-dual algorithm)
% * pdas.MTvN           :   result of minimized Tikhonov functional 
% * pdas.M1vN           :   result of first part of min. functional
% * pdas.M2vN           :   result of second part of min. functional
% * See also <recon.html>.
%
% * pdas.relLinDisInPda :   quotient: relLinDis = disLin / dis (vector of size seti.pdaN x 1)
% * pdas.disLinInPda    :   discrepancy of linearized problem for each inner iteration step (vector of size seti.pdaN x 1)
% * pdas.errInPda       :   relative error of the reconstructed contrast qROI (vector of size seti.pdaN x 1)
% * pdas.minf = minf    :   struct with parts of the minimization functional (Tikhonov functional)
% * See also <pda.html>.
%
% * pdas.disLin :   last relative discrepancy of linearized problem in inner iteration of pda
%                   for each outer iteration 
%                   (vector of size 1 x seti.nOut)
% * pdas.relDis :   quotient disLin/dis 
%                   for each outer iteration
%                   (rel. discrepancy of linearized problem / rel. discrepancy of the non-linearized problem)
%                   (vector of size 1 x seti.nOut)
% * pdas.pdaNv  :   number of inner iterations for each outer iteration
%                   (vector of size 1 x seti.nOut)
% * pdas.ThetaiOutV     :   inner tolerance for each outer iteration
%                           (in case of inner tolerance principle, see
%                           <minTolIn.html>.)
%
%% recon: minPda: Output Arguments
%
% * seti    :   structure array
% * seti.qROIcomp   :   new computed reconstruction of the contrast
%                       (complex vector of size seti.nROI^seti.dim x 1)
% * pdas    :   structure array for pda specific output, 
%               expanded for currently computed iteration step 
%               (see "Input Arguments" for a reference).
%

%%%%%

%% recon: minPda: minTolIn: Input Arguments
%
% Parameters for inner tolerance, see Section "More About".
% We give values of parameters, which worked fine.
%
% * seti.ThetaStart = 0.925
% * seti.ThetaMax = 0.95
% * seti.TolGamma = 0.90
%
% Other parameters in struct seti:
%
% * |seti.tau|    :   Discrepancy principle stops at parameter tau . delta.
% * |seti.delta|  :   relative noise level
% * |seti.FmeasDelta| :   scattered field at receivers positions (i.e. 'the data')
%                         (uScaRX = uTotRX - uIncRX, i.e. total field minus incident field)
%                         (complex matrix of size seti.measNb x seti.incNb)
%
%% recon: minPda: minTolIn: Output Arguments
%
% No fields of struct seti as output.

%%%%%

%% recon: minPda: minTolOut: Input Arguments
%
% * |seti.pdaNmax|      : Maximal iteration number, if |seti.useTolIn| or
%                         |seti.useTolOut| are 1 (default : 250).
% * |seti.relDisTol|    : Outer tolerance (default: 0.05).
%
%% recon: minPda: minTolOut: Output Arguments
%
% No fields of struct seti.
%

%%%%%

%% recon: minPda: pda: Input Arguments
%
% Most important input arguments (there are more...)
%
% * seti        :   struct seti
%
% *Some of the fields in struct seti*
%
% Note that several fields in struct seti are necessary to run pda.
%
% The routine |pda| is an internal one that needs the specific environment created in
% the package.
%
% * seti.dim    :   dimension of the problem (2 or 3)
% * seti.nROI   :   discretization points for each dimension
%                   of region of interest (ROI) (in samples)
%
% * seti.pdaN   :   number of inner iterations (PDA) (is called nPda)
% * seti.pdaStepsize    :   method to choose primal and dual stepsizes
%                           (|'fix'| (default) or |'adaptive'|), 
%                           see <pdaChoosingStepsizes.html>.
% * seti.vartheta       :   parameter in (0,1) used in case of adaptive stepsizes,
%                           see <pdaChoosingStepsizes.html>.
%
%
%% recon: minPda: pda: Output Arguments
%
% No fields of struct seti.

%%%%%


