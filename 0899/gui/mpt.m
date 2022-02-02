function varargout = mpt(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MPT M-file for mpt.fig
%      MPT, by itself, creates a new MPT or raises the existing
%      singleton*.
%
%      H = MPT returns the handle to a new MPT or the handle to
%      the existing singleton*.
%
%      MPT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPT.M with the given input arguments.
%
%      MPT('Property','Value',...) creates a new MPT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mpt_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mpt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mpt

% Last Modified by GUIDE v2.5 10-Jul-2008 15:13:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mpt_OpeningFcn, ...
                   'gui_OutputFcn',  @mpt_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



function mpt_OpeningFcn(hObject, eventdata, handles, varargin)        % --- Executes just before mpt is made visible.

  handles.output = hObject;                                                % Choose default command line output for signal
  guidata(hObject, handles);                                               % Update handles structure
  function varargout = mpt_OutputFcn(hObject, eventdata, handles)       % --- Outputs from this function are returned to the command line.
  varargout{1} = handles.output;                                           % Get default command line output from handles structure
  handles.isFunctionFirstTime = true;
  guidata(hObject,handles);


% --------------------------------------------------------------------
% ----- example functions, 1d ----------------------------------------
% --------------------------------------------------------------------

function examplesM_Callback(h, eventdata, ha)
function xFunction_Callback(h, eventdata, ha)
     exampleFunctionSetup1d(h, ha, 0);
function mStep_Callback(h, eventdata, ha)
     exampleFunctionSetup1d(h, ha, 1);
function combo2_Callback(h, eventdata, ha)
    exampleFunctionSetup1d(h, ha, 2);
function discontSin_Callback(h, eventdata, ha)
    exampleFunctionSetup1d(h, ha, 3);
function wFunction_Callback(h, eventdata, ha)   
    exampleFunctionSetup1d(h, ha, 4);
function shizgalFunc7_Callback(h, eventdata, ha)
  exampleFunctionSetup1d(h, ha, 5);
function smooth_Callback(h, eventdata, ha)
  exampleFunctionSetup1d(h, ha, 7);
function comboFunction3_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 8);
function sawToothWave_Callback(h, eventdata, ha)
  exampleFunctionSetup1d(h, ha, 9);
function sharpPeak_Callback(h, eventdata, ha)
  exampleFunctionSetup1d(h, ha, 10);
function centerStep_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 6);
function sheppLoganSlice_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 11);
function sinCos_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 12);
function discontinuousDerivative_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 13);
function discontDerivPeriodic_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 14);
function analyticPeriodic_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 15);
function discontFandFp_Callback(h, eventdata, ha)
   exampleFunctionSetup1d(h, ha, 16);

% --------------------------------------------------------------------
% ----- example functions, 2d ----------------------------------------
% --------------------------------------------------------------------

function examples2dmenu_Callback(h, eventdata, ha)

function sq2dMenu_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 1);
function circular_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 2);
function sheppLogan_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 3);
function circularConstant_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 4);
function circularNonCompact_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 5);
function periodicDiscontinuous2d_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 6);
function multipleJumps_Callback(h, eventdata, ha)
   exampleFunctionSetup2d(h, ha, 7);

% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------


function radiobutton1_Callback(hObject, eventdata, handles)
function radiobutton2_Callback(hObject, eventdata, handles)
function radiobutton3_Callback(hObject, eventdata, handles)
function radiobutton4_Callback(hObject, eventdata, handles)
function radiobutton5_Callback(hObject, eventdata, handles)
function ppTypeMenu_CreateFcn(hObject, eventdata, handles)
function ppTypeMenu_Callback(hObject, eventdata, handles)
function polyModeMenu_CreateFcn(hObject, eventdata, handles)
function polyModeMenu_Callback(hObject, eventdata, handles)
function processModeMenu_CreateFcn(hObject, eventdata, handles)
function processModeMenu_Callback(hObject, eventdata, handles)
function gridMenu_CreateFcn(hObject, eventdata, handles)
function gridMenu_Callback(hObject, eventdata, handles)
function filterMenu_CreateFcn(hObject, eventdata, handles)
function filterMenu_Callback(hObject, eventdata, handles)
function edgeDetectCriticalJ_CreateFcn(hObject, eventdata, handles)
function edgeDetectCriticalJ_Callback(hObject, eventdata, handles)
function Q_CreateFcn(hObject, eventdata, handles)
function Q_Callback(hObject, eventdata, handles)
function nonlinearEnhancementParam_CreateFcn(hObject, eventdata, handles)
function nonlinearEnhancementParam_Callback(hObject, eventdata, handles)
function filterAlpha_CreateFcn(hObject, eventdata, handles)
function filterAlpha_Callback(hObject, eventdata, handles)
function spectralFilterOrder_CreateFcn(hObject, eventdata, handles)
function spectralFilterOrder_Callback(hObject, eventdata, handles)
function edit7_CreateFcn(hObject, eventdata, handles)
function edit7_Callback(hObject, eventdata, handles)
function gegenbauerL_CreateFcn(hObject, eventdata, handles)
function gegenbauerL_Callback(hObject, eventdata, handles)
function M_CreateFcn(hObject, eventdata, handles)
function M_Callback(hObject, eventdata, handles)
function N_CreateFcn(hObject, eventdata, handles)
function N_Callback(hObject, eventdata, handles)
function edit19_CreateFcn(hObject, eventdata, handles)
function edit19_Callback(hObject, eventdata, handles)
function NP_CreateFcn(hObject, eventdata, handles)
function NP_Callback(hO, eventdata, handles)


function goButton_Callback(h, eventdata, ha)                    % green button

     if ha.is2d == false

         if get(ha.approximateButton,'Value')                   % 1d approximation

         elseif get(ha.edgeButton,'Value');                     % 1d edge detection
            edgeDetectionDriver(h,ha);
         elseif get(ha.postProcessButton,'Value');
            postProcessDriver(h,ha);      % 1d postprocessing
         end

      else
           postProcessDriver2d(h,ha);     % 2d postprocessing
end  % goButton


function padeM_CreateFcn(hObject, eventdata, handles)
function padeM_Callback(hObject, eventdata, handles)
function dtvIterations_CreateFcn(hObject, eventdata, handles)
function dtvIterations_Callback(hObject, eventdata, handles)
function mollGamma_CreateFcn(hObject, eventdata, handles)
function mollGamma_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function pdeExMenu_Callback(h, eventdata, ha)
   ha.is2d = false;
   ha.isPdeSolution = true;
   ha.isFunctionFirstTime = true;
   guidata(h,ha);


% --------------------------------------------------------------------

function B_CreateFcn(hObject, eventdata, handles)
function B_Callback(hObject, eventdata, handles)

function A_CreateFcn(hObject, eventdata, handles)
function A_Callback(hObject, eventdata, handles)

function gridBeta_CreateFcn(hObject, eventdata, handles)
function gridBeta_Callback(hObject, eventdata, handles)

function gridAlpha_CreateFcn(hObject, eventdata, handles)
function gridAlpha_Callback(hObject, eventdata, handles)
function optionsMenu_Callback(hObject, eventdata, handles)
function edgesEditText_CreateFcn(hObject, eventdata, handles)
function edgesEditText_Callback(hObject, eventdata, handles)


% -----------------------------------------------------------
% -------------- edge detection -----------------------------
% -----------------------------------------------------------

function manualEdges_Callback(h, eventdata, ha)
    manualEdges(h,ha);

function textEdges_Callback(h, eventdata, ha)
   textEdges(h,ha);


% ------------------------------------------------------------
% --------- export Data to Matlab workspace -------------------------------
% ------------------------------------------------------------

function savePP_Callback(h, eventdata, ha)                             
   answer = inputdlg('variable name for postprocessed:','postprocessed export',1,{'p'});
   assignin('base',answer{1},ha.fPost);

function exportFmenu_Callback(h, eventdata, ha)
    answer = inputdlg('variable name for spectral approximation:','approximation export',1,{'f'});
    assignin('base',answer{1},ha.fa);

function exportXmenu_Callback(h, eventdata, ha)
    answer = inputdlg('variable name for x:','grid export',1,{'x'});
    assignin('base',answer{1},ha.xa);

function exportExactMenu_Callback(h, eventdata, ha)
    answer = inputdlg('variable name for exact:','exact export',1,{'fe'});
    assignin('base',answer{1},ha.fExact);



% --------------------------------------------------------------------

function filterButtonGroup_SelectionChangeFcn(h, eventdata, ha)


function methodButtonGroup_SelectionChangeFcn(h, eventdata, ha)


    if get(ha.filterButton,'Value')
        set(ha.expFilterButton,'Enable','on');
        set(ha.erfLogFilterButton,'Enable','on');
        set(ha.vandevenFilterButton,'Enable','on');
        set(ha.spectralFilterOrder,'Enable','on');
     else
        set(ha.expFilterButton,'Enable','off');
        set(ha.erfLogFilterButton,'Enable','off');
        set(ha.vandevenFilterButton,'Enable','off');
        set(ha.spectralFilterOrder,'Enable','off');
     end
     
     if get(ha.grpButton,'Value') | get(ha.igrpButton,'Value') | get(ha.freudButton,'Value')
       set(ha.projectButton,'Enable','on');
       set(ha.reprojectButton,'Enable','on');
     else
       set(ha.projectButton,'Enable','off');
       set(ha.reprojectButton,'Enable','off');
     end
     
     if get(ha.dtvFilterButton,'Value')
       set(ha.dtvIterations,'Enable','on');
       set(ha.dtvLambda,'Enable','on');
       if ha.is2d
         set(ha.dtv8_radioButton,'Enable','on');
         set(ha.dtv4_radioButton,'Enable','on');
       end
     else
       set(ha.dtvIterations,'Enable','off');
       set(ha.dtvLambda,'Enable','off');
       if ha.is2d
         set(ha.dtv8_radioButton,'Enable','off');
         set(ha.dtv4_radioButton,'Enable','off');
       end
     end
     
     if get(ha.padeButton,'Value')
       set(ha.padeM,'Enable','on');
       set(ha.PadeNc,'Enable','on');
     else
       set(ha.padeM,'Enable','off');
       set(ha.PadeNc,'Enable','off');
     end
     
     if get(ha.grpButton,'Value') 
        set(ha.gegenbauerL,'Enable','on');
        set(ha.gegenbauerM,'Enable','on');
     else
        set(ha.gegenbauerL,'Enable','off');
        set(ha.gegenbauerM,'Enable','off');
     end
     
     if get(ha.igrpButton,'Value') | get(ha.grpButton,'Value')
       set(ha.gegenbauerM,'Enable','on');
     else
       set(ha.gegenbauerM,'Enable','off');
     end
     

function modeButtonGroup_SelectionChangeFcn(h, eventdata, ha)

    if get(ha.approximateButton,'Value')
      set(ha.N,'Enable','on');
      set(ha.M,'Enable','on');
      set(ha.examplesM,'Enable','on');
      set(ha.examples2dMenu,'Enable','on');
      set(ha.filterButton,'Enable','off');
      set(ha.dtvFilterButton,'Enable','off');
      set(ha.grpButton,'Enable','off');
      set(ha.igrpButton,'Enable','off');
      set(ha.padeButton,'Enable','off');
      set(ha.freudButton,'Enable','off');
      
      set(ha.expFilterButton,'Enable','off');
      set(ha.erfLogFilterButton,'Enable','off');
      set(ha.vandevenFilterButton,'Enable','off');
      set(ha.spectralFilterOrder,'Enable','off');
      
      set(ha.edgesEditText,'Enable','off');
      set(ha.edgeDetectCriticalJ,'Enable','off');
      set(ha.Q,'Enable','off');
      set(ha.nonlinearEnhancementParam,'Enable','off');
      
      set(ha.padeM,'Enable','off');
      set(ha.PadeNc,'Enable','off');
      
      set(ha.dtvIterations,'Enable','off');
      set(ha.dtvLambda,'Enable','off');
      set(ha.dtv8_radioButton,'Enable','off');
      set(ha.dtv4_radioButton,'Enable','off');
      
      set(ha.projectButton,'Enable','off');
      set(ha.reprojectButton,'Enable','off');
      
      set(ha.gegenbauerL,'Enable','off');
      set(ha.gegenbauerM,'Enable','off');
    end

    if get(ha.edgeButton,'Value') & ~ha.is2d
       set(ha.examplesM,'Enable','off');
       set(ha.examples2dMenu,'Enable','off');
       set(ha.N,'Enable','off');
       set(ha.M,'Enable','off');
       set(ha.filterButton,'Enable','off');
       set(ha.dtvFilterButton,'Enable','off');
       set(ha.grpButton,'Enable','off');
       set(ha.igrpButton,'Enable','off');
       set(ha.padeButton,'Enable','off');
       set(ha.freudButton,'Enable','off');
       set(ha.expFilterButton,'Enable','off');
       set(ha.erfLogFilterButton,'Enable','off');
       set(ha.vandevenFilterButton,'Enable','off');
       set(ha.spectralFilterOrder,'Enable','off');
       
       set(ha.edgesEditText,'Enable','on');
       set(ha.edgeDetectCriticalJ,'Enable','on');
       set(ha.Q,'Enable','on');
       set(ha.nonlinearEnhancementParam,'Enable','on');
       
       set(ha.padeM,'Enable','off');
       set(ha.PadeNc,'Enable','off');
       set(ha.dtvIterations,'Enable','off');
       set(ha.dtvLambda,'Enable','off');
       set(ha.dtv8_radioButton,'Enable','off');
       set(ha.dtv4_radioButton,'Enable','off');
       
       set(ha.projectButton,'Enable','off');
       set(ha.reprojectButton,'Enable','off');
       
       set(ha.gegenbauerL,'Enable','off');
       set(ha.gegenbauerM,'Enable','off');
    end

    if get(ha.postProcessButton,'Value')

       
       if get(ha.filterButton,'Value')
        set(ha.expFilterButton,'Enable','on');
        set(ha.erfLogFilterButton,'Enable','on');
        set(ha.vandevenFilterButton,'Enable','on');
        set(ha.spectralFilterOrder,'Enable','on');
       end
       if get(ha.grpButton,'Value') | get(ha.igrpButton,'Value') | get(ha.freudButton,'Value')
         set(ha.projectButton,'Enable','on');
         set(ha.reprojectButton,'Enable','on');
       end
       if get(ha.dtvFilterButton,'Value')
         set(ha.dtvIterations,'Enable','on');
         set(ha.dtvLambda,'Enable','on');
         if ha.is2d
           set(ha.dtv8_radioButton,'Enable','on');
           set(ha.dtv4_radioButton,'Enable','on');
         end
       end
       if get(ha.padeButton,'Value')
         set(ha.padeM,'Enable','on');
         set(ha.PadeNc,'Enable','on');
       end
       if get(ha.grpButton,'Value')
         set(ha.gegenbauerL,'Enable','on');
         set(ha.gegenbauerM,'Enable','on');
       end
       
       if get(ha.igrpButton,'Value') | get(ha.grpButton,'Value')
         set(ha.gegenbauerM,'Enable','on');
       end
    
       set(ha.examplesM,'Enable','off');
       set(ha.examples2dMenu,'Enable','off');
       set(ha.N,'Enable','off');
       set(ha.M,'Enable','off');
       set(ha.filterButton,'Enable','on');
       set(ha.dtvFilterButton,'Enable','on');
       
       if ~ha.is2d
         set(ha.grpButton,'Enable','on');
         if get(ha.fourierButton,'Value')
            set(ha.igrpButton,'Enable','on');
         end
         set(ha.padeButton,'Enable','on');
         set(ha.freudButton,'Enable','on');
       end
       
       set(ha.edgesEditText,'Enable','off');
       set(ha.edgeDetectCriticalJ,'Enable','off');
       set(ha.Q,'Enable','off');
       set(ha.nonlinearEnhancementParam,'Enable','off');
    end


% --------------------------------------------------------------------


function gegenbauerM_Callback(hObject, eventdata, handles)
function gegenbauerM_CreateFcn(hObject, eventdata, handles)


function graphicsOptionsMenu_Callback(hObject, eventdata, handles)


function surfaceMenuItem_Callback(hObject, eventdata, handles)

  if handles.is2d
    rotate3d on
    if get(handles.postProcessButton,'Value');
      surf(handles.axes1,handles.xa,handles.xa, handles.fPost);
      surf(handles.axes2,handles.xa,handles.xa, abs(handles.fPost-handles.fExact));
    else
      surf(handles.axes1,handles.xa,handles.xa, handles.fa);
      surf(handles.axes2,handles.xa,handles.xa, abs(handles.fa-handles.fExact));
    end

    if ~strcmp(get(handles.surfaceMenuItem, 'Checked'), 'on')
        set(handles.surfaceMenuItem, 'Checked', 'on');
        set(handles.countouMenuItem, 'Checked', 'off');
    end
  end

function countouMenuItem_Callback(hObject, eventdata, handles)

  if handles.is2d
    rotate3d off
    if get(handles.postProcessButton,'Value');
      contour(handles.axes1,handles.xa,handles.xa, handles.fPost);
      contour(handles.axes2,handles.xa,handles.xa, abs(handles.fPost-handles.fExact));
    else
      contour(handles.axes1,handles.xa,handles.xa, handles.fa);
      contour(handles.axes2,handles.xa,handles.xa, abs(handles.fa-handles.fExact));
    end

    if ~strcmp(get(handles.countouMenuItem, 'Checked'), 'on')
        set(handles.countouMenuItem, 'Checked', 'on');
        set(handles.surfaceMenuItem, 'Checked', 'off');
    end
  end



function displayLegends_Callback(h, eventdata, ha)
    if ~strcmp(get(ha.displayLegends, 'Checked'), 'on'), set(ha.displayLegends, 'Checked', 'on');
    else, set(ha.displayLegends, 'Checked', 'off');  end


function outputTextBox_Callback(hObject, eventdata, handles)
function outputTextBox_CreateFcn(hObject, eventdata, handles)

function PadeNc_Callback(hObject, eventdata, handles)
function PadeNc_CreateFcn(hObject, eventdata, handles)


% --------------------------------------------------------------------
% --------------- PDE examples ---------------------------------------
% --------------------------------------------------------------------


function fourierAdv_Callback(h, eventdata, ha)

function fourierAdvection_64_Callback(h, eventdata, ha)
   set(ha.fourierButton,'Value',1);
   setUpPdeExample(h,ha,0,0,'fourierAdvection64.txt','fourierAdvection64exact.txt');

function fourierAdvection128_Callback(h, eventdata, ha)
   set(ha.fourierButton,'Value',1);
   setUpPdeExample(h,ha,0,0,'fourierAdvection128.txt','fourierAdvection128exact.txt');

function fourierAdvection256_Callback(h, eventdata, ha)
    set(ha.fourierButton,'Value',1);
    setUpPdeExample(h,ha,0,0,'fourierAdvection256.txt','fourierAdvection256exact.txt');
    
function fourierAdvection_512_Callback(h, eventdata, ha)
     set(ha.fourierButton,'Value',1);
     setUpPdeExample(h,ha,0,0,'fourierAdvection512.txt','fourierAdvection512exact.txt');


function fourierBurgers_Callback(h, eventdata, ha)

function fourierBurgers_64_Callback(h, eventdata, ha)
    set(ha.fourierButton,'Value',1);
    setUpPdeExample(h,ha,0,0,'fourierBurgers64.txt','fourierBurgers64exact.txt');

function fourierBurgers128_Callback(h, eventdata, ha)
    set(ha.fourierButton,'Value',1);
    setUpPdeExample(h,ha,0,0,'fourierBurgers128.txt','fourierBurgers128exact.txt');

function fourierBurgers256_Callback(h, eventdata, ha)
    set(ha.fourierButton,'Value',1);
    setUpPdeExample(h,ha,0,0,'fourierBurgers256.txt','fourierBurgers256exact.txt');
    
function fourierBurgers_512_Callback(h, eventdata, ha)
     set(ha.fourierButton,'Value',1);
     setUpPdeExample(h,ha,0,0,'fourierBurgers512.txt','fourierBurgers512exact.txt');



function chebyshevAdvectionMenu_Callback(h, eventdata, ha)

function chebyshevAdvection_64_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.99,'chebyshevAdvection64.txt','chebyshevAdvection64exact.txt');

function chebyshevAdvection_128_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.99,'chebyshevAdvection128.txt','chebyshevAdvection128exact.txt');

function chebyshevAdvection_256_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.99,'chebyshevAdvection256.txt','chebyshevAdvection256exact.txt');

function chebyshevAdvection_512_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.999,'chebyshevAdvection512.txt','chebyshevAdvection512exact.txt');





function eulerEquationsMenu_Callback(hObject, eventdata, handles)

function chebyshevEuler_128_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.99,'chebyshevEuler128.txt','chebyshevEuler128exact.txt');

function eulerChebyshev_256_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.999,'chebyshevEuler256.txt','chebyshevEuler256exact.txt');


function chebyshevEuler_512_Callback(h, eventdata, ha)
    set(ha.chebyButton,'Value',1);
    setUpPdeExample(h,ha,1,0.999,'chebyshevEuler512.txt','chebyshevEuler512exact.txt');



% -----------------------------------------------------------------

% --------------------------------------------------------------------
function examples2dMenu_Callback(h, eventdata, ha)



