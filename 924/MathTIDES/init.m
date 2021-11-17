(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES: Taylor Integrator of Differential EquationS*)


(* ::Text:: *)
(*MathTIDES. Version 1.3.0.*)
(*This file is part of TIDES.*)
(*  *)


(* ::Text:: *)
(*Alberto Abad, Roberto Barrio, Fernando Blesa, Marcos Rodriguez*)
(*Grupo de Mec\[AAcute]nica Espacial.  IUMA.*)
(*University of Zaragoza*)
(*50009 Zaragoza. Spain.*)
(**)
(* http://gme.unizar.es/software/tides*)
(* Contact: <tides@unizar.es>*)
(**)


(* ::Title:: *)
(*init*)


(* ::Section:: *)
(*Iniciaci\[OAcute]n*)


(* ::Subsection::Closed:: *)
(*Fichero: LKFunctions*)


DeclarePackage["MathTIDES`LKFunctions`", 
	{   "LKF","LKFC", "ToLKF","ToTaylorLKF","ExtractConstants",
		"LKFPar","LKFCPar", "ToLKFPar","ToTaylorLKFPar","IterationLKF",
		"ToLKFC","ToTaylorLKFC","ToLKFCPar","ToTaylorLKFCPar",
		"RightIterationLKF","LeftIterationLKF","LieDer",
		"NumberOfVariables", "NumberOfParameters", "NumberOfFunctions",
		"NumberOfLinks", "LinksVariables", "LinksFunctions",
		"strFunC","strFunF","strFunCmp","StringCNumber","StringFNumber",
		"ListCTextConstants","ListFTextConstants","ListDoubleConstants",
		"TextMPConstants","TextDPConstants","Double", "Multiple"}];



(* ::Subsection::Closed:: *)
(*Fichero: Odes*)


DeclarePackage["MathTIDES`ODES`", 
	{   "FirstOrderODE$", "FirstOrderODE","HamiltonianToODE",
		"PotentialToODE","NthOrderODE", "Screen","Points", "Delta","Only","Until"}];


(* ::Subsection::Closed:: *)
(*Fichero: Iterations*)


DeclarePackage["MathTIDES`Iterations`", 
	{   "ListIndexFunDer","CountDerivatives","PrevousList",
		"PreviousIndexList",
		"IteratorsList","IteratorsListStar","ListToString",
		"CompleteIteratorsList","CompleteIteratorsListStar",
		"SortListDer", "DerOutputList","ListIndexFunLastOrder"}];



(* ::Subsection::Closed:: *)
(*Fichero: MinimalCode*)


DeclarePackage["MathTIDES`MinimalCode`", 
	{"MinCCText","MinCHText","DriverMinC","MinFFText","DriverMinFortran",
		"gmecopyrightC","gmecopyrightF", "gmecopyrightDrC", "gmecopyrightDrF"}];



(* ::Subsection::Closed:: *)
(*Fichero: StandardCode*)


DeclarePackage["MathTIDES`StandardCode`", 
	{"StandardCText", "StandardHText", "DriverStdC"}];



(* ::Subsection::Closed:: *)
(*Fichero : Texts*)


DeclarePackage["MathTIDES`Texts`",
	{ "mincgen", "minhgen","minfgen","headstdDP","headstdMP"}];


(* ::Subsection::Closed:: *)
(*Fichero: Codes*)


DeclarePackage["MathTIDES`Codes`",
	{ "TSMCodeFiles", "CodeFiles", "PrecisionDigits","MinTIDES", "Driver","OnlyDriver",
	"ParametersValue", "InitialConditions", "IntegrationPoints", 
	"Output", "OutputCoefficients","DataMatrix", "Factor1","Factor2",
	"Factor3","MaxStepRatio","MinStepRatio","MaxIterationsNumber",
	"OrderIncrement","MinOrder","MaxOrder",	"RelativeTolerance",
	"AbsoluteTolerance", "DefectErrorControl", "AddFunctions", "AddPartials",
	"Optimization"
}];



(* ::Subsection::Closed:: *)
(*Fichero : PWSTides*)


DeclarePackage["MathTIDES`PWSTides`",
	{ "PWS", "ZeroPWS","UnitPWS","OrderPWS","PWSeries", "TSMSolve"}];


(* ::Subsection::Closed:: *)
(*Mensaje Inicial y s\[IAcute]mbolos en el contexto general*)


mathTIDESVersion$ = "1.30"


mathTIDESInit$ = "    MathTIDES "<>mathTIDESVersion$<>"\n" <>
"    MathTIDES    is   part   of   the   TIDES   project.\n"<>
"    Abad, A.,  Barrio, R.,  Blesa, F.  and  Rodriguez, M.\n";


FrameBox[RowBox[{StyleBox[mathTIDESInit$,FontFamily->"Geneva"],
	StyleBox[Hyperlink["http://gme.unizar.es/software/tides"]]}],
	Background->LightYellow]//DisplayForm
