(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`MinimalCode:  C and Fortran Minimal Codes*)


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
(*MinimalCode*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`MinimalCode`",
				"MathTIDES`Iterations`",
						"MathTIDES`ODES`",
							"MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
MinCCText,
MinCHText,
DriverMinC,
MinFFText,
DriverMinFortran,
gmecopyrightC,
gmecopyrightF,
gmecopyrightDrC,
gmecopyrightDrF
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`MinimalCode`*"]
Clear @@ Names["MathTIDES`MinimalCode`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsubsection::Closed:: *)
(*General*)


$numberOfVariables = 0
$numberOfAddConstants = 0
$numberOfAddParams = 0
$listConstants = {}


$ParametersValue = Null 
$InitialConditions = Null 
$IntInt = Null 
$FileOutput = False 
$DataMatrix = False 
$Factor1 = Null 
$Factor2 = Null 
$Factor3 = Null 
$MaxStepRatio = Null 
$MinStepRatio = Null 
$MaxIterationsNumber = Null 
$ExcessOrder = Null 
$MinOrder = Null 
$MaxOrder = Null
$RelativeTolerance = Null
$AbsoluteTolerance = Null
$StepSizeEstimator = False


CleanMParams[]:=
(
	$ParametersValue = Null; 
	$InitialConditions = Null ;
	$IntInt = Null ;
	$FileOutput = False ;
	$DataMatrix = False ;
	$Factor1 = Null ;
	$Factor2 = Null ;
	$Factor3 = Null ;
	$MaxStepRatio = Null ;
	$MinStepRatio = Null ;
	$MaxIterationsNumber = Null ;
	$ExcessOrder = Null ;
	$MinOrder = Null; 
	$MaxOrder = Null;
	$RelativeTolerance = Null;
	$AbsoluteTolerance = Null;
	$StepSizeEstimator = False;
)


MinIP[Null]:= Null
MinIP[{x_,Delta[dt_],Points[n_?Positive]}]:={N[x, 16], N[x + n dt,16], N[dt, 16] }
MinIP[{x_,Points[n_?Positive], Delta[dt_]}]:={N[x, 16], N[x + n dt, 16], N[dt, 16] }
MinIP[{x_,y_,Delta[dt_]}]:={N[x, 16],N[y, 16],N[dt, 16]}
MinIP[{x_,y_,Points[n_?Positive]}]:={N[x, 16],N[y, 16],N[(y-x)/n, 16]}
MinIP[{x_,y_}]:= {N[x, 16],N[y, 16], N[y-x, 16]}
MinIP[{x_,y_,z__}]:= Null


ValueMParams[params_]:=
(
	$ParametersValue = params[[1]]; 
	$InitialConditions = params[[2]] ;
	$IntInt = MinIP[params[[3]]];
	$FileOutput = params[[4]] ;
	$DataMatrix = params[[5]] ;
	$Factor1 = params[[6]] ;
	$Factor2 = params[[7]] ;
	$Factor3 = params[[8]] ;
	$MaxStepRatio = params[[9]] ;
	$MinStepRatio = params[[10]] ;
	$MaxIterationsNumber = params[[11]] ;
	$ExcessOrder = params[[12]] ;
	$MinOrder = params[[13]]; 
	$MaxOrder = params[[14]];
	$RelativeTolerance = params[[15]];
	$AbsoluteTolerance = params[[16]];
	$StepSizeEstimator = params[[17]];

)


TextGen = "This file is part of TIDES"
TextDate[]:=
	Module[{dt, texto, mes,min, tmin},
	mes = {"January ", "February ", "March ", 
	"April ","May ", "June ", "July ", "August ", 
	"September ", "October ","November ", "December "};
	dt = Date[];
	texto = "It has been created by MathTIDES ("<> Global`mathTIDESVersion$ <> ") ";
	texto = texto <> mes[[dt[[2]]]] <>ToString[dt[[3]]] <> ", ";
	texto = texto <> ToString[dt[[1]]]<>", ";
	min = Round[dt[[5]]];
	tmin = If[min < 10, "0"<>ToString[min],  ToString[min]];
	texto = texto <> ToString[dt[[4]]] <> ":" <>tmin;
	texto
]


lineCComienzo = "/****************************************************************************\n";
lineCFinal = "\n*****************************************************************************/\n\n";
lineFComienzo = "C****************************************************************************\n";
lineFFinal = "\nC*****************************************************************************\n\n";


Cgmegnulic = 
	"\n\n\thttp://gme.unizar.es/software/tides" <>
	"\n\tContact: <tides@unizar.es>\n" 


Fgmegnulic =	
	"\nC       " <>
	"\nC     http://gme.unizar.es/software/tides" <>
	"\nC     Contact: <tides@unizar.es>" <>
	"\nC        "



gmecopyrightC[]:= lineCComienzo <> "\t"<> TextGen<>"\n\t"<> TextDate[] <> Cgmegnulic <> lineCFinal


gmecopyrightDrC["minc"]:= lineCComienzo <>
	"\tDriver file of the minc_tides program" <>
	"\n\t"<> TextGen<>"\n\t"<>TextDate[] <> Cgmegnulic <> lineCFinal

gmecopyrightDrC["dp"]:= lineCComienzo <>
	"\tDriver file of the dp_tides program" <>
	"\n\t"<> TextGen <>"\n\t"<>TextDate[] <> Cgmegnulic <> lineCFinal

gmecopyrightDrC["mp"]:= lineCComienzo <>
	"\tDriver file of the mp_tides program" <>
	"\n\t"<> TextGen <>"\n\t"<> TextDate[] <> Cgmegnulic <> lineCFinal


gmecopyrightF[]:= lineFComienzo <>
	"C     "<> TextGen <> 
	"\nC     "<> TextDate[] <> Fgmegnulic<> lineFFinal


gmecopyrightDrF[]:= lineFComienzo <>
	"C     Driver file of the minf_tides program" <>
	"\nC     "<> TextGen <> 
	"\nC     "<> TextDate[] <> Fgmegnulic<> lineFFinal


newFline[]:= "\n      "


separadorF[]:="C-------------------------------------------------------------------------------";
separadorC[]:="\n"


StringNumberC[x_]:=
	StringCNumber[N[x,16]]

StringNumberF[x_]:=
	StringFNumber[N[x,16]]


(* ::Subsubsection::Closed:: *)
(*C\[OAcute]digo C*)


strTayC[LKF$Plus]:= "add_mc";
strTayC[LKF$Times]:= "mul_mc";
strTayC[LKF$Divide]:= "div_mc";
strTayC[LKF$Power]:= "pow_mc";
strTayC[LKF$Sin]:= "sin_mc";
strTayC[LKF$Cos]:= "cos_mc";
strTayC[LKF$Tan]:= "tan_mc";
strTayC[LKF$Sinh]:= "sinh_mc";
strTayC[LKF$Cosh]:= "cosh_mc";
strTayC[LKF$Tanh]:= "tanh_mc";
strTayC[LKF$ArcSin]:= "asin_mc";
strTayC[LKF$ArcCos]:= "acos_mc";
strTayC[LKF$ArcTan]:= "atan_mc";
strTayC[LKF$ArcSinh]:= "asinh_mc";
strTayC[LKF$ArcCosh]:= "acosh_mc";
strTayC[LKF$ArcTanh]:= "atanh_mc";
strTayC[LKF$Log]:= "log_mc";


strOBJC[LKF$Var[n_]]:= ToString[n-1];
strOBJC[LKF$Par[n_]]:= 
	If[$numberOfAddParams > 0, 
		"pr[" <> ToString[n-1] <> "]",
		"p[" <> ToString[n-1] <> "]"];
strOBJC[LKF$Link[n_]]:= ToString[n+$numberOfVariables] ;
strOBJC[LKF$Const[n_]]:= $listConstants[[n]];
strOBJC[LKF$Constant[x_]]:= StringNumberC[x];
strOBJC[1]:= StringNumberC[1];
strOBJC[]:= ""


addParamsC[np_] := 
	"\n\tdouble pr["<> ToString[np]<> "];"<>
	"\n\tfor(i=0; i<PAR; i++) pr[i] = p[i];"

CtextPARAMS[np_,pit_,pr_ ]:= addParamsC[np] <>
	(StringJoin @@ Map[CstrPITER[pit[[#]],#,pr]&, Range[Length[pit]]]) <>
	separadorC[];


CstrPITER[LKF$Plus[x_,y_], i_, p_]:= 
	"\n\tpr[" <> ToString[i+p-1] <> "] = " <> strOBJC[x] <> " + " <> strOBJC[y] <> ";"

CstrPITER[LKF$Times[x_,y_], i_, p_]:= 	
	"\n\tpr[" <> ToString[i+p-1] <> "] = " <> strOBJC[x] <> " * " <> strOBJC[y] <> ";"
	
CstrPITER[LKF$Power[x_,y_], i_, p_]:= 
	"\n\tpr[" <> ToString[i+p-1] <> "] = pow(" <> strOBJC[x] <> "," <> strOBJC[y] <> ");"

CstrPITER[LKF$Divide[x_,y_], i_, p_]:= 
	"\n\tpr[" <> ToString[i+p-1] <> "] = " <> strOBJC[x] <> " / " <> strOBJC[y] <> ";"

CstrPITER[h_[x_], i_, p_]:= 
	"\n\tpr[" <> ToString[i+p-1] <> "] = " <> strFunC[h] <> "(" <> strOBJC[x] <> ");"


CtextITERS[it_]:=
	(StringJoin @@ Map[CstrLITER[it[[#]],#]&, Range[Length[it]]]);


CstrLITER[ LKF$Plus[x_,y_?numberADPQ], j_ ]:= 
	CstrLITER[ LKF$Plus[y,x], j ];

CstrLITER[ LKF$Minus[x_,y_?numberADPQ], j_ ]:= 
	Module[{texto},
		$numberOfAddConstants++;
		texto = "\n\t\tif(i==0) XX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
			"XX[" <> strOBJC[x] <> "][i]-"<> strOBJC[y] <> "; \n\t\telse " <>
			"    XX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
			"XX[" <> strOBJC[x] <> "][i];" ;
		texto]

CstrLITER[ LKF$Plus[x_?numberADPQ, y_], j_ ]:= 
	Module[{texto},
		$numberOfAddConstants++;
		texto = "\n\t\tif(i==0) XX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
			strOBJC[x] <> "+XX[" <> strOBJC[y] <> "][i]; \n\t\telse " <>
			"    XX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
			"XX[" <> strOBJC[y] <> "][i];" ;
		texto]

CstrLITER[ LKF$Minus[x_?numberADPQ, y_], j_ ]:= 
	Module[{texto},
		$numberOfAddConstants++;
		texto = "\n\t\tif(i==0) XX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
			strOBJC[x] <> "-XX[" <> strOBJC[y] <> "][i]; \n\t\telse " <>
			"    XX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
			"-XX[" <> strOBJC[y] <> "][i];" ;
		texto]

CstrLITER[ LKF$Plus[x_,y_], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
	"XX["  <> strOBJC[x] <> "][i]+XX[" <> strOBJC[y] <> "][i];";

CstrLITER[ LKF$Minus[x_,y_], j_ ]:=
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
	"XX["  <> strOBJC[x] <> "][i]-XX[" <> strOBJC[y] <> "][i];";

CstrLITER[ LKF$Times[x_,y_?numberADPQ], j_ ]:= 
	CstrLITER[ LKF$Times[y,x], j ];

CstrLITER[ LKF$Times[x_?numberADPQ, y_], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
	strOBJC[x] <> "*XX[" <> strOBJC[y] <> "][i];";

CstrLITER[ LKF$Times[x_, y_], j_ ]:=
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = mul_mc("<>
	"XX["  <> strOBJC[x] <> "],XX[" <> strOBJC[y] <> "],i);";

CstrLITER[ LKF$Divide[x_?numberADPQ, y_], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = inv_mc("<>
	strOBJC[x] <> ",XX[" <> strOBJC[y] <> "],XX[" <> 
	ToString[j+$numberOfVariables] <>"], i);";

CstrLITER[ LKF$Divide[x_, y_?numberADPQ], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = "<>
	"XX[" <> strOBJC[x] <> "][i]/" <> strOBJC[y]<>";";

CstrLITER[ LKF$Divide[x_,y_], j_ ]:=  
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = div_mc(XX["  <> 
	strOBJC[x]<>"],XX[" <> strOBJC[y] <>"],XX["<> 
	ToString[j+$numberOfVariables] <>"],i);";


CstrLITER[ LKF$Exp[y_], j_]:= 
	" \n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = exp_mc(XX["  <>  
	strOBJC[y] <>"],XX["<> ToString[j+$numberOfVariables]<> "],i);";

CstrLITER[ LKF$Power[x_, y_?numberADPQ], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = pow_mc_c(XX["  <> 
	strOBJC[x] <> "]," <> strOBJC[y] <>",XX["<> 
	ToString[j+$numberOfVariables] <>"],i);";

CstrLITER[ LKF$Log[x_LKF$Par], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = log("  <> 
	strOBJC[x] <> ");";

CstrLITER[ LKF$Log[x_], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = log_mc(XX["  <> 
	strOBJC[x] <> "],XX[" <> ToString[j+$numberOfVariables] <>"],i);";

CstrLITER[ LKF$Sin[x_, y_], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = sin_mc(XX["  <> 
	strOBJC[x] <> "],XX[" <> strOBJC[y] <>"],i);";

CstrLITER[ LKF$Cos[x_, y_], j_ ]:= 
	"\n\t\tXX[" <> ToString[j+$numberOfVariables] <> "][i] = cos_mc(XX["  <> 
	strOBJC[x] <> "],XX[" <> strOBJC[y] <>"],i);";


CstrLITER[ LKF$Const[y_], j_]:= 
	" \n\t\tXX[" <> ToString[j+$numberOfVariables] <> "] = "  <>  
	strOBJC[LKF$Const[y]] <>";";



bloqueDivC[]:=
	"\n\t\tinext = i + 1;"
CtextVARS[lf_]:=
	StringJoin @@ Map[CstrVITER[lf[[#]],#]&, Range[Length[lf]]];


CstrVITER[it_?numberADPQ,j_]:= ""

CstrVITER[it_,j_]:= 
	"\n\t\tXX[" <> ToString[j] <> "][inext] =  XX[" <> 
	strOBJC[it] <> "][i]/inext;"




CtextVARSC[lf_]:=
	Module[{fc, len, texto},
		texto = "";
		fc = Map[CVITERC[lf[[#]],#]&, Range[Length[lf]]];
		fc = Select[fc, (#!={})&];
		If[fc =!= {},
			fc = Transpose[fc];
			len = Length[fc[[1]]]; 
			texto = texto <> StringJoin @@ Map[CstrVITERC1[fc[[1,#]], fc[[2,#]]]&, Range[len]]];

		texto
		]


CVITERC[it_?numberADPQ,j_]:= {j, strOBJC[it]}
CVITERC[it_,j_]:= {}

CstrVITERC1[j_, str_]:= 
	"\n\tXX[" <> ToString[j] <> "][1] = " <> str <> ";"


bloqueDecl1CC[v_,p_,tt_]:= 
	"\n\n#include \"minc_tides.h\"\n"<>
	"\nvoid    mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO)" <> 
	"\n{\n\tint VAR,PAR,TT,i,j, inext;" <> 
	"\n\tVAR = " <> ToString[v] <> ";" <>
	"\n\tPAR = " <> ToString[p] <> ";" <>
	"\n\tTT = " <> ToString[tt] <> ";" <>
	"\n\tdouble XX[TT+1][MO+1];" 


bloqueDecl2CC[]:=
	"\n\tfor(j=0; j<=TT; j++)"<>
	"\n\t\tfor(i=0; i<=ORDER; i++)"<>
	"\n\t\t\tXX[j][i] = 0.e0;"<>
	"\n\tXX[0][0] = t;" <>
	"\n\tXX[0][1] = 1.e0;" <>
	"\n\tfor(i=1;i<=VAR;i++) {"<>
	"\n\t\tXX[i][0] = v[i-1];"<>
	"\n\t}"


bloqueEndCC[]:=
	"\n\tfor(j=0; j<=VAR; j++)"<>
	"\n\t\tfor(i=0; i<=ORDER; i++)"<>
	"\n\t\t\tXVAR[i][j] = XX[j][i];"<>
	"\n}\n\n"


MinCHText[name_?StringQ, adf_LKFPar]:=
	Module[{texto},
		texto = gmecopyrightC[];


		texto
	]


MinCCText[name_?StringQ, oldadf_LKFPar]:=
	Module[{v,p,np,tt,texto,adf},
		adf = ExtractConstants[oldadf];
		v = adf[[1]]-1;
		p = adf[[2]];
		$listConstants = ListCTextConstants[adf[[3]]];
		If[p == 0, p = 1];
		$numberOfAddParams = Length[adf[[4]]];
		$numberOfVariables = v;
		np = $numberOfAddParams+p;
		tt = Length[adf[[5]]]+v;
		texto = gmecopyrightC[];
		texto = texto <> bloqueDecl1CC[v,p,tt];
		texto = texto <> separadorC[];
		If[$numberOfAddParams > 0, 
			texto = texto <> CtextPARAMS[np, adf[[4]],p]];
		texto = texto <> bloqueDecl2CC[];
		texto = texto <> CtextVARSC[LinksVariables[adf]] ;
		texto = texto <> separadorC[];
		texto = texto <> "\n\tfor(i=0;i<ORDER;i++) {";
		texto = texto <> CtextITERS[adf[[5]]];
		texto = texto <> bloqueDivC[] ;
		texto = texto <> CtextVARS[LinksVariables[adf]] ;
		texto = texto <> "\n\t}";

		texto = texto <> separadorC[];
		texto = texto <> bloqueEndCC[];
		$numberOfVariables = 0;
		$numberOfAddConstants = 0;
		$numberOfAddParams = 0;
		$listConstants ={};
		texto
	]





(* ::Subsubsection::Closed:: *)
(*C\[OAcute]digo FORTRAN*)


strTayF[LKF$Plus]:= "add_mf";
strTayF[LKF$Minus]:= "sub_mf";
strTayF[LKF$Times]:= "mul_mf";
strTayF[LKF$Divide]:= "div_mf";
strTayF[LKF$Power]:= "pow_mf";
strTayF[LKF$Sin]:= "sin_mf";
strTayF[LKF$Cos]:= "cos_mf";
strTayF[LKF$Tan]:= "tan_mf";
strTayF[LKF$Sinh]:= "sinh_mf";
strTayF[LKF$Cosh]:= "cosh_mf";
strTayF[LKF$Tanh]:= "tanh_mf";
strTayF[LKF$ArcSin]:= "asin_mf";
strTayF[LKF$ArcCos]:= "acos_mf";
strTayF[LKF$ArcTan]:= "atan_mf";
strTayF[LKF$ArcSinh]:= "asinh_mf";
strTayF[LKF$ArcCosh]:= "acosh_mf";
strTayF[LKF$ArcTanh]:= "atanh_mf";
strTayF[LKF$Log]:= "log_mf\n";


strOBJF[LKF$Var[n_]]:= ToString[n-1];
strOBJF[LKF$Par[n_]]:= 
	If[$numberOfAddParams > 0, 
		"pr(" <> ToString[n] <> ")",
		"p(" <> ToString[n] <> ")"];
strOBJF[LKF$Link[n_]]:= ToString[n+$numberOfVariables] ;
strOBJF[LKF$Const[n_]]:= $listConstants[[n]];
strOBJF[LKF$Constant[x_]]:= StringNumberF[x];
strOBJF[1]:= StringNumberF[1];
strOBJF[]:= ""


addParamsF = newFline[]<>
	"DO j = 1, NPAR"<>newFline[]<>
	"    pr(j) = p(j)"<>newFline[]<>
	"END DO"

FtextPARAMS[pit_,pr_ ]:= addParamsF <>
	(StringJoin @@ Map[FstrPITER[pit[[#]],#,pr]&, Range[Length[pit]]]);


FstrPITER[LKF$Constant[x_], i_, p_]:= newFline[] <>
	"pr(" <> ToString[i+p] <> ") = " <> strOBJF[x]
	
FstrPITER[LKF$Plus[x_,y_], i_, p_]:= newFline[] <>
	"pr(" <> ToString[i+p] <> ") = " <> strOBJF[x] <> " + " <> strOBJF[y] 

FstrPITER[LKF$Times[x_,y_], i_, p_]:= newFline[] <>	
	"pr(" <> ToString[i+p] <> ") = " <> strOBJF[x] <> " * " <> strOBJF[y] 
	
FstrPITER[LKF$Power[x_,y_], i_, p_]:= newFline[] <>
	"pr(" <> ToString[i+p] <> ") = " <> strOBJF[x] <> " ** " <> strOBJF[y]

FstrPITER[LKF$Divide[x_,y_], i_, p_]:= newFline[] <>
	"pr(" <> ToString[i+p] <> ") = " <> strOBJF[x] <> " / " <> strOBJF[y] 

FstrPITER[h_[x_], i_, p_]:= newFline[] <>
	"pr(" <> ToString[i+p] <> ") = " <> strFunF[h] <> "(" <> strOBJF[x] <> ")"



FtextITERS[it_]:=
	(StringJoin @@ Map[FstrLITER[it[[#]],#]&, Range[Length[it]]]);


FstrLITER[ LKF$Plus[x_,y_?numberADPQ], j_ ]:= 
	FstrLITER[ LKF$Plus[y,x], j ];

FstrLITER[ LKF$Minus[x_,y_?numberADPQ], j_ ]:= 
	Module[{texto},
		$numberOfAddConstants++;
		texto = newFline[] <>
			"    IF (i .EQ. 0) THEN" <> newFline[] <>
			"        XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
			"XX(i," <> strOBJF[x] <> ")-"<> strOBJF[y] <> newFline[] <>
			"    ELSE" <> newFline[] <>
			"        XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
			"XX(i," <> strOBJF[x] <> ")"<> newFline[] <>
			"    END IF" ;
		texto]

FstrLITER[ LKF$Plus[x_?numberADPQ, y_], j_ ]:= 
	Module[{texto},
		$numberOfAddConstants++;
		texto = newFline[] <>
			"    IF (i .EQ. 0) THEN" <> newFline[] <>
			"        XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
			strOBJF[x] <> "+XX(i," <> strOBJF[y] <> ")"<> newFline[] <>
			"    ELSE" <> newFline[] <>
			"        XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
			"XX(i," <> strOBJF[y] <> ")"<> newFline[] <>
			"    END IF" ;
		texto]

FstrLITER[ LKF$Minus[x_?numberADPQ, y_], j_ ]:= 
	Module[{texto},
		$numberOfAddConstants++;
		texto = newFline[] <>
			"    IF (i .EQ. 0) THEN" <> newFline[] <>
			"        XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
			strOBJF[x] <> "-XX(i," <> strOBJF[y] <> ")"<> newFline[] <>
			"    ELSE" <> newFline[] <>
			"        XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
			"-XX(i," <> strOBJF[y] <> ")"<> newFline[] <>
			"    END IF" ;
		texto]

FstrLITER[ LKF$Plus[x_,y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
	"XX(i,"  <> strOBJF[x] <> ")+XX(i," <> strOBJF[y] <> ")";

FstrLITER[ LKF$Minus[x_,y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
	"XX(i,"  <> strOBJF[x] <> ")-XX(i," <> strOBJF[y] <> ")";

FstrLITER[ LKF$Times[x_,y_?numberADPQ], j_ ]:= 
	FstrLITER[ LKF$Times[y,x], j ];

FstrLITER[ LKF$Times[x_?numberADPQ, y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
	strOBJF[x] <>"* XX(i," <> strOBJF[y] <> ")";

FstrLITER[ LKF$Times[x_, y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = mul_mf("  <> 
	strOBJF[x] <> "," <> strOBJF[y] <> ",i,XX,TT,MO)";

FstrLITER[ LKF$Divide[x_?numberADPQ, y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = inv_mf("<>
	strOBJF[x] <>"," <> strOBJF[y] <> ","<> ToString[j+$numberOfVariables]<> 
	",i,XX,TT,MO)";

FstrLITER[ LKF$Divide[x_, y_?numberADPQ], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = "<>
	"XX(i," <> strOBJF[x] <> ")/"<> strOBJF[y];

FstrLITER[ LKF$Divide[x_,y_], j_ ]:=  newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = div_mf("  <> 
	strOBJF[x] <> "," <> strOBJF[y] <>","<> 
	ToString[j+$numberOfVariables] <>",i,XX,TT,MO)";


FstrLITER[ LKF$Exp[y_], j_]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = exp_mf("  <>  
	strOBJF[y] <>","<> ToString[j+$numberOfVariables]<> ",i,XX,TT,MO)";

FstrLITER[ LKF$Power[x_, y_?numberADPQ], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = pow_mf_c("  <> 
	strOBJF[x] <> "," <> strOBJF[y] <>","<> 
	ToString[j+$numberOfVariables] <>",i,XX,TT,MO)";

FstrLITER[ LKF$Log[x_LKF$Par], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = DLOG("  <> 
	strOBJF[x] <> ")";

FstrLITER[ LKF$Log[x_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = log_mf("  <> 
	strOBJF[x] <> "," <>ToString[j+$numberOfVariables] <>",i,XX,TT,MO)";

FstrLITER[ LKF$Sin[x_, y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = sin_mf("  <> 
	strOBJF[x] <> "," <> strOBJF[y] <>",i,XX,TT,MO)";

FstrLITER[ LKF$Cos[x_, y_], j_ ]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = cos_mf("  <> 
	strOBJF[x] <> "," <> strOBJF[y] <>",i,XX,TT,MO)";

CstrLITER[ LKF$Const[y_], j_]:= 
	" \n\t\tXX[" <> ToString[j+$numberOfVariables] <> "] = "  <>  
	strOBJC[y] <>";";

FstrLITER[ LKF$Const[y_], j_]:= newFline[] <>
	"    XX(i," <> ToString[j+$numberOfVariables] <> ") = "  <> strOBJF[LKF$Const[y]] ;



bloqueDiv[]:=newFline[] <>
	"    inext = i + 1"
FtextVARS[lf_]:=
	StringJoin @@ Map[FstrVITER[lf[[#]],#]&, Range[Length[lf]]];


FstrVITER[it_?numberADPQ,j_]:= ""

FstrVITER[it_,j_]:= newFline[] <>
	"    XX(inext," <> ToString[j] <> ") = " <> 
	"XX(i," <> strOBJF[it] <> ")/inext"



FtextVARSC[lf_]:=
	Module[{fc, len, texto},
		texto = "";
		fc = Map[FVITERC[lf[[#]],#]&, Range[Length[lf]]];
		fc = Select[fc, (#!={})&];
		If[fc =!= {},
			fc = Transpose[fc];
			len = Length[fc[[1]]]; 
			texto = texto <> StringJoin @@ Map[FstrVITERC1[fc[[1,#]], fc[[2,#]]]&, Range[len]]];
			texto = texto <> "\n"; 
		texto
		]


FVITERC[it_?numberADPQ,j_]:= {j, strOBJF[it]}
FVITERC[it_,j_]:= {}

FstrVITERC1[j_, str_]:= newFline[] <>
	"XX(1," <> ToString[j] <> ") = " <> str




bloqueDecl1FF[tt_]:= newFline[] <>
	"SUBROUTINE minfseries(t,v,NVAR,p,NPAR,XVAR,ORDER,MO)" <> newFline[]<> 
	"IMPLICIT   NONE" <> newFline[] <>
	"INTEGER    i,inext,j,ORDER,NVAR,NPAR,TT,MO" <> newFline[] <>
	"REAL*8     XVAR(0:MO,0:NVAR)" <> newFline[] <>
	"REAL*8     XX(0:MO,0:"<>ToString[tt]<>")" <> newFline[] <>
	"REAL*8     v(NVAR)" <> newFline[] <>
	"REAL*8     p(NPAR)" <> newFline[] <>
	"REAL*8     t\n" 


bloqueDecl2FF[]:= newFline[] <>
	"REAL*8     mul_mf" <>newFline[] <>
	"REAL*8     div_mf" <>newFline[] <>
	"REAL*8     inv_mf" <>newFline[] <>
	"REAL*8     exp_mf" <>newFline[] <>
	"REAL*8     pow_mf_c" <>newFline[] <>
	"REAL*8     sin_mf" <>newFline[] <>
	"REAL*8     cos_mf" <>newFline[] <>
	"REAL*8     log_mf\n"


bloqueInitFF[tt_]:= newFline[]<>
	"TT = "<>ToString[tt]<>newFline[]<>
	"DO j=0, TT"<>newFline[]<>
	"    DO i=0, ORDER"<>newFline[]<>
	"        XX(i,j) = 0.d0"<>newFline[]<>
	"    END DO"<>newFline[]<>
	"END DO"<>newFline[]<>
	"XX(0,0) = t"<>newFline[]<>
	"XX(1,0) = 1.d0"<>newFline[]<>
	"DO j = 1, NVAR"<>newFline[]<>
	"    XX(0,j) = v(j)"<>newFline[]<>
	"END DO"


bloqueDOAFF[]:=newFline[]<> "DO i=0, ORDER-1"
bloqueDOCFF[]:=newFline[]<> "END DO\n"


bloqueEnd1FF[]:= newFline[]<>
	"DO j=0, NVAR"<>newFline[]<>
	"    DO i=0, ORDER"<>newFline[]<>
	"        XVAR(i,j) = XX(i,j)"<>newFline[]<>
	"    END DO"<>newFline[]<>
	"END DO"

bloqueEnd2FF[]:= newFline[]<>
	"RETURN"<>newFline[]<>
	"END SUBROUTINE\n"



AddConstant[LKF$Minus[n_?numberADPQ,y_]]:= n
AddConstant[LKF$Minus[y_,n_?numberADPQ]]:= n
AddConstant[LKF$Plus[n_?numberADPQ,y_]]:= n
AddConstant[LKF$Plus[y_,n_?numberADPQ]]:= n
AddConstant[h_[y_,n_]]:= Null


AddConstants[lkf_LKFPar]:=
	Module[{iter,list},
		iter = lkf[[4]];
		list = Map[AddConstant[iter[[#]]]&, Range[Length[iter]]];
		Select[list,(#=!=Null)&]]


MinFFText[name_?StringQ, oldadf_LKFPar]:=
	Module[{v,p,np,tt,texto,adf},
		adf = ExtractConstants[oldadf];
		v = adf[[1]]-1;
		p = adf[[2]];
		$listConstants = ListFTextConstants[adf[[3]]];
		If[p == 0, p = 1];
		$numberOfAddParams = Length[adf[[4]]];
		$numberOfVariables = v;
		np = $numberOfAddParams+p;
		tt = Length[adf[[5]]]+v;
		texto = gmecopyrightF[];
		texto = texto <> bloqueDecl1FF[tt];
		If[$numberOfAddParams > 0, 
			texto = texto <> newFline[] <> 
				"REAL*8     pr("<>ToString[np]<>")\n"];
		texto = texto <> separadorF[];
		texto = texto <> bloqueDecl2FF[];
		texto = texto <> separadorF[];
		If[$numberOfAddParams > 0, 
			texto = texto <> FtextPARAMS[adf[[4]],p]];
		texto = texto <> bloqueInitFF[tt];
		texto = texto <> FtextVARSC[LinksVariables[adf]] ;
		texto = texto <> separadorF[];
		texto = texto <> bloqueDOAFF[];
		texto = texto <> FtextITERS[adf[[5]]];
		texto = texto <> bloqueDiv[] ;
		texto = texto <> FtextVARS[LinksVariables[adf]] ;
		texto = texto <> bloqueDOCFF[];
		texto = texto <> separadorF[];
		texto = texto <> bloqueEnd1FF[];
		texto = texto <> bloqueEnd2FF[];
		texto = texto <> separadorF[];
		$numberOfVariables = 0;
		$numberOfAddConstants = 0;
		$numberOfAddParams = 0;
		$listConstants = {};
		texto
	]



(* ::Subsubsection::Closed:: *)
(*Driver C*)


bloqueDRC2[]:= 
	"#include \"minc_tides.h\""<> 
	"\n\nint main() {" 

bloqueDRC3[v_,np_]:= "\n" <>
	"\n\tint  i, VARS, PARS; "<>
	"\n\tVARS = " <> ToString[v] <> ";" <> 
	"\n\tPARS = " <> ToString[np] <> ";" <> 
	"\n\tdouble tolrel, tolabs, tini, tend, dt; "<>
	"\n\tdouble v[VARS], p[PARS]; " <>
	If[$FileOutput === False, 
			"\n\textern int dense_output; \n",
			"\n\textern FILE     *fd; \n"];


bloqueDRC4[]:=
	Module[{texto="\n", separador, cond = False},
		separador ="/***********************************************************/";
		If[$Factor1 =!= Null || $Factor2 =!= Null || $Factor3 =!= Null ||
			$MaxStepRatio =!= Null || $MinStepRatio =!= Null ||
			$MaxIterationsNumber =!= Null ||$ExcessOrder =!= Null ||
			$MinOrder =!= Null ||$MaxOrder =!= Null, 
				texto = texto <> separador <> "\n";
				texto = texto <> separador <> "\n";
				texto = texto <> "/*       CONSTANTS OF THE METHOD                            */"<> "\n";
				texto = texto <> separador <> "\n";
				texto = texto <> separador <> "\n";
				cond = True];
		If[$Factor1 =!= Null,  texto = texto <> 
			"\n\textern double fac1;" <>
			"\n\tfac1 = "<>StringNumberC[$Factor1]<>";"];
		If[$Factor2 =!= Null,  texto = texto <> "\n\tfac2 = " 
			"\n\textern double fac2;"  <>
			StringNumberC[$Factor2]<>";"];
		If[$Factor3 =!= Null,  texto = texto  
			"\n\textern double fac3;"  <>
			"\n\tfac3 = "<>StringNumberC[$Factor3]<>";"];
		If[$MaxStepRatio =!= Null,  texto = texto <> 
			"\n\textern double rmaxstep;"  <>
			"\n\trmaxstep = "<>StringNumberC[$MaxStepRatio]<>";"];
		If[$MinStepRatio =!= Null,  texto = texto <> 
			"\n\textern double rminstep;"  <>
			"\n\trminstep = "<>StringNumberC[$MinStepRatio]<>";"];
		If[$MaxIterationsNumber =!= Null,  texto = texto <> 
			"\n\textern int nitermax;"  <>
			"\n\tnitermax = "<>ToString[$MaxIterationsNumber]<>";"];
		If[$ExcessOrder =!= Null,  texto  = texto <> 
			"\n\textern int nordinc;"  <>
			"\n\tnordinc = "<>ToString[$ExcessOrder]<>";"];
		If[$MinOrder =!= Null,  texto = texto <> 
			"\n\textern int minord;"  <>
			"\n\tminord = "<>ToString[$MinOrder]<>";"];
		If[$StepSizeEstimator === True,  texto = texto <> 
			"\n\textern int defect_error_control;"  <>
			"\n\tdefect_error_control = 1;"];
		If[cond, texto = texto <> "\n"];
		texto
		]



bloqueDRC5a[v_,p_]:= 
	Module[{texto="\n\n", separador},
		separador ="/************************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*      INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */"<> "\n";
		If[$ParametersValue === Null || $InitialConditions === Null || $IntInt === Null, 
			texto = texto <> "/*      Change *****  by numerical values if it is necesary */"<> "\n";];
		texto = texto <> separador <> "\n";
		texto = texto <> separador ;
		If[p > 0,
			texto = texto <> "\n\n/* --- PARAMETERS VALUE --- */";
			Map[(texto = texto <> "\n\tp[" <> ToString[#] <>
				 "] = "<> 
				If[$ParametersValue === Null, "*****", StringNumberC[$ParametersValue[[#+1]]]] <>
				" ; ")&, Range[0,p-1]]];
		texto = texto <> "\n\n/* --- INITIAL VALUES --- */";
		Map[(texto = texto <>  "\n\tv[" <> ToString[#] <>
				 "] = "<> 
				If[$InitialConditions === Null, "*****", StringNumberC[$InitialConditions[[#+1]]]] <>
				" ; ")&, Range[0,v-1]];
		texto
	]

bloqueDRC5b[]:= 
	Module[{texto="", dtstr},
		texto = texto <> "\n\n/* --- INITIAL INTEGRATION POINT --- */";
		texto = texto <> "\n\ttini = "<> 
		If[$IntInt === Null, "*****", StringNumberC[$IntInt[[1]]]] <> " ;";
		texto = texto <> "\n\n/* --- ENDPOINT OF INTEGRATION   --- */";
		texto = texto <> "\n\ttend = "<> 
		If[$IntInt === Null, "*****", StringNumberC[$IntInt[[2]]]] <> " ;";
		texto = texto <> "\n\n/* --- DELTA t FOR DENSE OUTPUT  --- */";
		texto = texto <> "\n\tdt   = "<> 
		If[$IntInt === Null, "*****", StringNumberC[$IntInt[[3]]]] <> " ;";
		texto]

bloqueDRC5c[]:= 
	Module[{texto="", trel, tabs},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			tabs = trel = 10^-16];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = $AbsoluteTolerance];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = $RelativeTolerance];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = $RelativeTolerance; tabs = $AbsoluteTolerance];
		texto = texto <> "\n\n/* --- REQUIRED TOLERANCES --- */";
		texto = texto <> "\n\ttolrel = "<> StringNumberC[trel]<>" ;";
		texto = texto <> "\n\ttolabs = "<> StringNumberC[tabs]<>" ;";
		texto
	]

bloqueDRC5d[]:= 
	Module[{texto="\n\n", separador},	
		separador = "/***********************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*             DENSE OUTPUT (file, screen or none)          */"<> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		If[$FileOutput === False, 
			texto = texto <> "\n\tdense_output = 0;"];
		If[$FileOutput === Screen, 
			texto = texto <> "\n\tfd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tfd = fopen(\""<>$FileOutput<>"\", \"w\");"];
		texto	
	]



bloqueDRC6[]:= 
	Module[{texto="\n\n", separador},
		separador ="/***********************************************************/";
		If[$FileOutput =!= Screen, 
			texto = texto <> separador <> "\n";
			texto = texto <> separador <> "\n";
			texto = texto <> "/*       SHOW INITIAL POINT ON THE SCREEN                   */"<> "\n";
			texto = texto <> separador <> "\n";
			texto = texto <> separador <> "\n";
			texto = texto <> "\n\tprintf(\"Initial time = %25.15le, x = \", tini);";
			texto = texto <> "\n\tfor(i = 0; i < VARS-1; i++) printf(\"%24.15le,\", v[i]);";
			texto = texto <> "\n\tprintf(\"%25.15le\\n\", v[VARS-1]);\n\n";];
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*       CALL THE INTEGRATOR                               */"<> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "\n\tminc_tides(v,VARS,p,PARS,tini,tend,dt,tolrel,tolabs);\n\n";
		If[$FileOutput =!= Screen, 
			texto = texto <> separador <> "\n";
			texto = texto <> separador <> "\n";
			texto = texto <> "/*       SHOW FINAL POINT ON THE SCREEN                    */"<> "\n";
			texto = texto <> "/*       SHOW ACCEPTED AND REJECTED STEPS ON THE SCREEN    */"<> "\n";
			texto = texto <> separador <> "\n";
			texto = texto <> separador <> "\n";
			texto = texto <> "\n\tprintf(\"  Final time = %25.15le, x = \", tend);";
			texto = texto <> "\n\tfor(i = 0; i < VARS-1; i++) printf(\"%24.15le,\", v[i]);";
			texto = texto <> "\n\tprintf(\"%25.15le\\n\", v[VARS-1]);\n";];
		
		If[Head[$FileOutput] === String,
			texto = texto  <>  "\n\tfclose(fd); \n"];
		texto = texto <> "\n\treturn 0;";
		texto = texto <> "\n}\n";
		texto]


bloqueDRC7[name_]:= 
	  "\n/***********************************************************/" <>
      "\n#include \"minc_tides.c\" "<>  
      "\n#include \""<> name <>".c\" \n" <>
      "/***********************************************************/" 


DriverMinC[name_?StringQ, adf_LKFPar, params_List]:=
	Module[{v, p, np, texto},
		ValueMParams[params]; 
		texto = gmecopyrightDrC["minc"];
		v = adf[[1]]-1;
		p = adf[[2]];
		If[p == 0, np = 1, np = p];
		texto = texto <> bloqueDRC2[];
		texto = texto <> bloqueDRC3[v,np];
		texto = texto <> bloqueDRC4[];  	
		texto = texto <> bloqueDRC5a[v,p];	
		texto = texto <> bloqueDRC5b[];	
		texto = texto <> bloqueDRC5c[];	
		texto = texto <> bloqueDRC5d[];	
		texto = texto <> bloqueDRC6[];
		texto = texto <> separadorC[];
		(*texto = texto <> bloqueDRC7[name] *)
		texto = texto <> "\n\n\n";
		CleanMParams[];
		texto
	]




(* ::Subsubsection::Closed:: *)
(*Driver Fortran*)


bloqueDRF2[name_]:= newFline[]<>
	"Program  dr_"<> name <> newFline[]<>
	"IMPLICIT NONE" <> newFline[] <>
	"INTEGER  i,j\n"  

bloqueDRF3[v_,np_]:= 
	"C --- NUMBER OF VARIABLES AND PARAMETERS" <> newFline[] <>
	"INTEGER  NVAR,NPAR"<>newFline[]<>
	"PARAMETER  (NVAR = " <> ToString[v] <> ")" <> newFline[] <>
	"PARAMETER  (NPAR = " <> ToString[np] <> ")\n"

bloqueDRF4[]:= 
	"C --- TOLERANCES" <> newFline[] <>
	"REAL*8 tolabs,tolrel\n" <> 
	"C --- TIMES: INITIAL, FINAL, INCREMENT" <> newFline[] <>
	"REAL*8 tini, tend, dt\n" <> 
	"C --- VARIABLES AND PARAMETERS" <> newFline[] <>
	"REAL*8 v(NVAR)" <> newFline[] <>
	"REAL*8 p(NPAR)\n"<>
	"C --- FILE NAME AND UNIT NUMBER OF DENSE OUTPUT" <> newFline[] <>
	"CHARACTER fname*20" <> newFline[] <> 
	"INTEGER   FL\n" <> 
	"C --- OPTIONS" <> newFline[] <>
	"LOGICAL dense_output, defect_error_control\n" <>  
	"C --- COUNTERS" <> newFline[] <>
	"INTEGER accepted_steps, rejected_steps\n" <>  
	"C --- CONSTANTS OF THE METHOD (safety factors, maximum order, ...)" <> newFline[] <>
	"REAL*8 fac1,fac2,fac3,rminstep,rmaxstep" <> newFline[] <>
	"INTEGER nitermax,nordinc,minord,maxord\n" <>  
	"C --- GLOBALS" <> newFline[] <>
	"COMMON /OPT/ dense_output, defect_error_control" <> newFline[] <>
	"COMMON /ARS/ accepted_steps, rejected_steps" <> newFline[] <>
	"COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep" <> newFline[] <>
	"COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord" <> newFline[] <>  
	"COMMON /FILE/ FL\n\n" 


bloqueDRF4b[]:=
	Module[{texto="\n", separador},
		separador ="/***********************************************************/";
		If[$Factor1 =!= Null || $Factor2 =!= Null || $Factor3 =!= Null ||
			$MaxStepRatio =!= Null || $MinStepRatio =!= Null ||
			$MaxIterationsNumber =!= Null ||$ExcessOrder =!= Null ||
			$MinOrder =!= Null ||$MaxOrder =!= Null || $StepSizeEstimator === True, 
				texto = texto <> separadorF[] <> "\n";
				texto = texto <> separadorF[] <> "\n";
				texto = texto <> "C       CONSTANTS OF THE METHOD" <> "\n";
				texto = texto <> separadorF[] <> "\n";
				texto = texto <> separadorF[] <> "\n";];
		If[$Factor1 =!= Null,  texto = texto <> newFline[] <>"fac1 = "<>StringNumberF[$Factor1]];
		If[$Factor2 =!= Null,  texto = texto <> newFline[] <>"fac2 = "<>StringNumberF[$Factor2]];
		If[$Factor3 =!= Null,  texto = texto <> newFline[] <>"fac3 = "<>StringNumberF[$Factor3]];
		If[$MaxStepRatio =!= Null,  texto = texto <> newFline[] <>"rmaxstep = "<>StringNumberF[$MaxStepRatio]];
		If[$MinStepRatio =!= Null,  texto = texto <> newFline[] <>"rminstep = "<>StringNumberF[$MinStepRatio]];
		If[$MaxIterationsNumber =!= Null,  texto = texto <> newFline[] <>"nitermax = "<>ToString[$MaxIterationsNumber]];
		If[$ExcessOrder =!= Null,  texto = texto <> newFline[] <>"nordinc = "<>ToString[$ExcessOrder]];
		If[$MinOrder =!= Null,  texto = texto <> newFline[] <>"minord = "<>ToString[$MinOrder]];
		If[$MaxOrder =!= Null,  texto = texto <> newFline[] <>"maxord = "<>ToString[$MaxOrder]];
		If[$StepSizeEstimator === True,  texto = texto <> newFline[] <> "defect_error_control = .TRUE."];
		texto
		]


bloqueDRF5a[v_,p_]:= 
	Module[{texto="\n\n"},
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> "C     INITIAL CONDITIONS,  INTEGRATION TIMES, TOLERANCES"<> "\n";
		If[$ParametersValue === Null || $InitialConditions === Null, 
			texto = texto <> "C     Change ***** by numerical values if it is necesary"<> "\n";];
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> separadorF[] <> "\n";
		If[p > 0,
			texto = texto <> "\nC --- PARAMETERS VALUE";
			Map[(texto = texto <> newFline[] <> "p(" <> ToString[#] <>
				 ") = "<> 
				If[$ParametersValue === Null, "*****", StringNumberF[$ParametersValue[[#]]]] 
			)&, Range[p]]];
		texto = texto <> "\n\nC --- INITIAL VALUES";
		Map[(texto = texto <> newFline[] <> "v(" <> ToString[#] <>
				 ") = "<>
				If[$InitialConditions === Null, "*****", StringNumberF[$InitialConditions[[#]]]] 
			)&, Range[v]];
		texto
	]

bloqueDRF5b[]:= 
	Module[{texto="", dtstr},
		texto = texto <> "\n\nC --- INITIAL INTEGRATION POINT"<>newFline[];
		texto = texto <> "tini = " <>
		If[$IntInt === Null, "*****", StringNumberF[$IntInt[[1]]]] <>"\n";
		texto = texto <> "\nC --- ENDPOINT OF INTEGRATION"<>newFline[];
		texto = texto <> "tend = "<>
		If[$IntInt === Null, "*****", StringNumberF[$IntInt[[2]]]] <>"\n";
		texto = texto <> "\nC --- DELTA t FOR DENSE OUTPUT"<>newFline[];
		texto = texto <> "dt   = "<>
		If[$IntInt === Null, "*****", StringNumberF[$IntInt[[3]]]] <>"\n";
		texto
	]

bloqueDRF5c[]:= 
	Module[{texto="", trel, tabs},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			tabs = trel = 10^-16];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = $AbsoluteTolerance];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = $RelativeTolerance];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = $RelativeTolerance; tabs = $AbsoluteTolerance];
		texto = texto <> "\n\nC --- REQUIRED TOLERANCES"<>newFline[];
		texto = texto <> "tolrel = "<> StringNumberF[trel]<>newFline[];
		texto = texto <> "tolabs = "<> StringNumberF[tabs];
		texto
	]

bloqueDRF5d[]:= 
	Module[{texto="\n\n"},
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> "C       DENSE OUTPUT (file , screen or none)"<> "\n";
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> separadorF[] <> "\n" <> newFline[];
		If[$FileOutput === False, 
			texto = texto <> "dense_output = .FALSE."<> newFline[]];
		If[Head[$FileOutput] === String,
			texto = texto <> "FL = 72"<> newFline[];
			texto = texto <> "OPEN (UNIT = FL, FILE = '"<>$FileOutput <>"', STATUS = 'UNKNOWN')"<> newFline[];];
		If[ $FileOutput === Screen, 
			texto = texto <> "FL = 6"<> newFline[];];
		texto	
	]


bloqueDRF6[]:= 
	Module[{texto="\n"},
		If[$FileOutput =!= Screen, 
			texto = texto <> separadorF[] <> "\n";
			texto = texto <> separadorF[] <> "\n";
			texto = texto <> "C       SHOW INITIAL POINT ON THE SCREEN"<> "\n";
			texto = texto <> separadorF[] <> "\n";
			texto = texto <> separadorF[] <> "\n"<> newFline[];
			texto = texto <> "WRITE (*,91) tini,(v(i),i=1,NVAR)"<>"\n\n";];
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> "C       CALL THE INTEGRATOR"<> "\n";
		texto = texto <> separadorF[] <> "\n";
		texto = texto <> separadorF[] <> "\n"<> newFline[];
		texto = texto <> "CALL minf_tides(v,NVAR,p,NPAR,tini,tend,dt,\n";
		texto = texto <> "     &   tolrel,tolabs)\n\n";
		If[$FileOutput =!= Screen, 
			texto = texto <> separadorF[] <> "\n";
			texto = texto <> separadorF[] <> "\n";
			texto = texto <> "C       SHOW FINAL POINT ON THE SCREEN"<> "\n";
			texto = texto <> "C       SHOW ACCEPTED AND REJECTED STEPS ON THE SCREEN"<> "\n";
			texto = texto <> separadorF[] <> "\n";
			texto = texto <> separadorF[] <> "\n"<> newFline[];
			texto = texto <> "WRITE (*,91) tend,(v(i),i=1,NVAR)\n";
			texto = texto <> "91    FORMAT(1X,'t =',E25.16,'    X =',90E25.16)"<> "\n";];
		
		If[Head[$FileOutput] === String,
			texto = texto  <> newFline[] <> "CLOSE(FL)\n"];

		texto = texto <> newFline[] <> "STOP"<> newFline[];
		texto = texto <> "END\n\n";
		texto]


bloqueDRF7[name_]:= newFline[]<>
      "INCLUDE 'minf_tides.for'  "<>  newFline[]<>
      "INCLUDE '"<> name <>".for'\n"


DriverMinFortran[name_?StringQ, adf_LKFPar, params_List]:=
	Module[{v, p, np, texto},
		ValueMParams[params];
		texto = gmecopyrightDrF[];
		v = adf[[1]]-1;
		p = adf[[2]];
		If[p == 0, np = 1, np = p];
		texto = texto <> bloqueDRF2[name];
		texto = texto <> bloqueDRF3[v,np];
		texto = texto <> bloqueDRF4[];
		texto = texto <> bloqueDRF4b[];
		texto = texto <> bloqueDRF5a[v,p];
		texto = texto <> bloqueDRF5b[];
		texto = texto <> bloqueDRF5c[];
		texto = texto <> bloqueDRF5d[];
		texto = texto <> bloqueDRF6[];
		(*texto = texto <> separadorF[];
		texto = texto <> bloqueDRF7[name];
		texto = texto <> separadorF[];*)
		texto = texto <> "\n\n\n";
		CleanMParams[];
		texto
	]


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`MinimalCode`"]

EndPackage[]

Null
