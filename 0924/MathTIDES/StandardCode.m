(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`StandardCode:  C Standard TIDES Code*)


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
(*StandardCode*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`StandardCode`",
				"MathTIDES`Iterations`",
					"MathTIDES`ODES`",
						"MathTIDES`Texts`",
							"MathTIDES`LKFunctions`",
								"MathTIDES`MinimalCode`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
StandardCText, 
StandardHText, 
DriverStdC, 
SinCosLKFList$, 
SinCoshLKFList$
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`StandardCode`*"]
Clear @@ Names["MathTIDES`StandardCode`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*Constantes e Iteraciones*)


ListOfPartials[fun_FirstOrderODE$]:=
	If[fun[[5]] =={}, {}, fun[[5,1]]]

NumberOfPartials[fun_FirstOrderODE$]:= 
	Length[ListOfPartials[fun]]


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


StdIPText[n_,x_List,dp_]:= 
	Module[{texto = ""},
		If[dp > 16, 
			Map[( texto = texto<>"\n\tmpfrts_set_str(&lt[" <>
			ToString[#] <>"], "<>  strNUMBERC[x[[#+1]],dp] <>
			" ); ")&, Range[0,n-1]],
			Map[(texto = texto<>"\n\tlt[" <> ToString[#] <>
			"] = "<> StringNumber[x[[#+1]]] <>
			" ; ")&, Range[0,n-1]]];
		texto]

StdIPTextDP[{t0_, dt_}]:=
	Module[{texto = ""},
		texto = texto <> "\n\tint i; ";
		texto = texto <> "\n\tfor(i=0; i<nipt; i++) ";
		texto = texto <>"\n\t\tlt[i]= " <> t0 <> " + i * " <> dt <> ";";
		texto]

StdIPTextMP[{t0_, dt_}]:=
	Module[{texto = ""},
		texto = texto <> "\n\tmpfr_t t0,dt,idt;";
		texto = texto <> "\n\tmpfrts_init(&t0);";
		texto = texto <> "\n\tmpfrts_init(&dt);";
		texto = texto <> "\n\tmpfrts_init(&idt);";
		texto = texto <> "\n\tmpfrts_set_str(&t0, "<>t0<>");";
		texto = texto <> "\n\tmpfrts_set_str(&dt, "<>dt<>");";
		texto = texto <> "\n\tfor(i=0; i<nipt; i++) {";
		texto = texto <> "\n\t\tmpfrts_mul_i(&idt,dt,i);";
		texto = texto <> "\n\t\tmpfrts_add(&lt[i],t0,idt);";
		texto = texto <> "\n\t}";
		texto]


StdIP[x_, 3999999999]:= Null

StdIP[Null,dp_]:= Null;

StdIP[{x_,y_},dp_]:= 
	{2,StdIPText[2,{x,y},dp]}
			
StdIP[{x_,Delta[dt_],Points[n_?Positive]},dp_]:=
	Module[{st0,std,texto},
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[dt,dp];
			texto =  StdIPTextMP[{st0,std}], 
			st0 = StringNumber[x];
			std =  StringNumber[dt];
			texto =  StdIPTextDP[{st0,std}]];
		{n+1, texto}
	]

StdIP[{x_,Points[n_?Positive],Delta[dt_]},dp_]:=
	Module[{st0,std,texto},
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[dt,dp];
			texto =  StdIPTextMP[{st0,std}], 
			st0 = StringNumber[x];
			std =  StringNumber[dt];
			texto =  StdIPTextDP[{st0,std}]];
		{n+1, texto}
	]

StdIP[{x_,y_,Points[n_?Positive]},dp_]:=
	Module[{st0,std,texto, dt},
		dt = N[(y-x)/n,dp]; 
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[dt,dp];
			texto =  StdIPTextMP[{st0,std}], 
			st0 = StringNumber[x];
			std =  StringNumber[dt];
			texto =  StdIPTextDP[{st0,std}]];
		{n+1, texto}
	]

StdIP[{x_,y_,Delta[dt_]}, dp_]:=
	Module[{st0,std,texto,n},
		n = Floor[Abs[(y-x)/dt]];
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[dt,dp];
			texto =  StdIPTextMP[{st0,std}], 
			st0 = StringNumber[x];
			std =  StringNumber[dt];
			texto =  StdIPTextDP[{st0,std}]];
		{n+1, texto}
	]
	

StdIP[x_List, dp_]:=
	{Length[x],StdIPText[Length[x],x,dp]}

StdIP[x_, dp_]:= Null



ValueMParams[params_, dp_]:=
(
	$ParametersValue = params[[1]]; 
	$InitialConditions = params[[2]] ;
	$IntInt = StdIP[params[[3]], dp] ;
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


(* ::Subsection::Closed:: *)
(*C\[OAcute]digo C*)


(* ::Subsubsection::Closed:: *)
(*Repeticiones codigo C*)


separador =
"/*********************************************************************************************/"


headFunction1C[name_?StringQ]:=
	"long  "<> name <> 
	"(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1])\n{\n"


headFunction2C[name_?StringQ, n_Integer]:=
	"long  "<> name <> 
	"_columns()\n{\n\t return "<>
	ToString[n+1]<> ";\n}\n\n";

headFunction3C[name_?StringQ, ldp_]:=
	"long  "<> name <> 
	"_pos_der(char *der)\n{\n"<>
	ListStringsToStringC["STR_DER", ldp]<>
	"\tlong i;\n"<>
	"\tfor(i=0; i < "<> ToString[Length[ldp]]<>"; i++)\n"<>
	"\t\tif(strcmp(der,STR_DER[i]) == 0) return i;\n"<>
	"\treturn -1;\n}\n\n";

headFunction4C[name_?StringQ]:=
	"long  "<> name <> 
	"_variable_column(int v, char *der)\n{\n\t return "<>
	"position_variable(v, "<> name <>"_pos_der, der);\n}\n\n"<>
	"long  "<> name <> 
	"_function_column(int f, char *der)\n{\n\t return "<>
	"position_function(f, "<> name <>"_pos_der, der);\n}\n\n\n"


numColumns[]:=
	"\tif(ORDER < 0) return NUM_COLUMNS;\n\n" 


setFirstC1 = 
	"\tif(NOT_INITIALIZED)\n\t{\n"<>
	"\t\tset_iterations();\n"<>
	"\t\tNOT_INITIALIZED = 0; \n\t}\n"<>
	"\tset_max_order(ORDER);\n\n"<>
	"\trealNUM var[VARIABLES+1][NUM_DERIVATIVES][ORDER+1];\n"<>
	"\trealNUM par[PARAMETERS][NUM_DERIVATIVES][ORDER+1];\n"<>
	"\trealNUM link[LINKS][NUM_DERIVATIVES][ORDER+1];\n"<>
	"\tvariables_init(var,v,t);\n"<>
	"\tparameters_init(par,p);\n"<>
	"\tlinks_init(link);\n"<>
	"\tderivatives_init(var,par,v);\n\n"<>
	"\tint i;\n"
	

setFirstC2 := 
	"\tfor(i=0;  i<=ORDER; i++) {";


setLast1C = 
	"\n\t}\n\n"<>
	"\twrite_solution(cvfd,var,link);\n"

setLast2C = 
	"\n\treturn NUM_COLUMNS;\n}\n\n";


includehC[name_?StringQ]:= "#include \""<> name <>".h\"\n"
includeSystemhC[name_?StringQ]:= "#include <"<> name <>">\n"


EndMessage[name_?StringQ]:= 
	"Files \"" <> name <> ".h\" and \"" <> name <> ".c\" written on directory \""<> Directory[]<> "\"."


(* ::Subsubsection::Closed:: *)
(*Listas y constantes en C*)


ConstantStringC[tipo_?StringQ, name_?StringQ, n_Integer]:=
	Module[{texto = "\tstatic "},
		texto = texto <> tipo <> " " <> name;
		texto = texto <>  " = " <> ToString[n]<> ";\n";
		texto]


ListIntegersToStringC[tipo_?StringQ, name_?StringQ,{}]:=
	Module[{texto="\tstatic "},
		texto = texto <> tipo <>" "; 
		texto = texto <> name <>"["<>ToString[1]<>"] = {0};\n";
		texto]


ListIntegersToStringC[tipo_?StringQ, name_?StringQ,mu_List]:=
	Module[{texto="\tstatic ",lmu},
		lmu = Length[mu];
		texto = texto <> tipo <>" "; 
		texto = texto <> name <>"["<>ToString[lmu]<>"] = {";
		Do[texto = texto <> ToString[mu[[i]]]<>",", {i,lmu-1}];
		texto=  texto <> ToString[mu[[lmu]]]<> "};\n";
		texto]


ListStringsToStringC[name_?StringQ,mu_List]:=
	Module[{texto="\tstatic char* ",lmu},
		lmu = Length[mu];
		texto = texto <> name <>"["<>ToString[lmu]<>"] = {";
		Do[texto = texto <> "\""<>mu[[i]]<> "\""<>",", {i,lmu-1}];
		texto= texto <> "\""<>mu[[lmu]]<> "\""<>"};\n";
		texto]


(* ::Subsubsection::Closed:: *)
(*Variables e iteraciones C*)


StringNumber[x_]:= StringCNumber[x]


strNUMBERC[x_,  3999999999]:= 
	"\"****\""

strNUMBERC[x_,  nd_]:= 
	"\""<> ToString[CForm[N[x,nd]]]<>"\""


strTayC[LKF$Plus]:= "add_t";
strTayC[LKF$Minus]:= "sub_t";
strTayC[LKF$Times]:= "mul_t";
strTayC[LKF$Divide]:= "divide_t";
strTayC[LKF$Power]:= "pow_t";
strTayC[LKF$Sin]:= "sin_t";
strTayC[LKF$Cos]:= "cos_t";
strTayC[LKF$Tan]:= "tan_t";
strTayC[LKF$Sinh]:= "sinh_t";
strTayC[LKF$Cosh]:= "cosh_t";
strTayC[LKF$Tanh]:= "tanh_t";
strTayC[LKF$ArcSin]:= "asin_t";
strTayC[LKF$ArcCos]:= "acos_t";
strTayC[LKF$ArcTan]:= "atan_t";
strTayC[LKF$ArcSinh]:= "asinh_t";
strTayC[LKF$ArcCosh]:= "acosh_t";
strTayC[LKF$ArcTanh]:= "atanh_t";
strTayC[LKF$Log]:= "log_t";
strTayC[LKF$Exp]:= "exp_t";
strTayC[LKF$Der]:= "der_t";


strOBJC[LKF$Var[n_]]:= "var[" <> ToString[n-1] <> "]";
strOBJC[LKF$Par[n_]]:= "par[" <> ToString[n-1] <> "]";
strOBJC[LKF$Link[n_]]:= "link[" <> ToString[n-1] <> "]";
strOBJC[LKF$Const[n_]]:= "ct[" <> ToString[n-1] <> "]";


restaDepthDer[i_]:= 
	Module[{dpt, text},
		dpt = depthDer[LKF$Link[i]];
		If[dpt === 0, 
			text = "",
			text = "-" <> ToString[dpt] ];
		text
	]
strLITERC[ LKF$Sin[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCos[i,j,x]
strLITERC[ LKF$Cos[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCos[j,i,x]
strLITERCSinCos[i_,j_,x_]:=
	Module[{prev, sal},
		prev = MemberQ[SinCosLKFList$, {i,j}];
		If[prev , 
			sal = ""; 
			Drop[SinCosLKFList$, Position[SinCosLKFList$,{i,j}][[1]]],
			sal = strLITERCSinCosPrint[i,j,x]; 
			AppendTo[SinCosLKFList$,{i,j}]];
		sal
	]
strLITERCSinCosPrint[i_,j_,x_]:=
	"\n\t\tsincos_t(" <> strOBJC[x] <> "," <> strOBJC[LKF$Link[i]] <>
	"," <> strOBJC[LKF$Link[j]] <> ",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Sinh[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCosh[i,j,x]
strLITERC[ LKF$Cosh[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCosh[j,i,x]
strLITERCSinCosh[i_,j_,x_]:=
	Module[{prev, sal},
		prev = MemberQ[SinCoshLKFList$, {i,j}];
		If[prev , 
			sal = ""; 
			Drop[SinCoshLKFList$, Position[SinCoshLKFList$,{i,j}][[1]]],
			sal = strLITERCSinCoshPrint[i,j,x]; 
			AppendTo[SinCoshLKFList$,{i,j}]];
		sal
	]
strLITERCSinCoshPrint[i_,j_,x_]:=
	"\n\t\tsincosh_t(" <> strOBJC[x] <> "," <> strOBJC[LKF$Link[i]] <>
	"," <> strOBJC[LKF$Link[j]] <> ",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Divide[1,z_], i_ ]:= 
	"\n\t\tinv_t(" <> strOBJC[z] <> "," <> strOBJC[LKF$Link[i]] <> 
	",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Divide[x_?numberADQ,z_], i_ ]:= 
	"\n\t\tdivide_t_cv(" <> strOBJC[x] <> "," <> strOBJC[z] <> ","<> strOBJC[LKF$Link[i]] <> 
	",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Divide[x_,z_?numberADQ], i_ ]:= 
	"\n\t\tdivide_t_vc(" <> strOBJC[x] <> "," <> strOBJC[z] <> ","<> strOBJC[LKF$Link[i]] <> 
	",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Power[E, y_], i_ ]:= 
	"\n\t\texp_t("  <>  StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ LKF$Power[x_, y_?numberADQ], i_ ]:= 
	"\n\t\t" <> strTayC[LKF$Power] <> "_cc" <> "("  <> strOBJC[x] <> "," <> 
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ h_[x_?numberADQ, y__], i_ ]:= 
	"\n\t\t" <> strTayC[h] <> "_cc" <> "("  <> strOBJC[x] <> "," <>
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ LKF$Minus[y_,x_?numberADQ], i_ ]:= 
	strLITERC[ LKF$Plus[y,-x], i ]

strLITERC[ h_[y_,x_?numberADQ, k___], i_ ]:= 
	"\n\t\t" <> strTayC[h] <> "_cc" <> "("  <> strOBJC[x] <> "," <>
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y,k}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ h_[y___], i_ ]:= 
	"\n\t\t" <> strTayC[h] <> "("  <> 
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";
	
strVITERC[it_?numberADQ,i_]:=
	"\n\t\tvar_t_cc(" <> strOBJC[it] <> ",var[" <> ToString[i] <> "], i);"
	
strVITERC[it_,i_]:=
	"\n\t\tvar_t(" <> strOBJC[it] <> ",var[" <> ToString[i] <> "], i);"


textVARSC[lf_]:=
	StringJoin @@ Map[strVITERC[lf[[#]],#]&, Range[Length[lf]]];
	
textITERSC[it_]:=
	StringJoin @@ Map[strLITERC[it[[#]],#]&, Range[Length[it]]];


(* ::Subsubsection::Closed:: *)
(*Cadena C*)


StandardCText[name_?StringQ, oldfun_, n_Integer]:=
	Module[{texto = "",lkf,lkv,iter,pp,oder,tdl,cil,cils, dgt, nder,ncol},
		$digitsNUMTIDES = n;
		SinCosLKFList$ = {};
		SinCoshLKFList$ = {};
		lkf = ToTaylorLKF[oldfun];
		lkf = ExtractConstants[lkf];
		iter =  IterationLKF[lkf[[4]]];
		lkv  =  LinksVariables[lkf];
		{pp,oder} = If[oldfun[[5]]=={}, {{},Until[ NumberOfVariables[lkf],0]}, oldfun[[5]]];
		tdl   =  DerOutputList[oder];
		cil   = CompleteIteratorsList[oder];
		cils = CompleteIteratorsListStar[oder];
		ChangeIndexDer[];
		Map[depthDer,iter];
		nder = tdl[[1]];
		ncol = nder (NumberOfVariables[lkf] + NumberOfFunctions[lkf]);


		texto = texto <> gmecopyrightC[];
		texto = texto <> includehC[name]<>"\n\n";
		texto = texto <> headFunction1C[name];

		texto = texto <> ConstantStringC["int  " , "VARIABLES       ", NumberOfVariables[lkf]];
		texto = texto <> ConstantStringC["int  " , "PARAMETERS      ", NumberOfParameters[lkf]];
		texto = texto <> ConstantStringC["int  " , "FUNCTIONS       ", NumberOfFunctions[lkf]];
		texto = texto <> ConstantStringC["int  " , "LINKS           ", NumberOfLinks[lkf]];
		texto = texto <> ConstantStringC["int  " , "PARTIALS_VARS   ", NumberOfPartials[oldfun]];
		texto = texto <> ConstantStringC["long ", "NUM_DERIVATIVES ", tdl[[1]]];
		texto = texto <> ConstantStringC["long ", "NUM_COLUMNS     ", ncol];
		texto = texto <> "\n";

		texto = texto <> ListIntegersToStringC["int  ", "POS__PARTIALS",pp];
		texto = texto <> ListIntegersToStringC["int  ", "POS_FUNCTIONS",LinksFunctions[lkf]];
		texto = texto <> "\n";

		texto = texto <> ListIntegersToStringC["long ", "POS_ACCUM", cil[[1]]];
		texto = texto <> ListIntegersToStringC["long ", "POS_COEFS", cil[[2]]];
		texto = texto <> ListIntegersToStringC["long ", "POS_PREVI", cil[[3]]];
		texto = texto <> ListIntegersToStringC["long ", "POS_PREIV", cil[[4]]];
		texto = texto <> "\n";

		texto = texto <> ListIntegersToStringC["long ", "POS_ACCUM_S", cils[[1]]];
		texto = texto <> ListIntegersToStringC["long ", "POS_COEFS_S", cils[[2]]];
		texto = texto <> ListIntegersToStringC["long ", "POS_PREVI_S", cils[[3]]];
		texto = texto <> ListIntegersToStringC["long ", "POS_PREIV_S", cils[[4]]];
		texto = texto <> "\n\n";
		texto = texto <> numColumns[];	
		texto = texto <> ConstantStringC["int ", "NOT_INITIALIZED", 1];

		texto =  texto <> setFirstC1;
		If[n > 16, 
			texto = texto <> TextMPConstants[lkf[[3]]]<>"\n",
			texto = texto <> TextDPConstants[lkf[[3]]]<>"\n"];

		texto =  texto <> setFirstC2;
		texto =  texto <> textVARSC[lkv] ;
		texto =  texto <> textITERSC[lkf[[4]]] ;

		texto =  texto <> setLast1C;
		If[n > 16, 
			texto = texto <> "\tclear(var,par,link);\n";
			texto = texto <> TextMPConstants[]<>"\n"];
		texto =  texto <> setLast2C;
	
		texto = texto <> headFunction2C[name, ncol];
		texto = texto <> headFunction3C[name, tdl[[2]]];
		texto = texto <> headFunction4C[name];
		texto		
	]


StandardHText[name_?StringQ, n_Integer]:=
	Module[{texto =""}, 
		texto = texto <> gmecopyrightC[];
		If[n == 16, 
			texto = texto <> headstdDP,
			texto = texto <> headstdMP];
		texto = texto <> "\n\n"<> separador <> "\n\n";
		
		texto = texto <> "#ifndef "<> name <>"_tides_h\n";
		texto = texto <> "#define "<> name <>"_tides_h\n\n";

		texto = texto <> "long  "<> name <> "_columns();\n";
		texto = texto <> "long  "<> name <> 
			"(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1]);\n";
		texto = texto <> "long  "<> name <> "_pos_der(char *der);\n";
		texto = texto <> "long  "<> name <> "_variable_column(int v, char *der);\n";
		texto = texto <> "long  "<> name <> "_function_column(int v, char *der);\n\n";
		texto = texto <> "#endif\n\n\n";
		texto
	]



(* ::Subsubsection:: *)
(*Drivers*)


TextDate[]:=
	Module[{dt, texto, mes,min, tmin},
	mes = {"January ", "February ", "March ", 
	"April ","May ", "June ", "July ", "August ", 
	"September ", "October ","November ", "December "};
	dt = Date[];
	texto = "This file has been created by MathTIDES ("<> Global`mathTIDESVersion$ <> ") ";
	texto = texto <> mes[[dt[[2]]]] <>ToString[dt[[3]]] <> ", ";
	texto = texto <> ToString[dt[[1]]]<>", ";
	min = Round[dt[[5]]];
	tmin = If[min < 10, "0"<>ToString[min],  ToString[min]];
	texto = texto <> ToString[dt[[4]]] <> ":" <>tmin;
	texto
]


NameDataMatrix[name_]:= 
	If[Head[$DataMatrix] === String, $DataMatrix, name<> "_DataMatrix"]


bloqueDBC2[name_]:= 
	"#include <stdio.h>"<> 
	"\n#include <stdlib.h>"<> 
	"\n#include \""<> name <> ".h\""<>
	"\n\nint main() {" 

bloqueDBC3[v_,p_,nfun_, name_]:= 
	Module[{texto = "\n"},
	texto = texto <> "\n\tint nvar = " <> ToString[v] <> ";" ;
	texto = texto <> "\n\tint npar = " <> If[p > 0, ToString[p], "0"] <> ";" ;
	texto = texto <> "\n\tint nfun = " <> ToString[nfun] <> ";" ; 
	If[$IntInt === Null, 
		texto = texto <> "\n\tint nipt = ***** ;" , 
		texto = texto <> "\n\tint nipt = " <> ToString[$IntInt[[1]]] <> ";" ]; 
	texto = texto <> "\n\tdouble v[nvar]";
	If[p > 0,  texto = texto <> ", p[npar]"] ;
	texto = texto <> ", lt[nipt];";
	texto = texto <> "\n\tdouble tolrel, tolabs; ";
	If[$DataMatrix =!= False, texto= texto <>"\n\tdouble** "<> NameDataMatrix[name] <> ";"];
	texto = texto <> "\n\tFILE   *fd;\n" ;
	texto]


bloqueMPC2[name_]:= 
	"#include <stdio.h>"<> 
	"\n#include <stdlib.h>"<> 
	"\n#include \"mpfr.h\""<> 
	"\n#include \""<> name <> ".h\""<>
	"\n\nint main() {" 

bloqueMPC3[v_,p_,nfun_, name_, n_]:= 
	Module[{texto = "\n"},
	If[n === 3999999999, 
		texto = texto <> "\n\tset_precision_digits(****);",
		texto = texto <> "\n\tset_precision_digits(" <> ToString[n] <> ");" ];
	texto = texto <> "\n\n\tint i;" ; 
	texto = texto <> "\n\tint nvar = " <> ToString[v] <> ";" ;
	texto = texto <> "\n\tint npar = " <> If[p > 0, ToString[p], "0"] <> ";" ;
	texto = texto <> "\n\tint nfun = " <> ToString[nfun] <> ";" ; 
	If[$IntInt === Null, 
		texto = texto <> "\n\tint nipt = ***** ;" , 
		texto = texto <> "\n\tint nipt = " <> ToString[$IntInt[[1]]] <> ";" ]; 
	texto = texto <> "\n\tmpfr_t v[nvar]";
	If[p > 0,  texto = texto <> ", p[npar]"] ;
	texto = texto <> ", lt[nipt];";
	texto = texto <> "\n\tmpfr_t tolrel, tolabs; ";
	If[$DataMatrix =!= False, texto= texto <>"\n\tmpfr_t** "<> NameDataMatrix[name] <> ";"];
	texto = texto <> "\n\tFILE   *fd;\n" ;
	texto]


bloqueDBC4[]:=
	Module[{texto="\n", separador, cond = False},
		separador ="/***********************************************************/";
		If[$Factor1 =!= Null || $Factor2 =!= Null || $Factor3 =!= Null ||
			$MaxStepRatio =!= Null || $MinStepRatio =!= Null ||
			$MaxIterationsNumber =!= Null ||$ExcessOrder =!= Null  ||
			$MinOrder =!= Null ||$MaxOrder =!= Null ||$StepSizeEstimator === True, 
				texto = texto <> separador <> "\n";
				texto = texto <> separador <> "\n";
				texto = texto <> "/*       CONSTANTS OF THE METHOD                            */"<> "\n";
				texto = texto <> separador <> "\n";
				texto = texto <> separador <> "\n";
				cond = True];
		If[$Factor1 =!= Null,  texto = texto <> 
			"\n\textern double fac1;" <>
			"\n\tfac1 = "<>StringNumber[$Factor1]<>";"];
		If[$Factor2 =!= Null,  texto = texto <> "\n\tfac2 = " 
			"\n\textern double fac2;"  <>
			StringNumber[$Factor2]<>";"];
		If[$Factor3 =!= Null,  texto = texto  
			"\n\textern double fac3;"  <>
			"\n\tfac3 = "<>StringNumber[$Factor3]<>";"];
		If[$MaxStepRatio =!= Null,  texto = texto <> 
			"\n\textern double rmaxstep;"  <>
			"\n\trmaxstep = "<>StringNumber[$MaxStepRatio]<>";"];
		If[$MinStepRatio =!= Null,  texto = texto <> 
			"\n\textern double rminstep;"  <>
			"\n\trminstep = "<>StringNumber[$MinStepRatio]<>";"];
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



bloqueDBC5a[v_,p_]:= 
	Module[{texto="\n", separador},
		separador ="/************************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */"<> "\n";
		If[$ParametersValue === Null || $InitialConditions === Null || $IntInt === Null, 
			texto = texto <> "/*     Change *****  by numerical values if it is necesary  */"<> "\n";];
		texto = texto <> separador <> "\n";
		texto = texto <> separador ;
		If[p > 0,
			texto = texto <> "\n\n/* --- PARAMETERS VALUE --- */";
			Map[(texto = texto <> "\n\tp[" <> ToString[#] <>
				 "] = "<> 
				If[$ParametersValue === Null, "*****", StringNumber[$ParametersValue[[#+1]]]] <>
				" ; ")&, Range[0,p-1]]];
		texto = texto <> "\n\n/* --- INITIAL VALUES --- */";
		Map[(texto = texto <>  "\n\tv[" <> ToString[#] <>
				 "] = "<> 
				If[$InitialConditions === Null, "*****", StringNumber[$InitialConditions[[#+1]]]] <>
				" ; ")&, Range[0,v-1]];
		texto
	]


bloqueMPC5a[v_,p_,n_]:= 
	Module[{texto="\n", separador},
		separador ="/************************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */"<> "\n";
		If[$ParametersValue === Null || $InitialConditions === Null || $IntInt === Null, 
			texto = texto <> "/*     Change *****  by numerical values if it is necesary */"<> "\n";];
		texto = texto <> separador <> "\n";
		texto = texto <> separador ;
		If[p > 0,
			texto = texto <> "\n\n/* --- PARAMETERS VALUE --- */";
			texto = texto <> "\n\tfor(i=0; i<npar; i++) mpfrts_init(&p[i]);";
			Map[(texto = texto <> "\n\tmpfrts_set_str(&p[" <> ToString[#] <> "], "<> 
				If[$ParametersValue === Null, "\"*****\"", strNUMBERC[$ParametersValue[[#+1]],n]] <>
				"); ")&, Range[0,p-1]]];
		texto = texto <> "\n\n/* --- INITIAL VALUES --- */";
		texto = texto <> "\n\tfor(i=0; i<nvar; i++) mpfrts_init(&v[i]);";
		Map[(texto = texto <>  "\n\tmpfrts_set_str(&v[" <> ToString[#] <> "], "<> 
				If[$InitialConditions === Null, "\"*****\"", strNUMBERC[$InitialConditions[[#+1]],n]] <>
				"); ")&, Range[0,v-1]];
		texto
	]


bloqueDBC5b[]:= 
	Module[{texto=""},
		texto = texto <> "\n\n/* ---     INTEGRATION POINTS    --- */";
		If[$IntInt === Null,
			texto = texto <>  "\n\tlt[0] = ***** ;";
			texto = texto <>  "\n\tlt[1] = ***** ;";
			texto = texto <>  "\n\t************* ;",
			texto = texto <> $IntInt[[2]]];
		texto
	]



bloqueMPC5b[n_]:= 
	Module[{texto=""},
		texto = texto <> "\n\n/* ---     INTEGRATION POINTS    --- */";
		texto = texto <> "\n\tfor(i=0; i<nipt; i++) mpfrts_init(&lt[i]);";
		If[$IntInt === Null,
			texto = texto <>  "\n\tmpfrts_set_str(&lt[0], \"*****\");";
			texto = texto <>  "\n\tmpfrts_set_str(&lt[1], \"*****\");";
			texto = texto <>  "\n\t****************************** ;",
			texto = texto <> $IntInt[[2]]];
		texto
	]



bloqueDBC5c[]:= 
	Module[{texto="", trel, tabs},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			tabs = trel = N[10^-16]];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = N[$AbsoluteTolerance,16]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = N[$RelativeTolerance,16]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = N[$RelativeTolerance,16]; tabs = N[$AbsoluteTolerance,16]];
		texto = texto <> "\n\n/* --- REQUIRED TOLERANCES --- */";
		texto = texto <> "\n\ttolrel = "<> StringNumber[trel]<>" ;";
		texto = texto <> "\n\ttolabs = "<> StringNumber[tabs]<>" ;";
		texto
	]


bloqueMPC5c[n_]:= 
	Module[{texto="", trel="", tabs=""},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			If[n === 3999999999, 
				tabs = trel = strNUMBERC[0,n],
				tabs = trel = strNUMBERC[10^-(n-1),n]]];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = strNUMBERC[$AbsoluteTolerance,n]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = strNUMBERC[$RelativeTolerance,n]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = strNUMBERC[$RelativeTolerance,n]; 
			tabs = strNUMBERC[$AbsoluteTolerance,n]];
		texto = texto <> "\n\n/* --- REQUIRED TOLERANCES --- */";
		texto = texto <> "\n\tmpfrts_init(&tolrel); "; 
		texto = texto <> "\n\tmpfrts_init(&tolabs); "; 
		texto = texto <> "\n\tmpfrts_set_str(&tolrel, "<> trel <>");";
		texto = texto <> "\n\tmpfrts_set_str(&tolabs, "<> tabs <>");";
		texto
	]


bloqueDBC5d[name_]:= 
	Module[{texto="\n\n", separador , donde},	
		separador =      "/***********************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*        OUTPUT: ";
		donde = If[$FileOutput === False, "        ", If[$FileOutput === Screen, "screen &", "file   &"]];
		donde = donde  <> If[$DataMatrix === False, "             ", "data matrix"];
		texto = texto <> donde <> "                     */\n";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		If[$FileOutput === Screen , 
			texto = texto <> "\n\tfd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tfd = fopen(\""<>$FileOutput<>"\", \"w\");"];
		If[$DataMatrix =!= False, 
			texto =texto <> "\n\tArray2DB_init(&"<> NameDataMatrix[name]<>", ";
			texto =texto <> If[$IntInt === Null, ToString[2],"nipt"];
			texto =texto <> ", "<> name <>"_columns());\n"];
		texto	
	]


bloqueMPC5d[name_]:= 
	Module[{texto="\n\n", separador , donde},	
		separador =      "/***********************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*        OUTPUT: ";
		donde = If[$FileOutput === False, "        ", If[$FileOutput === Screen, "screen &", "file   &"]];
		donde = donde  <> If[$DataMatrix === False, "             ", "data matrix"];
		texto = texto <> donde <> "                     */\n";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		If[$FileOutput === Screen , 
			texto = texto <> "\n\tfd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tfd = fopen(\""<>$FileOutput<>"\", \"w\");"];
		If[$DataMatrix =!= False, 
			texto =texto <> "\n\tArray2MP_init(&"<> NameDataMatrix[name]<>", ";
			texto =texto <> If[$IntInt === Null, ToString[2],"nipt"];
			texto =texto <> ", "<> name <>"_columns());\n"];
		texto	
	]


bloqueDBC6b[name_,p_]:= 
	Module[{texto="\n\n", separador},
		separador ="/***********************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*       CALL THE INTEGRATOR                               */"<> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "\n\tdp_tides(" <> name <>", nvar, npar, nfun, v, ";
		texto = texto <> If[p === 0, "NULL, ", "p, "] ;
		texto = texto <> "\n\t\t\tlt, nipt, tolrel, tolabs, "; 
		texto = texto <> 
			If[$DataMatrix === False, "NULL", NameDataMatrix[name]] <> ", ";
		texto = texto <> 
			If[$FileOutput === False, "NULL", "fd"] <> ");\n";
		texto]


bloqueMPC6b[name_,p_]:= 
	Module[{texto="\n", separador},
		separador ="/***********************************************************/";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "/*       CALL THE INTEGRATOR                               */"<> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> separador <> "\n";
		texto = texto <> "\n\tmp_tides(" <> name <>", nvar, npar, nfun, v, ";
		texto = texto <> If[p === 0, "NULL, ", "p, "] ;
		texto = texto <> "\n\t\t\tlt, nipt, tolrel, tolabs, "; 
		texto = texto <> 
			If[$DataMatrix === False, "NULL", NameDataMatrix[name]] <> ", ";
		texto = texto <> 
			If[$FileOutput === False, "NULL", "fd"] <> ");\n";
		texto]


bloqueDBC6d[]:= 
	Module[{texto="\n"},
		If[Head[$FileOutput] === String,
			texto = texto  <>  "\n\tfclose(fd); \n"];
		texto = texto <> "\n\treturn 0;";
		texto = texto <> "\n}\n\n\n";
		texto]


DriverDB[name_?StringQ, v_, p_, nf_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["dp"];
		texto = texto <> bloqueDBC2[name];
		texto = texto <> bloqueDBC3[v,p, nf,name];
		texto = texto <> bloqueDBC4[];
		texto = texto <> bloqueDBC5a[v,p];
		texto = texto <> bloqueDBC5b[];
		texto = texto <> bloqueDBC5c[];
		texto = texto <> bloqueDBC5d[name];
		texto = texto <> bloqueDBC6b[name,p];
		texto = texto <> bloqueDBC6d[];
		texto
	]


DriverMP[name_?StringQ, n_, v_, p_, nf_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["mp"];
		texto = texto <> bloqueMPC2[name];
		texto = texto <> bloqueMPC3[v,p, nf,name, n];
		texto = texto <> bloqueDBC4[];
		texto = texto <> bloqueMPC5a[v,p,n];
		texto = texto <> bloqueMPC5b[n];
		texto = texto <> bloqueMPC5c[n];
		texto = texto <> bloqueMPC5d[name];
		texto = texto <> bloqueMPC6b[name,p];
		texto = texto <> bloqueDBC6d[];
		texto
	]


DriverStdC[name_?StringQ, n_, fun_, params_List]:=
	Module[{v, p, np, nf, texto},
		ValueMParams[params,n]; 
		adf = ToTaylorLKF[fun];
		v = adf[[1]]-1;
		p = adf[[2]];
		nf = NumberOfFunctions[adf];
		texto = If[n >16, DriverMP[name,n,v,p, nf], DriverDB[name,v,p, nf]];
		CleanMParams[];
		texto
	]




(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`StandardCode`"]

EndPackage[]

Null
