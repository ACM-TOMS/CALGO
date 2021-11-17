(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`Codes:  C and Fortran Codes*)


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
(*Codes*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`Codes`",
				"MathTIDES`StandardCode`",
					"MathTIDES`MinimalCode`",
						"MathTIDES`Iterations`",
							"MathTIDES`ODES`",
								"MathTIDES`Texts`",
								    "MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
	TSMCodeFiles,
	CodeFiles,
	PrecisionDigits,
	MinTIDES,
	Driver,
	OnlyDriver, 
	ParametersValue,
	InitialConditions,
	IntegrationPoints,
	Output,
	DataMatrix,
	Factor1,
	Factor2,
	Factor3,
	MaxStepRatio,
	MinStepRatio,
	MaxIterationsNumber,
	OrderIncrement,
	MinOrder,
	MaxOrder,
	RelativeTolerance,
	AbsoluteTolerance,
	DefectErrorControl,
	AddFunctions,
	AddPartials,
	Optimization,
	Readmcc,
	Readmch,
	Readmff
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`Codes`*"]
Clear @@ Names["MathTIDES`Codes`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


TSMCodeFiles::"badopt"="Bad option"
TSMCodeFiles::"minFun"="Minimal TIDES version does not check functions"

TSMCodeFiles::usage = 
ToString[Style["TSMCodeFiles ",Bold],StandardForm] <>
"creates the C or Fortran code with the Taylor Series Integrator."<>
" The arguments of TSMCodeFiles[] are:"<>
ToString[Style["\n\nFirst argument: ",Bold,Italic],StandardForm]<>
"the first order differential equation. This is an expression with head FirstOrderODE$"<>
" created by one of the previously following expressions: "<>
"FirstOrderODE[], NthOrderODE[], HamiltonianToODE[] or PotentialToODE[] "<>
ToString[Style["\n\nSecond argument: ",Bold,Italic],StandardForm]<>
"an string that represents name of the files. "<>
"With this name MathTIDES writes several files "<>
"(depending on the options) with extension \".h\",\".c\" or \".f\""<>
ToString[Style["\n\nOptions: ",Bold,Italic],StandardForm]<>
"to modify the cade generated by TSMCodeFiles. The possible options of TSMCodeFiles are: "<>
"MinTIDES, PrecisionDigits, Driver, OnlyDriver, "<>
"InitialConditions, ParametersValue, IntegrationPoints, "<>
"RelativeTolerance, AbsoluteTolerance, Output, DataMatrix, "<>
"Factor1, Factor2, Factor3, MaxStepRatio, MinStepRatio, "<>
"MaxIterationsNumber, OrderIncrement, MinOrder and DefectErrorControl."

CodeFiles::usage = 
ToString[Style["CodeFiles ",Bold],StandardForm] <>
"creates the C or Fortran code with the Taylor Series Integrator."<>
" The arguments of CodeFiles[] are:"<>
ToString[Style["\n\nFirst argument: ",Bold,Italic],StandardForm]<>
"the first order differential equation. This is an expression with head FirstOrderODE$"<>
" created by one of the previously following expressions: "<>
"FirstOrderODE[], NthOrderODE[], HamiltonianToODE[] or PotentialToODE[] "<>
ToString[Style["\n\nSecond argument: ",Bold,Italic],StandardForm]<>
"an string that represents name of the files. "<>
"With this name MathTIDES writes several files "<>
"(depending on the options) with extension \".h\",\".c\" or \".f\""<>
ToString[Style["\n\nOptions: ",Bold,Italic],StandardForm]<>
"to modify the cade generated by CodeFiles. The possible options of CodeFiles are: "<>
"MinTIDES, PrecisionDigits, Driver, OnlyDriver, "<>
"InitialConditions, ParametersValue, IntegrationPoints, "<>
"RelativeTolerance, AbsoluteTolerance, Output, DataMatrix, "<>
"Factor1, Factor2, Factor3, MaxStepRatio, MinStepRatio, "<>
"MaxIterationsNumber, OrderIncrement, MinOrder and DefectErrorControl."

IntegrationPoints::usage = 
"With this option we declare, on the driver, the list of points in which the solution is computed."<>
" There are several versions of this option: "<>
ToString[Style["\n\nIntegrationPoints -> {t0, t1, ..., tf}",Bold],StandardForm] <>
"\n\t*) t0 is the initial integration point (where the initial conditions are given)."<>
"\n\t*) t1,...,tf are the points where we want to compute the solution."<>
"\n\t*) tf is the final integration point."<>
"\n\t*) This option is only valid for the standard version."<>
"\n\t*) In minimal version you can use this option, with the initial to and final tf point only, "<>
"for non-dense output."<>
"\n\t*) t0,t1,...,tf are in order (crescent or decrescent). They can be non-equidistant points."<>
ToString[Style["\n\nIntegrationPoints -> {t0, tf, Delta[dt]}",Bold],StandardForm] <>
"\n\t*) t0 is the initial integration point (where the initial conditions are given)."<>
"\n\t*) tf is the final integrarion point."<>
"\n\t*) dt is the interval between points in dense output."<>
"\n\t*) The solution is computed in {t0, t0+dt,...,t0+k*dt}, where k such us: t0+k*dt <= tf < t0 + (k + 1)*dt"<>
"\n\t*) Not always the last point of the dense output coincides with the end integration point tf."<>
ToString[Style["\n\nIntegrationPoints -> {t0, tf, Points[k]}",Bold],StandardForm] <>
"\n\t*) t0 is the initial integration point (where the initial conditions are given)."<>
"\n\t*) tf is the final integrarion point."<>
"\n\t*) k is an integer with the number of equidistant points in which the solution is computed. "<>
"\n\t*) dt for dense output is equal to (tf - t0)/k."<>
"\n\t*) The solution is computed in {t0, t0+dt,...,t0+k*dt}."<>
ToString[Style["\n\nIntegrationPoints -> {t0, Delta[dt], Points[k]}",Bold],StandardForm] <>
"\n\t*) t0 is the initial integration point (where the initial conditions are given)."<>
"\n\t*) dt is the interval between points in dense output."<>
"\n\t*) k is an integer with the number of equidistant points in which the solution is computed. "<>
"\n\t*) The solution is computed in {t0, t0+dt,...,t0+k*dt}."




MinTIDES::usage = 
ToString[Style["MinTIDES ",Bold],StandardForm] <>
"is used to create files to use with the minimum version of TIDES. Use "<>
ToString[Style["MinTIDES ->\"C\" ",Bold], StandardForm]<>
" to create the C minimum version minc-tides  and "<>
ToString[Style["MinTIDES ->\"Fortran\" ",Bold], StandardForm]<>
" to create the FORTRAN minimum version minf-tides."

PrecisionDigits::usage = 
"By default, when the option "<>
ToString[Style["MinTIDES ",Bold], StandardForm]<>
" is not used an standard version is created.  We choose between "<>
ToString[Style["dp-tides",Bold], StandardForm]<>" or "<>
ToString[Style["mp-tides",Bold], StandardForm]<>
" by means of the option "<>ToString[Style["PrecisionDigits",Bold], StandardForm]<>".
\nBy default this option has the value "<>
ToString[Style["PrecisionDigits -> 16 ",Bold], StandardForm]<>
" This means that the standard double precision version dp-tides is created."<>
" With a number greater than 16 this option declares the number of digits of "<>
"precision of the integrator and creates the multiple precision version "<>
ToString[Style["mp-tides",Bold], StandardForm]<>"." 

Driver::usage = 
"By default a driver with the main program is created. "<>
ToString[Style["CodeFiles",Bold], StandardForm]<>
" does not write a driver if we use the option "<>
ToString[Style["Driver -> False",Bold], StandardForm]<>
", but it writes the rest of the files."

OnlyDriver::usage = 
"The option OnlyDriver -> True creates only the driver with the main program, "<>
"and no other file. "
 
ParametersValue::usage = "With the option "<>
ToString[Style["ParametersValue -> {0.1, -2.3, ...} ",Bold], StandardForm]<>
" we change, on the driver, the value of the parameters. "<>
"The length of the list must be equal to the number of parameters."<>
" If we do not use this options stars, ******, instead of values appear on the driver. "

InitialConditions::usage = "With the option "<>
ToString[Style["InitialConditions -> {0.1, -2.3, ...} ",Bold], StandardForm]<>
" we change, on the driver,  the initial value of the vector of variables."<>
"  The length of the list must be equal to the number of variables."<>
" If we do not use this options stars, ******, instead of values appear on the driver."

Output::usage = 
"This options declares where the solution (dense or not) is written. "<>
"There are two posiblilities: "<>
ToString[Style["Output -> Screen",Bold], StandardForm]<>" or  "<>
ToString[Style["Output -> \"file\"",Bold], StandardForm]<>
". In the first case the solution is written on the screen, 
in the second case into a file named file. By default ("<>
ToString[Style["Output -> False",Bold], StandardForm]<>
") no output is written. In the Minimal version if the output is not 
sending into the screen the solution in t0 and"<>
" the solution in tf is written on the screen. "

DataMatrix::usage = 
"Option only for standard version. By default "<>
ToString[Style["DataMatrix->False",Bold], StandardForm]<>
", but there are two other posibilities: "<>
ToString[Style["DataMatrix->True",Bold], StandardForm]<>" or "<>
ToString[Style["DataMatrix->\"nameDM\"",Bold], StandardForm]<>
". DataMatrix declares a bidimensional array where the solution is stored."<> 
"The name is nameDM in the second case or the name of the file joined to \"_DataMatrix\" "<>
"in the first case. Each row corresponds to the solution in the point ti of the "<>
"integration interval. The first row represents the initial point. "<>
"The last row represents the final point. The number of columns "<>
"is sufficient to store ti, the  variables in ti, "<>
"the functions in ti, the partial derivatives of variables and functions in ti."

Factor1::usage = 
ToString[Style[Factor1,Bold], StandardForm]<>" changes the parameter fac1."

Factor2::usage = 
ToString[Style[Factor2,Bold], StandardForm]<>" changes the parameter fac2."

Factor3::usage =  
ToString[Style[Factor3,Bold], StandardForm]<>" changes the parameter fac3."

MaxStepRatio::usage = 
ToString[Style[MaxStepRatio,Bold], StandardForm]<>" changes the parameter rmaxstep."

MinStepRatio::usage = 
ToString[Style[MinStepRatio,Bold], StandardForm]<>" changes the parameter rminstep."

MaxIterationsNumber::usage = 
ToString[Style[MaxIterationsNumber,Bold], StandardForm]<>" changes the parameter nitermax."

OrderIncrement::usage = 
ToString[Style[OrderIncrement,Bold], StandardForm]<>" changes the parameter nordinc."

MinOrder::usage = ToString[Style[MinOrder,Bold], StandardForm]<> " changes the parameter minord."

RelativeTolerance::usage = 
"Declares the value of the relative tolerance in the application of  the method. "<>
ToString[Style["RelativeTolerance -> rtol",Bold], StandardForm]<>
", where rtol is a real numbers. The default value is 10^-p, where "<>  
"p  is the value of the option PrecisionDigits. "<> 
"If only this tolerance is declared both (relative and absolute) are taken equal."

AbsoluteTolerance::usage = 
"Declares the value of the absolute tolerance in the application of  the method. "<>
ToString[Style["AbsoluteTolerance -> atol",Bold], StandardForm]<>
", where atol is a real numbers. The default value is 10^-p, where "<>
"p  is the value of the option PrecisionDigits. "<> 
"If only this tolerance is declared both (relative and absolute) are taken equal."

DefectErrorControl::usage = ToString[Style[DefectErrorControl -> True,Bold], StandardForm]<>
" selects the use of the defect error control in the integrator. By default "<> 
ToString[Style[DefectErrorControl -> False,Bold], StandardForm]<> "."



AddFunctions::usage =
"The integration of the system "<>
"\!\(\*FormBox[
RowBox[{FractionBox[
RowBox[{\"\[DifferentialD]\", \" \", OverscriptBox[\"x\", \"_\"], \" \"}], 
RowBox[{\"\[DifferentialD]\", \"t\"}]], \"=\", \" \", 
RowBox[{OverscriptBox[\"F\", \"_\"], \"(\", 
RowBox[{\"t\", \",\", \" \", OverscriptBox[\"x\", \"_\"], \",\", \" \", OverscriptBox[\"p\", \"_\"]}], \" \", \")\"}]}],
TraditionalForm]\) , t\[Element]\[DoubleStruckCapitalR], \!\(\*FormBox[
RowBox[{
RowBox[{OverscriptBox[\"x\", \"_\"], \"\[Element]\", SuperscriptBox[\"\[DoubleStruckCapitalR]\", \"n\"]}], \" \", \",\", \" \", 
RowBox[{OverscriptBox[\"p\", \"_\"], \" \", \"\[Element]\", SuperscriptBox[\"\[DoubleStruckCapitalR]\", \"m\"], \" \"}]}],
TraditionalForm]\)"<>
"gives the function "<> "\!\(\*OverscriptBox[\"x\", \"_\"]\)\!\(\*
StyleBox[\"(\",\nFontFamily->\"Courier\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"t\",\nFontFamily->\"Courier\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontFamily->\"Courier\",\nFontSlant->\"Italic\"]\)"<> 
"i.e. the evolution over time of the variables. "<>
"Sometimes we are interested  in the evolution, along "<>
"the solution of the system, of a dynamical variable defined by a function" <> "\!\(\*FormBox[
RowBox[{\"G\", \"(\", 
RowBox[{\"t\", \",\", \" \", OverscriptBox[\"x\", \"_\"], \",\", \" \", OverscriptBox[\"p\", \"_\"]}], \" \", \")\"}],
TraditionalForm]\)"<>
", i.e. the funtion" <> "G(t) = G\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[OverscriptBox[\"x\", \"_\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"(\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[OverscriptBox[\"p\", \"_\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\")\",\nFontSlant->\"Italic\"]\)." <>
"Writting the option "<>
ToString[Style["AddFunctions -> {G1, G2, ...} ",Bold],StandardForm]<>
"we redefine the differential equation to extend the application of the "<> 
"Taylor method to find the time evolution of the functions G1,G2,...\n\n"<>
ToString[Style["Example: ",Bold],StandardForm]<>"With " <>
ToString[Style["FirstOrderODE[{y,-x},t,{x,y},AddFunctions->{x/y}] ",Bold],StandardForm]<> 
"we compute the evolution over the time of x, y and x/y (sin, cos and tan)"


AddPartials::usage =
"Together with the time evolution of the variables and functions we may compute the evolution of the partials of the variables ( and partials of the functions) with respect to the initial conditions and with respec to to the parametes.The option to do that has four possible  formats: \n"<>
ToString[Style["\nAddPartials ->{{u,v,...}, n} ",Bold],StandardForm]<>
ToString[Style["\nAddPartials ->{{u,v,..}, n, Until}",Bold],StandardForm] <>
ToString[Style["\nAddPartials ->{{u,v,..}, n, Only}",Bold],StandardForm]<>
ToString[Style["\nAddPartials ->{{u,v,..}, listOfOrders}",Bold],StandardForm]<> "\n\n"<>
"The list "<>
ToString[Style["{u,v,...}",Bold],StandardForm]<>
" represents the symbols of the elements with respect to we want the derivatives. The symbols of this list are symbols of the variables or symbols of the parameters. If the symbol corresponds to a variable the partials with respect to the initial value of this variables computed. If the symbol correspond to a parameter the partial with respect to the parameter is computed. "<>
ToString[Style["n",Bold],StandardForm]<>
" represents the total maximum order  of the partials to compute. If no third argument appear (or the third argument is the symbol "<>
ToString[Style["Until",Bold],StandardForm]<>
" , all the partials until total order  n are computed. If the third argument is the symbol "<>
ToString[Style["Only",Bold],StandardForm]<>
", only the partial derivatives of total order n are computed. If the second argument, "<>
ToString[Style["listOfOrders",Bold],StandardForm]<>
", is a list, only the partials of the orders in the list are computed."<>
ToString[Style["\n\nExample: ",Bold],StandardForm]<>
"Let's assume a differential equation with two variables x,y,z and three parameters a,b: "<>
ToString[Style["\n\nAddPartials ->{{y,a}, 2}",Bold],StandardForm]<>
"\ncomputes:  \!\(\*FormBox[FractionBox[
RowBox[{\"\[PartialD]\", \"x\"}], 
RowBox[{\"\[PartialD]\", SubscriptBox[\"y\", \"o\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{\"\[PartialD]\", \"y\"}], 
RowBox[{\"\[PartialD]\", SubscriptBox[\"y\", \"o\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{\"\[PartialD]\", \"z\"}], 
RowBox[{\"\[PartialD]\", SubscriptBox[\"y\", \"o\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{\"\[PartialD]\", \"x\"}], 
RowBox[{\"\[PartialD]\", \"a\"}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{\"\[PartialD]\", \"y\"}], 
RowBox[{\"\[PartialD]\", \"a\"}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{\"\[PartialD]\", \"z\"}], 
RowBox[{\"\[PartialD]\", \"a\"}]],
TraditionalForm]\), \!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"x\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"x\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"x\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", \"a\"}]}]],
TraditionalForm]\), \!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"y\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"y\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"y\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", \"a\"}]}]],
TraditionalForm]\), \!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"z\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"z\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]],
TraditionalForm]\),\!\(\*FormBox[FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"z\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", \"a\"}]}]],
TraditionalForm]\)"<>
ToString[Style["\nAddPartials ->{{y,a}, 2, Only}",Bold],StandardForm]<>
"\ncomputes: \!\(\*FormBox[
RowBox[{FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"x\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]], \",\", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"x\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]], \",\", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"x\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", \"a\"}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"y\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]], \",\", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"y\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]], \",\", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"y\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", \"a\"}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"z\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]], \",\", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"z\"}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]], \",\", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"2\"], \"z\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", \"a\"}]}]]}],
TraditionalForm]\)"<>
ToString[Style["\nAddPartials ->{{y,a}, 2, {{2,3},{1,2}}}",Bold],StandardForm]<>
"\ncomputes: \!\(\*FormBox[
RowBox[{FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"5\"], \"x\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"3\"]}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"3\"], \"x\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"5\"], \"y\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"3\"]}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"3\"], \"y\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"5\"], \"z\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SuperscriptBox[SubscriptBox[\"y\", \"o\"], \"2\"]}]}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"3\"]}]}]], \",\", \" \", FractionBox[
RowBox[{SuperscriptBox[\"\[PartialD]\", \"3\"], \"z\"}], 
RowBox[{
RowBox[{\"\[PartialD]\", 
RowBox[{\"\[InvisiblePrefixScriptBase]\", SubscriptBox[\"y\", \"o\"]}]}], 
RowBox[{\"\[PartialD]\", SuperscriptBox[\"a\", \"2\"]}]}]]}],
TraditionalForm]\)"


Optimization::usage =
"With the default option "<>ToString[Style["Optimization -> 1 ",Bold],StandardForm]<>
"mathTIDES uses the expression Simplify[] to simplify "<>
"the linked function used to apply the Taylor method to the ODE. With "<>
ToString[Style["Optimization -> 2 ",Bold],StandardForm]<>
"mathTIDES uses the FullSimplify[], and with "<>
ToString[Style["Optimization -> 0 ",Bold],StandardForm]<>
"no simplification is made. The option "<>
ToString[Style["Optimization -> 2 ",Bold],StandardForm]<>
"not ensure a drastic simplification  with respect the "<>
"default but,sometimes,it takes a very long time of computation."


End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*Mensajes*)


EndMessageC1[name_?StringQ, dr_, odr_]:= 
	Module[{texto},
		If[odr, texto  = "File ", texto  = "Files "];
		If[odr || dr, texto = texto <> "\"dr_" <> name <> ".c\""];
		If[!odr, If[dr, texto =texto <> ", "];
			texto = texto <> "\"" <> name <> ".h\" and ";
			texto = texto <> name <> ".c\""];
	    texto = texto <> " written on directory \""<> Directory[]<> "\".";
		texto]


EndMessageC2[name_?StringQ, dr_, odr_]:= 
	Module[{texto},
		If[odr, texto  = "File ", texto  = "Files "];
		If[odr || dr, texto = texto <> "\"dr_" <> name <> ".c\""];
		If[!odr, If[dr, texto =texto <> ", "];
			texto = texto <> "\"" <> name <> ".c\" , ";
			texto = texto <> "\"minc_tides.c\" and \"minc_tides.h\""];
	    texto = texto <> " written on directory \""<> Directory[]<> "\".";
		texto]


EndMessageF[name_?StringQ,  v_, p_, dr_,odr_]:=
	Module[{texto},
		If[odr, texto  = "File ", texto  = "Files "];
		If[odr || dr, texto = texto <> "\"dr_" <> name <> ".f\""];
		If[!odr, If[dr, texto =texto <> ", "];
			texto = texto <> "\"" <> name <> ".f\" and ";
			texto = texto <> "\"minc_tides.f\""];
	    texto = texto <> " written on directory \""<> Directory[]<> "\".";
		texto]


(* ::Subsection::Closed:: *)
(*Codigo*)


(* ::Subsubsection::Closed:: *)
(*TSMCodeFiles y CodeFiles*)


CodeFiles = TSMCodeFiles


(* ::Subsubsection::Closed:: *)
(*Auxiliares listas derivadas parciales*)


PositionVARPAR[var_List, par_List, der_]:=
	Module[{pv,pp, num},
		pv = Position[var, der];
		pp = Position[par,der];
		num = If[pv == {}, pp[[1,1]]+Length[var],pv[[1,1]]];
		num] 

DerivativesList[var_List, par_List, der_List]:=
	Module[{dv,dp, dx, ndv,ndp},
		dv = Intersection[var,der];
		dp =  Intersection[par,der];
		dx =Complement[der, Join[dv,dp]];
		If[Length[dx] != 0, Message[FirstOrderODE$::"badder"]; Abort[]];
		Map[PositionVARPAR[var,par,#]&,der]]


OnlyUntilList[{}]={}
OnlyUntilList[{a_List, n_, s_Symbol}]:= {a, s[Length[a], n]}
OnlyUntilList[{a_List, n_Integer}]:= {a, Until[Length[a],n]}
OnlyUntilList[{a_List, n_List}]:= {a, Only[Length[a], n]}


NewDerivativesList[var_List, par_List, pwrt_List]:=
	Module[{ldir, sal, dim}, 
		ldir = If[pwrt =!= {}, 
					DerivativesList[var,par,pwrt[[1]]]];
		sal = OnlyUntilList[If[pwrt == {}, {}, 
					pwrt/.{pwrt[[1]]->ldir}]];
		If[sal != {} && Head[sal[[2,2]]] == List, 
			dim = Union[Union[Map[Length,sal[[2,2]]]]];
			If[Length[dim] != 1 || sal[[2,1]] =!= dim[[1]], 
				Message[FirstOrderODE$::"badderord"]; Abort[]]];
		sal
	]


(* ::Subsubsection::Closed:: *)
(*Modificacion de ecuaciones de primer orden *)


ModifyFirstOrderODE$[x_, nfun_, dpar_, optz_]:=
  Module[{ffun, pwrt}, 
	ffun = Join[x[[1]],nfun]; 
	ffun = Switch[optz,
				0, ffun,
				1, Simplify[ffun],
				2, FullSimplify[ffun],
				_, ffun]; 
	pwrt = NewDerivativesList[x[[3]],x[[4]], dpar];
	FirstOrderODE$[ffun,x[[2]],x[[3]], x[[4]], pwrt]
]


(* ::Subsubsection::Closed:: *)
(*TSMCodeFiles*)


TSMCodeFiles/:
	Options[TSMCodeFiles]= {
		Precision -> Double,
		PrecisionDigits -> 16,
		MinTIDES -> False,
		Driver -> True,
		OnlyDriver -> False,
		ParametersValue -> Null,
		InitialConditions -> Null,
		IntegrationPoints -> Null,
		Output -> False,
		DataMatrix -> False,
		Factor1 -> Null,
		Factor2 -> Null,
		Factor3 -> Null,
		MaxStepRatio->Null,
		MinStepRatio->Null,
		MaxIterationsNumber->Null,
		OrderIncrement->Null,
		MinOrder->Null,
		MaxOrder->Null,
		RelativeTolerance->Null,
		AbsoluteTolerance->Null,
		DefectErrorControl->False,
		AddFunctions->{},
		AddPartials->{},
		Optimization->1
};


FirstOrderODE$ /: 
	TSMCodeFiles[x_FirstOrderODE$, name_?StringQ, opt___Rule]:=
		Module[{nx, pre, pd, vers, dr, odr, nvers = 100, nf, 
					sigo = True, params, nfun, dpar, optz, nd},
			{pre, pd, vers, dr, odr, params, nfun, dpar, optz}=
				{Precision, PrecisionDigits, MinTIDES, Driver, OnlyDriver,
				 {ParametersValue, InitialConditions, IntegrationPoints,
					Output, DataMatrix,Factor1,Factor2,Factor3,
						MaxStepRatio, MinStepRatio, MaxIterationsNumber, 
							OrderIncrement, MinOrder, MaxOrder,	
								RelativeTolerance,AbsoluteTolerance,DefectErrorControl},
									AddFunctions, AddPartials, Optimization
				}/.{opt}/.Options[TSMCodeFiles];
			If[pd < 16, pd = 16];
			If[pd === 16 && pre === Multiple, pd = 3999999999];
			If[pd === 16 && Head[pre] === Multiple, pd = pre[[1]]];
			If[vers==False, nvers = 0];
			If[vers=="C", nvers= 1];
			If[vers=="FORTRAN" || vers=="Fortran", nvers= 2];
			If[nvers == 100, Message[TSMCodeFiles::"badopt"]; sigo = False];
			nx = ModifyFirstOrderODE$[x, nfun, dpar, optz];
			If[nvers > 0 && Length[x[[1]]] > Length[nx[[3]]], 
					Message[TSMCodeFiles::"minFun"]; sigo = False ];
			If[sigo, 
				Switch[nvers,
					0, StdFile$[nx, name, pd, dr, odr, params],
					1, CMinFile$[nx, name, dr, odr, params],
					2, FMinFile$[nx, name, dr, odr, params]]]
		]


StdFile$[fun_, name_, n_, dr_, odr_, params_]:=
		Module[{namet, nameh, namec, namedr},
			namet = name <> ".txt";
			namec = name <> ".c";
			nameh = name <> ".h";
			namedr = "dr_" <> name <> ".c";
			Off[DeleteFile::nffil];
			If[dr, 
				Export[namet, DriverStdC[name,n,fun, params]];
				DeleteFile[namedr];
				RenameFile[namet,namedr];];
			If[!odr,
				Export[namet, StandardCText[name,fun,n]];
				DeleteFile[namec];
				RenameFile[namet,namec];
				Export[namet, StandardHText[name,n]];
				DeleteFile[nameh];
				RenameFile[namet,nameh];];
			On[DeleteFile::nffil];
			EndMessageC1[name, dr,odr]
		]


CMinFile$[fun_, name_,dr_, odr_, params_]:=
		Module[{namet, namec, namedr, adf},
			namet = name <> ".txt";
			namec = name <> ".c";
			namedr = "dr_" <> name <> ".c";
			adf = ToTaylorLKFPar[fun];
			Off[DeleteFile::nffil];
			If[dr, 
				Export[namet, DriverMinC[name,adf, params]];
				DeleteFile[namedr];
				RenameFile[namet,namedr];];
			If[!odr,
				Export[namet, MinCCText[name,adf]];
				DeleteFile[namec];
				RenameFile[namet,namec];
				Readmcc["minc_tides.c"];
				Readmch["minc_tides.h"];];
			On[DeleteFile::nffil];
			EndMessageC2[name,dr,odr]
		]	


FMinFile$[fun_, name_,dr_, odr_, params_]:=
		Module[{namet, namec, namedr, adf},
			namet = name <> ".txt";
			namec = name <> ".f";
			namedr = "dr_" <> name <> ".f";
			adf = ToTaylorLKFPar[fun];
			Off[DeleteFile::nffil];
			If[dr, 
				Export[namet, DriverMinFortran[name,adf, params]];
				DeleteFile[namedr];
				RenameFile[namet,namedr];];
			If[!odr,
				Export[namet, MinFFText[name,adf]];
				DeleteFile[namec];
				RenameFile[namet,namec];
				Readmff["minf_tides.f"];];
			On[DeleteFile::nffil];
			EndMessageF[name, adf[[1]]-1,adf[[2]],dr, odr]
		]	




Readmcc[file_String]:=
	Module[{},
		DeleteFile[file];
		Export[file, mincgen, "Text"];]

Readmch[file_String]:=
	Module[{},
		DeleteFile[file];
		Export[file, minhgen, "Text"];]

Readmff[file_String]:=
	Module[{},
		DeleteFile[file];
		Export[file, minfgen, "Text"];]



(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`Codes`"]

EndPackage[]

Null
