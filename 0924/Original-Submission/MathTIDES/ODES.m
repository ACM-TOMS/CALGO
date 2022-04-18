(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`ODES: Ordinary Differential Equations*)


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
(*ODES*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`ODES`",
				"MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
	FirstOrderODE$,
	FirstOrderODE,
	HamiltonianToODE,
	PotentialToODE,
	NthOrderODE,
	Screen,
	Delta,
	Points,
	Only,
	Until
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`ODES`*"]
Clear @@ Names["MathTIDES`ODES`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


FirstOrderODE$::"badvarg"="Bad declaration of variables"
FirstOrderODE$::"badparg"="Bad declaration of parameters"
FirstOrderODE$::"badFunLength"="Number of equations different than number of variables"
FirstOrderODE$::"badDVarLength"="Bad Length in List of AddVariables"
FirstOrderODE$::"badCM"="Number of coordinates different than number of momenta"
FirstOrderODE$::"badHam"="Hamiltonian must be scalar"
FirstOrderODE$::"badPot"="Potential must be scalar"
FirstOrderODE$::"badnthODE"="Bad nth ODE definition"
FirstOrderODE$::"badopt"="Option AddVariables not implemented in NthOrderODE"
FirstOrderODE$::"badt"="Bad independent variable"
FirstOrderODE$::"badder"="Bad list of variables in partial derivatives"
FirstOrderODE$::"badderord"="Incompatible orders of AddPartials"


FirstOrderODE$::usage = "..."


FirstOrderODE::usage =ToString[Style["FirstOrderODE[] ",Bold],StandardForm]<>
"declares a first order ODE to be used with the expression CodeFiles[] "<>"
to create the Taylor Integration Code. "<>"The arguments of FirstOrderODE[] are:"<>
ToString[Style["\n\nFirst argument: ",Bold,Italic],StandardForm]<>
"The list with the expressions of the derivatives of the variables. "<>
"The number (n) of elements of the list must be equal to the number of variables. "<>
"If n = 1 the argument is not a list."<>
ToString[Style["\nSecond argument: ",Bold,Italic],StandardForm]<>
"The symbol that represents the independent variable. "<>
"This symbol may appear explicitly or not in the first argument."<>
ToString[Style["\nThird argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the variables. "<>
"It has the same number of elements than the first argument. "<>
"If there is only one variable the argument is not a list."<>
ToString[Style["\nFourth argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the parameters. "<>
"If the number of parameters is equal to 1 the argument is not a list. "<>
"If there is no parameter this argument may be avoided."<>
ToString[Style["\nOptions: ",Bold,Italic],StandardForm]<>
"To modify the declaration of the ODE. The possible options are: "<>
ToString[Style["ReplaceSymbol, Optimization,  AddFunctions & AddPartials. ",Bold],StandardForm]<>
ToString[Style["\n\nExample: ",Bold,Italic],StandardForm]<>
"Let's suppose the first order differential system "<>
ToString[{Dt[x,t]==y,Dt[y,t]==-x},TraditionalForm]<>
". To represent the previous system in "<>
ToString[Style["MathTIDES",Bold],StandardForm]<>
" we will write "<>ToString[Style["FirstOrderODE[{y, -x}, t, {x,y}]",Bold],StandardForm]


NthOrderODE::usage = 
ToString[Style["NthOrderODE[] ",Bold],StandardForm]<>"declares a n-th order ODE to be used with the expression CodeFiles[] "<>
"to create the Taylor Integration Code. "<>
"The system is transformed into a first order systems by extending the variables including the derivatives until order (n-1). "<>
"The arguments of NthOrderODE[] are:"<>
ToString[Style["\n\nFirst argument: ",Bold,Italic],StandardForm]<>
"The list with the system of equations with the format defined by the following rules: "<>
"\na) The derivatives of a variable of symbol x must be represented by quotes: x, x', x'', x''', ..."<>
"\nb) The equations are represented by means of the symbol =="<>
"\nc) The number of equations is equal to the number of variables"<>
"\nd) If the number of variables is equal to one, the first and the third arguments are not lists"<>
"\ne) In the system must appear the derivatives of grater order of all the variables."<>
ToString[Style["\nSecond argument: ",Bold,Italic],StandardForm]<>
"The symbol that represents the independent variable. "<>
"This symbol may appear explicitly or not in the first argument."<>
ToString[Style["\nThird argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the variables. "<>
"It has the same number of elements than the first argument. "<>
"If there is only one variable the argument is not a list."<>
ToString[Style["\nFourth argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the parameters. "<>
"If the number of parameters is equal to 1 the argument is not a list. "<>
"If there is no parameter this argument may be avoided."<>
ToString[Style["\nOptions: ",Bold,Italic],StandardForm]<>
"To modify the declaration of the ODE. The possible options are: "<>
ToString[Style["ReplaceSymbol, Optimization,  AddFunctions & AddPartials. ",Bold],StandardForm]<>
ToString[Style["\n\n Example: ",Bold,Italic],StandardForm]<>
"Let's suppose the third order differential system "<>
ToString[{








\!\(\*SuperscriptBox["x", 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)-2 








\!\(\*SuperscriptBox["y", "\[Prime]\[Prime]",
MultilineFunction->None]\)+ 








\!\(\*SuperscriptBox["x", "\[Prime]",
MultilineFunction->None]\)==2 x^2-y, 4 








\!\(\*SuperscriptBox["y", 
TagBox[
RowBox[{"(", "3", ")"}],
Derivative],
MultilineFunction->None]\)-2 








\!\(\*SuperscriptBox["x", "\[Prime]\[Prime]",
MultilineFunction->None]\) 








\!\(\*SuperscriptBox["y", "\[Prime]",
MultilineFunction->None]\)==2 x+y^2},TraditionalForm]<>
". To represent the previous system in "<>ToString[Style["MathTIDES",Bold],StandardForm]<>" we will write "<>ToString[Style["NthOrderODE[{x''' + x '- 2 y'' \[Equal] 2 x^2 - y, 4 y''' - 2 x '' y' \[Equal] 2 x + y^2},t, {x,y}]",Bold],StandardForm]<>
". This system is transformed into a first order system with the variables "<>" {x, y, x', y', x'', y''}"


PotentialToODE::usage =
ToString[Style["PotentialToODE ",Bold],StandardForm]<>"declares a first order ODE, equivalent to the Newton's equation "<>
ToString[








\!\(\*SuperscriptBox[
OverscriptBox["x", "_"], "\[Prime]\[Prime]",
MultilineFunction->None]\)==-\[Del]V,TraditionalForm]<>
"t, where V is the potential function, and  creates the Taylor Integration Code. "<>
"The arguments of"<> 
ToString[Style["PotentialToODE[] ",Bold],StandardForm] <>"are:"<>
ToString[Style["\n\nFirst argument: ",Bold,Italic],StandardForm]<>
"The expression of the potential V. This expression is never a list."<>
ToString[Style["\nSecond argument: ",Bold,Italic],StandardForm]<>
"The symbol that represents the independent variable. "<>
ToString[Style["Third argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the variables. "<>
"If there is only one variable the argument is not a list."<>
ToString[Style["\nFourth argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the parameters. "<>
"If the number of parameters is equal to 1 the argument is not a list. "<>
"If there is no parameter this argument may be avoided."<>
ToString[Style["\nOptions: ",Bold,Italic],StandardForm]<>
"To modify the declaration of the ODE. The possible options are: "<>
ToString[Style["ReplaceSymbol, Optimization,  AddFunctions & AddPartials. ",Bold],StandardForm]<>
ToString[Style["\n\nExample: ",Bold,Italic],StandardForm]<>
"Let's suppose the potential "<>ToString[V==\[Mu]/Sqrt[x^2+y^2],TraditionalForm]<>
", with the parameter \[Mu]. Then the Newton' s equations will be represented in "<>
ToString[Style["MathTIDES",Bold],StandardForm]<>" by "<>
ToString[Style[" PotentialToODE[\[Mu]/Sqrt[x^2+y^2] , t,{x,y},\[Mu]]",Bold],StandardForm]<>
". This system is transformed into a first order system with the variables "<>" {x, y, x', y'}"


HamiltonianToODE::usage = ToString[Style["HamiltonianToODE[] ",Bold],StandardForm]<>
"declares a first order ODE, equivalent to the Hamilton equations "<>
"of a dynamical system with Hamiltonian H, and  creates the Taylor Integration Code. "<>
"The arguments of HamiltonianToODE[] are:"<>
ToString[Style["\n\nFirst argument: ",Bold,Italic],StandardForm]<>
"The expression of the Hamiltonian H. This expression is never a list."<>
ToString[Style["\nSecond argument: ",Bold,Italic],StandardForm]<>
"The symbol that represents the independent variable. "<>
ToString[Style["\nThird argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the variables and momenta {x,X}. "<>
"The length of this list is always an even number. "<>
"The order of the momenta corresponds with the order of the associated variables."<>
ToString[Style["\nFourth argument: ",Bold,Italic],StandardForm]<>
"The list of symbols of the parameters. "<>
"If the number of parameters is equal to 1 the argument is not a list. "<>
"If there is no parameter this argument may be avoided."<>
ToString[Style["\nOptions: ",Bold,Italic],StandardForm]<>
"To modify the declaration of the ODE. The possible options are: "<>
ToString[Style["ReplaceSymbol, Optimization,  AddFunctions & AddPartials. ",Bold],StandardForm]<>
ToString[Style["\n\nExample: ",Bold,Italic],StandardForm]<>
"Let's suppose the Hamiltonian "<>
ToString[\[ScriptCapitalH]==1/2 (X^2+Y^2)-1/Sqrt[x^2+y^2],TraditionalForm]<>
" in the variables (x,y) and the momenta (X,Y). Then the Hamilton equations will be represented in "<>
ToString[Style["MathTIDES",Bold],StandardForm]<>
" by "<>ToString[Style["HamiltonianToODE[(X^2 + Y^2) /2 - 1/Sqrt[x^2+y^2], t, {x,y,X,Y}]",Bold],StandardForm]


End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*Ecuaciones de primer orden *)


FirstOrderODE[fun_List,t_, var_List, par_List]:=
	 Module[{ nfun, nfund, nvar, cfun},
		If[Head[t]=!=Symbol,Message[FirstOrderODE$::"badt"];Abort[]];
		If[!varListQ[var],Message[FirstOrderODE$::"badvarg"];Print[var];Abort[]];
		If[!varListQ[par],Message[FirstOrderODE$::"badparg"];Print[par];Abort[]];
		If[(Length[fun] != Length[var]),Message[FirstOrderODE$::"badFunLength"];Abort[]];
		FirstOrderODE$[fun,t, var,par]
     ] 

FirstOrderODE[fun_,t_,var_,par_List]:=
	If[((Head[fun]=!=List) && (Head[var]=!= List)), 
			FirstOrderODE[{fun},t,{var},par],
			Message[FirstOrderODE$::"badFunLength"];Abort[];]

FirstOrderODE[fun_,t_,var_]:=
	FirstOrderODE[fun,t,var,{}]

FirstOrderODE[fun_,t_,var_,par_]:=
	FirstOrderODE[fun,t,var,{par}]



(* ::Subsection::Closed:: *)
(*Hamiltonianos y ecuaciones de Hamilton*)


DefDer = If[#1 =!= #2, #1 /: Dt[#1, #2] = 0] & ;
DefDerD[x_, dvd_, dvv_] := 
  (#1 /: Dt[#1, x[[#2]]] = dvd[[Position[dvv, #1][[1,1]],#2]]) & 


HamEquations[nham_, t_, x_, par_ ]:=
	 Module[{nvar,var,mom,derH},
		nvar = Length[x]/2;
		var = Take[x,nvar];
		mom = Drop[x,nvar];
		Map[(t/: Dt[t,#] = 0)&, x];
		Outer[DefDer,x,x];
		Outer[DefDer,par,x];
		derH = Join[Map[Dt[nham,#]&, mom],
				Map[-Dt[nham,#]&, var]];
		Map[Clear,x];
		Map[Clear,par];
		derH ]


HamiltonianToODE[ham_,t_, var_List,par_List]:=
	  Module[{fun},
		If[Head[ham]===List,Message[FirstOrderODE$::"badHam"]; Abort[]];
		If[Head[t]=!=Symbol,Message[FirstOrderODE$::"badt"];Abort[]];
		If[!varListQ[var],Message[FirstOrderODE$::"badvarg"];Print[var];Abort[]];
		If[!varListQ[par],Message[FirstOrderODE$::"badparg"];Print[par];Abort[]];
		If[OddQ[Length[var]],Message[FirstOrderODE$::"badCM"];Abort[]];
		fun = HamEquations[ham, t, var, par];
		FirstOrderODE$[fun,t, var,par]
]

HamiltonianToODE[ham_,t_, var_List]:=
	HamiltonianToODE[ham,t,var,{}]

HamiltonianToODE[ham_,t_, var_List, par_]:=
	HamiltonianToODE[ham,t,var,{par}]


(* ::Subsection::Closed:: *)
(*Potenciales*)


add$dn[x_Symbol,n_Integer]:=
	ToExpression[StringJoin[ToString[x],"$d",ToString[n]]]


PotToForce[npot_, t_, x_, par_]:=
	 Module[{derH},
		Map[(t/: Dt[t,#] = 0)&, x];
		Outer[DefDer,x,x];
		Outer[DefDer,par,x];
		derH = Map[-Dt[npot,#]&,x];
		Map[Clear,x];
		Map[Clear,par];
		derH ]


PotentialToODE[pot_,t_, var_List, par_List]:=
	  Module[{fun, der,nvar,nfun},
		If[Head[pot]===List,Message[FirstOrderODE$::"badPot"]; Abort[]];
		If[Head[t]=!=Symbol,Message[FirstOrderODE$::"badt"];Abort[]];
		If[!varListQ[var],Message[FirstOrderODE$::"badvarg"];Print[var];Abort[]];
		If[!varListQ[par],Message[FirstOrderODE$::"badparg"];Print[par];Abort[]];
		fun = PotToForce[pot,t,var, par];
		der=Map[add$dn[#,1]&,var]; 
		nvar = Join[var, der];
		nfun = Join[der,fun];
		FirstOrderODE$[nfun,t, nvar, par]
]

PotentialToODE[pot_,t_,var_]:=
	PotentialToODE[pot,t,var,{}]

PotentialToODE[pot_,t_, var_, par_List]:=
	PotentialToODE[pot,t,{var}, par ]

PotentialToODE[pot_,t_, var_, par_]:=
	PotentialToODE[pot,t,var,{par}]


(* ::Subsection::Closed:: *)
(*Ecuaciones de orden n*)


nth$ODE[nfun_,var_,par_]:=
	Module[{ ecu, listader, ordmax, ldomax, sol,funder,der,dern,nvar,funfin},
		ecu=Map[(#[[1]]-#[[2]]==0)&,nfun];
		listader=Cases[ecu,Derivative[_][_],Infinity];
		ordmax=Max[Map[#[[0,1]]&,listader]];
		ldomax=Map[Derivative[ordmax][#]&,var];
		If[!(And@@Map[MemberQ[listader,#]&,ldomax]),Message[FirstOrderODE$::"badnthODE"];];
		sol=Solve[ecu,ldomax];
		If[Length[sol]!=1,Message[FirstOrderODE$::"badode"]];
		funder=ldomax/.sol[[1]]/.{Derivative[n_][s_]:>add$dn[s,n]};
		der[n_]:=Map[add$dn[#,n]&,var];
		dern=Map[der[#]&,Range[ordmax-1]];
		nvar=Join[var,Flatten[dern]];
		funfin=Join[Flatten[dern],funder];
		{funfin,nvar}]


NthOrderODE[fun_List,t_, var_List,par_List]:=
	 Module[{ nfun,   nvar},
		If[Head[t]=!=Symbol,Message[FirstOrderODE$::"badt"];Abort[]];
		If[!varListQ[var],Message[FirstOrderODE$::"badvarg"];Print[var];Abort[]];
		If[!varListQ[par],Message[FirstOrderODE$::"badparg"];Print[par];Abort[]];
		If[(Length[fun] != Length[var]),Message[FirstOrderODE$::"badFunLength"];Abort[]];
		{nfun,nvar} = nth$ODE[fun, var,par];
		FirstOrderODE$[nfun,t, nvar,par]
     ] 

NthOrderODE[fun_,t_,var_,par_List]:=
	If[((Head[fun]=!=List) && (Head[var]=!= List)), 
			NthOrderODE[{fun},t,{var},par],
			Message[FirstOrderODE$::"badFunLength"];Abort[];]

NthOrderODE[fun_,t_,var_]:=
	NthOrderODE[fun,t,var,{}]

NthOrderODE[fun_,t_,var_,par_]:=
	NthOrderODE[fun,t,var,{par}]


(* ::Subsection::Closed:: *)
(*Interacci\[OAcute]n con LKF*)


FirstOrderODE$/: ToTaylorLKF[FirstOrderODE$[fun_, t_,{y__}, par_, opt___]]:=
					ToTaylorLKF[fun, {t,y}, par]

FirstOrderODE$/: ToLKF[FirstOrderODE$[fun_, t_,{y__}, par_, opt___]]:=
			Module[{ly, nfun}, 
				ly = Length[{y}];
				nfun = Take[fun,ly]; 
				ToLKF[nfun, {t,y}, par]
				]


FirstOrderODE$/: ToTaylorLKFPar[FirstOrderODE$[fun_, t_,{y__}, par_, opt___]]:=
					ToTaylorLKFPar[fun, {t,y}, par]

FirstOrderODE$/: ToLKFPar[FirstOrderODE$[fun_, t_,{y__}, par_, opt___]]:=
			Module[{ly, nfun}, 
				ly = Length[{y}];
				nfun = Take[fun,ly]; 
				ToLKFPar[nfun, {t,y}, par]
				]


(* ::Subsection::Closed:: *)
(*Simplificaciones*)


FirstOrderODE$/: Simplify[FirstOrderODE$[fun_, restode___], restsymp___ ]:=
	FirstOrderODE$[Simplify[fun,restsymp], restode]

FirstOrderODE$/: FullSimplify[FirstOrderODE$[fun_, restode___], restsymp___ ]:=
	FirstOrderODE$[FullSimplify[fun,restsymp], restode]


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`ODES`"]

EndPackage[]

Null
