(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`LKFunctions: Linked Functions*)


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
(*LKFunctions*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{ 
LKF, LKFC, ToLKF, ToTaylorLKF, ExtractConstants,
ToLKFC, ToTaylorLKFC, LKFPar, LKFCPar,
ToLKFPar, ToTaylorLKFPar, ToLKFCPar,
ToTaylorLKFCPar, IterationLKF, RightIterationLKF,
LeftIterationLKF, NumberOfVariables, 
NumberOfParameters, NumberOfFunctions,
NumberOfLinks, LinksVariables, LinksFunctions,
LieDer, TMinus, strFunC, strFunF, strFunCmp,
StringCNumber, StringFNumber, ListDoubleConstants,
ListCTextConstants, ListFTextConstants,
TextDPConstants, TextMPConstants, 
Double, Multiple
}
	
{
LKF$Link, LKF$Var,LKF$Par,LKF$Constant,LKF$Const,
LKF$Plus, LKF$Minus,LKF$Times, LKF$Power, LKF$Divide,LKF$Sin,
LKF$Cos,LKF$Tan,LKF$Sinh,LKF$Cosh,LKF$Tanh,LKF$ArcSin,
LKF$ArcCos,LKF$ArcTan,LKF$ArcSinh,LKF$ArcCosh,LKF$ArcTanh,
LKF$Csc,LKF$Sec,LKF$Cot,LKF$ArcCsc,LKF$ArcSec,LKF$ArcCot,
LKF$Csch,LKF$Sech,LKF$Coth,LKF$ArcCsch,LKF$ArcSech,
LKF$ArcCoth,LKF$Log,LKF$Exp,LKF$Constant,LKF$Der,
LKF$List, LKF$Names, 
fromLKF$Link,depthDer
}

{
varListQ, ChangeIndexDer, depthDer, numberADQ, numberADPQ,
SymbolsList,$digitsNUMTIDES,numberSQ
}



(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`LKFunctions`*"]
Clear @@ Names["MathTIDES`LKFunctions`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


ToLKF ::usage = "..."
ToTaylorLKF ::usage = "..."
LKF ::usage = "..."

IterationLKF ::usage = "..."
RightIterationLKF ::usage = "..."
LeftIterationLKF ::usage = "..."
LieDer ::usage = "..."

ToLKF::"badpar"="Bad declaration of arguments. The output is wrong"
ToLKF::"badvarg"="Bad declaration of variables"
ToLKF::"badparg"="Bad declaration of parameters"
ToLKF::"badopts"="Incompatible options"



End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*Cabeceras de funciones*)


LKF$Rules = {
	Plus == LKF$Plus, 
	TMinus == LKF$Minus, 
	Times == LKF$Times, 
	Power == LKF$Power, 
	Divide == LKF$Divide,
	Sin == LKF$Sin,
	Cos == LKF$Cos,
	Tan == LKF$Tan,
	Sinh == LKF$Sinh,
	Cosh == LKF$Cosh,
	Tanh == LKF$Tanh,
	Csc == LKF$Csc,
	Sec == LKF$Sec,
	Cot == LKF$Cot,
	Csch == LKF$Csch,
	Sech == LKF$Sech,
	Coth == LKF$Coth,
	ArcSin == LKF$ArcSin,
	ArcCos == LKF$ArcCos,
	ArcTan == LKF$ArcTan,
	ArcSinh == LKF$ArcSinh,
	ArcCosh == LKF$ArcCosh,
	ArcTanh == LKF$ArcTanh,
	ArcCsc == LKF$ArcCsc,
	ArcSec == LKF$ArcSec,
	ArcCot == LKF$ArcCot,
	ArcCsch == LKF$ArcCsch,
	ArcSech == LKF$ArcSech,
	ArcCoth == LKF$ArcCoth,
	Log == LKF$Log,
	Exp == LKF$Exp,
	LieDer == LKF$Der
}


ToLKF$Link = Map[(#[[1]]->#[[2]])&, LKF$Rules]

fromLKF$Link = Map[(#[[2]]->#[[1]])&, LKF$Rules]

twoToOneArg$ = 
	Map[(#[x_,y_]->#[x])&,
		{LKF$Sin, LKF$Cos, 
		LKF$Sinh, LKF$Cosh, 
		LKF$ArcSin, LKF$ArcCos, LKF$ArcTan, 
		LKF$ArcSinh, LKF$ArcCosh, LKF$ArcTanh}] 

SetAttributes[LKF$Times, Orderless];
SetAttributes[LKF$Plus, Orderless];

TMinus[x_,y_]:= x-y;

LKF$Plus /: LKF$Plus[x_]:= x;
LKF$Plus /: LKF$Plus[x_?numberSQ, y_?numberSQ]:= x + y;
LKF$Plus /: LKF$Plus[x_, y_, z__]:= Fold[LKF$Plus[#2,#1]&, x,{y,z}];
LKF$Plus /: LKF$Plus[0, y_]:= y;
LKF$Plus /: LKF$Plus[x_?numberSQ, y_]:= LKF$Plus[LKF$Constant[x] , y];
LKF$Plus /: LKF$Plus[x_, y_?numberSQ]:= LKF$Plus[LKF$Constant[y] , x];
LKF$Plus /: LKF$Plus[x_, LKF$Times[-1,y_]]:= LKF$Minus[x,y];
LKF$Plus /: LKF$Plus[LKF$Times[-1,y_],x_]:= LKF$Minus[x,y];
LKF$Plus /: LKF$Plus[LKF$Times[-1,x_], LKF$Times[-1,y_]]:= 
	LKF$Times[-1,LKF$Plus[x,y]];

LKF$Times /: LKF$Times[x_]:= x;
LKF$Times /: LKF$Times[x_?numberSQ, y_?numberSQ]:= x  y;
LKF$Times /: LKF$Times[x_, y_, z__]:= Fold[LKF$Times[#2,#1]&, x,{y,z}];
LKF$Times /: LKF$Times[1, y_]:= y;
LKF$Times /: LKF$Times[0, y_]:= 0;
LKF$Times /: LKF$Times[x_?numberSQ, y_]:= LKF$Times[LKF$Constant[x] , y];
LKF$Times /: LKF$Times[x_, y_?numberSQ]:= LKF$Times[LKF$Constant[y] , x];

LKF$Divide /: LKF$Divide[1, LKF$Divide[1,x_]]:= x;
LKF$Divide /: LKF$Divide[1, x_?numberSQ]:= 1/x;
LKF$Divide /: LKF$Times[x_, LKF$Divide[1,y_]]:= LKF$Divide[x, y]/; !numberSQ[x]
LKF$Divide /: LKF$Divide[x_?numberSQ, y_?numberSQ]:= x/y;
LKF$Divide /: LKF$Divide[x_?numberSQ/;x =!=1, y_]:= LKF$Divide[LKF$Constant[x] , y];
LKF$Divide /: LKF$Divide[x_, y_?numberSQ]:= LKF$Divide[x,LKF$Constant[y]];

LKF$Power /: LKF$Power[x_,  0]:= 1;
LKF$Power /: LKF$Power[x_,  1]:= x;
LKF$Power /: LKF$Power[x_, -1]:= LKF$Divide[1,x];
LKF$Power /: LKF$Power[x_,  2]:= LKF$Times[x, x];
LKF$Power /: LKF$Power[x_, -2]:= LKF$Divide[1,LKF$Times[x, x]];
LKF$Power /: LKF$Power[x_,  3]:= LKF$Times[x, LKF$Times[x, x]];

LKF$Power /: LKF$Power[x_, - y_]:= 
	LKF$Divide[1,LKF$Power[x,y]]/;!numberSQ[y]
LKF$Power /: LKF$Power[E, x_]:= LKF$Exp[x];
LKF$Power /: LKF$Power[a_?numberADPQ, x_]:= 
	LKF$Exp[LKF$Times[LKF$Log[a], x]];

LKF$Power /: LKF$Power[x_, a_?numberSQ]:= 
	LKF$Power[x , LKF$Constant[a]]


LKF$Log   /: LKF$Log[a_?numberSQ]:= LKF$Constant[Log,a]; 
LKF$Exp   /: LKF$Exp[a_?numberSQ]:= LKF$Constant[Exp,a]; 

LKF$Sin   /: LKF$Sin[a_?numberSQ]:= LKF$Constant[Sin,a]; 
LKF$Cos   /: LKF$Cos[a_?numberSQ]:= LKF$Constant[Cos,a]; 
LKF$Tan   /: LKF$Tan[a_?numberSQ]:= LKF$Constant[Tan,a]; 

LKF$Sinh   /: LKF$Sinh[a_?numberSQ]:= LKF$Constant[Sinh,a]; 
LKF$Cosh   /: LKF$Cosh[a_?numberSQ]:= LKF$Constant[Cosh,a]; 
LKF$Tanh  /: LKF$Tanh[a_?numberSQ]:= LKF$Constant[Tanh,a]; 

LKF$Csc  /: LKF$Csc[a_?numberSQ]:= LKF$Constant[Csc,a]; 
LKF$Sec  /: LKF$Sec[a_?numberSQ]:= LKF$Constant[Sec,a]; 
LKF$Cot  /: LKF$Cot[a_?numberSQ]:= LKF$Constant[Cot,a]; 

LKF$Csch  /: LKF$Csch[a_?numberSQ]:= LKF$Constant[Csch,a]; 
LKF$Sech  /: LKF$Sech[a_?numberSQ]:= LKF$Constant[Sech,a]; 
LKF$Coth  /: LKF$Coth[a_?numberSQ]:= LKF$Constant[Coth,a]; 

LKF$ArcSin  /: LKF$ArcSin[a_?numberSQ]:= LKF$Constant[ArcSin,a]; 
LKF$ArcCos  /: LKF$ArcCos[a_?numberSQ]:= LKF$Constant[ArcCos,a]; 
LKF$ArcTan  /: LKF$ArcTan[a_?numberSQ]:= LKF$Constant[ArcTan,a]; 

LKF$ArcCsc  /: LKF$ArcCsc[a_?numberSQ]:= LKF$Constant[ArcCsc,a]; 
LKF$ArcSec  /: LKF$ArcSec[a_?numberSQ]:= LKF$Constant[ArcSec,a]; 
LKF$ArcCot  /: LKF$ArcCot[a_?numberSQ]:= LKF$Constant[ArcCot,a]; 

LKF$ArcSinh  /: LKF$ArcSinh[a_?numberSQ]:= LKF$Constant[ArcSinh,a]; 
LKF$ArcCosh  /: LKF$ArcCosh[a_?numberSQ]:= LKF$Constant[ArcCosh,a]; 
LKF$ArcTanh  /: LKF$ArcTanh[a_?numberSQ]:= LKF$Constant[ArcTanh,a]; 

LKF$ArcCsch  /: LKF$ArcCsch[a_?numberSQ]:= LKF$Constant[ArcCsch,a]; 
LKF$ArcSech  /: LKF$ArcSech[a_?numberSQ]:= LKF$Constant[ArcSech,a]; 
LKF$ArcCoth  /: LKF$ArcCoth[a_?numberSQ]:= LKF$Constant[ArcCoth,a]; 


LKF$Tan  /: LKF$Tan[x_]:= LKF$Divide[LKF$Sin[x],LKF$Cos[x]];
LKF$Tanh /: LKF$Tanh[x_]:= LKF$Divide[LKF$Sinh[x],LKF$Cosh[x]];

LKF$Csc /: LKF$Csc[x_]:= LKF$Divide[1,LKF$Sin[x]];
LKF$Sec /: LKF$Sec[x_]:= LKF$Divide[1,LKF$Cos[x]];
LKF$Cot /: LKF$Cot[x_]:= LKF$Divide[LKF$Cos[x],LKF$Sin[x]];
LKF$Csch /: LKF$Csch[x_]:= LKF$Divide[1,LKF$Sinh[x]];
LKF$Sech /: LKF$Sech[x_]:= LKF$Divide[1,LKF$Cosh[x]];
LKF$Coth /: LKF$Coth[x_]:= LKF$Divide[LKF$Cosh[x],LKF$Sinh[x]];

LKF$ArcCsc /: LKF$ArcCsc[x_]:= LKF$ArcSin[LKF$Divide[1,x]];
LKF$ArcSec /: LKF$ArcSec[x_]:= LKF$ArcCos[LKF$Divide[1,x]];
LKF$ArcCot /: LKF$ArcCot[x_]:= LKF$ArcTan[LKF$Divide[1,x]];
LKF$ArcCsch /: LKF$ArcCsch[x_]:= LKF$ArcSinh[LKF$Divide[1,x]];
LKF$ArcSech /: LKF$ArcSech[x_]:= LKF$ArcCosh[LKF$Divide[1,x]];
LKF$ArcCoth /: LKF$ArcCoth[x_]:= LKF$ArcTanh[LKF$Divide[1,x]];


LKF$Der	/: LKF$Der[n_Integer][x_]:= LKF$Der[LKF$Der[n-1][x]];
LKF$Der	/: LKF$Der[1][x_]:= LKF$Der[x];


LKF$Names = {
	LKF$Plus,LKF$Minus,LKF$Times,LKF$Power,
	LKF$Divide,LKF$Sin,LKF$Cos,
	LKF$Tan,LKF$Sinh,LKF$Cosh,LKF$Tanh,
	LKF$Csc,LKF$Sec,LKF$Cot,LKF$Csch,
	LKF$Sech,LKF$Coth,LKF$ArcSin,LKF$ArcCos,
	LKF$ArcTan,LKF$ArcSinh,LKF$ArcCosh,
	LKF$ArcTanh,LKF$ArcCsc,LKF$ArcSec,
	LKF$ArcCot,LKF$ArcCsch,LKF$ArcSech,
	LKF$ArcCoth,LKF$Log,LKF$Exp,LKF$Der,
	LKF$Constant, LKF$List, E}


(* ::Subsection::Closed:: *)
(*Test para comprobaciones y otras funciones auxiliares*)


$digitsNUMTIDES = 16 


varADQ[x_LKF$Var]:= True
varADQ[x_LKF$Link]:= True
varADQ[x_LKF$Par]:= True
varADQ[x_LKF$Constant]:= False
varADQ[x_?NumberQ]:= False
varADQ[x_]:= False


varListQ[{}]:= True
varListQ[x_List]:= And @@ Map[(Head[x[[#]]] === Symbol)&, Range[Length[x]]]
varListQ[x_]:= False /; Head[x] != List
varListQ[x__]:= False


numberSQ[x_?NumberQ]:= True
numberSQ[E]:= True
numberSQ[Pi]:= True
numberSQ[_]:= False


(*Para trabajar con parametros como series*)
numberADQ[x_?numberSQ]:= True
numberADQ[x_LKF$Constant]:= True
numberADQ[x_LKF$Const]:= True
numberADQ[x_LKF$Par]:= False
numberADQ[x_LKF$Var]:= False
numberADQ[x_LKF$Link]:= False
numberADQ[x_]:= False


(*Para trabajar con parametros sin series*)
numberADPQ[x_?numberSQ]:= True
numberADPQ[x_LKF$Constant]:= True
numberADPQ[x_LKF$Const]:= True
numberADPQ[x_LKF$Par]:= True
numberADPQ[x_LKF$Var]:= False
numberADPQ[x_LKF$Link]:= False
numberADPQ[x_]:= False


(* ::Subsection::Closed:: *)
(*Busqueda y sustituci\[OAcute]n de Parametros y constantes*)


findPar[LKF$Plus[x_LKF$Par, y_LKF$Par]] := True
findPar[LKF$Plus[x_LKF$Constant, y_LKF$Par]] := True
findPar[LKF$Plus[x_LKF$Par, y_LKF$Constant]] := True
findPar[LKF$Minus[x_LKF$Par, y_LKF$Par]] := True
findPar[LKF$Minus[x_LKF$Constant, y_LKF$Par]] := True
findPar[LKF$Minus[x_LKF$Par, y_LKF$Constant]] := True
findPar[LKF$Times[x_LKF$Par, y_LKF$Par]] := True
findPar[LKF$Times[x_LKF$Constant, y_LKF$Par]] := True
findPar[LKF$Times[x_LKF$Par, y_LKF$Constant]] := True
findPar[LKF$Power[x_LKF$Par, y_LKF$Par]] := True
findPar[LKF$Power[x_LKF$Constant, y_LKF$Par]] := True
findPar[LKF$Power[x_LKF$Par, y_LKF$Constant]] := True
findPar[LKF$Divide[x_LKF$Par, y_LKF$Par]] := True
findPar[LKF$Divide[x_LKF$Constant, y_LKF$Par]] := True
findPar[LKF$Divide[x_LKF$Par, y_LKF$Constant]] := True
findPar[LKF$Sin[x_LKF$Par]] := True
findPar[LKF$Cos[x_LKF$Par]] := True
findPar[LKF$Tan[x_LKF$Par]] := True
findPar[LKF$Sinh[x_LKF$Par]] := True
findPar[LKF$Cosh[x_LKF$Par]] := True
findPar[LKF$Tanh[x_LKF$Par]] := True
findPar[LKF$ArcSin[x_LKF$Par]] := True
findPar[LKF$ArcCos[x_LKF$Par]] := True
findPar[LKF$ArcTan[x_LKF$Par]] := True
findPar[LKF$ArcSinh[x_LKF$Par]] := True
findPar[LKF$ArcCosh[x_LKF$Par]] := True
findPar[LKF$ArcTanh[x_LKF$Par]] := True
findPar[LKF$Log[x_LKF$Par]] := True
findPar[LKF$Exp[x_LKF$Par]] := True
findPar[_] := False


findConst[LKF$Plus[x_LKF$Const, y_LKF$Const]] := True
findConst[LKF$Minus[x_LKF$Const, y_LKF$Const]] := True
findConst[LKF$Times[x_LKF$Const, y_LKF$Const]] := True
findConst[LKF$Power[x_LKF$Const, y_LKF$Const]] := True
findConst[LKF$Divide[x_LKF$Const,y_LKF$Const]] := True
findConst[LKF$Sin[x_LKF$Const]] := True
findConst[LKF$Cos[x_LKF$Const]] := True
findConst[LKF$Tan[x_LKF$Const]] := True
findConst[LKF$Sinh[x_LKF$Const]] := True
findConst[LKF$Cosh[x_LKF$Const]] := True
findConst[LKF$Tanh[x_LKF$Const]] := True
findConst[LKF$ArcSin[x_LKF$Const]] := True
findConst[LKF$ArcCos[x_LKF$Const]] := True
findConst[LKF$ArcTan[x_LKF$Const]] := True
findConst[LKF$ArcSinh[x_LKF$Const]] := True
findConst[LKF$ArcCosh[x_LKF$Const]] := True
findConst[LKF$ArcTanh[x_LKF$Const]] := True
findConst[LKF$Log[x_LKF$Const]] := True
findConst[LKF$Exp[x_LKF$Const]] := True
findConst[_] := False


(* ::Subsection::Closed:: *)
(*Creaci\[OAcute]n funci\[OAcute]n de iteraci\[OAcute]n*)


SymbolsList[nexp_]:=
	Union[Cases[nexp, _Symbol , 
		Infinity, Heads->True]]


ToLKF$[fun_, var_, par_]:=
	Module[{vpl, nexp,slist, vn,pn,rulv,rulp,indice , 
			links = {}, binfun, nfun,ii, npar,norder,
			listiter ={}, listpar = {}, iteritem, nivel},
		nfun = If[Head[fun]===List,LKF$List@@fun, LKF$List[fun]];
		vpl = Union[var,par]; 
		nexp = nfun/.ToLKF$Link;
		slist = SymbolsList[nexp/.LKF$Constant[f_,a__]->LKF$Constant[1]];
		slist = Select[slist,!MemberQ[LKF$Names,#]&];
		slist = Select[slist,!MemberQ[vpl,#]&];
		If[Length[slist] =!= 0,  Message[ToLKF::"badpar"]; Print[slist]; Abort[]];
		vn = Length[var];
		pn = Length[par];
		rulv = Map[(var[[#]]->LKF$Var[#])&, Range[vn]];
		rulp = Map[(par[[#]]->LKF$Par[#])&, Range[pn]];
		nexp = nexp/.rulv/.rulp;
		
		indice = 0;
		binfun = Union[Level[nexp,{-3}]];
		While[Length[binfun] > 0 && Head[binfun[[1]]] =!= LKF$List,
			listiter = Flatten[Append[listiter,
					Map[(binfun[[#]] -> LKF$Link[#+indice])&,
						 Range[Length[binfun]]]]];
			links = Flatten[Append[links,binfun]];
			indice = Length[listiter];
			nexp = nexp //.listiter;
			binfun = Union[Level[nexp,{-3}]];
			]; 
		nexp = nexp/.LKF$List->List;
		nexp = Map[If[numberSQ[#], LKF$Constant[#],#]&, nexp];

		ii = 1;
		npar = 0;
		listpar = Map[#[[1]]&,listpar];
		listiter = Map[#[[1]]&,listiter];
		
		LKF[vn,pn,listiter,nexp]
	]


ToLKFPar$[fun_, var_, par_]:=
	Module[{vpl, nexp,slist, vn,pn,rulv,rulp,indice , 
			links = {}, binfun, nfun,ii, npar,norder,
			listiter ={}, listpar = {}, iteritem, nivel},
		nfun = If[Head[fun]===List,LKF$List@@fun, LKF$List[fun]];
		vpl = Union[var,par]; 
		nexp = nfun/.ToLKF$Link;
		slist = SymbolsList[nexp/.LKF$Constant[f_,a__]->LKF$Constant[1]];
		slist = Select[slist,!MemberQ[LKF$Names,#]&];
		slist = Select[slist,!MemberQ[vpl,#]&];
		If[Length[slist] =!= 0,  Message[ToLKF::"badpar"]; Print[slist]; Abort[]];
		vn = Length[var];
		pn = Length[par];
		rulv = Map[(var[[#]]->LKF$Var[#])&, Range[vn]];
		rulp = Map[(par[[#]]->LKF$Par[#])&, Range[pn]];
		nexp = nexp/.rulv/.rulp;
		
		indice = 0;
		binfun = Union[Level[nexp,{-3}]];
		While[Length[binfun] > 0 && Head[binfun[[1]]] =!= LKF$List,
			listiter = Flatten[Append[listiter,
					Map[(binfun[[#]] -> LKF$Link[#+indice])&,
						 Range[Length[binfun]]]]];
			links = Flatten[Append[links,binfun]];
			indice = Length[listiter];
			nexp = nexp //.listiter;
			binfun = Union[Level[nexp,{-3}]];
			]; 
		nexp = nexp/.LKF$List->List;
		nexp = Map[If[numberSQ[#], LKF$Constant[#],#]&, nexp];

		ii = 1;
		listpar = {};
		npar = 0;
		While[ii <= Length[listiter],
  			iteritem = listiter[[ii]];
  			If[findPar[iteritem[[1]]],
    				norder = iteritem[[2,1]];
    				listiter = listiter/.LKF$Link[norder]->LKF$Par[pn+npar+1];
    				AppendTo[listpar,listiter[[ii]][[1]]];
    				listiter = Drop[listiter,{ii}]/.changeIndex[norder,-1];
    				nexp = nexp/.LKF$Link[norder]->LKF$Par[pn+npar+1];
    				nexp = nexp/.changeIndex[norder,-1];
    				ii--; npar++;];
  			ii++];
		listiter = Map[#[[1]]&,listiter];

		LKFPar[vn,pn, listpar,listiter,nexp]
	]


(* ::Subsubsection:: *)
(*Funci\[OAcute]n del usuario*)


ToLKF[fun_,var_,par_]:=
  Module[{sal},
	If[!varListQ[var],Message[ToLKF::"badvarg"]; Print[var]; Abort[]];
	If[!varListQ[par],Message[ToLKF::"badparg"]; Print[par]; Abort[]];
	If[Head[fun] === List, 
		sal =  ToLKF$[fun,var,par],
		sal =  ToLKF$[{fun},var,par]];
	sal
     ]

ToLKF[fun_,var_]:=
  ToLKF[fun, var, {}]


ToLKFPar[fun_,var_,par_]:=
  Module[{sal},
	If[!varListQ[var],Message[ToLKF::"badvarg"]; Print[var]; Abort[]];
	If[!varListQ[par],Message[ToLKF::"badparg"]; Print[par]; Abort[]];
	If[Head[fun] === List, 
		sal =  ToLKFPar$[fun,var,par],
		sal =  ToLKFPar$[{fun},var,par]];
	sal
     ]

ToLKFPar[fun_,var_]:=
  ToLKFPar[fun, var, {}]


(* ::Subsection::Closed:: *)
(*Extraccion de constantes*)


ExtractConstants[lkf_LKF]:=
	Module[{consts, nconsts,newlkf,iter,itef,item},
		consts = Union[Cases[lkf, LKF$Constant[__], Infinity]];
		nconsts = Length[consts];
		newlkf = lkf/.Map[(consts[[#]]-> LKF$Const[#])&, Range[nconsts]];
		iter = newlkf[[3]];
		itef = newlkf[[4]];
		item = FindFirstConstant[iter];
		While[item =!= False,
			consts = Append[consts,item];
			nconsts = Length[consts];
			iter = iter/.{item ->LKF$Const[Position[consts,item][[1,1]]]};
			{iter,itef} = ExtractIsolateConstants[iter,itef];
			item = FindFirstConstant[iter];
		];
		LKFC[newlkf[[1]],newlkf[[2]],consts,iter,itef]
	]


ExtractConstants[lkf_LKFPar]:=
	Module[{consts, nconsts,newlkf,iter,iterp,itef,item},
		consts = Union[Cases[lkf, LKF$Constant[__], Infinity]];
		nconsts = Length[consts];
		newlkf = lkf/.Map[(consts[[#]]-> LKF$Const[#])&, Range[nconsts]];
		iterp = newlkf[[3]];
		iter = newlkf[[4]];
		itef = newlkf[[5]];
		item = FindFirstConstant[iter];
		While[item=!= False,
			consts = Append[consts,item];
			nconsts = Length[consts];
			iter = iter/.{item ->LKF$Const[Position[consts,item][[1,1]]]};
			iterp = iterp/.{item ->LKF$Const[Position[consts,item][[1,1]]]};
			{iter,itef} = ExtractIsolateConstants[iter,itef];
			item =  FindFirstConstant[iter];
		];
		LKFCPar[newlkf[[1]],newlkf[[2]], consts, iterp, iter,itef]
	]


ExtractIsolateConstants[x_,xf_]:=
	Module[{pos, nx, nxf },
		nx = x;
		nxf= xf;
		pos =Flatten[Position[x, LKF$Const[_],1]] ;
		If[pos =!= {}, 
			nx = nx/.Map[(LKF$Link[#]->x[[#]])&, pos];
			Map[(nx = nx/.LKF$Link[n_/;n>#]->LKF$Link[n-1];)&, Reverse[pos]];
			Map[(nxf = nxf/.LKF$Link[n_/;n>#]->LKF$Link[n-1];)&, Reverse[pos]];
			nx = DeleteCases[nx,LKF$Const[_]]];
		{nx, nxf}
	]


FindFirstConstant[iter_]:=
	Module[{cases},
		cases = Select[iter, findConst];
		If[cases === {}, False, cases[[1]]]
	]


ToLKFC[x__]:=
	ExtractConstants[ToLKF[x]]

ToLKFCPar[x__]:=
	ExtractConstants[ToLKFPar[x]]

ToTaylorLKFC[x__]:=
	ExtractConstants[ToTaylorLKF[x]]

ToTaylorLKFCPar[x__]:=
	ExtractConstants[ToTaylorLKFPar[x]]


(* ::Subsection::Closed:: *)
(*Informaci\[OAcute]n funci\[OAcute]n de iteraci\[OAcute]n *)


IterationLKF[iter_List]:=
	Module[{listder, listizq, np},
		np = Length[iter];
		listizq = Map[LKF$Link, Range[Length[iter]]];
		listder = iter;
		Map[(listizq[[#]] == listder[[#]])&, Range[np]]
	]


IterationLKF[LKF[v_Integer,c_Integer,iter_List,fun_]]:=
	Module[{listder, listizq, np},
		np = Length[iter];
		listizq = Map[LKF$Link, Range[Length[iter]]];
		listder = iter;
		Map[(listizq[[#]] == listder[[#]])&, Range[np]]
	]

IterationLKF[LKFPar[v_Integer,c_Integer,itpar_List, iter_List,fun_]]:=
	Module[{listder, listizq, np},
		np = Length[itpar]+Length[iter];
		listizq = Join[Map[LKF$Par, Range[c+1,Length[itpar]+c]], 
			Map[LKF$Link, Range[Length[iter]]]];
		listder = Join[itpar, iter];
		Map[(listizq[[#]] == listder[[#]])&, Range[np]]
	]

RightIterationLKF[adf_LKF]:=
	Map[(#[[1]]->#[[2]])&, IterationLKF[adf]]
	
LeftIterationLKF[adf_LKF]:=
	Map[(#[[2]]->#[[1]])&, IterationLKF[adf]]
	
RightIterationLKF[adf_LKFPar]:=
	Map[(#[[1]]->#[[2]])&, IterationLKF[adf]]
	
LeftIterationLKF[adf_LKFPar]:=
	Map[(#[[2]]->#[[1]])&, IterationLKF[adf]]



(* ::Subsection::Closed:: *)
(*Modificaci\[OAcute]n necesaria para aumentar los elementos de iteraci\[OAcute]n para series de Taylor*)


changeIndex[k_,i_]  := {LKF$Link[(n_)?(#1 > k & )]  :> LKF$Link[n +i]}
changeIndexE[k_,i_] := {LKF$Link[(n_)?(#1 >= k & )] :> LKF$Link[n +i]}


addIterT[{itg_,fl_}, LKF$Link[i_]==LKF$Sin[x_]]:=
  	Module[{sel,ix,nitg,nfl},
		sel = Select[itg,(#[[2]]==LKF$Cos[x])&];
		nitg = itg;
		nfl = fl;
		If[sel == {},
			nitg = nitg/.changeIndex[i,1];
			nfl = nfl/.changeIndex[i,1];
			nitg[[i]] = (LKF$Link[i]== LKF$Sin[x,LKF$Link[i+1]]);
			nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Cos[x,LKF$Link[i]]),i+1],
			
			ix =sel[[1,1,1]];
			nitg[[i]] = (LKF$Link[i]== LKF$Sin[x,LKF$Link[ix]]);
			nitg[[ix]] = (LKF$Link[ix]== LKF$Cos[x,LKF$Link[i]]);
		];
		{nitg,nfl}
	]
	
addIterT[{itg_,fl_}, LKF$Link[i_]==LKF$Cos[x_]]:=
  	Module[{sel,ix, nitg,nfl},
		sel = Select[itg,(#[[2]]==LKF$Sin[x])&];
		nitg = itg;
		nfl = fl;
		If[sel == {},
			nitg = nitg/.changeIndex[i,1];
			nfl = nfl/.changeIndex[i,1];
			nitg[[i]] = (LKF$Link[i]== LKF$Cos[x,LKF$Link[i+1]]);
			nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Sin[x,LKF$Link[i]]),i+1],
			
			ix =sel[[1,1,1]];
			nitg[[i]] = (LKF$Link[i]== LKF$Cos[x,LKF$Link[ix]]);
			nitg[[ix]] = (LKF$Link[ix]== LKF$Sin[x,LKF$Link[i]]);
		];
		{nitg,nfl}
	]
	
addIterT[{itg_,fl_}, LKF$Link[i_]==LKF$Sinh[x_]]:=
  	Module[{sel,ix,nitg,nfl},
		sel = Select[itg,(#[[2]]==LKF$Cosh[x])&];
		nitg = itg;
		nfl = fl;
		If[sel == {},
			nitg = nitg/.changeIndex[i,1];
			nfl = nfl/.changeIndex[i,1];
			nitg[[i]] = (LKF$Link[i]== LKF$Sinh[x,LKF$Link[i+1]]);
			nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Cosh[x,LKF$Link[i]]),i+1],
			
			ix =sel[[1,1,1]];
			nitg[[i]] = (LKF$Link[i]== LKF$Sinh[x,LKF$Link[ix]]);
			nitg[[ix]] = (LKF$Link[ix]== LKF$Cosh[x,LKF$Link[i]]);
		];
		{nitg,nfl}
	]
	
addIterT[{itg_,fl_}, LKF$Link[i_]==LKF$Cosh[x_]]:=
  	Module[{sel,ix, nitg,nfl},
		sel = Select[itg,(#[[2]]==LKF$Sinh[x])&];
		nitg = itg;
		nfl = fl;
		If[sel == {},
			nitg = nitg/.changeIndex[i,1];
			nfl = nfl/.changeIndex[i,1];
			nitg[[i]] = (LKF$Link[i]== LKF$Cosh[x,LKF$Link[i+1]]);
			nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Sinh[x,LKF$Link[i]]),i+1],
			
			ix =sel[[1,1,1]];
			nitg[[i]] = (LKF$Link[i]== LKF$Cosh[x,LKF$Link[ix]]);
			nitg[[ix]] = (LKF$Link[ix]== LKF$Sinh[x,LKF$Link[i]]);
		];
		{nitg,nfl}
	]
	

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$ArcSin[x_]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,4];
		nfl = fl/.changeIndexE[i,4]; 
		nitg[[i]] = (LKF$Link[i]== LKF$Power[x,2]);
		nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Times[-1,LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]== LKF$Plus[1,LKF$Link[i+1]]),i+2];
		nitg = Insert[nitg, (LKF$Link[i+3]== LKF$Power[LKF$Link[i+2],Rational[1,2]]),i+3];
		nitg = Insert[nitg, (LKF$Link[i+4]== LKF$ArcSin[x,LKF$Link[i+3]]),i+4];
		{nitg,nfl}
	]

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$ArcCos[x_]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,5];
		nfl = fl/.changeIndexE[i,5];
		nitg[[i]] = (LKF$Link[i]== LKF$Power[x,2]);
		nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Times[-1,LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]== LKF$Plus[1,LKF$Link[i+1]]),i+2];
		nitg = Insert[nitg, (LKF$Link[i+3]== LKF$Power[LKF$Link[i+2],Rational[1,2]]),i+3];
		nitg = Insert[nitg, (LKF$Link[i+4]== LKF$Times[-1,LKF$Link[i+3]]),i+4];
		nitg = Insert[nitg, (LKF$Link[i+5]== LKF$ArcCos[x,LKF$Link[i+4]]),i+5];
		{nitg,nfl}
	]

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$ArcTan[x_]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,2];
		nfl = fl/.changeIndexE[i,2];
		nitg[[i]] = (LKF$Link[i]== LKF$Power[x,2]);
		nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Plus[1,LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]== LKF$ArcTan[x,LKF$Link[i+1]]),i+2];
		{nitg,nfl}
	]

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$ArcSinh[x_]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,3];
		nfl = fl/.changeIndexE[i,3];
		nitg[[i]] = (LKF$Link[i]== LKF$Power[x,2]);
		nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Plus[1,LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]== LKF$Power[LKF$Link[i+1], Rational[1,2]]),i+2];
		nitg = Insert[nitg, (LKF$Link[i+3]== LKF$ArcSinh[x,LKF$Link[i+2]]),i+3];
		{nitg,nfl}
	]

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$ArcCosh[x_]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,3];
		nfl = fl/.changeIndexE[i,3];
		nitg[[i]] = (LKF$Link[i]== LKF$Power[x,2]);
		nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Plus[-1,LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]== LKF$Power[LKF$Link[i+1], Rational[1,2]]),i+2];
		nitg = Insert[nitg, (LKF$Link[i+3]== LKF$ArcCosh[x,LKF$Link[i+2]]),i+3];
		{nitg,nfl}
	]

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$ArcTanh[x_]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,3];
		nfl = fl/.changeIndexE[i,3];
		nitg[[i]] = (LKF$Link[i]== LKF$Power[x,2]);
		nitg = Insert[nitg, (LKF$Link[i+1]== LKF$Times[-1,LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]== LKF$Plus[1,LKF$Link[i+1]]),i+2];
		nitg = Insert[nitg, (LKF$Link[i+3]== LKF$ArcTanh[x,LKF$Link[i+2]]),i+3];
		{nitg,nfl}
	]

addIterT[{itg_,fl_}, LKF$Link[i_]== LKF$Power[x_?varADQ,y_?varADQ]]:=
  	Module[{nitg,nfl},
		nitg = itg/.changeIndexE[i,2];
		nfl = fl/.changeIndexE[i,2];
		nitg[[i]] = (LKF$Link[i] ==  LKF$Log[x]);
		nitg = Insert[nitg, (LKF$Link[i+1]==  LKF$Times[y, LKF$Link[i]]),i+1];
		nitg = Insert[nitg, (LKF$Link[i+2]==  LKF$Exp[LKF$Link[i+1]]),i+2];
		{nitg,nfl}
	]


addIterT[{itg_,fl_},_]:= {itg,fl}



addTaylorIter[itg_List, fl_List]:=
	Module[{nitg,nfl,ii},
		ii = 1;
		nitg = itg;
		nfl = fl;
		While[ii <= Length[nitg],
			{nitg,nfl} = addIterT[{nitg,nfl},nitg[[ii]]];
    		ii++];	
		{nitg,nfl}
	]



ToTaylorLKF[adf_LKF]:=
	Module[{auxg, auxf, it, itg, v, c,  nit},
		v = adf[[1]];
		c = adf[[2]];
		it = IterationLKF[adf];
		itg = Select[it,(Head[#[[1]]]=== LKF$Link)&];
		{auxg, auxf} = addTaylorIter[itg, adf[[4]]];
		auxg = Map[#[[2]]&,auxg];
		LKF[v,c,auxg, auxf]
	]

ToTaylorLKF[args__]:=
	ToTaylorLKF[ToLKF[args]]



ToTaylorLKFPar[adf_LKFPar]:=
	Module[{auxg, auxf,it,itp,itg, v, c,  nit},
		v = adf[[1]];
		c = adf[[2]];
		it = IterationLKF[adf];
		itp = adf[[3]];
		itg = Select[it,(Head[#[[1]]]=== LKF$Link)&];
		{auxg, auxf} = addTaylorIter[itg, adf[[5]]];
		auxg = Map[#[[2]]&,auxg];
		LKFPar[v,c,itp, auxg, auxf]
	]

ToTaylorLKFPar[args__]:=
	ToTaylorLKFPar[ToLKFPar[args]]




(* ::Subsection::Closed:: *)
(*Avance de \[IAcute]ndices para evaluar derivadas*)


ChangeIndexDer[]:= Module[{},
	Clear[depthDer];
	
	depthDer[E]:= 0;
	depthDer[LKF$Par[_]]:= 0;
	depthDer[LKF$Var[_]]:= 0;
	depthDer[LKF$Constant[__]]:= 0;
	depthDer[LKF$Const[__]]:= 0;
	depthDer[x_?numberSQ]:= 0;
	depthDer[LKF$Link[_]]:= 0;

	depthDer[LKF$Link[n_] == LKF$Der[x_]]:= 
		depthDer[LKF$Link[n]] = depthDer[x]+1;
	
	depthDer[LKF$Link[n_] == LKF$h_[x__]]:=
		depthDer[LKF$Link[n]] = Max @@ Map[depthDer, {x}];

	depthDer[LKF$Par[n_] == _]:= Null;
]




(* ::Subsection::Closed:: *)
(*Evaluaci\[OAcute]n de funciones de iteraci\[OAcute]n*)


LKF[v_Integer,c_Integer,iter_List,fun_][x___]:=
	Module[{totalvc, totaliter,  xval, niter,  vc, 
				rulvc,rulpar,rulgen, lmax,final},
		xval = {x};
		totalvc   = v + c;
		vc = Join[Map[LKF$Var, Range[v]], Map[LKF$Par, Range[c]]];
		rulvc = Inner[Rule[#1,#2]&, vc,xval,List];
		rulpar = {LKF$Constant[ct_]->ct, LKF$Constant[f_,ct__]->f[ct]};		
		niter = iter/.twoToOneArg$;
		rulgen = Inner[Rule[#1,#2]&, Map[LKF$Link, Range[Length[niter]]],niter,List];
		final = fun//.rulgen//.rulpar/.rulvc/.fromLKF$Link;
		If[Length[final] === 1, final[[1]], final]
	]/;  Length[{x}] == v + c
	



LKFPar[v_Integer,c_Integer,itpar_List, iter_List,fun_][x___]:=
	Module[{totalvc, totaliter,  xval, niter,  vc, 
				rulvc,rulpar,rulgen, lmax,final},
		xval = {x};
		totalvc   = v + c;
		vc = Join[Map[LKF$Var, Range[v]], Map[LKF$Par, Range[c]]];
		rulvc = Inner[Rule[#1,#2]&, vc,xval,List];
		rulpar = Inner[Rule[#1,#2]&, 
			Map[LKF$Par, Range[c+1, c+Length[itpar]]],itpar,List]/.
					rulvc/.{LKF$Constant[ct_]->ct, LKF$Constant[f_,ct__]->f[ct]};
		niter = iter/.twoToOneArg$;
		rulgen = Inner[Rule[#1,#2]&, Map[LKF$Link, Range[Length[niter]]],niter,List];
		final = fun//.rulgen//.rulpar/.rulvc/.fromLKF$Link;
		If[Length[final] === 1, final[[1]], final]
	]/;  Length[{x}] == v + c



(* ::Subsection::Closed:: *)
(*Constantes e Iteraciones*)


NumberOfVariables[fun_LKF]:= fun[[1]]-1
NumberOfParameters[fun_LKF]:= fun[[2]]
NumberOfFunctions[fun_LKF]:=  Length[fun[[4]]]-fun[[1]]+1
NumberOfLinks[fun_LKF]:=  Length[IterationLKF[fun]]

NumberOfVariables[fun_LKFC]:= fun[[1]]-1
NumberOfParameters[fun_LKFC]:= fun[[2]]
NumberOfFunctions[fun_LKFC]:=  Length[fun[[5]]]-fun[[1]]+1
NumberOfLinks[fun_LKFC]:=  Length[IterationLKF[fun[[4]]]]


NumberOfVariables[fun_LKFPar]:= fun[[1]]-1
NumberOfParameters[fun_LKFPar]:= fun[[2]]
NumberOfFunctions[fun_LKFPar]:=  Length[fun[[5]]]-fun[[1]]+1
NumberOfLinks[fun_LKFPar]:=  Length[IterationLKF[fun]]

NumberOfVariables[fun_LKFCPar]:= fun[[1]]-1
NumberOfParameters[fun_LKFCPar]:= fun[[2]]
NumberOfFunctions[fun_LKFCPar]:=  Length[fun[[6]]]-fun[[1]]+1
NumberOfLinks[fun_LKFCPar]:=  Length[IterationLKF[fun[[5]]]]


LinksVariables[fun_LKF]:=
	Drop[fun[[4]],-NumberOfFunctions[fun]]

LinksVariables[fun_LKFC]:=
	Drop[fun[[5]],-NumberOfFunctions[fun]]

LinksVariables[fun_LKFPar]:=
	Drop[fun[[5]],-NumberOfFunctions[fun]]

LinksVariables[fun_LKFCPar]:=
	Drop[fun[[6]],-NumberOfFunctions[fun]]


PosFunction[LKF$Link[n_]]:=  n-1;
PosFunction[LKF$Var[n_]]:= -(n-1);

LinksFunctions[fun_LKF]:=
	Map[PosFunction,Take[fun[[4]],-NumberOfFunctions[fun]]]

LinksFunctions[fun_LKFC]:=
	Map[PosFunction,Take[fun[[5]],-NumberOfFunctions[fun]]]

LinksFunctions[fun_LKFPar]:=
	Map[PosFunction,Take[fun[[5]],-NumberOfFunctions[fun]]]

LinksFunctions[fun_LKFCPar]:=
	Map[PosFunction,Take[fun[[6]],-NumberOfFunctions[fun]]]


(* ::Subsection::Closed:: *)
(*Textos comunes*)


strFunC[LKF$Power]= "pow";
strFunC[LKF$Exp]= "exp";
strFunC[LKF$Sin]= "sin";
strFunC[LKF$Cos]= "cos";
strFunC[LKF$Tan]= "tan";
strFunC[LKF$Sinh]= "sinh";
strFunC[LKF$Cosh]= "cosh";
strFunC[LKF$Tanh]= "tanh";
strFunC[LKF$ArcSin]= "asin";
strFunC[LKF$ArcCos]= "acos";
strFunC[LKF$ArcTan]= "atan";
strFunC[LKF$ArcSinh]= "asinh";
strFunC[LKF$ArcCosh]= "acosh";
strFunC[LKF$ArcTanh]= "atanh";
strFunC[LKF$Log]= "log";


strFunF[LKF$Sin]= "DSIN";
strFunF[LKF$Cos]= "DCOS";
strFunF[LKF$Tan]= "DTAN";
strFunF[LKF$Sinh]= "DSINH";
strFunF[LKF$Cosh]= "DCOSH";
strFunF[LKF$Tanh]= "DTANH";
strFunF[LKF$ArcSin]= "DASIN";
strFunF[LKF$ArcCos]= "DACOS";
strFunF[LKF$ArcTan]= "DATAN";
strFunF[LKF$ArcSinh]= "DASINH";
strFunF[LKF$ArcCosh]= "DACOSH";
strFunF[LKF$ArcTanh]= "DATANH";
strFunF[LKF$Log]= "DLOG";
strFunF[LKF$Exp]= "DEXP";


strFunCmp[LKF$Power]= "mpfrts_pow";
strFunCmp[LKF$Sin]= "mpfrts_sin";
strFunCmp[LKF$Cos]= "mpfrts_cos";
strFunCmp[LKF$Tan]= "mpfrts_tan";
strFunCmp[LKF$Sinh]= "mpfrts_sinh";
strFunCmp[LKF$Cosh]= "mpfrts_cosh";
strFunCmp[LKF$Tanh]= "mpfrts_tanh";
strFunCmp[LKF$ArcSin]= "mpfrts_asin";
strFunCmp[LKF$ArcCos]= "mpfrts_acos";
strFunCmp[LKF$ArcTan]= "mpfrts_atan";
strFunCmp[LKF$ArcSinh]= "mpfrts_asinh";
strFunCmp[LKF$ArcCosh]= "mpfrts_acosh";
strFunCmp[LKF$ArcTanh]= "mpfrts_atanh";
strFunCmp[LKF$Log]= "mpfrts_log";
strFunCmp[LKF$Exp]= "mpfrts_exp";

strFunCmp[Power]= "mpfrts_pow";
strFunCmp[Sin]= "mpfrts_sin";
strFunCmp[Cos]= "mpfrts_cos";
strFunCmp[Tan]= "mpfrts_tan";
strFunCmp[Sinh]= "mpfrts_sinh";
strFunCmp[Cosh]= "mpfrts_cosh";
strFunCmp[Tanh]= "mpfrts_tanh";
strFunCmp[ArcSin]= "mpfrts_asin";
strFunCmp[ArcCos]= "mpfrts_acos";
strFunCmp[ArcTan]= "mpfrts_atan";
strFunCmp[ArcSinh]= "mpfrts_asinh";
strFunCmp[ArcCosh]= "mpfrts_acosh";
strFunCmp[ArcTanh]= "mpfrts_atanh";
strFunCmp[Log]= "mpfrts_log";
strFunCmp[Exp]= "mpfrts_exp";


ConstantToReal[LKF$Constant[x_]]:= N[x,16]
ConstantToReal[LKF$Constant[f_,x__]]:=N[f[x],16]
ConstantToReal[x_Real]:=N[x,16]


StringCNumber[0]:= "0.0" 
StringCNumber[x_Real]:=
	ToString[CForm[N[x,16]]];
StringCNumber[x_]:=
	StringCNumber[N[x,16]];

StringFNumber[0]:= "0.d0" 
StringFNumber[x_Real]:=  
	Module[{name,es},
		name = ToString[FortranForm[x]];
		es = StringCases[name, "e"];
		name = If[es === {}, name <> "d0", StringReplace[name, "e"->"d"]];
		name] 
StringFNumber[x_]:=  
	StringFNumber[N[x,16]]


ListDoubleConstants[x_]:=
	Module[{ct,nct,cnum,mod, first},
		ct =Cases[x,LKF$Constant[__]];
		mod =Cases[x,Except[LKF$Constant[__]]];
		nct = Length[ct];
		cnum   =Map[ConstantToReal,ct];
		While[Length[mod] > 0,
			first = mod[[1]];
			first = first/. Map[(LKF$Const[#]->cnum[[#]])&, Range[nct]];
			AppendTo[cnum,ConstantToReal[first/.fromLKF$Link]];
		nct = Length[cnum];
		mod = Drop[mod,1]];
		cnum
	]


ListCTextConstants[x_]:=
		Map[StringCNumber,ListDoubleConstants[x]]

ListFTextConstants[x_]:=
		Map[StringFNumber,ListDoubleConstants[x]]


TextDPConstants[{}]:= ""

TextDPConstants[x_]:= 
	Module[{len,texto="", const},
		len = Length[x];
		const = ListCTextConstants[x];
		texto = texto <>"\tdouble ct[] = {"<> const[[1]];
		Map[(texto = texto <> ", " <> const[[#]])&, Range[2,len]];
		texto = texto <>"};\n";
		texto] 


CstrConstmpAux[LKF$Constant[E], var_String]:= 
	"\tmpfrts_set_i(&"<> var <> ", 1);\n"<>
	"\tmpfrts_exp(&"<> var <> ","<> var <>");\n";
CstrConstmpAux[LKF$Constant[Pi], var_String]:= 
	"\tmpfr_const_pi("<> var <> ", GMP_RND);\n";
CstrConstmpAux[LKF$Constant[x_Integer], var_String]:= 
	"\tmpfrts_set_i(&"<> var <> ","<> ToString[x] <>");\n";
CstrConstmpAux[LKF$Constant[x_Rational], var_String]:= 
	"\tmpfrts_set_i(&"<> var <> ","<> ToString[Numerator[x]] <>");\n" <>
	"\tmpfrts_div_i(&"<> var <> ","<> var <> ","<> ToString[Denominator[x]] <>");\n" ;
CstrConstmpAux[LKF$Constant[x_Real], var_String]:= 
	"\tmpfrts_set_str(&"<> var <> ","<> ToString[N[x, Precision[x]]] <>");\n";



CstrConstmp[LKF$Const[n_]]:= "ct[" <> ToString[n-1] <> "]";

CstrConstmp[LKF$Constant[f_,x_], i_]:= 
	Module[{par,texto},
		par = "ct["<>ToString[i-1] <> "]";
		texto = CstrConstmpAux[LKF$Constant[x], par ];
		texto = texto <> "\t"<>strFunCmp[f]<>"(&" <> par <> ", " <> par <> ");\n";
		texto]
CstrConstmp[x_LKF$Constant, i_]:= 
	CstrConstmpAux[x, "ct["<> ToString[i-1] <> "]" ]

CstrConstmp[LKF$Plus[x_,y_],i_]:= 
	"\tmpfrts_add(&ct["<> ToString[i-1] <> "], "<> CstrConstmp[x] <> "," <> CstrConstmp[y] <>");\n";
CstrConstmp[LKF$Times[x_,y_], i_]:= 	
	"\tmpfrts_mul(&ct["<> ToString[i-1] <> "], "<> CstrConstmp[x] <> "," <> CstrConstmp[y] <>");\n";
CstrConstmp[LKF$Power[x_,y_], i_]:= 
	"\tmpfrts_pow(&ct["<> ToString[i-1] <> "], "<> CstrConstmp[x] <> "," <> CstrConstmp[y] <>");\n";
CstrConstmp[LKF$Divide[x_,y_], i_]:= 
	"\tmpfrts_div(&ct["<> ToString[i-1] <> "], "<> CstrConstmp[x] <> "," <> CstrConstmp[y] <>");\n";
CstrConstmp[h_[x_], i_]:= 
	"\t"<>strFunCmp[h]<>"(&ct["<>ToString[i-1] <> "],"<> CstrConstmp[x] <> ");\n";



TextMPConstants[x_]:= 
	Module[{len,texto=""},
		len = Length[x];
		texto = texto <>"\tint NCONST = "<>ToString[len]<>";\n";
		texto = texto <>"\tmpfr_t ct[NCONST];\n";
		texto = texto <> "\tfor(i = 0; i < NCONST ; i++ ) { ";
		texto = texto<>"\n\t\tmpfrts_init(&ct[i]);\n\t}\n";
		texto = texto<>(StringJoin @@ Map[CstrConstmp[x[[#]],#]&, Range[len]]);
		texto] 

TextMPConstants[]:= 
	Module[{len,texto=""},
		texto = texto <> "\tfor(i = 0; i < NCONST ; i++ ) { ";
		texto = texto<>"\n\t\tmpfr_clear(ct[i]);\n\t}\n";
		texto] 


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`LKFunctions`"]

EndPackage[]

Null
