(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`Iterations: Iterations for partial derivatives*)


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
(*Iterations*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`Iterations`",
				"MathTIDES`ODES`",
					"MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
ListIndexFunDer,
CountDerivatives,
PreviousList,
PreviousIndexList,
IteratorsList,
IteratorsListStar,
CompleteIteratorsList,
CompleteIteratorsListStar,
ListToString,
SortListDer,
DerOutputList,
ListIndexFunLastOrder
}


{NTuples, NTuplesLastOrder, OrderedNTuples, Kji, TIndex, VarOrder}


{MinusOnePosition,MinusOnePosition,MinusOneOptimum, ComputeBinomials, PreviousList}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`Iterations`*"]
Clear @@ Names["MathTIDES`Iterations`*"]


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
(*n - tuplas*)


NTuples[nvar_,mord_]:=
	Select[Tuples[Range[0,mord], nvar],
		(Plus@@# <= mord)&]


NTuplesLastOrder[nvar_,mord_]:=
	Select[Tuples[Range[0,mord], nvar],
		(Plus@@# == mord)&]


OrderedNTuples[n_,mord_, Index_]:=
	Sort[ NTuples[n,mord], (Index[#1] <Index[#2])&]


Kji[j_, mu_List]:= 
	Sum[mu[[Length[mu]-i]],{i,0,j-1}]-1

TIndex[mu_List]:= 
	Sum[Binomial[Kji[j,mu]+j,j],{j,1,Length[mu]}]


ListIndexFunDer[nvar_,mord_]:=
	OrderedNTuples[nvar,mord, TIndex]


CountDerivatives[nvar_,mord_]:= 
	Binomial[mord+nvar,nvar]


ListToString[mu_List]:=
	StringJoin@@Map[ToString,mu]


ListIndexFunLastOrder[nvar_,mord_]:=
	Prepend[
		SortListDer[
			NTuplesLastOrder[nvar,mord]],Table[0,{nvar}]]


VarOrder[{mu__List}]:= 
	Module[{nv,  mord},
		nv = Map[Length,{mu}][[1]];
		mord = Max[Map[(Plus @@ #)&, {mu}]];
		{nv,mord}
	]


SortListDer[{mu__List}]:= 
	Module[{lt},
		lt = ListIndexFunDer[Sequence@@VarOrder[{mu}]];
		Sort[{mu}, (Position[lt,#1][[1,1]]<Position[lt,#2][[1,1]])&]
]/; Equal@@Map[Length,{mu}]


(* ::Subsubsection::Closed:: *)
(*Binomials (Duda \[DownQuestion]hacerlo en real con precision?*)


BinomialsInt[ord_]:= 
	Map[(BinTaylor[ord,#] = Binomial[ord,#];)&, 
		Range[0,ord]];


ComputeBinomials[max_]:=  
	( Map[BinomialsInt, Range[0,max]];)


DeleteBinomials[]:= Clear[BinTaylor]


BinomialProduct[li_List, lv_List]:=
	Inner[BinTaylor,li,lv,Times]


(* ::Subsubsection::Closed:: *)
(*Listas de anteriores*)


PreviousList[{mu__List}]:= 
	SortListDer[
		Union[
			Flatten[
				Map[PreviousList,{mu}],1]]]

PreviousList[mu_List]:= 
	Module[{nmu},
		nmu = Union[
			Flatten[
				Outer[List,Sequence @@ Map[Range[0,#]&, mu]],
					Length[mu]-1]];
		SortListDer[nmu]
	]


PreviousIndexList[mu_]:=
	Map[TIndex, PreviousList[mu]]


MinusOnePosition[mu_, pos_]:= 
	If[mu[[pos]] == 0 , {}, 
		{ReplacePart[mu,pos->mu[[pos]]-1]}][[1]]

FirstNonZeroMinimum[mu_List]:=
	Module[{munc, min},
		munc = Select[mu,(#!=0)&];
		min = Min[munc];
		Position[mu, min][[1]]
	]

MinusOneOptimum[mu_]:= 
	MinusOnePosition[mu,FirstNonZeroMinimum[mu][[1]]]


IteratorsList[mu_]:=
	Module[{la , lc, pb},
		la = Union[Flatten[Outer[List,
			Sequence @@ Map[Range[0,#]&, mu]],Length[mu]-1]];
		lc = Map[(mu-#)&, la];
		pb = Map[BinomialProduct [mu, #]&, la];
		{pb,Map[TIndex,la],Map[TIndex, lc]}]


IteratorsListStar[mu_]:=
	{{1},{0},{0}}/; mu === Table[0,{Length[mu]}]

IteratorsListStar[mu_]:=
	Module[{me, la , lc, pb},
		me = MinusOneOptimum[mu];
		la = Union[Flatten[Outer[List,
			Sequence @@ Map[Range[0,#]&, me]],Length[mu]-1]];
		lc = Map[(mu-#)&, la];
		pb = Map[BinomialProduct [me, #]&, la];
		{pb,Map[TIndex,la],Map[TIndex, lc]}]


IteratorsList[nvar_, ordmax_]:=
	Map[IteratorsList, ListIndexFunDer[nvar,ordmax]]

IteratorsList[{mu__List}]:= 
	Map[IteratorsList, {mu}]


IteratorsListStar[nvar_, ordmax_]:=
	Map[IteratorsListStar, ListIndexFunDer[nvar,ordmax]]

IteratorsListStar[{mu__List}]:= 
	Map[IteratorsListStar, {mu}]


(* ::Subsection::Closed:: *)
(*Iteraciones generales*)


DerOutputList[Until[nvar_Integer, mord_Integer]]:=
	DerOutputList[nvar, mord]

DerOutputList[Only[nvar_Integer, mord_Integer]]:=
	DerOutputList[ListIndexFunLastOrder[nvar, mord]]

DerOutputList[Until[nvar_Integer, ord_List]]:=
	DerOutputList[ord, All]

DerOutputList[Only[nvar_Integer, ord_List]]:=
	DerOutputList[ord]


DerOutputList[nvar_, ordmax_]:=
	Module[{ ls,lo, cd},
		cd = CountDerivatives[nvar,ordmax];
		ls = Map[ListToString,ListIndexFunDer[nvar,ordmax]];
		lo = Map[TIndex, ListIndexFunDer[nvar,ordmax]];
		{cd,ls, lo}
	 ]

DerOutputList[{mu__List}, All]:=
	Module[{lmu, nmu, smu, omu},
		nmu = PreviousList[{mu}];
		lmu = Length[nmu];
		smu = Map[ListToString,nmu];
		omu = Range[0,lmu-1];
		{lmu, smu, omu}  ]

DerOutputList[{mu__List}]:=
	Module[{lmu, nmu, xmu, pmu,  smu, omu},
		nmu = PreviousList[{mu}];
		lmu = Length[nmu];
		xmu = {mu};
		PrependTo[xmu,Table[0,{Length[xmu[[1]]]}]];
		xmu = Union[xmu];
		pmu = SortListDer[xmu];
		smu = Map[ListToString,pmu];
		omu = Map[Position[nmu, #][[1,1]]&,pmu];
		{lmu, smu, omu-1}  ]


CompleteIteratorsList[Until[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsList[nvar, mord]

CompleteIteratorsList[Only[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsList[ListIndexFunLastOrder[nvar, mord]]

CompleteIteratorsList[Until[nvar_Integer, ord_List]]:=
	CompleteIteratorsList[ord]

CompleteIteratorsList[Only[nvar_Integer, ord_List]]:=
	CompleteIteratorsList[ord]


CompleteIteratorsListStar[Until[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsListStar[nvar, mord]

CompleteIteratorsListStar[Only[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsListStar[ListIndexFunLastOrder[nvar, mord]]

CompleteIteratorsListStar[Until[nvar_Integer, ord_List]]:=
	CompleteIteratorsListStar[ord]

CompleteIteratorsListStar[Only[nvar_Integer, ord_List]]:=
	CompleteIteratorsListStar[ord]


CompleteIteratorsList[nvar_, ordmax_]:=
	Module[{ todos, acum, coef, lv,lc, cd},
		cd = CountDerivatives[nvar,ordmax];
		ComputeBinomials[ordmax];
		todos = IteratorsList[nvar,ordmax];
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];
		coef = Flatten[Map[#[[1]]&, todos]];
		lv = Flatten[Map[#[[2]]&, todos]];
		lc = Flatten[Map[#[[3]]&, todos]];
		{acum, coef, lv, lc}
	 ]


CompleteIteratorsListStar[nvar_, ordmax_]:=
	Module[{ todos, acum, coef, lv,lc, cd},
		cd = CountDerivatives[nvar,ordmax];
		ComputeBinomials[ordmax];
		todos = IteratorsListStar[nvar,ordmax]; 
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];  
		coef = Flatten[Map[#[[1]]&, todos]]; 
		lv =  Flatten[Map[#[[2]]&, todos]];
		lc =  Flatten[Map[#[[3]]&, todos]]; 
		{acum, coef, lv, lc}
	 ]


CompleteIteratorsList[{mu__List}]:=
	Module[{nvar,ordmax, todos, acum, coef, lv,lc, cd, pl,mpl},
		{nvar,ordmax} = VarOrder[{mu}];
		pl = PreviousList[{mu}];
		cd = Length[pl];
		ComputeBinomials[ordmax];
		todos = IteratorsList[pl];
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];
		coef = Flatten[Map[#[[1]]&, todos]];
		mpl = Map[TIndex, pl];
		rules = Inner[Rule, mpl, Range[0,Length[mpl]-1],List];
		lv = Flatten[Map[#[[2]]&, todos]]/.rules;
		lc = Flatten[Map[#[[3]]&, todos]]/.rules;
		{acum, coef, lv, lc}
	 ]


CompleteIteratorsListStar[{mu__List}]:=
	Module[{nvar,ordmax, todos, acum, coef, lv,lc, cd, pl,mpl},
		{nvar,ordmax} = VarOrder[{mu}];
		pl = PreviousList[{mu}];
		cd = Length[pl];
		ComputeBinomials[ordmax];
		todos = IteratorsListStar[pl];
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];
		coef = Flatten[Map[#[[1]]&, todos]];
		mpl = Map[TIndex, pl];
		rules = Inner[Rule, mpl, Range[0,Length[mpl]-1],List];
		lv = Flatten[Map[#[[2]]&, todos]]/.rules;
		lc = Flatten[Map[#[[3]]&, todos]]/.rules;
		{acum, coef, lv, lc}
	 ]


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`Iterations`"]

EndPackage[]

Null
