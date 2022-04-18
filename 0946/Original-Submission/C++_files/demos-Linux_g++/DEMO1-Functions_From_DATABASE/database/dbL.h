/*	==================================================================
	RELIADIFF Laplace Transforms/Inverses database
	==================================================================
 
	89 Laplace Transforms are defined in file dbLaplace.c, to test RELIADIFF software 
	and compare its results with the 89 Laplace Inverse Transforms 
	contained in the file dbInvLaplace.c
	
	Each Transform function is of the kind
	T<double> fzXX(T<double> z)
	where XX is the number of the function in the database
	
	Each Inverse function is of the kind
	double gzXX(double)
	where XX is the number of the function in the database.
	
	==================================================================
	AUTHORS
	==================================================================

		Luisa D'Amore - University of Naples, Federico II
		Rosanna Campagna - University of Naples, Federico II
		Valeria Mele - University of Naples, Federico II
		Almerico Murli - CMCC and SPACI

	==================================================================	
 
*/
#ifndef DBL_H
	#define DBL_H	
	#include "../../utility/Util.h"
		
	/*DB LAPLACE Transforms*/
	/*
	===========================================================================
	ARGUMENTS
	===========================================================================
	
	each function requires in input: 
	z: 		(TADIFF) double precision: the evaluation point of the Transform
	
	===========================================================================
	RETURN VALUE
	===========================================================================
	
	each function returns:
	(TADIFF) double precision: the evaluation of the Transform in z
	
	===========================================================================	
	*/
	T<double> fz1( T<double> );
	T<double> fz2( T<double> );
	T<double> fz3( T<double> );
	T<double> fz4( T<double> );
	T<double> fz5( T<double> );
	T<double> fz6( T<double> );
	T<double> fz7( T<double> );
	T<double> fz8( T<double> );
	T<double> fz9( T<double> );
	T<double> fz10(T<double> );
	T<double> fz11(T<double> );
	T<double> fz12(T<double> );
	T<double> fz13(T<double> );
	T<double> fz14(T<double> );
	T<double> fz15(T<double> );
	T<double> fz16(T<double> );
	T<double> fz17(T<double> );
	T<double> fz18(T<double> );
	T<double> fz19(T<double> );
	T<double> fz20(T<double> );
	T<double> fz21(T<double> );
	T<double> fz22(T<double> );
	T<double> fz23(T<double> );
	T<double> fz24(T<double> );
	T<double> fz25(T<double> );
	T<double> fz26(T<double> );
	T<double> fz27(T<double> );
	T<double> fz28(T<double> );
	T<double> fz29(T<double> );
	T<double> fz30(T<double> );
	T<double> fz31(T<double> );
	T<double> fz32(T<double> );
	T<double> fz33(T<double> );
	
	T<double> fz34(T<double> );
	T<double> fz35( T<double> );
	T<double> fz36( T<double> );
	T<double> fz37( T<double> );
	T<double> fz38( T<double> );
	T<double> fz39( T<double> );
	T<double> fz40( T<double> );
	T<double> fz41( T<double> );
	T<double> fz42( T<double> );
	T<double> fz43( T<double> );
	T<double> fz44(T<double> );
	T<double> fz45(T<double> );
	T<double> fz46(T<double> );
	T<double> fz47(T<double> );
	T<double> fz48(T<double> );
	T<double> fz49(T<double> );
	T<double> fz50(T<double> );
	T<double> fz51(T<double> );
	T<double> fz52(T<double> );
	T<double> fz53(T<double> );
	T<double> fz54(T<double> );
	T<double> fz55(T<double> );
	T<double> fz56(T<double> );
	T<double> fz57(T<double> );
	T<double> fz58(T<double> );
	T<double> fz59(T<double> );
	T<double> fz60(T<double> );
	T<double> fz61(T<double> );
	T<double> fz62(T<double> );
	T<double> fz63(T<double> );
	T<double> fz64(T<double> );
	T<double> fz65(T<double> );
	T<double> fz66(T<double> );
	T<double> fz67(T<double> );
	T<double> fz68(T<double> );
	T<double> fz69( T<double> );
	T<double> fz70( T<double> );
	T<double> fz71( T<double> );
	T<double> fz72( T<double> );
	T<double> fz73( T<double> );
	T<double> fz74( T<double> );
	T<double> fz75( T<double> );
	T<double> fz76( T<double> );
	T<double> fz77( T<double> );
	T<double> fz78(T<double> );
	T<double> fz79(T<double> );
	T<double> fz80(T<double> );
	T<double> fz81(T<double> );
	T<double> fz82(T<double> );
	T<double> fz83(T<double> );
	T<double> fz84(T<double> );
	
	T<double> fz85(T<double> );
	T<double> fz86(T<double> );
	T<double> fz87(T<double> );
	T<double> fz88(T<double> );
	T<double> fz89(T<double> );
	T<double> fz90(T<double> );
	T<double> fz91(T<double> );
	T<double> fz92(T<double> );


	

	/*DB LAPLACE Inverses*/
	/*
	===========================================================================
	ARGUMENTS
	===========================================================================
	
	each function requires in input: 
	t: 		double precision: the evaluation point of the Inverse function
	
	===========================================================================
	RETURN VALUE
	===========================================================================
	
	each function returns:
	double precision: the evaluation of the Inverse function in t
	
	===========================================================================	
	*/
	
	double gz1( double );
	double gz2( double );
	double gz3( double );
	double gz4( double );
	double gz5( double );
	double gz6( double );
	double gz7( double );
	double gz8( double );
	double gz9( double );
	double gz10(double );
	double gz11(double );
	double gz12(double );
	double gz13(double );
	double gz14(double );
	double gz15(double );
	double gz16(double );
	double gz17(double );
	double gz18(double );
	double gz19(double );
	double gz20(double );
	double gz21(double );
	double gz22(double );
	double gz23(double );
	double gz24(double );
	double gz25(double );
	double gz26(double );
	double gz27(double );
	double gz28(double );
	double gz29(double );
	double gz30(double );
	double gz31(double );
	double gz32(double );
	double gz33(double );
	double gz34(double );
	
	double gz35( double );
	double gz36( double );
	double gz37( double );
	double gz38( double );
	double gz39( double );
	double gz40( double );
	double gz41( double );
	double gz42( double );
	double gz43( double );
	double gz44(double );
	double gz45(double );
	double gz46(double );
	double gz47(double );
	double gz48(double );
	double gz49(double );
	double gz50(double );
	double gz51(double );
	double gz52(double );
	double gz53(double );
	double gz54(double );
	double gz55(double );
	double gz56(double );
	double gz57(double );
	double gz58(double );
	double gz59(double );
	double gz60(double );
	double gz61(double );
	double gz62(double );
	double gz63(double );
	double gz64(double );
	double gz65(double );
	double gz66(double );
	double gz67(double );
	double gz68(double );
	
	double gz69( double );
	double gz70( double );
	double gz71( double );
	double gz72( double );
	double gz73( double );
	double gz74( double );
	double gz75( double );
	double gz76( double );
	double gz77( double );
	double gz78(double );
	double gz79(double );
	double gz80(double );
	double gz81(double );
	double gz82(double );
	double gz83(double );
	double gz84(double );
	
	double gz85(double );
	double gz86(double );
	double gz87(double );
	double gz88(double );
	double gz89(double );
	double gz90(double );
	double gz91(double );
	double gz92(double );

#endif
