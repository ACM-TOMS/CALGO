% Generated by ADiMat 0.6.0-4870
% Copyright 2009-2013 Johannes Willkomm, Fachgebiet Scientific Computing,
% TU Darmstadt, 64289 Darmstadt, Germany
% Copyright 2001-2008 Andre Vehreschild, Institute for Scientific Computing,
% RWTH Aachen University, 52056 Aachen, Germany.
% Visit us on the web at http://www.adimat.de
% Report bugs to adimat-users@lists.sc.informatik.tu-darmstadt.de
%
%
%                             DISCLAIMER
%
% ADiMat was prepared as part of an employment at the Institute
% for Scientific Computing, RWTH Aachen University, Germany and is
% provided AS IS. NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL
% REPUBLIC OF GERMANY NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
% FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR
% PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
% PRIVATELY OWNED RIGHTS.
%
% Global flags were:
% FORWARDMODE -- Apply the forward mode to the files.
% NOOPEROPTIM -- Do not use optimized operators. I.e.:
%		 g_a*b*g_c -/-> mtimes3(g_a, b, g_c)
% NOLOCALCSE  -- Do not use local common subexpression elimination when
%		 canonicalizing the code.
% NOGLOBALCSE -- Prevents the application of global common subexpression
%		 elimination after canonicalizing the code.
% NOPRESCALARFOLDING -- Switch off folding of scalar constants before
%		 augmentation.
% NOPOSTSCALARFOLDING -- Switch off folding of scalar constants after
%		 augmentation.
% NOCONSTFOLDMULT0 -- Switch off folding of product with one factor
%		 being zero: b*0=0.
% FUNCMODE    -- Inputfile is a function (This flag can not be set explicitly).
% NOTMPCLEAR  -- Suppress generation of clear g_* instructions.
% UNBOUND_XML  -- Write list of unbound identifiers in XML format.
% DEPENDENCIES_XML  -- Write list of functions in XML format.
% UNBOUND_ERROR	-- Stop with error if unbound identifiers found (default).
% FUNCTION_LIST_XML	-- Write list of functions to XML file.
% VERBOSITYLEVEL=5
% AD_IVARS= x
% AD_DVARS= F

function [g_F, F]= g_MinpackHHD_F(g_x, x, Prob)
   % MinpackHHD_F: Computes function for the Human Heart Dipole (HHD) problem
   %               from the MINPACK-2 collection.  
   % USE:
   %         F=MinpackHHD_F(x,Prob)
   % where
   %   x(8,1) : is an arbitrary vector x
   %   Prob   : structure returned by MinpackHHD_Prob
   %
   % Returns
   %   F        : function value
   
   % AUTHOR: S.A.Forth & K. Lenton
   % DATE: 26/06/06
   % Copyright 2006-2009: S.A. Forth, Cranfield University
   % REVISIONS:
   % DATE  WHO   WHAT
   % 30/6/09 SAF tidied some coding
   
   % Original Fortran Header comments follow
   %      character*(*) task, prob
   %      integer n, ldfjac
   %   originally
   %      double precision x(n), fvec(n), fjac(ldfjac,n)
   %   changed to
   %      double precision x(n), fvec(n), fjac(n,n)
   % **********
   %
   % Subroutine dhhdfj
   %
   % This subroutine computes the function and the Jacobian matrix of
   % the human heart dipole problem.
   %
   % The subroutine statement is
   %
   %   subroutine dhhdfj(n,x,fvec,fjac,ldfjac,task,prob)
   %
   % where
   %
   %   n is an integer variable.
   %     On entry n is the number of variables. n = 8.
   %     On exit n is unchanged.
   %
   %   x is a double precision array of dimension n.
   %     On entry x specifies the vector x if task = 'F', 'J', or 'FJ'.
   %        Otherwise x need not be specified.
   %     On exit x is unchanged if task = 'F', 'J', or 'FJ'. Otherwise
   %        x is set according to task.
   %
   %   fvec is a double precision array of dimension n.
   %     On entry fvec need not be specified.
   %     On exit fvec contains the function evaluated at x if
   %        task = 'F' or 'FJ'.
   %
   %   fjac is a double precision array of dimension (ldfjac,n).
   %     On entry fjac need not be specified.
   %     On exit fjac contains the Jacobian matrix evaluated at x if
   %        task = 'J' or 'FJ'.
   %
   %   ldfjac is an integer variable.
   %      On entry ldfjac is the leading dimension of fjac.
   %      On exit ldfjac is unchanged.
   %
   %   task is a character variable.
   %     On entry task specifies the action of the subroutine:
   %
   %        task               action
   %        ----               ------
   %         'F'     Evaluate the function at x.
   %         'J'     Evaluate the Jacobian matrix at x.
   %         'FJ'    Evaluate the function and the Jacobian at x.
   %         'XS'    Set x to the standard starting point xs.
   %         'XL'    Set x to the lower bound xl.
   %         'XU'    Set x to the upper bound xu.
   %
   %     On exit task is unchanged.
   %
   %   prob is a character*5 variable.
   %     On entry prob specifies the version of the problem. The
   %        experiment label is the same as in Dennis, Gay, and Vu.
   %
   %               prob             experiment
   %               ----             ---------
   %              'DHHD1'             791129
   %              'DHHD2'             791226
   %              'DHHD3'              0121a
   %              'DHHD4'              0121b
   %              'DHHD5'              0121c
   %
   %     On exit prob is unchanged.
   %
   % MINPACK-2 Project. November 1993.
   % Argonne National Laboratory and University of Minnesota.
   % Brett M. Averick.
   %
   % **********
   n= length(x); 
   if n~= 8
      error(['length(x) =', n, ' is illegal, must = 8'])
   end
   % unpack constants
   summx= Prob.user.summx; 
   summy= Prob.user.summy; 
   suma= Prob.user.suma; 
   sumb= Prob.user.sumb; 
   sumc= Prob.user.sumc; 
   sumd= Prob.user.sumd; 
   sume= Prob.user.sume; 
   sumf= Prob.user.sumf; 
   % function and Jacobian value
   g_tmp_x_00000= g_x(1);
   tmp_x_00000= x(1);
   g_a= g_tmp_x_00000;
   a= tmp_x_00000; 
   g_tmp_x_00001= g_x(2);
   tmp_x_00001= x(2);
   g_b= g_tmp_x_00001;
   b= tmp_x_00001; 
   g_tmp_x_00002= g_x(3);
   tmp_x_00002= x(3);
   g_c= g_tmp_x_00002;
   c= tmp_x_00002; 
   g_tmp_x_00003= g_x(4);
   tmp_x_00003= x(4);
   g_d= g_tmp_x_00003;
   d= tmp_x_00003; 
   g_tmp_x_00004= g_x(5);
   tmp_x_00004= x(5);
   g_t= g_tmp_x_00004;
   t= tmp_x_00004; 
   g_tmp_x_00005= g_x(6);
   tmp_x_00005= x(6);
   g_u= g_tmp_x_00005;
   u= tmp_x_00005; 
   g_tmp_x_00006= g_x(7);
   tmp_x_00006= x(7);
   g_v= g_tmp_x_00006;
   v= tmp_x_00006; 
   g_tmp_x_00007= g_x(8);
   tmp_x_00007= x(8);
   g_w= g_tmp_x_00007;
   w= tmp_x_00007; 
   g_tv= g_t* v+ t* g_v;
   tv= t* v; 
   g_tt= g_t* t+ t* g_t;
   tt= t* t; 
   g_vv= g_v* v+ v* g_v;
   vv= v* v; 
   g_tsvs= g_tt- g_vv;
   tsvs= tt- vv; 
   g_tmp_MinpackHHD_F_00000= 3* g_vv;
   tmp_MinpackHHD_F_00000= 3* vv;
   g_ts3vs= g_tt- g_tmp_MinpackHHD_F_00000;
   ts3vs= tt- tmp_MinpackHHD_F_00000; 
   g_tmp_MinpackHHD_F_00001= 3* g_tt;
   tmp_MinpackHHD_F_00001= 3* tt;
   g_vs3ts= g_vv- g_tmp_MinpackHHD_F_00001;
   vs3ts= vv- tmp_MinpackHHD_F_00001; 
   g_uw= g_u* w+ u* g_w;
   uw= u* w; 
   g_uu= g_u* u+ u* g_u;
   uu= u* u; 
   g_ww= g_w* w+ w* g_w;
   ww= w* w; 
   g_usws= g_uu- g_ww;
   usws= uu- ww; 
   g_tmp_MinpackHHD_F_00002= 3* g_ww;
   tmp_MinpackHHD_F_00002= 3* ww;
   g_us3ws= g_uu- g_tmp_MinpackHHD_F_00002;
   us3ws= uu- tmp_MinpackHHD_F_00002; 
   g_tmp_MinpackHHD_F_00003= 3* g_uu;
   tmp_MinpackHHD_F_00003= 3* uu;
   g_ws3us= g_ww- g_tmp_MinpackHHD_F_00003;
   ws3us= ww- tmp_MinpackHHD_F_00003; 
   % % Evaluate the function
   g_tmp_MinpackHHD_F_00004= g_a+ g_b+ g_zeros(size(summx));
   tmp_MinpackHHD_F_00004= a+ b- summx;
   g_tmp_MinpackHHD_F_00005= g_c+ g_d+ g_zeros(size(summy));
   tmp_MinpackHHD_F_00005= c+ d- summy;
   g_tmp_MinpackHHD_F_00006= g_t* a+ t* g_a;
   tmp_MinpackHHD_F_00006= t* a;
   g_tmp_MinpackHHD_F_00007= g_u* b+ u* g_b;
   tmp_MinpackHHD_F_00007= u* b;
   g_tmp_MinpackHHD_F_00008= g_v* c+ v* g_c;
   tmp_MinpackHHD_F_00008= v* c;
   g_tmp_MinpackHHD_F_00009= g_w* d+ w* g_d;
   tmp_MinpackHHD_F_00009= w* d;
   g_tmp_MinpackHHD_F_00010= g_tmp_MinpackHHD_F_00006+ g_tmp_MinpackHHD_F_00007- g_tmp_MinpackHHD_F_00008- g_tmp_MinpackHHD_F_00009+ g_zeros(size(suma));
   tmp_MinpackHHD_F_00010= tmp_MinpackHHD_F_00006+ tmp_MinpackHHD_F_00007- tmp_MinpackHHD_F_00008- tmp_MinpackHHD_F_00009- suma;
   g_tmp_MinpackHHD_F_00011= g_v* a+ v* g_a;
   tmp_MinpackHHD_F_00011= v* a;
   g_tmp_MinpackHHD_F_00012= g_w* b+ w* g_b;
   tmp_MinpackHHD_F_00012= w* b;
   g_tmp_MinpackHHD_F_00013= g_t* c+ t* g_c;
   tmp_MinpackHHD_F_00013= t* c;
   g_tmp_MinpackHHD_F_00014= g_u* d+ u* g_d;
   tmp_MinpackHHD_F_00014= u* d;
   g_tmp_MinpackHHD_F_00015= g_tmp_MinpackHHD_F_00011+ g_tmp_MinpackHHD_F_00012+ g_tmp_MinpackHHD_F_00013+ g_tmp_MinpackHHD_F_00014+ g_zeros(size(sumb));
   tmp_MinpackHHD_F_00015= tmp_MinpackHHD_F_00011+ tmp_MinpackHHD_F_00012+ tmp_MinpackHHD_F_00013+ tmp_MinpackHHD_F_00014- sumb;
   g_tmp_MinpackHHD_F_00016= g_a* tsvs+ a* g_tsvs;
   tmp_MinpackHHD_F_00016= a* tsvs;
   g_tmp_MinpackHHD_F_00017= 2* g_c* tv+ 2* c* g_tv;
   tmp_MinpackHHD_F_00017= 2* c* tv;
   g_tmp_MinpackHHD_F_00018= g_b* usws+ b* g_usws;
   tmp_MinpackHHD_F_00018= b* usws;
   g_tmp_MinpackHHD_F_00019= 2* g_d* uw+ 2* d* g_uw;
   tmp_MinpackHHD_F_00019= 2* d* uw;
   g_tmp_MinpackHHD_F_00020= g_tmp_MinpackHHD_F_00016- g_tmp_MinpackHHD_F_00017+ g_tmp_MinpackHHD_F_00018- g_tmp_MinpackHHD_F_00019+ g_zeros(size(sumc));
   tmp_MinpackHHD_F_00020= tmp_MinpackHHD_F_00016- tmp_MinpackHHD_F_00017+ tmp_MinpackHHD_F_00018- tmp_MinpackHHD_F_00019- sumc;
   g_tmp_MinpackHHD_F_00021= g_c* tsvs+ c* g_tsvs;
   tmp_MinpackHHD_F_00021= c* tsvs;
   g_tmp_MinpackHHD_F_00022= 2* g_a* tv+ 2* a* g_tv;
   tmp_MinpackHHD_F_00022= 2* a* tv;
   g_tmp_MinpackHHD_F_00023= g_d* usws+ d* g_usws;
   tmp_MinpackHHD_F_00023= d* usws;
   g_tmp_MinpackHHD_F_00024= 2* g_b* uw+ 2* b* g_uw;
   tmp_MinpackHHD_F_00024= 2* b* uw;
   g_tmp_MinpackHHD_F_00025= g_tmp_MinpackHHD_F_00021+ g_tmp_MinpackHHD_F_00022+ g_tmp_MinpackHHD_F_00023+ g_tmp_MinpackHHD_F_00024+ g_zeros(size(sumd));
   tmp_MinpackHHD_F_00025= tmp_MinpackHHD_F_00021+ tmp_MinpackHHD_F_00022+ tmp_MinpackHHD_F_00023+ tmp_MinpackHHD_F_00024- sumd;
   g_tmp_MinpackHHD_F_00026= g_a* t* ts3vs+ a* g_t* ts3vs+ a* t* g_ts3vs;
   tmp_MinpackHHD_F_00026= a* t* ts3vs;
   g_tmp_MinpackHHD_F_00027= g_c* v* vs3ts+ c* g_v* vs3ts+ c* v* g_vs3ts;
   tmp_MinpackHHD_F_00027= c* v* vs3ts;
   g_tmp_MinpackHHD_F_00028= g_b* u* us3ws+ b* g_u* us3ws+ b* u* g_us3ws;
   tmp_MinpackHHD_F_00028= b* u* us3ws;
   g_tmp_MinpackHHD_F_00029= g_d* w* ws3us+ d* g_w* ws3us+ d* w* g_ws3us;
   tmp_MinpackHHD_F_00029= d* w* ws3us;
   g_tmp_MinpackHHD_F_00030= g_tmp_MinpackHHD_F_00026+ g_tmp_MinpackHHD_F_00027+ g_tmp_MinpackHHD_F_00028+ g_tmp_MinpackHHD_F_00029+ g_zeros(size(sume));
   tmp_MinpackHHD_F_00030= tmp_MinpackHHD_F_00026+ tmp_MinpackHHD_F_00027+ tmp_MinpackHHD_F_00028+ tmp_MinpackHHD_F_00029- sume;
   g_tmp_MinpackHHD_F_00031= g_c* t* ts3vs+ c* g_t* ts3vs+ c* t* g_ts3vs;
   tmp_MinpackHHD_F_00031= c* t* ts3vs;
   g_tmp_MinpackHHD_F_00032= g_a* v* vs3ts+ a* g_v* vs3ts+ a* v* g_vs3ts;
   tmp_MinpackHHD_F_00032= a* v* vs3ts;
   g_tmp_MinpackHHD_F_00033= g_d* u* us3ws+ d* g_u* us3ws+ d* u* g_us3ws;
   tmp_MinpackHHD_F_00033= d* u* us3ws;
   g_tmp_MinpackHHD_F_00034= g_b* w* ws3us+ b* g_w* ws3us+ b* w* g_ws3us;
   tmp_MinpackHHD_F_00034= b* w* ws3us;
   g_tmp_MinpackHHD_F_00035= g_tmp_MinpackHHD_F_00031- g_tmp_MinpackHHD_F_00032+ g_tmp_MinpackHHD_F_00033- g_tmp_MinpackHHD_F_00034+ g_zeros(size(sumf));
   tmp_MinpackHHD_F_00035= tmp_MinpackHHD_F_00031- tmp_MinpackHHD_F_00032+ tmp_MinpackHHD_F_00033- tmp_MinpackHHD_F_00034- sumf;
   g_F= [g_tmp_MinpackHHD_F_00004;
      ;
      g_tmp_MinpackHHD_F_00005;
      ;
      g_tmp_MinpackHHD_F_00010;
      ;
      g_tmp_MinpackHHD_F_00015;
      ;
      g_tmp_MinpackHHD_F_00020;
      ;
      g_tmp_MinpackHHD_F_00025;
      ;
      g_tmp_MinpackHHD_F_00030;
      ;
      g_tmp_MinpackHHD_F_00035];
   F= [tmp_MinpackHHD_F_00004;
      ;
      tmp_MinpackHHD_F_00005;
      ;
      tmp_MinpackHHD_F_00010;
      ;
      tmp_MinpackHHD_F_00015;
      ;
      tmp_MinpackHHD_F_00020;
      ;
      tmp_MinpackHHD_F_00025;
      ;
      tmp_MinpackHHD_F_00030;
      ;
      tmp_MinpackHHD_F_00035]; 
end
