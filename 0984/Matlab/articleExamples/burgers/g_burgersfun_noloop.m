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
% AD_IVARS= y
% AD_DVARS= out

function [g_out, out]= g_burgersfun_noloop(t, g_y, y, N)
   %   Jacek Kierzenka and Lawrence F. Shampine
   %   Copyright 1984-2014 The MathWorks, Inc.
   % Derivative function.
   tmp_burgersfun_noloop_00000= 1: N;
   g_tmp_y_00000= g_y(tmp_burgersfun_noloop_00000);
   tmp_y_00000= y(tmp_burgersfun_noloop_00000);
   g_a= g_tmp_y_00000;
   a= tmp_y_00000; 
   g_tmp_a_00000= g_a(: );
   tmp_a_00000= a(: );
   g_a= g_tmp_a_00000;
   a= tmp_a_00000; 
   g_tmp_y_00001= g_y(N+ 1: end);
   tmp_y_00001= y(N+ 1: end);
   g_x= g_tmp_y_00001;
   x= tmp_y_00001; 
   g_tmp_x_00000= g_x(: );
   tmp_x_00000= x(: );
   g_x= g_tmp_x_00000;
   x= tmp_x_00000; 
   x0= 0; 
   a0= 0; 
   xNP1= 1; 
   aNP1= 0; 
   g_tmp_y_00002= g_y(1);
   tmp_y_00002= y(1);
   g_y1= g_tmp_y_00002;
   y1= tmp_y_00002; 
   g_g= zeros(2* N, 1).* g_y1;
   g= zeros(2* N, 1).* y1; 
   i= 2: N- 1; 
   %for i = 2:N-1
   tmp_burgersfun_noloop_00001= i+ 1;
   g_tmp_x_00001= g_x(tmp_burgersfun_noloop_00001);
   tmp_x_00001= x(tmp_burgersfun_noloop_00001);
   tmp_burgersfun_noloop_00002= i- 1;
   g_tmp_x_00002= g_x(tmp_burgersfun_noloop_00002);
   tmp_x_00002= x(tmp_burgersfun_noloop_00002);
   g_delx= g_tmp_x_00001- g_tmp_x_00002;
   delx= tmp_x_00001- tmp_x_00002; 
   tmp_burgersfun_noloop_00003= i+ 1;
   g_tmp_a_00001= g_a(tmp_burgersfun_noloop_00003);
   tmp_a_00001= a(tmp_burgersfun_noloop_00003);
   g_tmp_a_00002= g_a(i);
   tmp_a_00002= a(i);
   g_tmp_burgersfun_noloop_00004= g_tmp_a_00001- g_tmp_a_00002;
   tmp_burgersfun_noloop_00004= tmp_a_00001- tmp_a_00002;
   tmp_burgersfun_noloop_00005= i+ 1;
   g_tmp_x_00003= g_x(tmp_burgersfun_noloop_00005);
   tmp_x_00003= x(tmp_burgersfun_noloop_00005);
   g_tmp_x_00004= g_x(i);
   tmp_x_00004= x(i);
   g_tmp_burgersfun_noloop_00006= g_tmp_x_00003- g_tmp_x_00004;
   tmp_burgersfun_noloop_00006= tmp_x_00003- tmp_x_00004;
   g_tmp_burgersfun_noloop_00007= (g_tmp_burgersfun_noloop_00004.* tmp_burgersfun_noloop_00006- tmp_burgersfun_noloop_00004.* g_tmp_burgersfun_noloop_00006)./ tmp_burgersfun_noloop_00006.^ 2;
   tmp_burgersfun_noloop_00007= tmp_burgersfun_noloop_00004./ tmp_burgersfun_noloop_00006;
   g_tmp_a_00003= g_a(i);
   tmp_a_00003= a(i);
   tmp_burgersfun_noloop_00008= i- 1;
   g_tmp_a_00004= g_a(tmp_burgersfun_noloop_00008);
   tmp_a_00004= a(tmp_burgersfun_noloop_00008);
   g_tmp_burgersfun_noloop_00009= g_tmp_a_00003- g_tmp_a_00004;
   tmp_burgersfun_noloop_00009= tmp_a_00003- tmp_a_00004;
   g_tmp_x_00005= g_x(i);
   tmp_x_00005= x(i);
   tmp_burgersfun_noloop_00010= i- 1;
   g_tmp_x_00006= g_x(tmp_burgersfun_noloop_00010);
   tmp_x_00006= x(tmp_burgersfun_noloop_00010);
   g_tmp_burgersfun_noloop_00011= g_tmp_x_00005- g_tmp_x_00006;
   tmp_burgersfun_noloop_00011= tmp_x_00005- tmp_x_00006;
   g_tmp_burgersfun_noloop_00012= (g_tmp_burgersfun_noloop_00009.* tmp_burgersfun_noloop_00011- tmp_burgersfun_noloop_00009.* g_tmp_burgersfun_noloop_00011)./ tmp_burgersfun_noloop_00011.^ 2;
   tmp_burgersfun_noloop_00012= tmp_burgersfun_noloop_00009./ tmp_burgersfun_noloop_00011;
   g_tmp_burgersfun_noloop_00013= g_tmp_burgersfun_noloop_00007- g_tmp_burgersfun_noloop_00012;
   tmp_burgersfun_noloop_00013= tmp_burgersfun_noloop_00007- tmp_burgersfun_noloop_00012;
   g_tmp_burgersfun_noloop_00014= 1e-4* g_tmp_burgersfun_noloop_00013;
   tmp_burgersfun_noloop_00014= 1e-4* tmp_burgersfun_noloop_00013;
   g_tmp_burgersfun_noloop_00015= 0.5* g_delx;
   tmp_burgersfun_noloop_00015= 0.5* delx;
   g_tmp_burgersfun_noloop_00016= (g_tmp_burgersfun_noloop_00014.* tmp_burgersfun_noloop_00015- tmp_burgersfun_noloop_00014.* g_tmp_burgersfun_noloop_00015)./ tmp_burgersfun_noloop_00015.^ 2;
   tmp_burgersfun_noloop_00016= tmp_burgersfun_noloop_00014./ tmp_burgersfun_noloop_00015;
   tmp_burgersfun_noloop_00017= i+ 1;
   g_tmp_a_00005= g_a(tmp_burgersfun_noloop_00017);
   tmp_a_00005= a(tmp_burgersfun_noloop_00017);
   g_tmp_burgersfun_noloop_00018= 2.* tmp_a_00005.^ (2- 1).* g_tmp_a_00005;
   tmp_burgersfun_noloop_00018= tmp_a_00005.^ 2;
   tmp_burgersfun_noloop_00019= i- 1;
   g_tmp_a_00006= g_a(tmp_burgersfun_noloop_00019);
   tmp_a_00006= a(tmp_burgersfun_noloop_00019);
   g_tmp_burgersfun_noloop_00020= 2.* tmp_a_00006.^ (2- 1).* g_tmp_a_00006;
   tmp_burgersfun_noloop_00020= tmp_a_00006.^ 2;
   g_tmp_burgersfun_noloop_00021= g_tmp_burgersfun_noloop_00018- g_tmp_burgersfun_noloop_00020;
   tmp_burgersfun_noloop_00021= tmp_burgersfun_noloop_00018- tmp_burgersfun_noloop_00020;
   g_tmp_burgersfun_noloop_00022= 0.5* g_tmp_burgersfun_noloop_00021;
   tmp_burgersfun_noloop_00022= 0.5* tmp_burgersfun_noloop_00021;
   g_tmp_burgersfun_noloop_00023= (g_tmp_burgersfun_noloop_00022.* delx- tmp_burgersfun_noloop_00022.* g_delx)./ delx.^ 2;
   tmp_burgersfun_noloop_00023= tmp_burgersfun_noloop_00022./ delx;
   g_g(i)= g_tmp_burgersfun_noloop_00016- g_tmp_burgersfun_noloop_00023;
   g(i)= tmp_burgersfun_noloop_00016- tmp_burgersfun_noloop_00023; 
   %end
   g_tmp_x_00007= g_x(2);
   tmp_x_00007= x(2);
   g_delx= g_tmp_x_00007+ g_zeros(size(x0));
   delx= tmp_x_00007- x0; 
   g_tmp_a_00007= g_a(2);
   tmp_a_00007= a(2);
   g_tmp_a_00008= g_a(1);
   tmp_a_00008= a(1);
   g_tmp_burgersfun_noloop_00024= g_tmp_a_00007- g_tmp_a_00008;
   tmp_burgersfun_noloop_00024= tmp_a_00007- tmp_a_00008;
   g_tmp_x_00008= g_x(2);
   tmp_x_00008= x(2);
   g_tmp_x_00009= g_x(1);
   tmp_x_00009= x(1);
   g_tmp_burgersfun_noloop_00025= g_tmp_x_00008- g_tmp_x_00009;
   tmp_burgersfun_noloop_00025= tmp_x_00008- tmp_x_00009;
   g_tmp_burgersfun_noloop_00026= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00024), (tmp_burgersfun_noloop_00024), (g_tmp_burgersfun_noloop_00025), (tmp_burgersfun_noloop_00025));
   tmp_burgersfun_noloop_00026= tmp_burgersfun_noloop_00024/ tmp_burgersfun_noloop_00025;
   g_tmp_a_00009= g_a(1);
   tmp_a_00009= a(1);
   g_tmp_burgersfun_noloop_00027= g_tmp_a_00009+ g_zeros(size(a0));
   tmp_burgersfun_noloop_00027= tmp_a_00009- a0;
   g_tmp_x_00010= g_x(1);
   tmp_x_00010= x(1);
   g_tmp_burgersfun_noloop_00028= g_tmp_x_00010+ g_zeros(size(x0));
   tmp_burgersfun_noloop_00028= tmp_x_00010- x0;
   g_tmp_burgersfun_noloop_00029= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00027), (tmp_burgersfun_noloop_00027), (g_tmp_burgersfun_noloop_00028), (tmp_burgersfun_noloop_00028));
   tmp_burgersfun_noloop_00029= tmp_burgersfun_noloop_00027/ tmp_burgersfun_noloop_00028;
   g_tmp_burgersfun_noloop_00030= g_tmp_burgersfun_noloop_00026- g_tmp_burgersfun_noloop_00029;
   tmp_burgersfun_noloop_00030= tmp_burgersfun_noloop_00026- tmp_burgersfun_noloop_00029;
   g_tmp_burgersfun_noloop_00031= 1e-4* g_tmp_burgersfun_noloop_00030;
   tmp_burgersfun_noloop_00031= 1e-4* tmp_burgersfun_noloop_00030;
   g_tmp_burgersfun_noloop_00032= 0.5* g_delx;
   tmp_burgersfun_noloop_00032= 0.5* delx;
   g_tmp_burgersfun_noloop_00033= adimat_g_mrdivide(g_tmp_burgersfun_noloop_00031, tmp_burgersfun_noloop_00031, (g_tmp_burgersfun_noloop_00032), (tmp_burgersfun_noloop_00032));
   tmp_burgersfun_noloop_00033= tmp_burgersfun_noloop_00031/ tmp_burgersfun_noloop_00032;
   g_tmp_a_00010= g_a(2);
   tmp_a_00010= a(2);
   g_tmp_burgersfun_noloop_00034= adimat_g_pow_left(g_tmp_a_00010, tmp_a_00010, 2);
   tmp_burgersfun_noloop_00034= tmp_a_00010^ 2;
   tmp_burgersfun_noloop_00035= a0^ 2;
   g_tmp_burgersfun_noloop_00036= g_tmp_burgersfun_noloop_00034+ g_zeros(size(tmp_burgersfun_noloop_00035));
   tmp_burgersfun_noloop_00036= tmp_burgersfun_noloop_00034- tmp_burgersfun_noloop_00035;
   g_tmp_burgersfun_noloop_00037= 0.5* g_tmp_burgersfun_noloop_00036;
   tmp_burgersfun_noloop_00037= 0.5* tmp_burgersfun_noloop_00036;
   g_tmp_burgersfun_noloop_00038= adimat_g_mrdivide(g_tmp_burgersfun_noloop_00037, tmp_burgersfun_noloop_00037, g_delx, delx);
   tmp_burgersfun_noloop_00038= tmp_burgersfun_noloop_00037/ delx;
   g_g(1)= g_tmp_burgersfun_noloop_00033- g_tmp_burgersfun_noloop_00038;
   g(1)= tmp_burgersfun_noloop_00033- tmp_burgersfun_noloop_00038; 
   tmp_burgersfun_noloop_00039= N- 1;
   g_tmp_x_00011= g_x(tmp_burgersfun_noloop_00039);
   tmp_x_00011= x(tmp_burgersfun_noloop_00039);
   g_delx= -g_tmp_x_00011+ g_zeros(size(xNP1));
   delx= xNP1- tmp_x_00011; 
   g_tmp_a_00011= g_a(N);
   tmp_a_00011= a(N);
   g_tmp_burgersfun_noloop_00040= -g_tmp_a_00011+ g_zeros(size(aNP1));
   tmp_burgersfun_noloop_00040= aNP1- tmp_a_00011;
   g_tmp_x_00012= g_x(N);
   tmp_x_00012= x(N);
   g_tmp_burgersfun_noloop_00041= -g_tmp_x_00012+ g_zeros(size(xNP1));
   tmp_burgersfun_noloop_00041= xNP1- tmp_x_00012;
   g_tmp_burgersfun_noloop_00042= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00040), (tmp_burgersfun_noloop_00040), (g_tmp_burgersfun_noloop_00041), (tmp_burgersfun_noloop_00041));
   tmp_burgersfun_noloop_00042= tmp_burgersfun_noloop_00040/ tmp_burgersfun_noloop_00041;
   g_tmp_a_00012= g_a(N);
   tmp_a_00012= a(N);
   tmp_burgersfun_noloop_00043= N- 1;
   g_tmp_a_00013= g_a(tmp_burgersfun_noloop_00043);
   tmp_a_00013= a(tmp_burgersfun_noloop_00043);
   g_tmp_burgersfun_noloop_00044= g_tmp_a_00012- g_tmp_a_00013;
   tmp_burgersfun_noloop_00044= tmp_a_00012- tmp_a_00013;
   g_tmp_x_00013= g_x(N);
   tmp_x_00013= x(N);
   tmp_burgersfun_noloop_00045= N- 1;
   g_tmp_x_00014= g_x(tmp_burgersfun_noloop_00045);
   tmp_x_00014= x(tmp_burgersfun_noloop_00045);
   g_tmp_burgersfun_noloop_00046= g_tmp_x_00013- g_tmp_x_00014;
   tmp_burgersfun_noloop_00046= tmp_x_00013- tmp_x_00014;
   g_tmp_burgersfun_noloop_00047= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00044), (tmp_burgersfun_noloop_00044), (g_tmp_burgersfun_noloop_00046), (tmp_burgersfun_noloop_00046));
   tmp_burgersfun_noloop_00047= tmp_burgersfun_noloop_00044/ tmp_burgersfun_noloop_00046;
   g_tmp_burgersfun_noloop_00048= g_tmp_burgersfun_noloop_00042- g_tmp_burgersfun_noloop_00047;
   tmp_burgersfun_noloop_00048= tmp_burgersfun_noloop_00042- tmp_burgersfun_noloop_00047;
   g_tmp_burgersfun_noloop_00049= 1e-4* g_tmp_burgersfun_noloop_00048;
   tmp_burgersfun_noloop_00049= 1e-4* tmp_burgersfun_noloop_00048;
   g_tmp_burgersfun_noloop_00050= adimat_g_mrdivide(g_tmp_burgersfun_noloop_00049, tmp_burgersfun_noloop_00049, g_delx, delx);
   tmp_burgersfun_noloop_00050= tmp_burgersfun_noloop_00049/ delx;
   tmp_burgersfun_noloop_00051= aNP1^ 2;
   tmp_burgersfun_noloop_00052= N- 1;
   g_tmp_a_00014= g_a(tmp_burgersfun_noloop_00052);
   tmp_a_00014= a(tmp_burgersfun_noloop_00052);
   g_tmp_burgersfun_noloop_00053= adimat_g_pow_left(g_tmp_a_00014, tmp_a_00014, 2);
   tmp_burgersfun_noloop_00053= tmp_a_00014^ 2;
   g_tmp_burgersfun_noloop_00054= -g_tmp_burgersfun_noloop_00053+ g_zeros(size(tmp_burgersfun_noloop_00051));
   tmp_burgersfun_noloop_00054= tmp_burgersfun_noloop_00051- tmp_burgersfun_noloop_00053;
   g_tmp_burgersfun_noloop_00055= 0.5* g_tmp_burgersfun_noloop_00054;
   tmp_burgersfun_noloop_00055= 0.5* tmp_burgersfun_noloop_00054;
   g_tmp_burgersfun_noloop_00056= adimat_g_mrdivide(g_tmp_burgersfun_noloop_00055, tmp_burgersfun_noloop_00055, g_delx, delx);
   tmp_burgersfun_noloop_00056= tmp_burgersfun_noloop_00055/ delx;
   g_g(N)= g_tmp_burgersfun_noloop_00050- g_tmp_burgersfun_noloop_00056;
   g(N)= tmp_burgersfun_noloop_00050- tmp_burgersfun_noloop_00056; 
   % Evaluate monitor function.
   g_M= zeros(N, 1).* g_y1;
   M= zeros(N, 1).* y1; 
   i= 2: N- 1; 
   %for i = 2:N-1
   tmp_burgersfun_noloop_00057= i+ 1;
   g_tmp_a_00015= g_a(tmp_burgersfun_noloop_00057);
   tmp_a_00015= a(tmp_burgersfun_noloop_00057);
   tmp_burgersfun_noloop_00058= i- 1;
   g_tmp_a_00016= g_a(tmp_burgersfun_noloop_00058);
   tmp_a_00016= a(tmp_burgersfun_noloop_00058);
   g_tmp_burgersfun_noloop_00059= g_tmp_a_00015- g_tmp_a_00016;
   tmp_burgersfun_noloop_00059= tmp_a_00015- tmp_a_00016;
   tmp_burgersfun_noloop_00060= i+ 1;
   g_tmp_x_00015= g_x(tmp_burgersfun_noloop_00060);
   tmp_x_00015= x(tmp_burgersfun_noloop_00060);
   tmp_burgersfun_noloop_00061= i- 1;
   g_tmp_x_00016= g_x(tmp_burgersfun_noloop_00061);
   tmp_x_00016= x(tmp_burgersfun_noloop_00061);
   g_tmp_burgersfun_noloop_00062= g_tmp_x_00015- g_tmp_x_00016;
   tmp_burgersfun_noloop_00062= tmp_x_00015- tmp_x_00016;
   g_tmp_burgersfun_noloop_00063= (g_tmp_burgersfun_noloop_00059.* tmp_burgersfun_noloop_00062- tmp_burgersfun_noloop_00059.* g_tmp_burgersfun_noloop_00062)./ tmp_burgersfun_noloop_00062.^ 2;
   tmp_burgersfun_noloop_00063= tmp_burgersfun_noloop_00059./ tmp_burgersfun_noloop_00062;
   g_tmp_burgersfun_noloop_00064= 2.* tmp_burgersfun_noloop_00063.^ (2- 1).* g_tmp_burgersfun_noloop_00063;
   tmp_burgersfun_noloop_00064= tmp_burgersfun_noloop_00063.^ 2;
   g_tmp_burgersfun_noloop_00065= g_tmp_burgersfun_noloop_00064+ g_zeros(1);
   tmp_burgersfun_noloop_00065= 1+ tmp_burgersfun_noloop_00064;
   M(i)= sqrt(tmp_burgersfun_noloop_00065); 
   %end
   g_M(i)= g_tmp_burgersfun_noloop_00065./ (2.* M(i));
   g_tmp_a_00017= g_a(1);
   tmp_a_00017= a(1);
   g_tmp_burgersfun_noloop_00066= g_tmp_a_00017+ g_zeros(size(a0));
   tmp_burgersfun_noloop_00066= tmp_a_00017- a0;
   g_tmp_x_00017= g_x(1);
   tmp_x_00017= x(1);
   g_tmp_burgersfun_noloop_00067= g_tmp_x_00017+ g_zeros(size(x0));
   tmp_burgersfun_noloop_00067= tmp_x_00017- x0;
   g_tmp_burgersfun_noloop_00068= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00066), (tmp_burgersfun_noloop_00066), (g_tmp_burgersfun_noloop_00067), (tmp_burgersfun_noloop_00067));
   tmp_burgersfun_noloop_00068= tmp_burgersfun_noloop_00066/ tmp_burgersfun_noloop_00067;
   g_tmp_burgersfun_noloop_00069= adimat_g_pow_left((g_tmp_burgersfun_noloop_00068), (tmp_burgersfun_noloop_00068), 2);
   tmp_burgersfun_noloop_00069= tmp_burgersfun_noloop_00068^ 2;
   g_tmp_burgersfun_noloop_00070= g_tmp_burgersfun_noloop_00069+ g_zeros(1);
   tmp_burgersfun_noloop_00070= 1+ tmp_burgersfun_noloop_00069;
   M0= sqrt(tmp_burgersfun_noloop_00070); 
   g_M0= g_tmp_burgersfun_noloop_00070./ (2.* M0);
   g_tmp_a_00018= g_a(2);
   tmp_a_00018= a(2);
   g_tmp_burgersfun_noloop_00071= g_tmp_a_00018+ g_zeros(size(a0));
   tmp_burgersfun_noloop_00071= tmp_a_00018- a0;
   g_tmp_x_00018= g_x(2);
   tmp_x_00018= x(2);
   g_tmp_burgersfun_noloop_00072= g_tmp_x_00018+ g_zeros(size(x0));
   tmp_burgersfun_noloop_00072= tmp_x_00018- x0;
   g_tmp_burgersfun_noloop_00073= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00071), (tmp_burgersfun_noloop_00071), (g_tmp_burgersfun_noloop_00072), (tmp_burgersfun_noloop_00072));
   tmp_burgersfun_noloop_00073= tmp_burgersfun_noloop_00071/ tmp_burgersfun_noloop_00072;
   g_tmp_burgersfun_noloop_00074= adimat_g_pow_left((g_tmp_burgersfun_noloop_00073), (tmp_burgersfun_noloop_00073), 2);
   tmp_burgersfun_noloop_00074= tmp_burgersfun_noloop_00073^ 2;
   g_tmp_burgersfun_noloop_00075= g_tmp_burgersfun_noloop_00074+ g_zeros(1);
   tmp_burgersfun_noloop_00075= 1+ tmp_burgersfun_noloop_00074;
   M(1)= sqrt(tmp_burgersfun_noloop_00075); 
   g_M(1)= g_tmp_burgersfun_noloop_00075./ (2.* M(1));
   tmp_burgersfun_noloop_00076= N- 1;
   g_tmp_a_00019= g_a(tmp_burgersfun_noloop_00076);
   tmp_a_00019= a(tmp_burgersfun_noloop_00076);
   g_tmp_burgersfun_noloop_00077= -g_tmp_a_00019+ g_zeros(size(aNP1));
   tmp_burgersfun_noloop_00077= aNP1- tmp_a_00019;
   tmp_burgersfun_noloop_00078= N- 1;
   g_tmp_x_00019= g_x(tmp_burgersfun_noloop_00078);
   tmp_x_00019= x(tmp_burgersfun_noloop_00078);
   g_tmp_burgersfun_noloop_00079= -g_tmp_x_00019+ g_zeros(size(xNP1));
   tmp_burgersfun_noloop_00079= xNP1- tmp_x_00019;
   g_tmp_burgersfun_noloop_00080= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00077), (tmp_burgersfun_noloop_00077), (g_tmp_burgersfun_noloop_00079), (tmp_burgersfun_noloop_00079));
   tmp_burgersfun_noloop_00080= tmp_burgersfun_noloop_00077/ tmp_burgersfun_noloop_00079;
   g_tmp_burgersfun_noloop_00081= adimat_g_pow_left((g_tmp_burgersfun_noloop_00080), (tmp_burgersfun_noloop_00080), 2);
   tmp_burgersfun_noloop_00081= tmp_burgersfun_noloop_00080^ 2;
   g_tmp_burgersfun_noloop_00082= g_tmp_burgersfun_noloop_00081+ g_zeros(1);
   tmp_burgersfun_noloop_00082= 1+ tmp_burgersfun_noloop_00081;
   M(N)= sqrt(tmp_burgersfun_noloop_00082); 
   g_M(N)= g_tmp_burgersfun_noloop_00082./ (2.* M(N));
   g_tmp_a_00020= g_a(N);
   tmp_a_00020= a(N);
   g_tmp_burgersfun_noloop_00083= -g_tmp_a_00020+ g_zeros(size(aNP1));
   tmp_burgersfun_noloop_00083= aNP1- tmp_a_00020;
   g_tmp_x_00020= g_x(N);
   tmp_x_00020= x(N);
   g_tmp_burgersfun_noloop_00084= -g_tmp_x_00020+ g_zeros(size(xNP1));
   tmp_burgersfun_noloop_00084= xNP1- tmp_x_00020;
   g_tmp_burgersfun_noloop_00085= adimat_g_mrdivide((g_tmp_burgersfun_noloop_00083), (tmp_burgersfun_noloop_00083), (g_tmp_burgersfun_noloop_00084), (tmp_burgersfun_noloop_00084));
   tmp_burgersfun_noloop_00085= tmp_burgersfun_noloop_00083/ tmp_burgersfun_noloop_00084;
   g_tmp_burgersfun_noloop_00086= adimat_g_pow_left((g_tmp_burgersfun_noloop_00085), (tmp_burgersfun_noloop_00085), 2);
   tmp_burgersfun_noloop_00086= tmp_burgersfun_noloop_00085^ 2;
   g_tmp_burgersfun_noloop_00087= g_tmp_burgersfun_noloop_00086+ g_zeros(1);
   tmp_burgersfun_noloop_00087= 1+ tmp_burgersfun_noloop_00086;
   MNP1= sqrt(tmp_burgersfun_noloop_00087); 
   % Spatial smoothing with gamma = 2, p = 2.
   g_MNP1= g_tmp_burgersfun_noloop_00087./ (2.* MNP1);
   g_SM= zeros(N, 1).* g_y1;
   SM= zeros(N, 1).* y1; 
   i= 3: N- 2; 
   %for i = 3:N-2
   tmp_burgersfun_noloop_00088= i- 2;
   g_tmp_M_00000= g_M(tmp_burgersfun_noloop_00088);
   tmp_M_00000= M(tmp_burgersfun_noloop_00088);
   g_tmp_burgersfun_noloop_00089= 2.* tmp_M_00000.^ (2- 1).* g_tmp_M_00000;
   tmp_burgersfun_noloop_00089= tmp_M_00000.^ 2;
   g_tmp_burgersfun_noloop_00090= 4* g_tmp_burgersfun_noloop_00089;
   tmp_burgersfun_noloop_00090= 4* tmp_burgersfun_noloop_00089;
   tmp_burgersfun_noloop_00091= i- 1;
   g_tmp_M_00001= g_M(tmp_burgersfun_noloop_00091);
   tmp_M_00001= M(tmp_burgersfun_noloop_00091);
   g_tmp_burgersfun_noloop_00092= 2.* tmp_M_00001.^ (2- 1).* g_tmp_M_00001;
   tmp_burgersfun_noloop_00092= tmp_M_00001.^ 2;
   g_tmp_burgersfun_noloop_00093= 6* g_tmp_burgersfun_noloop_00092;
   tmp_burgersfun_noloop_00093= 6* tmp_burgersfun_noloop_00092;
   g_tmp_M_00002= g_M(i);
   tmp_M_00002= M(i);
   g_tmp_burgersfun_noloop_00094= 2.* tmp_M_00002.^ (2- 1).* g_tmp_M_00002;
   tmp_burgersfun_noloop_00094= tmp_M_00002.^ 2;
   g_tmp_burgersfun_noloop_00095= 9* g_tmp_burgersfun_noloop_00094;
   tmp_burgersfun_noloop_00095= 9* tmp_burgersfun_noloop_00094;
   tmp_burgersfun_noloop_00096= i+ 1;
   g_tmp_M_00003= g_M(tmp_burgersfun_noloop_00096);
   tmp_M_00003= M(tmp_burgersfun_noloop_00096);
   g_tmp_burgersfun_noloop_00097= 2.* tmp_M_00003.^ (2- 1).* g_tmp_M_00003;
   tmp_burgersfun_noloop_00097= tmp_M_00003.^ 2;
   g_tmp_burgersfun_noloop_00098= 6* g_tmp_burgersfun_noloop_00097;
   tmp_burgersfun_noloop_00098= 6* tmp_burgersfun_noloop_00097;
   tmp_burgersfun_noloop_00099= i+ 2;
   g_tmp_M_00004= g_M(tmp_burgersfun_noloop_00099);
   tmp_M_00004= M(tmp_burgersfun_noloop_00099);
   g_tmp_burgersfun_noloop_00100= 2.* tmp_M_00004.^ (2- 1).* g_tmp_M_00004;
   tmp_burgersfun_noloop_00100= tmp_M_00004.^ 2;
   g_tmp_burgersfun_noloop_00101= 4* g_tmp_burgersfun_noloop_00100;
   tmp_burgersfun_noloop_00101= 4* tmp_burgersfun_noloop_00100;
   g_tmp_burgersfun_noloop_00102= g_tmp_burgersfun_noloop_00090+ g_tmp_burgersfun_noloop_00093+ g_tmp_burgersfun_noloop_00095+ g_tmp_burgersfun_noloop_00098+ g_tmp_burgersfun_noloop_00101;
   tmp_burgersfun_noloop_00102= tmp_burgersfun_noloop_00090+ tmp_burgersfun_noloop_00093+ tmp_burgersfun_noloop_00095+ tmp_burgersfun_noloop_00098+ tmp_burgersfun_noloop_00101;
   g_tmp_burgersfun_noloop_00103= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00102), (tmp_burgersfun_noloop_00102), 29);
   tmp_burgersfun_noloop_00103= tmp_burgersfun_noloop_00102/ 29;
   SM(i)= sqrt(tmp_burgersfun_noloop_00103); 
   %end
   g_SM(i)= g_tmp_burgersfun_noloop_00103./ (2.* SM(i));
   g_tmp_burgersfun_noloop_00104= adimat_g_pow_left(g_M0, M0, 2);
   tmp_burgersfun_noloop_00104= M0^ 2;
   g_tmp_burgersfun_noloop_00105= 9* g_tmp_burgersfun_noloop_00104;
   tmp_burgersfun_noloop_00105= 9* tmp_burgersfun_noloop_00104;
   g_tmp_M_00005= g_M(1);
   tmp_M_00005= M(1);
   g_tmp_burgersfun_noloop_00106= adimat_g_pow_left(g_tmp_M_00005, tmp_M_00005, 2);
   tmp_burgersfun_noloop_00106= tmp_M_00005^ 2;
   g_tmp_burgersfun_noloop_00107= 6* g_tmp_burgersfun_noloop_00106;
   tmp_burgersfun_noloop_00107= 6* tmp_burgersfun_noloop_00106;
   g_tmp_M_00006= g_M(2);
   tmp_M_00006= M(2);
   g_tmp_burgersfun_noloop_00108= adimat_g_pow_left(g_tmp_M_00006, tmp_M_00006, 2);
   tmp_burgersfun_noloop_00108= tmp_M_00006^ 2;
   g_tmp_burgersfun_noloop_00109= 4* g_tmp_burgersfun_noloop_00108;
   tmp_burgersfun_noloop_00109= 4* tmp_burgersfun_noloop_00108;
   g_tmp_burgersfun_noloop_00110= g_tmp_burgersfun_noloop_00105+ g_tmp_burgersfun_noloop_00107+ g_tmp_burgersfun_noloop_00109;
   tmp_burgersfun_noloop_00110= tmp_burgersfun_noloop_00105+ tmp_burgersfun_noloop_00107+ tmp_burgersfun_noloop_00109;
   g_tmp_burgersfun_noloop_00111= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00110), (tmp_burgersfun_noloop_00110), 19);
   tmp_burgersfun_noloop_00111= tmp_burgersfun_noloop_00110/ 19;
   SM0= sqrt(tmp_burgersfun_noloop_00111); 
   g_SM0= g_tmp_burgersfun_noloop_00111./ (2.* SM0);
   g_tmp_burgersfun_noloop_00112= adimat_g_pow_left(g_M0, M0, 2);
   tmp_burgersfun_noloop_00112= M0^ 2;
   g_tmp_burgersfun_noloop_00113= 6* g_tmp_burgersfun_noloop_00112;
   tmp_burgersfun_noloop_00113= 6* tmp_burgersfun_noloop_00112;
   g_tmp_M_00007= g_M(1);
   tmp_M_00007= M(1);
   g_tmp_burgersfun_noloop_00114= adimat_g_pow_left(g_tmp_M_00007, tmp_M_00007, 2);
   tmp_burgersfun_noloop_00114= tmp_M_00007^ 2;
   g_tmp_burgersfun_noloop_00115= 9* g_tmp_burgersfun_noloop_00114;
   tmp_burgersfun_noloop_00115= 9* tmp_burgersfun_noloop_00114;
   g_tmp_M_00008= g_M(2);
   tmp_M_00008= M(2);
   g_tmp_burgersfun_noloop_00116= adimat_g_pow_left(g_tmp_M_00008, tmp_M_00008, 2);
   tmp_burgersfun_noloop_00116= tmp_M_00008^ 2;
   g_tmp_burgersfun_noloop_00117= 6* g_tmp_burgersfun_noloop_00116;
   tmp_burgersfun_noloop_00117= 6* tmp_burgersfun_noloop_00116;
   g_tmp_M_00009= g_M(3);
   tmp_M_00009= M(3);
   g_tmp_burgersfun_noloop_00118= adimat_g_pow_left(g_tmp_M_00009, tmp_M_00009, 2);
   tmp_burgersfun_noloop_00118= tmp_M_00009^ 2;
   g_tmp_burgersfun_noloop_00119= 4* g_tmp_burgersfun_noloop_00118;
   tmp_burgersfun_noloop_00119= 4* tmp_burgersfun_noloop_00118;
   g_tmp_burgersfun_noloop_00120= g_tmp_burgersfun_noloop_00113+ g_tmp_burgersfun_noloop_00115+ g_tmp_burgersfun_noloop_00117+ g_tmp_burgersfun_noloop_00119;
   tmp_burgersfun_noloop_00120= tmp_burgersfun_noloop_00113+ tmp_burgersfun_noloop_00115+ tmp_burgersfun_noloop_00117+ tmp_burgersfun_noloop_00119;
   g_tmp_burgersfun_noloop_00121= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00120), (tmp_burgersfun_noloop_00120), 25);
   tmp_burgersfun_noloop_00121= tmp_burgersfun_noloop_00120/ 25;
   SM(1)= sqrt(tmp_burgersfun_noloop_00121); 
   g_SM(1)= g_tmp_burgersfun_noloop_00121./ (2.* SM(1));
   g_tmp_burgersfun_noloop_00122= adimat_g_pow_left(g_M0, M0, 2);
   tmp_burgersfun_noloop_00122= M0^ 2;
   g_tmp_burgersfun_noloop_00123= 4* g_tmp_burgersfun_noloop_00122;
   tmp_burgersfun_noloop_00123= 4* tmp_burgersfun_noloop_00122;
   g_tmp_M_00010= g_M(1);
   tmp_M_00010= M(1);
   g_tmp_burgersfun_noloop_00124= adimat_g_pow_left(g_tmp_M_00010, tmp_M_00010, 2);
   tmp_burgersfun_noloop_00124= tmp_M_00010^ 2;
   g_tmp_burgersfun_noloop_00125= 6* g_tmp_burgersfun_noloop_00124;
   tmp_burgersfun_noloop_00125= 6* tmp_burgersfun_noloop_00124;
   g_tmp_M_00011= g_M(2);
   tmp_M_00011= M(2);
   g_tmp_burgersfun_noloop_00126= adimat_g_pow_left(g_tmp_M_00011, tmp_M_00011, 2);
   tmp_burgersfun_noloop_00126= tmp_M_00011^ 2;
   g_tmp_burgersfun_noloop_00127= 9* g_tmp_burgersfun_noloop_00126;
   tmp_burgersfun_noloop_00127= 9* tmp_burgersfun_noloop_00126;
   g_tmp_M_00012= g_M(3);
   tmp_M_00012= M(3);
   g_tmp_burgersfun_noloop_00128= adimat_g_pow_left(g_tmp_M_00012, tmp_M_00012, 2);
   tmp_burgersfun_noloop_00128= tmp_M_00012^ 2;
   g_tmp_burgersfun_noloop_00129= 6* g_tmp_burgersfun_noloop_00128;
   tmp_burgersfun_noloop_00129= 6* tmp_burgersfun_noloop_00128;
   g_tmp_M_00013= g_M(4);
   tmp_M_00013= M(4);
   g_tmp_burgersfun_noloop_00130= adimat_g_pow_left(g_tmp_M_00013, tmp_M_00013, 2);
   tmp_burgersfun_noloop_00130= tmp_M_00013^ 2;
   g_tmp_burgersfun_noloop_00131= 4* g_tmp_burgersfun_noloop_00130;
   tmp_burgersfun_noloop_00131= 4* tmp_burgersfun_noloop_00130;
   g_tmp_burgersfun_noloop_00132= g_tmp_burgersfun_noloop_00123+ g_tmp_burgersfun_noloop_00125+ g_tmp_burgersfun_noloop_00127+ g_tmp_burgersfun_noloop_00129+ g_tmp_burgersfun_noloop_00131;
   tmp_burgersfun_noloop_00132= tmp_burgersfun_noloop_00123+ tmp_burgersfun_noloop_00125+ tmp_burgersfun_noloop_00127+ tmp_burgersfun_noloop_00129+ tmp_burgersfun_noloop_00131;
   g_tmp_burgersfun_noloop_00133= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00132), (tmp_burgersfun_noloop_00132), 29);
   tmp_burgersfun_noloop_00133= tmp_burgersfun_noloop_00132/ 29;
   SM(2)= sqrt(tmp_burgersfun_noloop_00133); 
   g_SM(2)= g_tmp_burgersfun_noloop_00133./ (2.* SM(2));
   tmp_burgersfun_noloop_00134= N- 3;
   g_tmp_M_00014= g_M(tmp_burgersfun_noloop_00134);
   tmp_M_00014= M(tmp_burgersfun_noloop_00134);
   g_tmp_burgersfun_noloop_00135= adimat_g_pow_left(g_tmp_M_00014, tmp_M_00014, 2);
   tmp_burgersfun_noloop_00135= tmp_M_00014^ 2;
   g_tmp_burgersfun_noloop_00136= 4* g_tmp_burgersfun_noloop_00135;
   tmp_burgersfun_noloop_00136= 4* tmp_burgersfun_noloop_00135;
   tmp_burgersfun_noloop_00137= N- 2;
   g_tmp_M_00015= g_M(tmp_burgersfun_noloop_00137);
   tmp_M_00015= M(tmp_burgersfun_noloop_00137);
   g_tmp_burgersfun_noloop_00138= adimat_g_pow_left(g_tmp_M_00015, tmp_M_00015, 2);
   tmp_burgersfun_noloop_00138= tmp_M_00015^ 2;
   g_tmp_burgersfun_noloop_00139= 6* g_tmp_burgersfun_noloop_00138;
   tmp_burgersfun_noloop_00139= 6* tmp_burgersfun_noloop_00138;
   tmp_burgersfun_noloop_00140= N- 1;
   g_tmp_M_00016= g_M(tmp_burgersfun_noloop_00140);
   tmp_M_00016= M(tmp_burgersfun_noloop_00140);
   g_tmp_burgersfun_noloop_00141= adimat_g_pow_left(g_tmp_M_00016, tmp_M_00016, 2);
   tmp_burgersfun_noloop_00141= tmp_M_00016^ 2;
   g_tmp_burgersfun_noloop_00142= 9* g_tmp_burgersfun_noloop_00141;
   tmp_burgersfun_noloop_00142= 9* tmp_burgersfun_noloop_00141;
   g_tmp_M_00017= g_M(N);
   tmp_M_00017= M(N);
   g_tmp_burgersfun_noloop_00143= adimat_g_pow_left(g_tmp_M_00017, tmp_M_00017, 2);
   tmp_burgersfun_noloop_00143= tmp_M_00017^ 2;
   g_tmp_burgersfun_noloop_00144= 6* g_tmp_burgersfun_noloop_00143;
   tmp_burgersfun_noloop_00144= 6* tmp_burgersfun_noloop_00143;
   g_tmp_burgersfun_noloop_00145= adimat_g_pow_left(g_MNP1, MNP1, 2);
   tmp_burgersfun_noloop_00145= MNP1^ 2;
   g_tmp_burgersfun_noloop_00146= 4* g_tmp_burgersfun_noloop_00145;
   tmp_burgersfun_noloop_00146= 4* tmp_burgersfun_noloop_00145;
   g_tmp_burgersfun_noloop_00147= g_tmp_burgersfun_noloop_00136+ g_tmp_burgersfun_noloop_00139+ g_tmp_burgersfun_noloop_00142+ g_tmp_burgersfun_noloop_00144+ g_tmp_burgersfun_noloop_00146;
   tmp_burgersfun_noloop_00147= tmp_burgersfun_noloop_00136+ tmp_burgersfun_noloop_00139+ tmp_burgersfun_noloop_00142+ tmp_burgersfun_noloop_00144+ tmp_burgersfun_noloop_00146;
   g_tmp_burgersfun_noloop_00148= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00147), (tmp_burgersfun_noloop_00147), 29);
   tmp_burgersfun_noloop_00148= tmp_burgersfun_noloop_00147/ 29;
   SM(N- 1)= sqrt(tmp_burgersfun_noloop_00148); 
   g_SM(N- 1)= g_tmp_burgersfun_noloop_00148./ (2.* SM(N- 1));
   tmp_burgersfun_noloop_00149= N- 2;
   g_tmp_M_00018= g_M(tmp_burgersfun_noloop_00149);
   tmp_M_00018= M(tmp_burgersfun_noloop_00149);
   g_tmp_burgersfun_noloop_00150= adimat_g_pow_left(g_tmp_M_00018, tmp_M_00018, 2);
   tmp_burgersfun_noloop_00150= tmp_M_00018^ 2;
   g_tmp_burgersfun_noloop_00151= 4* g_tmp_burgersfun_noloop_00150;
   tmp_burgersfun_noloop_00151= 4* tmp_burgersfun_noloop_00150;
   tmp_burgersfun_noloop_00152= N- 1;
   g_tmp_M_00019= g_M(tmp_burgersfun_noloop_00152);
   tmp_M_00019= M(tmp_burgersfun_noloop_00152);
   g_tmp_burgersfun_noloop_00153= adimat_g_pow_left(g_tmp_M_00019, tmp_M_00019, 2);
   tmp_burgersfun_noloop_00153= tmp_M_00019^ 2;
   g_tmp_burgersfun_noloop_00154= 6* g_tmp_burgersfun_noloop_00153;
   tmp_burgersfun_noloop_00154= 6* tmp_burgersfun_noloop_00153;
   g_tmp_M_00020= g_M(N);
   tmp_M_00020= M(N);
   g_tmp_burgersfun_noloop_00155= adimat_g_pow_left(g_tmp_M_00020, tmp_M_00020, 2);
   tmp_burgersfun_noloop_00155= tmp_M_00020^ 2;
   g_tmp_burgersfun_noloop_00156= 9* g_tmp_burgersfun_noloop_00155;
   tmp_burgersfun_noloop_00156= 9* tmp_burgersfun_noloop_00155;
   g_tmp_burgersfun_noloop_00157= adimat_g_pow_left(g_MNP1, MNP1, 2);
   tmp_burgersfun_noloop_00157= MNP1^ 2;
   g_tmp_burgersfun_noloop_00158= 6* g_tmp_burgersfun_noloop_00157;
   tmp_burgersfun_noloop_00158= 6* tmp_burgersfun_noloop_00157;
   g_tmp_burgersfun_noloop_00159= g_tmp_burgersfun_noloop_00151+ g_tmp_burgersfun_noloop_00154+ g_tmp_burgersfun_noloop_00156+ g_tmp_burgersfun_noloop_00158;
   tmp_burgersfun_noloop_00159= tmp_burgersfun_noloop_00151+ tmp_burgersfun_noloop_00154+ tmp_burgersfun_noloop_00156+ tmp_burgersfun_noloop_00158;
   g_tmp_burgersfun_noloop_00160= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00159), (tmp_burgersfun_noloop_00159), 25);
   tmp_burgersfun_noloop_00160= tmp_burgersfun_noloop_00159/ 25;
   SM(N)= sqrt(tmp_burgersfun_noloop_00160); 
   g_SM(N)= g_tmp_burgersfun_noloop_00160./ (2.* SM(N));
   tmp_burgersfun_noloop_00161= N- 1;
   g_tmp_M_00021= g_M(tmp_burgersfun_noloop_00161);
   tmp_M_00021= M(tmp_burgersfun_noloop_00161);
   g_tmp_burgersfun_noloop_00162= adimat_g_pow_left(g_tmp_M_00021, tmp_M_00021, 2);
   tmp_burgersfun_noloop_00162= tmp_M_00021^ 2;
   g_tmp_burgersfun_noloop_00163= 4* g_tmp_burgersfun_noloop_00162;
   tmp_burgersfun_noloop_00163= 4* tmp_burgersfun_noloop_00162;
   g_tmp_M_00022= g_M(N);
   tmp_M_00022= M(N);
   g_tmp_burgersfun_noloop_00164= adimat_g_pow_left(g_tmp_M_00022, tmp_M_00022, 2);
   tmp_burgersfun_noloop_00164= tmp_M_00022^ 2;
   g_tmp_burgersfun_noloop_00165= 6* g_tmp_burgersfun_noloop_00164;
   tmp_burgersfun_noloop_00165= 6* tmp_burgersfun_noloop_00164;
   g_tmp_burgersfun_noloop_00166= adimat_g_pow_left(g_MNP1, MNP1, 2);
   tmp_burgersfun_noloop_00166= MNP1^ 2;
   g_tmp_burgersfun_noloop_00167= 9* g_tmp_burgersfun_noloop_00166;
   tmp_burgersfun_noloop_00167= 9* tmp_burgersfun_noloop_00166;
   g_tmp_burgersfun_noloop_00168= g_tmp_burgersfun_noloop_00163+ g_tmp_burgersfun_noloop_00165+ g_tmp_burgersfun_noloop_00167;
   tmp_burgersfun_noloop_00168= tmp_burgersfun_noloop_00163+ tmp_burgersfun_noloop_00165+ tmp_burgersfun_noloop_00167;
   g_tmp_burgersfun_noloop_00169= adimat_g_mrdivide1((g_tmp_burgersfun_noloop_00168), (tmp_burgersfun_noloop_00168), 19);
   tmp_burgersfun_noloop_00169= tmp_burgersfun_noloop_00168/ 19;
   SMNP1= sqrt(tmp_burgersfun_noloop_00169); 
   g_SMNP1= g_tmp_burgersfun_noloop_00169./ (2.* SMNP1);
   i= 2: N- 1; 
   %for i = 2:N-1
   tmp_burgersfun_noloop_00170= i+ 1;
   g_tmp_SM_00000= g_SM(tmp_burgersfun_noloop_00170);
   tmp_SM_00000= SM(tmp_burgersfun_noloop_00170);
   g_tmp_SM_00001= g_SM(i);
   tmp_SM_00001= SM(i);
   g_tmp_burgersfun_noloop_00171= g_tmp_SM_00000+ g_tmp_SM_00001;
   tmp_burgersfun_noloop_00171= tmp_SM_00000+ tmp_SM_00001;
   tmp_burgersfun_noloop_00172= i+ 1;
   g_tmp_x_00021= g_x(tmp_burgersfun_noloop_00172);
   tmp_x_00021= x(tmp_burgersfun_noloop_00172);
   g_tmp_x_00022= g_x(i);
   tmp_x_00022= x(i);
   g_tmp_burgersfun_noloop_00173= g_tmp_x_00021- g_tmp_x_00022;
   tmp_burgersfun_noloop_00173= tmp_x_00021- tmp_x_00022;
   g_tmp_burgersfun_noloop_00174= g_tmp_burgersfun_noloop_00171.* tmp_burgersfun_noloop_00173+ tmp_burgersfun_noloop_00171.* g_tmp_burgersfun_noloop_00173;
   tmp_burgersfun_noloop_00174= tmp_burgersfun_noloop_00171.* tmp_burgersfun_noloop_00173;
   g_tmp_SM_00002= g_SM(i);
   tmp_SM_00002= SM(i);
   tmp_burgersfun_noloop_00175= i- 1;
   g_tmp_SM_00003= g_SM(tmp_burgersfun_noloop_00175);
   tmp_SM_00003= SM(tmp_burgersfun_noloop_00175);
   g_tmp_burgersfun_noloop_00176= g_tmp_SM_00002+ g_tmp_SM_00003;
   tmp_burgersfun_noloop_00176= tmp_SM_00002+ tmp_SM_00003;
   g_tmp_x_00023= g_x(i);
   tmp_x_00023= x(i);
   tmp_burgersfun_noloop_00177= i- 1;
   g_tmp_x_00024= g_x(tmp_burgersfun_noloop_00177);
   tmp_x_00024= x(tmp_burgersfun_noloop_00177);
   g_tmp_burgersfun_noloop_00178= g_tmp_x_00023- g_tmp_x_00024;
   tmp_burgersfun_noloop_00178= tmp_x_00023- tmp_x_00024;
   g_tmp_burgersfun_noloop_00179= g_tmp_burgersfun_noloop_00176.* tmp_burgersfun_noloop_00178+ tmp_burgersfun_noloop_00176.* g_tmp_burgersfun_noloop_00178;
   tmp_burgersfun_noloop_00179= tmp_burgersfun_noloop_00176.* tmp_burgersfun_noloop_00178;
   g_g(i+ N)= g_tmp_burgersfun_noloop_00174- g_tmp_burgersfun_noloop_00179;
   g(i+ N)= tmp_burgersfun_noloop_00174- tmp_burgersfun_noloop_00179; 
   %end
   g_tmp_SM_00004= g_SM(2);
   tmp_SM_00004= SM(2);
   g_tmp_SM_00005= g_SM(1);
   tmp_SM_00005= SM(1);
   g_tmp_burgersfun_noloop_00180= g_tmp_SM_00004+ g_tmp_SM_00005;
   tmp_burgersfun_noloop_00180= tmp_SM_00004+ tmp_SM_00005;
   g_tmp_x_00025= g_x(2);
   tmp_x_00025= x(2);
   g_tmp_x_00026= g_x(1);
   tmp_x_00026= x(1);
   g_tmp_burgersfun_noloop_00181= g_tmp_x_00025- g_tmp_x_00026;
   tmp_burgersfun_noloop_00181= tmp_x_00025- tmp_x_00026;
   g_tmp_burgersfun_noloop_00182= g_tmp_burgersfun_noloop_00180* tmp_burgersfun_noloop_00181+ tmp_burgersfun_noloop_00180* g_tmp_burgersfun_noloop_00181;
   tmp_burgersfun_noloop_00182= tmp_burgersfun_noloop_00180* tmp_burgersfun_noloop_00181;
   g_tmp_SM_00006= g_SM(1);
   tmp_SM_00006= SM(1);
   g_tmp_burgersfun_noloop_00183= g_tmp_SM_00006+ g_SM0;
   tmp_burgersfun_noloop_00183= tmp_SM_00006+ SM0;
   g_tmp_x_00027= g_x(1);
   tmp_x_00027= x(1);
   g_tmp_burgersfun_noloop_00184= g_tmp_x_00027+ g_zeros(size(x0));
   tmp_burgersfun_noloop_00184= tmp_x_00027- x0;
   g_tmp_burgersfun_noloop_00185= g_tmp_burgersfun_noloop_00183* tmp_burgersfun_noloop_00184+ tmp_burgersfun_noloop_00183* g_tmp_burgersfun_noloop_00184;
   tmp_burgersfun_noloop_00185= tmp_burgersfun_noloop_00183* tmp_burgersfun_noloop_00184;
   g_g(1+ N)= g_tmp_burgersfun_noloop_00182- g_tmp_burgersfun_noloop_00185;
   g(1+ N)= tmp_burgersfun_noloop_00182- tmp_burgersfun_noloop_00185; 
   g_tmp_SM_00007= g_SM(N);
   tmp_SM_00007= SM(N);
   g_tmp_burgersfun_noloop_00186= g_SMNP1+ g_tmp_SM_00007;
   tmp_burgersfun_noloop_00186= SMNP1+ tmp_SM_00007;
   g_tmp_x_00028= g_x(N);
   tmp_x_00028= x(N);
   g_tmp_burgersfun_noloop_00187= -g_tmp_x_00028+ g_zeros(size(xNP1));
   tmp_burgersfun_noloop_00187= xNP1- tmp_x_00028;
   g_tmp_burgersfun_noloop_00188= g_tmp_burgersfun_noloop_00186* tmp_burgersfun_noloop_00187+ tmp_burgersfun_noloop_00186* g_tmp_burgersfun_noloop_00187;
   tmp_burgersfun_noloop_00188= tmp_burgersfun_noloop_00186* tmp_burgersfun_noloop_00187;
   g_tmp_SM_00008= g_SM(N);
   tmp_SM_00008= SM(N);
   tmp_burgersfun_noloop_00189= N- 1;
   g_tmp_SM_00009= g_SM(tmp_burgersfun_noloop_00189);
   tmp_SM_00009= SM(tmp_burgersfun_noloop_00189);
   g_tmp_burgersfun_noloop_00190= g_tmp_SM_00008+ g_tmp_SM_00009;
   tmp_burgersfun_noloop_00190= tmp_SM_00008+ tmp_SM_00009;
   g_tmp_x_00029= g_x(N);
   tmp_x_00029= x(N);
   tmp_burgersfun_noloop_00191= N- 1;
   g_tmp_x_00030= g_x(tmp_burgersfun_noloop_00191);
   tmp_x_00030= x(tmp_burgersfun_noloop_00191);
   g_tmp_burgersfun_noloop_00192= g_tmp_x_00029- g_tmp_x_00030;
   tmp_burgersfun_noloop_00192= tmp_x_00029- tmp_x_00030;
   g_tmp_burgersfun_noloop_00193= g_tmp_burgersfun_noloop_00190* tmp_burgersfun_noloop_00192+ tmp_burgersfun_noloop_00190* g_tmp_burgersfun_noloop_00192;
   tmp_burgersfun_noloop_00193= tmp_burgersfun_noloop_00190* tmp_burgersfun_noloop_00192;
   g_g(N+ N)= g_tmp_burgersfun_noloop_00188- g_tmp_burgersfun_noloop_00193;
   g(N+ N)= tmp_burgersfun_noloop_00188- tmp_burgersfun_noloop_00193; 
   tau= 1e-3; 
   g_tmp_g_00000= g_g(1+ N: end);
   tmp_g_00000= g(1+ N: end);
   g_tmp_burgersfun_noloop_00194= -g_tmp_g_00000;
   tmp_burgersfun_noloop_00194= -tmp_g_00000;
   tmp_burgersfun_noloop_00195= 2* tau;
   g_g(1+ N: end)= adimat_g_mrdivide1(g_tmp_burgersfun_noloop_00194, tmp_burgersfun_noloop_00194, (tmp_burgersfun_noloop_00195));
   g(1+ N: end)= tmp_burgersfun_noloop_00194/ tmp_burgersfun_noloop_00195; 
   g_out= g_g;
   out= g; 
end