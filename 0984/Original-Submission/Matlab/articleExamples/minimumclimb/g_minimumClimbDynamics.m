% Generated by ADiMat 0.6.0-4975
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
% AD_DVARS= ydot

function [g_ydot, ydot]= g_minimumClimbDynamics(g_x, x, auxdata)
   % Dynamics file for minimum time to climb of supersonic aircraft
   
   CoF= auxdata.CoF; 
   
   % Extract state and control
   g_tmp_x_00000= g_x(: , 1);
   tmp_x_00000= x(: , 1);
   g_h= g_tmp_x_00000;
   h= tmp_x_00000; 
   g_tmp_x_00001= g_x(: , 2);
   tmp_x_00001= x(: , 2);
   g_v= g_tmp_x_00001;
   v= tmp_x_00001; 
   g_tmp_x_00002= g_x(: , 3);
   tmp_x_00002= x(: , 3);
   g_fpa= g_tmp_x_00002;
   fpa= tmp_x_00002; 
   g_tmp_x_00003= g_x(: , 4);
   tmp_x_00003= x(: , 4);
   g_u= g_tmp_x_00003;
   u= tmp_x_00003; 
   
   % Constants for scaled problem
   c1= 392.4; c2= 16818; c3= 86.138; 
   c4= 288.14; c5= 6.49; c6= 4.0519e9; 
   c7= 288.08; c8= 5.256; c9= 216.64; 
   c10= 9.06e8; c11= 1.73; c12= 0.157; 
   c13= 6e-5; c14= 4.00936; c15= 2.2; 
   
   % NOAA atmospheric model for pressure, temperature, density
   g_Nzeros= g_h.* 0;
   Nzeros= h.* 0; 
   g_T= g_Nzeros;
   T= Nzeros; 
   g_p= g_Nzeros;
   p= Nzeros; 
   
   ihlow= h< 11000; 
   g_tmp_h_00000= g_h(ihlow);
   tmp_h_00000= h(ihlow);
   g_tmp_minimumClimbDynamics_00000= c5* g_tmp_h_00000;
   tmp_minimumClimbDynamics_00000= c5* tmp_h_00000;
   g_T(ihlow)= -g_tmp_minimumClimbDynamics_00000+ g_zeros(size(c4));
   T(ihlow)= c4- tmp_minimumClimbDynamics_00000; 
   g_tmp_T_00000= g_T(ihlow);
   tmp_T_00000= T(ihlow);
   g_tmp_minimumClimbDynamics_00001= g_tmp_T_00000./ c7;
   tmp_minimumClimbDynamics_00001= tmp_T_00000./ c7;
   tmp_minimumClimbDynamics_00035= c8.* tmp_minimumClimbDynamics_00001.^ (c8- 1);
   tmp_minimumClimbDynamics_00035(tmp_minimumClimbDynamics_00001== 0& c8== 0)= 0; % Ensure, that derivative of 0.^0 is 0 and not NaN.
   g_tmp_minimumClimbDynamics_00002= tmp_minimumClimbDynamics_00035.* g_tmp_minimumClimbDynamics_00001;
   tmp_minimumClimbDynamics_00002= tmp_minimumClimbDynamics_00001.^ c8;
   g_p(ihlow)= c6* g_tmp_minimumClimbDynamics_00002;
   p(ihlow)= c6* tmp_minimumClimbDynamics_00002; 
   
   ihhigh= ~ihlow; 
   T(ihhigh)= c9; 
   g_T(ihhigh)= g_zeros(size(T(ihhigh)));
   g_tmp_h_00001= g_h(ihhigh);
   tmp_h_00001= h(ihhigh);
   g_tmp_minimumClimbDynamics_00003= c12* g_tmp_h_00001;
   tmp_minimumClimbDynamics_00003= c12* tmp_h_00001;
   g_tmp_minimumClimbDynamics_00004= -g_tmp_minimumClimbDynamics_00003+ g_zeros(size(c11));
   tmp_minimumClimbDynamics_00004= c11- tmp_minimumClimbDynamics_00003;
   g_tmp_exp_00000= g_tmp_minimumClimbDynamics_00004.* exp(tmp_minimumClimbDynamics_00004);
   tmp_exp_00000= exp(tmp_minimumClimbDynamics_00004);
   g_p(ihhigh)= c10* g_tmp_exp_00000;
   p(ihhigh)= c10* tmp_exp_00000; 
   
   g_tmp_minimumClimbDynamics_00005= c3* g_p;
   tmp_minimumClimbDynamics_00005= c3* p;
   g_rho= (g_tmp_minimumClimbDynamics_00005.* T- tmp_minimumClimbDynamics_00005.* g_T)./ T.^ 2;
   rho= tmp_minimumClimbDynamics_00005./ T; 
   g_q= 0.5.* g_rho.* v.* v.* c13+ 0.5.* rho.* g_v.* v.* c13+ 0.5.* rho.* v.* g_v.* c13;
   q= 0.5.* rho.* v.* v.* c13; 
   
   % Collect mach number powers
   tmp_sqrt_00000= sqrt(T);
   g_tmp_sqrt_00000= g_T./ (2.* tmp_sqrt_00000);
   g_a= c14.* g_tmp_sqrt_00000;
   a= c14.* tmp_sqrt_00000; g_M= (g_v.* a- v.* g_a)./ a.^ 2;
   M= v./ a; 
   Mp= cell(1, 6); 
   g_Mp= cell(size(Mp));
   for i= 1: 6
      tmp_minimumClimbDynamics_00006= i- 1;
      tmp_minimumClimbDynamics_00036= tmp_minimumClimbDynamics_00006.* M.^ (tmp_minimumClimbDynamics_00006- 1);
      tmp_minimumClimbDynamics_00036(M== 0& tmp_minimumClimbDynamics_00006== 0)= 0; % Ensure, that derivative of 0.^0 is 0 and not NaN.
      g_Mp{i}= tmp_minimumClimbDynamics_00036.* g_M;
      Mp{i}= M.^ tmp_minimumClimbDynamics_00006; 
   end
   
   % Use table data for drag coefficient (from Darby formulation)
   g_numeratorCD0= g_Nzeros;
   numeratorCD0= Nzeros; g_denominatorCD0= g_Nzeros;
   denominatorCD0= Nzeros; 
   g_numeratorK= g_Nzeros;
   numeratorK= Nzeros; g_denominatorK= g_Nzeros;
   denominatorK= Nzeros; 
   for i= 1: 6
      g_tmp_Mp_00000= g_Mp{i};
      tmp_Mp_00000= Mp{i};
      g_Mpi= g_tmp_Mp_00000;
      Mpi= tmp_Mp_00000; 
      if i< 6
         g_tmp_minimumClimbDynamics_00007= CoF(1, i).* g_Mpi;
         tmp_minimumClimbDynamics_00007= CoF(1, i).* Mpi;
         g_tmp_minimumClimbDynamics_00029= g_numeratorCD0+ g_tmp_minimumClimbDynamics_00007;
         tmp_minimumClimbDynamics_00029= numeratorCD0+ tmp_minimumClimbDynamics_00007; 
         % Update detected: numeratorCD0= some_expression(numeratorCD0,...)
         g_numeratorCD0= g_tmp_minimumClimbDynamics_00029;
         numeratorCD0= tmp_minimumClimbDynamics_00029;
         g_tmp_minimumClimbDynamics_00008= CoF(2, i).* g_Mpi;
         tmp_minimumClimbDynamics_00008= CoF(2, i).* Mpi;
         g_tmp_minimumClimbDynamics_00030= g_denominatorCD0+ g_tmp_minimumClimbDynamics_00008;
         tmp_minimumClimbDynamics_00030= denominatorCD0+ tmp_minimumClimbDynamics_00008; 
         % Update detected: denominatorCD0= some_expression(denominatorCD0,...)
         g_denominatorCD0= g_tmp_minimumClimbDynamics_00030;
         denominatorCD0= tmp_minimumClimbDynamics_00030;
         g_tmp_minimumClimbDynamics_00009= CoF(3, i).* g_Mpi;
         tmp_minimumClimbDynamics_00009= CoF(3, i).* Mpi;
         g_tmp_minimumClimbDynamics_00031= g_numeratorK+ g_tmp_minimumClimbDynamics_00009;
         tmp_minimumClimbDynamics_00031= numeratorK+ tmp_minimumClimbDynamics_00009; 
         % Update detected: numeratorK= some_expression(numeratorK,...)
         g_numeratorK= g_tmp_minimumClimbDynamics_00031;
         numeratorK= tmp_minimumClimbDynamics_00031;
      end
      g_tmp_minimumClimbDynamics_00010= CoF(4, i).* g_Mpi;
      tmp_minimumClimbDynamics_00010= CoF(4, i).* Mpi;
      g_tmp_minimumClimbDynamics_00032= g_denominatorK+ g_tmp_minimumClimbDynamics_00010;
      tmp_minimumClimbDynamics_00032= denominatorK+ tmp_minimumClimbDynamics_00010; 
      % Update detected: denominatorK= some_expression(denominatorK,...)
      g_denominatorK= g_tmp_minimumClimbDynamics_00032;
      denominatorK= tmp_minimumClimbDynamics_00032;
   end
   g_Cd0= (g_numeratorCD0.* denominatorCD0- numeratorCD0.* g_denominatorCD0)./ denominatorCD0.^ 2;
   Cd0= numeratorCD0./ denominatorCD0; 
   g_K= (g_numeratorK.* denominatorK- numeratorK.* g_denominatorK)./ denominatorK.^ 2;
   K= numeratorK./ denominatorK; 
   tmp_minimumClimbDynamics_00011= c2^ 2;
   tmp_minimumClimbDynamics_00012= c1^ 2;
   tmp_minimumClimbDynamics_00013= tmp_minimumClimbDynamics_00011.* tmp_minimumClimbDynamics_00012;
   g_tmp_minimumClimbDynamics_00014= 2.* q.^ (2- 1).* g_q;
   tmp_minimumClimbDynamics_00014= q.^ 2;
   g_tmp_minimumClimbDynamics_00015= (-tmp_minimumClimbDynamics_00013.* g_tmp_minimumClimbDynamics_00014)./ tmp_minimumClimbDynamics_00014.^ 2;
   tmp_minimumClimbDynamics_00015= tmp_minimumClimbDynamics_00013./ tmp_minimumClimbDynamics_00014;
   g_tmp_minimumClimbDynamics_00016= 2.* u.^ (2- 1).* g_u;
   tmp_minimumClimbDynamics_00016= u.^ 2;
   g_tmp_minimumClimbDynamics_00017= g_K.* tmp_minimumClimbDynamics_00015.* tmp_minimumClimbDynamics_00016+ K.* g_tmp_minimumClimbDynamics_00015.* tmp_minimumClimbDynamics_00016+ K.* tmp_minimumClimbDynamics_00015.* g_tmp_minimumClimbDynamics_00016;
   tmp_minimumClimbDynamics_00017= K.* tmp_minimumClimbDynamics_00015.* tmp_minimumClimbDynamics_00016;
   g_tmp_minimumClimbDynamics_00018= g_Cd0+ g_tmp_minimumClimbDynamics_00017;
   tmp_minimumClimbDynamics_00018= Cd0+ tmp_minimumClimbDynamics_00017;
   g_FD= g_q.* tmp_minimumClimbDynamics_00018+ q.* g_tmp_minimumClimbDynamics_00018;
   FD= q.* tmp_minimumClimbDynamics_00018; 
   
   % Use table data for thrust coefficient (from Darby formulation)
   g_FT= g_Nzeros;
   FT= Nzeros; 
   for i= 1: 6
      g_ei= g_Nzeros;
      ei= Nzeros; 
      for j= 1: 6
         g_tmp_Mp_00001= g_Mp{j};
         tmp_Mp_00001= Mp{j};
         g_tmp_minimumClimbDynamics_00019= CoF(4+ j, i).* g_tmp_Mp_00001;
         tmp_minimumClimbDynamics_00019= CoF(4+ j, i).* tmp_Mp_00001;
         g_tmp_minimumClimbDynamics_00033= g_ei+ g_tmp_minimumClimbDynamics_00019;
         tmp_minimumClimbDynamics_00033= ei+ tmp_minimumClimbDynamics_00019; 
         % Update detected: ei= some_expression(ei,...)
         g_ei= g_tmp_minimumClimbDynamics_00033;
         ei= tmp_minimumClimbDynamics_00033;
      end
      tmp_minimumClimbDynamics_00020= i- 1;
      tmp_minimumClimbDynamics_00037= tmp_minimumClimbDynamics_00020.* h.^ (tmp_minimumClimbDynamics_00020- 1);
      tmp_minimumClimbDynamics_00037(h== 0& tmp_minimumClimbDynamics_00020== 0)= 0; % Ensure, that derivative of 0.^0 is 0 and not NaN.
      g_tmp_minimumClimbDynamics_00021= tmp_minimumClimbDynamics_00037.* g_h;
      tmp_minimumClimbDynamics_00021= h.^ tmp_minimumClimbDynamics_00020;
      g_tmp_minimumClimbDynamics_00022= g_ei.* tmp_minimumClimbDynamics_00021+ ei.* g_tmp_minimumClimbDynamics_00021;
      tmp_minimumClimbDynamics_00022= ei.* tmp_minimumClimbDynamics_00021;
      g_tmp_minimumClimbDynamics_00034= g_FT+ g_tmp_minimumClimbDynamics_00022;
      tmp_minimumClimbDynamics_00034= FT+ tmp_minimumClimbDynamics_00022; 
      % Update detected: FT= some_expression(FT,...)
      g_FT= g_tmp_minimumClimbDynamics_00034;
      FT= tmp_minimumClimbDynamics_00034;
   end
   g_tmp_minimumClimbDynamics_00023= g_FT.* c1;
   tmp_minimumClimbDynamics_00023= FT.* c1;
   g_FT= adimat_g_mrdivide1(g_tmp_minimumClimbDynamics_00023, tmp_minimumClimbDynamics_00023, c15);
   FT= tmp_minimumClimbDynamics_00023/ c15; 
   
   % State dynamics
   g_tmp_sin_00000= g_fpa.* cos(fpa);
   tmp_sin_00000= sin(fpa);
   g_hdot= g_v.* tmp_sin_00000+ v.* g_tmp_sin_00000;
   hdot= v.* tmp_sin_00000; 
   g_tmp_minimumClimbDynamics_00024= g_FT- g_FD;
   tmp_minimumClimbDynamics_00024= FT- FD;
   g_tmp_minimumClimbDynamics_00025= g_tmp_minimumClimbDynamics_00024./ c2;
   tmp_minimumClimbDynamics_00025= tmp_minimumClimbDynamics_00024./ c2;
   g_tmp_sin_00001= g_fpa.* cos(fpa);
   tmp_sin_00001= sin(fpa);
   g_tmp_minimumClimbDynamics_00026= c1.* g_tmp_sin_00001;
   tmp_minimumClimbDynamics_00026= c1.* tmp_sin_00001;
   g_vdot= g_tmp_minimumClimbDynamics_00025- g_tmp_minimumClimbDynamics_00026;
   vdot= tmp_minimumClimbDynamics_00025- tmp_minimumClimbDynamics_00026; 
   g_tmp_cos_00000= g_fpa.* (-sin(fpa));
   tmp_cos_00000= cos(fpa);
   g_tmp_minimumClimbDynamics_00027= g_u- g_tmp_cos_00000;
   tmp_minimumClimbDynamics_00027= u- tmp_cos_00000;
   g_tmp_minimumClimbDynamics_00028= c1.* g_tmp_minimumClimbDynamics_00027;
   tmp_minimumClimbDynamics_00028= c1.* tmp_minimumClimbDynamics_00027;
   g_fpadot= (g_tmp_minimumClimbDynamics_00028.* v- tmp_minimumClimbDynamics_00028.* g_v)./ v.^ 2;
   fpadot= tmp_minimumClimbDynamics_00028./ v; 
   
   g_ydot= [g_hdot, g_vdot, g_fpadot];
   ydot= [hdot, vdot, fpadot]; 
end
