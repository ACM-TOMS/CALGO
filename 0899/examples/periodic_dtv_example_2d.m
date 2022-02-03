
% periodic_dtv_example_2d

  clear, home, close all
  
  N = 64;
  M = 98;
  
   h = @(X,Y) 0.*((X.^2 + Y.^2)>0.25) + 1.*((X.^2 + Y.^2)<=0.25);
  
        x = -1 + 2*(0:N-1)/N;
       xp = -1 + 2*(0:M-1)/M;  
    
     
       [X,Y] = meshgrid(x,x);
       [XP,YP] = meshgrid(xp,xp); 
       
         f = h(X,Y);
        fp = h(XP,YP);
   
      clear X Y XP YP 
       
      [uc,ak] = fourierInterpolation2d(f,M,M);
   
       clear f
       

       mesh(xp,xp,uc);
       
       lambda = 10;
       iterations = 100;
       
       tic
%      up = digitalTotalVariationFilterPeriodic_2d(uc,lambda,iterations);
       up = digitalTotalVariationFilterPeriodic_2d_8(uc,lambda,iterations);
       toc
       
       figure
       mesh(xp,xp,up);
