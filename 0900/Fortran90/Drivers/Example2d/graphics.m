%%% clear
aux = sprintf(' ') ; disp( aux )
aux = sprintf('cleaning ... ') ; disp( aux )
close all;
clear;

%%% load info
aux = sprintf(' ') ; disp( aux )
aux = sprintf('loading ... ') ; disp( aux )
load truth.txt truth
load model.txt model
load ekf.txt ekf
load rrsqrtkf.txt rrsqrtkf
load enkf.txt enkf
load rrsqrtenkf.txt rrsqrtenkf

%%% constants
aux = sprintf(' ') ; disp( aux )
aux = sprintf('constants ... ') ; disp( aux )
ni = 1;
nt = 1110;
nx = 21;
ny = 21;

%%% reconstructing t , x , y , truth , model , ekf , rrsqrtkf , enkf , rrsqrtenkf
aux = sprintf(' ') ; disp( aux )
aux = sprintf('building fields ... ') ; disp( aux )
counter = 1;
for l = 1 : nt
  for i = 1 : nx
    for j = 1 : ny
     t( l ) = truth( counter , 1 );
     x( i , j ) = truth( counter , 2 );
     y( i , j ) = truth( counter , 3 );
     truth1( l , i , j ) = truth( counter , 4 );
     model1( l , i , j ) = model( counter , 4 );
     ekf1( l , i , j ) = ekf( counter , 4 );
     rrsqrtkf1( l , i , j ) = rrsqrtkf( counter , 4 );
     enkf1( l , i , j ) = enkf( counter , 4 );
     rrsqrtenkf1( l , i , j ) = rrsqrtenkf( counter , 4 );
     counter = counter + 1;
    end
  end
end 
clear truth;
clear model;
clear ekf;
clear rrsqrtkf;
clear enkf;
clear rrsqrtenkf;

%%% printing on screen
aux = sprintf(' ') ; disp( aux )
aux = sprintf('absolute value of the percentage relative error') ; disp( aux )
for l = ni : nt
  for i = 2 : nx-1
    for j = 2 : ny-1
      xx( i , j ) = x( i , j );
      yy( i , j ) = y( i , j );
      truthl( i , j ) = truth1( l , i , j );
      modeltruthl( i , j ) = abs( ( model1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
      ekftruthl( i , j ) = abs( ( ekf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
      rrsqrtkftruthl( i , j ) = abs( ( rrsqrtkf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
      enkftruthl( i , j ) = abs( ( enkf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
      rrsqrtenkftruthl( i , j ) = abs( ( rrsqrtenkf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
    end
  end
  aux = sprintf('l = %d|model-truth=%f|ekf-truth=%f|rrsqrtkf-truth=%f|enkf-truth=%f|rrsqrtenkf-truth=%f' , l , max( max( modeltruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) , max( max( ekftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) , max( max( rrsqrtkftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) , max( max( enkftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) , max( max( rrsqrtenkftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) ) ; disp( aux )
end
clear xx;
clear yy;
clear truthl;
clear modeltruthl;
clear ekftruthl;
clear rrsqrtkftruthl;
clear enkftruthl;
clear rrsqrtenkftruthl;

%%% making movies 
aux = sprintf(' ') ; disp( aux )
aux = sprintf('making movies') ; disp( aux )

%%% model and truth
aux = sprintf(' ') ; disp( aux )
aux = sprintf('making movies model-truth') ; disp( aux )
fig = figure( 1 );
winsize = get( fig , 'Position' );
winsize( 1 : 2 ) = [ 0 0 ];
movie = moviein( nt - ni + 1 , fig , winsize );
set( fig , 'NextPlot' , 'replacechildren' );
for l = ni : nt
  for i = 2 : nx-1
    for j = 2 : ny-1
      xx( i , j ) = x( i , j );
      yy( i , j ) = y( i , j );
      truthl( i , j ) = truth1( l , i , j );
      modeltruthl( i , j ) = abs( ( model1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
    end
  end
  aux = sprintf('l = %d|model-truth=%f' , l , max( max( modeltruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) ) ; disp( aux )
  h = get( gcf , 'children' );
  set( h , 'fontsize' , 12 );
  xlabel('X-AXIS');
  ylabel('Y-AXIS');
  zlabel(' RELATIVE ERROR (%): |MODEL - TRUTH|');
  set( gca , 'NextPlot' , 'replacechildren' )
  axis( [ -1 1   -1 1   0 100 ] )
  surf( xx , yy , modeltruthl );
  grid on
  aux = sprintf( 'MODEL AND TRUTH -- TIME STEP L = %d' , l );
  title( aux );
  movie( : , l - ni + 1 ) = getframe( fig , winsize );
end
movie2avi( movie , 'movie_model_truth.avi' , 'fps' , 25 , 'quality' , 25 );
print -depsc example2d_model_truth.eps
close all;
clear xx;
clear yy;
clear truthl;
clear modeltruthl;
clear movie;

%%% ekf and truth
aux = sprintf(' ') ; disp( aux )
aux = sprintf('making movies ekf-truth') ; disp( aux )
fig = figure( 1 );
winsize = get( fig , 'Position' );
winsize( 1 : 2 ) = [ 0 0 ];
movie = moviein( nt - ni + 1 , fig , winsize );
set( fig , 'NextPlot' , 'replacechildren' );
for l = ni : nt
  for i = 2 : nx-1
    for j = 2 : ny-1
      xx( i , j ) = x( i , j );
      yy( i , j ) = y( i , j );
      truthl( i , j ) = truth1( l , i , j );
      ekftruthl( i , j ) = abs( ( ekf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
    end
  end
  aux = sprintf('l = %d|ekf-truth=%f' , l , max( max( ekftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) ) ; disp( aux )
  h = get( gcf , 'children' );
  set( h , 'fontsize' , 12 );
  xlabel('X-AXIS');
  ylabel('Y-AXIS');
  zlabel(' RELATIVE ERROR (%): |EKF - TRUTH|');
  set( gca , 'NextPlot' , 'replacechildren' )
  axis( [ -1 1   -1 1   0 100 ] )
  surf( xx , yy , ekftruthl );
  grid on
  aux = sprintf( 'EKF AND TRUTH -- TIME STEP L = %d' , l );
  title( aux );
  movie( : , l - ni + 1 ) = getframe( fig , winsize );
end
movie2avi( movie , 'movie_ekf_truth.avi' , 'fps' , 25 , 'quality' , 25 );
print -depsc example2d_ekf_truth.eps
close all;
clear xx;
clear yy;
clear truthl;
clear ekftruthl;
clear movie;

%%% rrsqrtkf and truth
aux = sprintf(' ') ; disp( aux )
aux = sprintf('making movies rrsqrtkf-truth') ; disp( aux )
fig = figure( 1 );
winsize = get( fig , 'Position' );
winsize( 1 : 2 ) = [ 0 0 ];
movie = moviein( nt - ni + 1 , fig , winsize );
set( fig , 'NextPlot' , 'replacechildren' );
for l = ni : nt
  for i = 2 : nx-1
    for j = 2 : ny-1
      xx( i , j ) = x( i , j );
      yy( i , j ) = y( i , j );
      truthl( i , j ) = truth1( l , i , j );
      rrsqrtkftruthl( i , j ) = abs( ( rrsqrtkf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
    end
  end
  aux = sprintf('l = %d|rrsqrtkf-truth=%f' , l , max( max( rrsqrtkftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) ) ; disp( aux )
  h = get( gcf , 'children' );
  set( h , 'fontsize' , 12 );
  xlabel('X-AXIS');
  ylabel('Y-AXIS');
  zlabel(' RELATIVE ERROR (%): |RRSQRTKF - TRUTH|');
  set( gca , 'NextPlot' , 'replacechildren' )
  axis( [ -1 1   -1 1   0 100 ] )
  surf( xx , yy , rrsqrtkftruthl );
  grid on
  aux = sprintf( 'RRSQRTKF AND TRUTH -- TIME STEP L = %d' , l );
  title( aux );
  movie( : , l - ni + 1 ) = getframe( fig , winsize );
end
movie2avi( movie , 'movie_rrsqrtkf_truth.avi' , 'fps' , 25 , 'quality' , 25 );
print -depsc example2d_rrsqrtkf_truth.eps
close all;
clear xx;
clear yy;
clear truthl;
clear rrsqrttruthl;
clear movie;

%%% enkf and truth
aux = sprintf(' ') ; disp( aux )
aux = sprintf('making movies enkf-truth') ; disp( aux )
fig = figure( 1 );
winsize = get( fig , 'Position' );
winsize( 1 : 2 ) = [ 0 0 ];
movie = moviein( nt - ni + 1 , fig , winsize );
set( fig , 'NextPlot' , 'replacechildren' );
for l = ni : nt
  for i = 2 : nx-1
    for j = 2 : ny-1
      xx( i , j ) = x( i , j );
      yy( i , j ) = y( i , j );
      truthl( i , j ) = truth1( l , i , j );
      enkftruthl( i , j ) = abs( ( enkf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
    end
  end
  aux = sprintf('l = %d|enkf-truth=%f' , l , max( max( enkftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) ) ; disp( aux )
  h = get( gcf , 'children' );
  set( h , 'fontsize' , 12 );
  xlabel('X-AXIS');
  ylabel('Y-AXIS');
  zlabel(' RELATIVE ERROR (%): |ENKF - TRUTH|');
  set( gca , 'NextPlot' , 'replacechildren' )
  axis( [ -1 1   -1 1   0 100 ] )
  surf( xx , yy , enkftruthl );
  grid on
  aux = sprintf( 'ENKF AND TRUTH -- TIME STEP L = %d' , l );
  title( aux );
  movie( : , l - ni + 1 ) = getframe( fig , winsize );
end
movie2avi( movie , 'movie_enkf_truth.avi' , 'fps' , 25 , 'quality' , 25 );
print -depsc example2d_enkf_truth.eps
close all;
clear xx;
clear yy;
clear truthl;
clear enkftruthl;
clear movie;

%%% rrsqrtenkf and truth
aux = sprintf(' ') ; disp( aux )
aux = sprintf('making movies rrsqrtenkf-truth') ; disp( aux )
fig = figure( 1 );
winsize = get( fig , 'Position' );
winsize( 1 : 2 ) = [ 0 0 ];
movie = moviein( nt - ni + 1 , fig , winsize );
set( fig , 'NextPlot' , 'replacechildren' );
for l = ni : nt
  for i = 2 : nx-1
    for j = 2 : ny-1
      xx( i , j ) = x( i , j );
      yy( i , j ) = y( i , j );
      truthl( i , j ) = truth1( l , i , j );
      rrsqrtenkftruthl( i , j ) = abs( ( rrsqrtenkf1( l , i , j ) - truth1( l , i , j ) ) / ( truth1( l , i , j ) ) * 100.0 );
    end
  end
  aux = sprintf('l = %d|rrsqrtenkf-truth=%f' , l , max( max( rrsqrtenkftruthl( 2 : nx - 1 , 2 : ny - 1 ) ) ) ) ; disp( aux )
  h = get( gcf , 'children' );
  set( h , 'fontsize' , 12 );
  xlabel('X-AXIS');
  ylabel('Y-AXIS');
  zlabel(' RELATIVE ERROR (%): |RRSQRTENKF - TRUTH|');
  set( gca , 'NextPlot' , 'replacechildren' )
  axis( [ -1 1   -1 1   0 100 ] )
  surf( xx , yy , rrsqrtenkftruthl );
  grid on
  aux = sprintf( 'RRSQRTENKF AND TRUTH -- TIME STEP L = %d' , l );
  title( aux );
  movie( : , l - ni + 1 ) = getframe( fig , winsize );
end
movie2avi( movie , 'movie_rrsqrtenkf_truth.avi' , 'fps' , 25 , 'quality' , 25 );
print -depsc example2d_rrsqrtenkf_truth.eps
close all;
clear xx;
clear yy;
clear truthl;
clear rrsqrtenkftruthl;
clear movie;

