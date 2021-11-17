% test tmiscR2017b

%#ok<*NASGU>
%#ok<*CTPCT>
%#ok<*ASGLU>

INIT

syms zz;
syms xx; 
assume( xx, 'real' );

% preconditions

%% test liminf
EXPECT_EQ( liminf(sym([10 1 9 2 8 3 7 4 6])), sym([1 1 2 2 3 3 4 4 6]) );


%% test limsup
EXPECT_EQ( limsup(sym([10 1 9 2 8 3 7 4 6])), sym([10 9 9 8 8 7 7 6 6]) );