%%************************************************************************
%%************************************************************************
   function [fname,fd] = problem3
   
%%
%% theta: random
%% idx1 = [1:14];   
%%
   fname{1,1} = 'theta4';      fname{1,2}=-49; fname{1,3}=-58; 
   fname{2,1} = 'theta42';     fname{2,2}=-23; fname{2,3}=-26; 
   fname{3,1} = 'theta6';      fname{3,2}=-62; fname{3,3}=-66; 
   fname{4,1} = 'theta62';     fname{4,2}=-29; fname{4,3}=-32; 
   fname{5,1} = 'theta8';      fname{5,2}=-73; fname{5,3}=-76; 
   fname{6,1} = 'theta82';    
   fname{7,1} = 'theta83';    
   fname{8,1} = 'theta10';    
   fname{9,1} = 'theta102';   
   fname{10,1} = 'theta103';  
   fname{11,1} = 'theta104';  
   fname{12,1} = 'theta12';   
   fname{13,1} = 'theta123';     
   fname{14,1} = 'theta162'; 
   fd(1:20) = ones(1,20);  
%%
%% theta: Dimacs 
%% idx2 = [21:38, 41:45, 51:54, 61:72];
%%
   fname{21,1} = 'MANN-a27';      
   fname{22,1} = 'johnson8-4-4';  
   fname{23,1} = 'johnson16-2-4'; 
   fname{24,1} = 'san200-0.7-1';  
   fname{25,1} = 'sanr200-0.7'; 
   fname{26,1} = 'c-fat200-1';    
   fname{27,1} = 'hamming-6-4';   
   fname{28,1} = 'hamming-8-4';   
   fname{29,1} = 'hamming-9-8';   
   fname{30,1} = 'hamming-10-2';  
   fname{31,1} = 'hamming-7-5-6'; 
   fname{32,1} = 'hamming-8-3-4'; 
   fname{33,1} = 'hamming-9-5-6'; 
   fname{34,1} = 'brock200-1';  
   fname{35,1} = 'brock200-4';  
   fname{36,1} = 'brock400-1';  
   fname{37,1} = 'keller4';     
   fname{38,1} = 'p-hat300-1';  

   fname{41,1} = 'G43';  
   fname{42,1} = 'G44';  
   fname{43,1} = 'G45';  
   fname{44,1} = 'G46';  
   fname{45,1} = 'G47';  
   fname{51,1} = []; 'G51'; 
   fname{52,1} = []; 'G52'; 
   fname{53,1} = []; 'G53'; 
   fname{54,1} = []; 'G54'; 

   fd(21:60) = 1.2*ones(1,40); 
%%
%% theta: Borchers
%% idx3 = [101:128];
%%
   fname{101,1} = '1dc.64';      fname{101,2} = 10; fname{101,3} = -13; 
   fname{102,1} = '1et.64';      fname{102,2} = 18; fname{102,3} = -20;   
   fname{103,1} = '1tc.64';      fname{103,2} = 20; fname{103,3} = -23;   
   fname{104,1} = '1dc.128';     fname{104,2} = 16; fname{104,3} = -19;   
   fname{105,1} = '1et.128';     fname{105,2} = 28; fname{105,3} = -31;   
   fname{106,1} = '1tc.128';     fname{106,2} = 38; fname{106,3} = -42;    
   fname{107,1} = '1zc.128';     fname{107,2} = 18; fname{107,3} = -21;   
   fname{108,1} = []; '2dc.128';     fname{108,2} = 5;    
   fname{109,1} = '1dc.256';     fname{109,2} = 30; fname{109,3} = -33;   
   fname{110,1} = '1et.256';     fname{110,2} = 50; fname{110,3} = -56;    
   fname{111,1} = '1tc.256';     fname{111,2} = 63; fname{111,3} = -65;   
   fname{112,1} = '1zc.256';     fname{112,2} = 36; fname{112,3} = -40;     
   fname{113,1} = []; '2dc.256';     fname{113,2} = 7;    
   fname{114,1} = '1dc.512';     fname{114,2} = 52;  fname{114,3} = -54;   
   fname{115,1} = '1et.512';     fname{115,2} = 100; fname{115,3} = -106;   
   fname{116,1} = '1tc.512';     fname{116,2} = 110; fname{116,3} = -116;  
   fname{117,1} = '1zc.512';     fname{117,2} = 62;  fname{117,3} = -70;   
   fname{118,1} = []; '2dc.512';     fname{118,2} = 11;  fname{118,3} = -14;     
   fname{119,1} = '1dc.1024';    fname{119,2} = 94;  fname{119,3} = -98;  
   fname{120,1} = '1et.1024';    fname{120,2} = 171; fname{120,3} = -190;    
   fname{121,1} = '1tc.1024';    fname{121,2} = 196; fname{121,3} = -208;   
   fname{122,1} = '1zc.1024';    fname{122,2} = 112; fname{122,3} = -132; fname{122,4} = -1;     
   fname{123,1} = []; '2dc.1024';    fname{123,2} = 16; fname{123,4} = -1;    
   fname{124,1} = '1dc.2048';
   fname{125,1} = '1et.2048';
   fname{126,1} = '1tc.2048';
   fname{127,1} = '1zc.2048';
   fname{128,1} = '2dc.2048';
   fname{129,1} = []; '1zc.4096';
   fd(101:130) = 1.3*ones(1,30); 
%%
%% QAP
%%
   fname{201,1} = 'bur26a'; fname{201,2} = 5426670; fname{201,3} = 5.2e6;
   fname{202,1} = 'bur26b'; fname{202,2} = 3817852; fname{202,3} = 3.5e6; 
   fname{203,1} = 'bur26c'; fname{203,2} = 5426795; fname{203,3} = 5.2e6;
   fname{204,1} = 'bur26d'; fname{204,2} = 3821225; fname{204,3} = 3.5e6;
   fname{205,1} = 'bur26e'; fname{205,2} = 5386879; fname{205,3} = 5.2e6;
   fname{206,1} = 'bur26f'; fname{206,2} = 3782044; fname{206,3} = 3.5e6;
   fname{207,1} = 'bur26g'; fname{207,2} = 10117172;fname{207,3} = 9.5e6;
   fname{208,1} = 'bur26h'; fname{208,2} = 7098658; fname{208,3} = 6.5e6;
   fname{209,1} = 'chr12a'; fname{209,2} = 9552;  fname{209,3} = 9.2e3;
   fname{210,1} = 'chr12b'; fname{210,2} = 9742;  fname{210,3} = 9.2e3;
   fname{211,1} = 'chr12c'; fname{211,2} = 11156; fname{211,3} = 1.05e4;
   fname{212,1} = 'chr15a'; fname{212,2} = 9896;  fname{212,3} = 9.0e3;
   fname{213,1} = 'chr15b'; fname{213,2} = 7990;  fname{213,3} = 7.5e3;
   fname{214,1} = 'chr15c'; fname{214,2} = 9504;  fname{214,3} = 9.2e3;
   fname{215,1} = 'chr18a'; fname{215,2} = 11098; fname{215,3} = 1.05e4;
   fname{216,1} = 'chr18b'; fname{216,2} = 1534;  fname{216,3} = 1.4e3;
   fname{217,1} = 'chr20a'; fname{217,2} = 2192;  fname{217,3} = 2.1e3;
   fname{218,1} = 'chr20b'; fname{218,2} = 2298;  fname{218,3} = 2.2e3;
   fname{219,1} = 'chr20c'; fname{219,2} = 14142; fname{219,3} = 1.35e4;
   fname{220,1} = 'chr22a'; fname{220,2} = 6156;  fname{220,3} = 6.0e3;
   fname{221,1} = 'chr22b'; fname{221,2} = 6194;  fname{221,3} = 6.1e3;
   fname{222,1} = 'chr25a'; fname{222,2} = 3796;  fname{222,3} = 3.6e3;
   fname{223,1} = 'els19';  fname{223,2} = 17212548; fname{223,3} = 1.6e7; 
   fname{224,1} = 'esc16a'; fname{224,2} = 68;
   fname{225,1} = 'esc16b'; fname{225,2} = 292;
   fname{226,1} = 'esc16c'; fname{226,2} = 160;
   fname{227,1} = 'esc16d'; fname{227,2} = 16;  fname{227,3} = 12;
   fname{228,1} = 'esc16e'; fname{228,2} = 28;
   fname{229,1} = []; % 'esc16f'; fname{229,2} = 0;
   fname{230,1} = 'esc16g'; fname{230,2} = 26;
   fname{231,1} = 'esc16h'; fname{231,2} = 996;
   fname{232,1} = 'esc16i'; fname{232,2} = 14;
   fname{233,1} = 'esc16j'; fname{233,2} = 8;
   fname{234,1} = 'esc32a'; fname{234,2} = 130; fname{234,4} = -1; 
   fname{235,1} = 'esc32b'; fname{235,2} = 168; fname{235,4} = -1; 
   fname{236,1} = 'esc32c'; fname{236,2} = 642; fname{236,4} = -1; 
   fname{237,1} = 'esc32d'; fname{237,2} = 200; fname{237,4} = -1; 
   fname{238,1} = 'esc32e'; fname{238,2} = 2;
   fname{239,1} = 'esc128'; fname{239,2} = 64;
   fname{240,1} = 'esc32g'; fname{240,2} = 6;
   fname{241,1} = 'esc32h'; fname{241,2} = 438; fname{241,4} = -1; 
   fname{242,1} = 'had12';  fname{242,2} = 1652;
   fname{243,1} = 'had14';  fname{243,2} = 2724;
   fname{244,1} = 'had16';  fname{244,2} = 3720;
   fname{245,1} = 'had18';  fname{245,2} = 5358;
   fname{246,1} = 'had20';  fname{246,2} = 6922; fname{246,3} = 6850;
   fname{247,1} = 'kra30a'; fname{247,2} = 88900; 
   fname{248,1} = 'kra30b'; fname{248,2} = 91420;
   fname{249,1} = 'kra32';  fname{249,2} = 88900;
   fname{250,1} = 'lipa20a'; fname{250,2} = 3683; fname{250,3} = 3500;
   fname{251,1} = 'lipa20b'; fname{251,2} = 27076;
   fname{252,1} = 'lipa30a'; fname{252,2} = 13178;
   fname{253,1} = 'lipa30b'; fname{253,2} = 151426;
   fname{254,1} = 'lipa40a'; fname{254,2} = 31538;
   fname{255,1} = 'lipa40b'; fname{255,2} = 476581;
   fname{256,1} = 'nug12';  fname{256,2} = 578;  fname{256,3} = 550;
   fname{257,1} = 'nug14';  fname{257,2} = 1014; fname{257,3} = 950;
   fname{258,1} = 'nug15';  fname{258,2} = 1150; fname{258,3} = 1100;
   fname{259,1} = 'nug16a'; fname{259,2} = 1610; fname{259,3} = 1550;
   fname{260,1} = 'nug16b'; fname{260,2} = 1240; fname{260,3} = 1200;
   fname{261,1} = 'nug17';  fname{261,2} = 1732; fname{261,3} = 1700;
   fname{262,1} = 'nug18';  fname{262,2} = 1930; fname{262,3} = 1900;
   fname{263,1} = 'nug20';  fname{263,2} = 2570; fname{263,3} = 2500;
   fname{264,1} = 'nug21';  fname{264,2} = 2438; fname{264,3} = 2350;
   fname{265,1} = 'nug22';  fname{265,2} = 3596; fname{265,3} = 3500;
   fname{266,1} = 'nug24';  fname{266,2} = 3488; fname{266,3} = 3400;
   fname{267,1} = 'nug25';  fname{267,2} = 3744; fname{267,3} = 3600;
   fname{268,1} = 'nug27';  fname{268,2} = 5234; fname{268,3} = 5100;
   fname{269,1} = 'nug28';  fname{269,2} = 5166; fname{269,3} = 4950;
   fname{270,1} = 'nug30';  fname{270,2} = 6124; fname{270,3} = 5900;
   fname{271,1} = 'rou12';  fname{271,2} = 235528;
   fname{272,1} = 'rou15';  fname{272,2} = 354210;
   fname{273,1} = 'rou20';  fname{273,2} = 725522;
   fname{274,1} = 'scr12';  fname{274,2} = 31410;
   fname{275,1} = 'scr15';  fname{275,2} = 51140;
   fname{276,1} = 'scr20';  fname{276,2} = 110030;
   fname{277,1} = 'ste36a'; fname{277,2} = 9526;
   fname{278,1} = 'ste36b'; fname{278,2} = 15852;
   fname{279,1} = 'ste36c'; fname{279,2} = 8239110; 
   fname{280,1} = 'tai12a'; fname{280,2} = 224416; 
   fname{281,1} = 'tai12b'; fname{281,2} = 39464925;
   fname{282,1} = 'tai15a'; fname{282,2} = 388214;
   fname{283,1} = 'tai15b'; fname{283,2} = 51765268;
   fname{284,1} = 'tai17a'; fname{284,2} = 491812;
   fname{285,1} = 'tai20a'; fname{285,2} = 703482;
   fname{286,1} = 'tai20b'; fname{286,2} = 122455319;
   fname{287,1} = 'tai25a'; fname{287,2} = 1167256;
   fname{288,1} = 'tai25b'; fname{288,2} = 344355646;
   fname{289,1} = 'tai30a'; fname{289,2} = 1818146;  fname{289,3} = 1.7e6; fname{289,4} = -1;  
   fname{290,1} = 'tai30b'; fname{290,2} = 637117113;fname{290,3} = 5.9e8; 
   fname{291,1} = 'tai35a'; fname{291,2} = 2422002;  fname{291,3} = 2.1e6; fname{291,4} = -1;  
   fname{292,1} = 'tai35b'; fname{292,2} = 283315445;fname{292,3} = 2.6e8;
   fname{293,1} = 'tai40a'; fname{293,2} = 3139370;  fname{293,3} = 2.8e6; fname{293,4} = -1;
   fname{294,1} = 'tai40b'; fname{294,2} = 637250948;fname{294,3} = 6.0e8;
   fname{295,1} = 'tho30';  fname{295,2} = 149936; fname{295,3} = 1.3e5;
   fname{296,1} = 'tho40';  fname{296,2} = 240516; fname{296,3} = 2.2e5; fname{296,4} = -1;
    %%% larger QAPs
   fname{297,1} = 'lipa50a'; fname{297,2} = 62093; 
   fname{298,1} = 'lipa50b'; fname{298,2} = 1210244; 
   fname{299,1} = 'lipa60a'; fname{299,2} = 107218; 
   fname{300,1} = 'lipa60b'; fname{300,2} = 2520135;
   fname{301,1} = 'lipa70a'; fname{301,2} = 169755; 
   fname{302,1} = 'lipa70b'; fname{302,2} = 4603200;   
   fname{303,1} = 'lipa80a'; fname{303,2} = 253195; 
   fname{304,1} = 'lipa80b'; fname{304,2} = 7763962;
   fname{305,1} = 'lipa90a'; fname{305,2} = 360630; 
   fname{306,1} = 'lipa90b'; fname{306,2} = 12490441;   
   fname{307,1} = 'sko42'; fname{307,2} = 15812 ;
   fname{308,1} = 'sko49'; fname{308,2} = 23386;
   fname{309,1} = 'sko56'; fname{309,2} = 34458;
   fname{310,1} = 'sko64'; fname{310,2} = 48498;
   fname{311,1} = 'sko72'; fname{311,2} = 66256;
   fname{312,1} = 'sko81'; fname{312,2} = 90998;
   fname{313,1} = 'sko90'; fname{313,2} = 115534;
   fname{314,1} = 'sko100a'; fname{314,2} = 152002;
   fname{315,1} = 'sko100b'; fname{315,2} = 153890;
   fname{316,1} = 'sko100c'; fname{316,2} = 147862;
   fname{317,1} = 'sko100d'; fname{317,2} = 149576;
   fname{318,1} = 'sko100e'; fname{318,2} = 149150;
   fname{319,1} = 'sko100f'; fname{319,2} = 149036;   
   fname{320,1} = 'tai50a'; fname{320,2} = 4938796;
   fname{321,1} = 'tai50b'; fname{321,2} = 458821517;
   fname{322,1} = 'tai60a'; fname{322,2} = 7205962;
   fname{323,1} = 'tai60b'; fname{323,2} = 608215054;
   fname{324,1} = 'tai64c'; fname{324,2} = 1855928;
   fname{325,1} = 'tai80a'; fname{325,2} = 13499184;
   fname{326,1} = 'tai80b'; fname{326,2} = 818415043;
   fname{327,1} = 'tai100a'; fname{327,2} = 21052466;
   fname{328,1} = 'tai100b'; fname{328,2} = 1185996137;
   fname{329,1} = 'tai150b'; fname{329,2} = 498896643;
   fname{330,1} = 'tai256c'; fname{330,2} = 44759294;
   fname{331,1} = 'tho150'; fname{331,2} = 8133398;
   fname{332,1} = 'wil50'; fname{332,2} = 48816;
   fname{333,1} = 'wil100'; fname{333,2} = 273038;
   
   fd(201:340) = 2*ones(1,140);    
%%
%% binary quadratic programming
%%
    fname{601,1} = 'be100.1';   fname{601,2} = -19412;
    fname{602,1} = 'be100.2';  fname{602,2} = -17290;
    fname{603,1} = 'be100.3';  fname{603,2} = -17565;   
    fname{604,1} = 'be100.4';  fname{604,2} = -19125;   
    fname{605,1} = 'be100.5';  fname{605,2} = -15868;
    fname{606,1} = 'be100.6';  fname{606,2} = -17368;   
    fname{607,1} = 'be100.7';  fname{607,2} = -18629;   
    fname{608,1} = 'be100.8';  fname{608,2} = -18649;
    fname{609,1} = 'be100.9';  fname{609,2} = -13294;   
    fname{610,1} = 'be100.10'; fname{610,2} = -15352;
    fname{611,1} = 'be120.3.1'; fname{611,2} = -13067;
    fname{612,1} = 'be120.3.2'; fname{612,2} = -13046;  
    fname{613,1} = 'be120.3.3'; fname{613,2} = -12418;    
    fname{614,1} = 'be120.3.4'; fname{614,2} = -13867;  
    fname{615,1} = 'be120.3.5'; fname{615,2} = -11403;    
    fname{616,1} = 'be120.3.6'; fname{616,2} = -12915;    
    fname{617,1} = 'be120.3.7'; fname{617,2} = -14068;  
    fname{618,1} = 'be120.3.8'; fname{618,2} = -14701;    
    fname{619,1} = 'be120.3.9'; fname{619,2} = -10458;  
    fname{620,1} = 'be120.3.10'; fname{620,2} = -12201;      
    fname{621,1} = 'be120.8.1';  fname{621,2} = -18691;    
    fname{622,1} = 'be120.8.2';  fname{622,2} = -18827;     
    fname{623,1} = 'be120.8.3';  fname{623,2} = -19302;    
    fname{624,1} = 'be120.8.4';  fname{624,2} = -20765;     
    fname{625,1} = 'be120.8.5';  fname{625,2} = -20417;     
    fname{626,1} = 'be120.8.6';  fname{626,2} = -18482;    
    fname{627,1} = 'be120.8.7';  fname{627,2} = -22194;     
    fname{628,1} = 'be120.8.8';  fname{628,2} = -19534;     
    fname{629,1} = 'be120.8.9';  fname{629,2} = -18195;    
    fname{630,1} = 'be120.8.10'; fname{630,2} = -19049;    
    fname{631,1} = 'be150.3.1';  fname{631,2} = -18889;    
    fname{632,1} = 'be150.3.2';  fname{632,2} = -17816;    
    fname{633,1} = 'be150.3.3';  fname{633,2} = -17314;    
    fname{634,1} = 'be150.3.4';  fname{634,2} = -19884;     
    fname{635,1} = 'be150.3.5';  fname{635,2} = -16817;    
    fname{636,1} = 'be150.3.6';  fname{636,2} = -16780;      
    fname{637,1} = 'be150.3.7';  fname{637,2} = -18001;    
    fname{638,1} = 'be150.3.8';  fname{638,2} = -18303;    
    fname{639,1} = 'be150.3.9';  fname{639,2} = -12838;    
    fname{640,1} = 'be150.3.10'; fname{640,2} = -17963;        
    fname{641,1} = 'be150.8.1';  fname{641,2} = -27089;        
    fname{642,1} = 'be150.8.2';  fname{642,2} = -26779;         
    fname{643,1} = 'be150.8.3';  fname{643,2} = -29438;         
    fname{644,1} = 'be150.8.4';  fname{644,2} = -26911;        
    fname{645,1} = 'be150.8.5';  fname{645,2} = -28017;         
    fname{646,1} = 'be150.8.6';  fname{646,2} = -29221;         
    fname{647,1} = 'be150.8.7';  fname{647,2} = -31209;        
    fname{648,1} = 'be150.8.8';  fname{648,2} = -29730;         
    fname{649,1} = 'be150.8.9';  fname{649,2} = -25388;         
    fname{650,1} = 'be150.8.10'; fname{650,2} = -28374;         
    fname{651,1} = 'be200.3.1';  fname{651,2} = -25453;   
    fname{652,1} = 'be200.3.2';  fname{652,2} = -25027;   
    fname{653,1} = 'be200.3.3';  fname{653,2} = -28023;   
    fname{654,1} = 'be200.3.4';  fname{654,2} = -27434;   
    fname{655,1} = 'be200.3.5';  fname{655,2} = -26355;   
    fname{656,1} = 'be200.3.6';  fname{656,2} = -26146;   
    fname{657,1} = 'be200.3.7';  fname{657,2} = -30483;   
    fname{658,1} = 'be200.3.8';  fname{658,2} = -27355;   
    fname{659,1} = 'be200.3.9';  fname{659,2} = -24683;   
    fname{660,1} = 'be200.3.10'; fname{660,2} = -23842;   
    fname{661,1} = 'be200.8.1';  fname{661,2} = -48534;    
    fname{662,1} = 'be200.8.2';  fname{662,2} = -40821;   
    fname{663,1} = 'be200.8.3';  fname{663,2} = -43207;   
    fname{664,1} = 'be200.8.4';  fname{664,2} = -43757;   
    fname{665,1} = 'be200.8.5';  fname{665,2} = -41482;   
    fname{666,1} = 'be200.8.6';  fname{666,2} = -49492;   
    fname{667,1} = 'be200.8.7';  fname{667,2} = -46828;   
    fname{668,1} = 'be200.8.8';  fname{668,2} = -44502;   
    fname{669,1} = 'be200.8.9';  fname{669,2} = -43241;   
    fname{670,1} = 'be200.8.10'; fname{670,2} = -42832;   
    fname{671,1} = 'be250.1';    fname{671,2} = -24076;   
    fname{672,1} = 'be250.2';    fname{672,2} = -22540;   
    fname{673,1} = 'be250.3';    fname{673,2} = -22923;   
    fname{674,1} = 'be250.4';    fname{674,2} = -24649;   
    fname{675,1} = 'be250.5';    fname{675,2} = -21057;   
    fname{676,1} = 'be250.6';    fname{676,2} = -22735;   
    fname{677,1} = 'be250.7';    fname{677,2} = -24095;   
    fname{678,1} = 'be250.8';    fname{678,2} = -23801;   
    fname{679,1} = 'be250.9';    fname{679,2} = -20051;   
    fname{680,1} = 'be250.10';   fname{680,2} = -23159;   
    fname{681,1} = 'bqp50-1';    fname{681,2} = -2098;        
    fname{682,1} = 'bqp50-2';    fname{682,2} = -3702;                
    fname{683,1} = 'bqp50-3';    fname{683,2} = -4626;                  
    fname{684,1} = 'bqp50-4';    fname{684,2} = -3544;                
    fname{685,1} = 'bqp50-5';    fname{685,2} = -4012;                  
    fname{686,1} = 'bqp50-6';    fname{686,2} = -3693;                 
    fname{687,1} = 'bqp50-7';    fname{687,2} = -4520;             
    fname{688,1} = 'bqp50-8';    fname{688,2} = -4216;                 
    fname{689,1} = 'bqp50-9';    fname{689,2} = -3780;        
    fname{690,1} = 'bqp50-10';   fname{690,2} = -3507;                      
    fname{691,1} = 'bqp100-1';   fname{691,2} = -7970;  fname{691,3} = -8.2e3;                       
    fname{692,1} = 'bqp100-2';   fname{692,2} = -11036; fname{692,3} = -1.2e4;       
    fname{693,1} = 'bqp100-3';   fname{693,2} = -12723; fname{693,3} = -1.3e4;           
    fname{694,1} = 'bqp100-4';   fname{694,2} = -10368; fname{694,3} = -1.1e4;         
    fname{695,1} = 'bqp100-5';   fname{695,2} = -9083;  fname{695,3} = -9.5e3;           
    fname{696,1} = 'bqp100-6';   fname{696,2} = -10210; fname{696,3} = -1.1e4;         
    fname{697,1} = 'bqp100-7';   fname{697,2} = -10125; fname{697,3} = -1.1e4;      
    fname{698,1} = 'bqp100-8';   fname{698,2} = -11435; fname{698,3} = -1.2e4;         
    fname{699,1} = 'bqp100-9';   fname{699,2} = -11455; fname{699,3} = -1.2e4;          
    fname{700,1} = 'bqp100-10';  fname{700,2} = -12565; fname{700,3} = -1.3e4;          
    fname{701,1} = 'bqp250-1';   fname{701,2} = -45607; fname{701,3} = -4.7e4;        
    fname{702,1} = 'bqp250-2';   fname{702,2} = -44810; fname{702,3} = -4.6e4;             
    fname{703,1} = 'bqp250-3';   fname{703,2} = -49037; fname{703,3} = -5.1e4;        
    fname{704,1} = 'bqp250-4';   fname{704,2} = -41274; fname{704,3} = -4.3e4;          
    fname{705,1} = 'bqp250-5';   fname{705,2} = -47961; fname{705,3} = -4.9e4;            
    fname{706,1} = 'bqp250-6';   fname{706,2} = -41014; fname{706,3} = -4.3e4;           
    fname{707,1} = 'bqp250-7';   fname{707,2} = -46757; fname{707,3} = -4.8e4;            
    fname{708,1} = 'bqp250-8';   fname{708,2} = -35726; fname{708,3} = -3.8e4;          
    fname{709,1} = 'bqp250-9';   fname{709,2} = -48916; fname{709,3} = -5.0e4;       
    fname{710,1} = 'bqp250-10';  fname{710,2} = -40442; fname{710,3} = -4.2e4;          
    fname{711,1} = 'bqp500-1';   fname{711,2} = -116586;                 
    fname{712,1} = 'bqp500-2';   fname{712,2} = -128223;                      
    fname{713,1} = 'bqp500-3';   fname{713,2} = -130812;                       
    fname{714,1} = 'bqp500-4';   fname{714,2} = -130097;               
    fname{715,1} = 'bqp500-5';   fname{715,2} = -125487;                       
    fname{716,1} = 'bqp500-6';   fname{716,2} = -121772;           
    fname{717,1} = 'bqp500-7';   fname{717,2} = -122201;  
    fname{718,1} = 'bqp500-8';   fname{718,2} = -123559;          
    fname{719,1} = 'bqp500-9';   fname{719,2} = -120798;           
    fname{720,1} = 'bqp500-10';  fname{720,2} = -130619;      
    fname{721,1} = 'gka1a';      fname{721,2} = -3414;      
    fname{722,1} = 'gka2a';      fname{722,2} = -6063;      
    fname{723,1} = 'gka3a';      fname{723,2} = -6037;      
    fname{724,1} = 'gka4a';      fname{724,2} = -8598;       
    fname{725,1} = 'gka5a';      fname{725,2} = -5737;      
    fname{726,1} = 'gka6a';      fname{726,2} = -3980;      
    fname{727,1} = 'gka7a';      fname{727,2} = -4541;      
    fname{728,1} = 'gka8a';      fname{728,2} = -11109;       
    fname{729,1} = 'gka1b';  fname{729,2} = -133;         
    fname{730,1} = 'gka2b';  fname{730,2} = -121;         
    fname{731,1} = 'gka3b';  fname{731,2} = -118;              
    fname{732,1} = 'gka4b';  fname{732,2} = -129;             
    fname{733,1} = 'gka5b';  fname{733,2} = -150;             
    fname{734,1} = 'gka6b';  fname{734,2} = -146;             
    fname{735,1} = 'gka7b';  fname{735,2} = -160;             
    fname{736,1} = 'gka8b';  fname{736,2} = -145;             
    fname{737,1} = 'gka9b';  fname{737,2} = -137;             
    fname{738,1} = 'gka10b'; fname{738,2} = -154;              
    fname{739,1} = 'gka1c';  fname{739,2} = -5058;         
    fname{740,1} = 'gka2c';  fname{740,2} = -6213;               
    fname{741,1} = 'gka3c';  fname{741,2} = -6665;             
    fname{742,1} = 'gka4c';  fname{742,2} = -7398;                 
    fname{743,1} = 'gka5c';  fname{743,2} = -7362;                 
    fname{744,1} = 'gka6c';  fname{744,2} = -5824;             
    fname{745,1} = 'gka7c';  fname{745,2} = -7225;                 
    fname{746,1} = 'gka1d';  fname{746,2} = -6333;                  
    fname{747,1} = 'gka2d';  fname{747,2} = -6579;                      
    fname{748,1} = 'gka3d';  fname{748,2} = -9261;                      
    fname{749,1} = 'gka4d';  fname{749,2} = -10727;                      
    fname{750,1} = 'gka5d';  fname{750,2} = -11626;                  
    fname{751,1} = 'gka6d';  fname{751,2} = -14207;                       
    fname{752,1} = 'gka7d';  fname{752,2} = -14476;                  
    fname{753,1} = 'gka8d';  fname{753,2} = -16352;                      
    fname{754,1} = 'gka9d';  fname{754,2} = -15656;                  
    fname{755,1} = 'gka10d'; fname{755,2} = -19102;                      
    fname{756,1} = 'gka1e';  fname{756,2} = -16464;                      
    fname{757,1} = 'gka2e';  fname{757,2} = -23395;                      
    fname{758,1} = 'gka3e';  fname{758,2} = -25243;                      
    fname{759,1} = 'gka4e';  fname{759,2} = -35594;                  
    fname{760,1} = 'gka5e';  fname{760,2} = -35154;                       
    fname{761,1} = 'gka1f';  fname{761,2} = -61194;   %% non-optimal 
    fname{762,1} = 'gka2f';  fname{762,2} = -100161;  %% non-optimal 
    fname{763,1} = 'gka3f';  fname{763,2} = -138035;  %% non-optimal
    fname{764,1} = 'gka4f';  fname{764,2} = -172771;  %% non-optimal  
    fname{765,1} = 'gka5f';  fname{765,2} = -190507;  %% non-optimal 

    fd(601:780) = 3*ones(1,180); 
%%
%% Qaudratic multiple knapsack problem
%%
    fname{801,1} = 'bqp100-1';                       
    fname{802,1} = 'bqp100-2';          
    fname{803,1} = 'bqp100-3';              
    fname{804,1} = 'bqp100-4';           
    fname{805,1} = 'bqp100-5';   
    fname{806,1} = 'bqp250-1';                       
    fname{807,1} = 'bqp250-2';          
    fname{808,1} = 'bqp250-3';              
    fname{809,1} = 'bqp250-4';           
    fname{810,1} = 'bqp250-5';
    %%%%% Kojima 2015/12/20 ---> 
    fname{811,1} = 'bqp500-1';                       
    fname{812,1} = 'bqp500-2';          
    fname{813,1} = 'bqp500-3';              
    fname{814,1} = 'bqp500-4';           
    fname{815,1} = 'bqp500-5';   
       
    fd(801:815) = 3.1*ones(1,15);
    %%%%% <--- Kojima 2015/12/20
%%
%% sparse qaudratic problem
%%
%     fname{901,1} = 'bqp100-1';                       
%     fname{902,1} = 'bqp100-2';          
%     fname{903,1} = 'bqp100-3';              
%     fname{904,1} = 'bqp100-4';           
%     fname{905,1} = 'bqp100-5';   
    
%     fd(901:910) = 3.2*ones(1,10);     
%%************************************************************************
