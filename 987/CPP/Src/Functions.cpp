#include "Functions.h"

//Functions of user's

long double example1(long double x) {
//Example 1
//The formula of the function is f(x)=sin(x^3)cos(x^2)-0.02
	return sin( pow(x, 3) ) * cos( pow(x, 2) ) - 0.02;
}

long double example2(long double x) {
//Example2
//type of function

	long double f = 1;
	double fct = 2;
	double sqrx = x * x; 
	int index = 2;
	double factor = 0;

	for ( int i = 1; i < 100; i++ ) {
		factor = 0;
		factor = sqrx / fct;
		if ( i % 2 ) {
			f -= factor;
		} else {
			f += factor;
		}
		index++;
		fct *= index;
		index++;
		fct *= index;
		sqrx *= ( x * x );
	}
	return f;
}


long double fun100(long double x) {
//Formula of Function
//with 100 root int the interval (0,1)

	return (x-0.00836362)*(x-0.0124448)*(x-0.0153484)*(x-0.0219223)*(x-0.0219523)*(x-0.0234919)*(x-0.0359076)*(x-0.0592197)*(x-0.0602863)*(x-0.0622228)*(x-0.0773369)*(x-0.083651)*(x-0.085323)*(x-0.0859828)*(x-0.0980329)*(x-0.10431)*(x-0.110675)*(x-0.147371)*(x-0.160367)*(x-0.181091)*(x-0.189041)*(x-0.190519)*(x-0.193729)*(x-0.20121)*(x-0.201905)*(x-0.222775)*(x-0.224292)*(x-0.248438)*(x-0.256589)*(x-0.257202)*(x-0.292521)*(x-0.292908)*(x-0.300599)*(x-0.321267)*(x-0.322693)*(x-0.329578)*(x-0.3582)*(x-0.381427)*(x-0.395688)*(x-0.405212)*(x-0.406746)*(x-0.4236)*(x-0.431124)*(x-0.436732)*(x-0.461748)*(x-0.467849)*(x-0.474095)*(x-0.475738)*(x-0.476985)*(x-0.498473)*(x-0.508556)*(x-0.527749)*(x-0.540636)*(x-0.54202)*(x-0.546126)*(x-0.548098)*(x-0.558173)*(x-0.559319)*(x-0.573197)*(x-0.577809)*(x-0.597179)*(x-0.600967)*(x-0.626183)*(x-0.665296)*(x-0.677172)*(x-0.690724)*(x-0.691258)*(x-0.701131)*(x-0.708842)*(x-0.712783)*(x-0.714582)*(x-0.717417)*(x-0.718308)*(x-0.744896)*(x-0.767514)*(x-0.778769)*(x-0.778769)*(x-0.791663)*(x-0.792761)*(x-0.801815)*(x-0.811681)*(x-0.8156)*(x-0.819628)*(x-0.821511)*(x-0.850246)*(x-0.854559)*(x-0.855139)*(x-0.871151)*(x-0.872632)*(x-0.880412)*(x-0.886919)*(x-0.917004)*(x-0.920016)*(x-0.947788)*(x-0.949985)*(x-0.959393)*(x-0.967495)*(x-0.98144)*(x-0.983395)*(x-0.985129);
}

long double fun1000(long double x) {
//Formula of Function
//with 100 root int the interval (0,1)

	return (x-0.000562666)*(x-0.00112346)*(x-0.00183367)*(x-0.0031963)*(x-0.00330965)*(x-0.00420889)*(x-0.0044626)*(x-0.00473886)*(x-0.00634148)*(x-0.00661171)*(x-0.00674534)*(x-0.00686931)*(x-0.00743619)*(x-0.00782999)*(x-0.0110646)*(x-0.0129715)*(x-0.015084)*(x-0.0172049)*(x-0.0178795)*(x-0.0184252)*(x-0.0257268)*(x-0.0262861)*(x-0.0266387)*(x-0.0275102)*(x-0.0304213)*(x-0.0313447)*(x-0.0313576)*(x-0.0313852)*(x-0.0326041)*(x-0.0332574)*(x-0.0339492)*(x-0.0341016)*(x-0.0343451)*(x-0.0345377)*(x-0.035091)*(x-0.0352085)*(x-0.0354316)*(x-0.0355405)*(x-0.0358288)*(x-0.0358545)*(x-0.0359323)*(x-0.0359596)*(x-0.0372596)*(x-0.0398126)*(x-0.0402963)*(x-0.0411317)*(x-0.0420591)*(x-0.0422122)*(x-0.0428267)*(x-0.0439453)*(x-0.0441143)*(x-0.0445568)*(x-0.0448521)*(x-0.0449988)*(x-0.0464765)*(x-0.0465497)*(x-0.0474731)*(x-0.0475527)*(x-0.0489576)*(x-0.0492647)*(x-0.0497726)*(x-0.0500797)*(x-0.0503482)*(x-0.0507847)*(x-0.053418)*(x-0.0535759)*(x-0.0549441)*(x-0.0550228)*(x-0.0561256)*(x-0.0569376)*(x-0.0578467)*(x-0.0578973)*(x-0.0582562)*(x-0.0585659)*(x-0.0590223)*(x-0.059351)*(x-0.0607373)*(x-0.0623308)*(x-0.0633453)*(x-0.0661731)*(x-0.0666139)*(x-0.0671766)*(x-0.067209)*(x-0.0688801)*(x-0.0696436)*(x-0.0718418)*(x-0.0718773)*(x-0.0745501)*(x-0.075355)*(x-0.0780202)*(x-0.0798689)*(x-0.0815193)*(x-0.0822409)*(x-0.0832546)*(x-0.0841485)*(x-0.0854632)*(x-0.0855403)*(x-0.0866496)*(x-0.0878774)*(x-0.0898861)*(x-0.091817)*(x-0.091862)*(x-0.0929765)*(x-0.0941634)*(x-0.0949234)*(x-0.0960301)*(x-0.0966953)*(x-0.0974045)*(x-0.0980475)*(x-0.0992965)*(x-0.0997719)*(x-0.101081)*(x-0.102089)*(x-0.105387)*(x-0.105487)*(x-0.108122)*(x-0.109804)*(x-0.115666)*(x-0.116202)*(x-0.116637)*(x-0.116915)*(x-0.117887)*(x-0.119458)*(x-0.120053)*(x-0.120232)*(x-0.120825)*(x-0.121757)*(x-0.122252)*(x-0.123456)*(x-0.123464)*(x-0.124248)*(x-0.124401)*(x-0.124953)*(x-0.125126)*(x-0.125671)*(x-0.126203)*(x-0.12747)*(x-0.127996)*(x-0.130977)*(x-0.130995)*(x-0.134277)*(x-0.135231)*(x-0.13546)*(x-0.136409)*(x-0.136464)*(x-0.136877)*(x-0.13776)*(x-0.138958)*(x-0.139291)*(x-0.139932)*(x-0.140271)*(x-0.14111)*(x-0.143813)*(x-0.143881)*(x-0.144186)*(x-0.1443)*(x-0.144398)*(x-0.147348)*(x-0.14759)*(x-0.148406)*(x-0.151174)*(x-0.15186)*(x-0.152173)*(x-0.152183)*(x-0.152679)*(x-0.152744)*(x-0.152935)*(x-0.153036)*(x-0.154343)*(x-0.154984)*(x-0.155162)*(x-0.156029)*(x-0.157167)*(x-0.157222)*(x-0.158116)*(x-0.158215)*(x-0.158987)*(x-0.162246)*(x-0.163014)*(x-0.164754)*(x-0.165635)*(x-0.165725)*(x-0.16669)*(x-0.166867)*(x-0.168564)*(x-0.169105)*(x-0.169185)*(x-0.169564)*(x-0.17094)*(x-0.172979)*(x-0.173697)*(x-0.174133)*(x-0.175161)*(x-0.175485)*(x-0.176636)*(x-0.176637)*(x-0.179984)*(x-0.180578)*(x-0.185588)*(x-0.185664)*(x-0.188298)*(x-0.189362)*(x-0.190403)*(x-0.190627)*(x-0.191488)*(x-0.193037)*(x-0.193524)*(x-0.193592)*(x-0.195138)*(x-0.197111)*(x-0.19762)*(x-0.198496)*(x-0.203001)*(x-0.203608)*(x-0.203937)*(x-0.203974)*(x-0.205061)*(x-0.205904)*(x-0.207394)*(x-0.207579)*(x-0.209016)*(x-0.211607)*(x-0.213356)*(x-0.213427)*(x-0.214611)*(x-0.214852)*(x-0.215109)*(x-0.215273)*(x-0.218184)*(x-0.21873)*(x-0.218927)*(x-0.219209)*(x-0.220043)*(x-0.221252)*(x-0.221489)*(x-0.225653)*(x-0.228418)*(x-0.229412)*(x-0.229608)*(x-0.23054)*(x-0.230905)*(x-0.231188)*(x-0.231894)*(x-0.233951)*(x-0.234321)*(x-0.235636)*(x-0.236552)*(x-0.238023)*(x-0.238241)*(x-0.238933)*(x-0.239045)*(x-0.240269)*(x-0.240627)*(x-0.240669)*(x-0.240808)*(x-0.240987)*(x-0.241609)*(x-0.242684)*(x-0.244304)*(x-0.24821)*(x-0.24822)*(x-0.250866)*(x-0.251119)*(x-0.251848)*(x-0.256413)*(x-0.257071)*(x-0.257469)*(x-0.261725)*(x-0.261728)*(x-0.262017)*(x-0.264937)*(x-0.265052)*(x-0.265219)*(x-0.265688)*(x-0.266666)*(x-0.267146)*(x-0.267237)*(x-0.269224)*(x-0.270051)*(x-0.27061)*(x-0.271877)*(x-0.272836)*(x-0.274612)*(x-0.276291)*(x-0.276315)*(x-0.276429)*(x-0.276469)*(x-0.276819)*(x-0.280026)*(x-0.280393)*(x-0.280931)*(x-0.282056)*(x-0.286017)*(x-0.286025)*(x-0.288828)*(x-0.29045)*(x-0.291104)*(x-0.292227)*(x-0.292303)*(x-0.292334)*(x-0.292826)*(x-0.292867)*(x-0.293305)*(x-0.294652)*(x-0.297005)*(x-0.297353)*(x-0.298231)*(x-0.299702)*(x-0.300651)*(x-0.300952)*(x-0.301035)*(x-0.303617)*(x-0.303683)*(x-0.304597)*(x-0.305046)*(x-0.306293)*(x-0.30656)*(x-0.306879)*(x-0.307149)*(x-0.307674)*(x-0.308381)*(x-0.309033)*(x-0.309078)*(x-0.309505)*(x-0.310165)*(x-0.312791)*(x-0.312893)*(x-0.313485)*(x-0.313824)*(x-0.315642)*(x-0.316163)*(x-0.316383)*(x-0.316582)*(x-0.317455)*(x-0.317818)*(x-0.319539)*(x-0.319579)*(x-0.320053)*(x-0.321034)*(x-0.321646)*(x-0.323443)*(x-0.323985)*(x-0.324047)*(x-0.324192)*(x-0.324391)*(x-0.324429)*(x-0.325044)*(x-0.326545)*(x-0.3283)*(x-0.3285)*(x-0.330262)*(x-0.331341)*(x-0.33161)*(x-0.331843)*(x-0.331874)*(x-0.333163)*(x-0.333867)*(x-0.334115)*(x-0.336282)*(x-0.336952)*(x-0.337787)*(x-0.338075)*(x-0.338408)*(x-0.339764)*(x-0.341316)*(x-0.344699)*(x-0.345045)*(x-0.345853)*(x-0.346022)*(x-0.346838)*(x-0.346881)*(x-0.34733)*(x-0.348546)*(x-0.348618)*(x-0.349729)*(x-0.350838)*(x-0.351524)*(x-0.354009)*(x-0.355431)*(x-0.355853)*(x-0.356038)*(x-0.356107)*(x-0.356829)*(x-0.357404)*(x-0.359498)*(x-0.36078)*(x-0.361316)*(x-0.361389)*(x-0.364824)*(x-0.366786)*(x-0.36686)*(x-0.367273)*(x-0.368634)*(x-0.36943)*(x-0.370294)*(x-0.370425)*(x-0.370435)*(x-0.373472)*(x-0.374542)*(x-0.375299)*(x-0.376206)*(x-0.3764)*(x-0.378294)*(x-0.378511)*(x-0.37903)*(x-0.379507)*(x-0.380188)*(x-0.380219)*(x-0.381911)*(x-0.383011)*(x-0.38307)*(x-0.383991)*(x-0.384356)*(x-0.38543)*(x-0.385481)*(x-0.385886)*(x-0.387468)*(x-0.387829)*(x-0.388565)*(x-0.389177)*(x-0.391421)*(x-0.391947)*(x-0.392319)*(x-0.392578)*(x-0.39298)*(x-0.393069)*(x-0.39386)*(x-0.395189)*(x-0.395722)*(x-0.396278)*(x-0.396459)*(x-0.397224)*(x-0.397749)*(x-0.399554)*(x-0.400212)*(x-0.400905)*(x-0.401732)*(x-0.404927)*(x-0.406573)*(x-0.407177)*(x-0.409986)*(x-0.410138)*(x-0.411959)*(x-0.413664)*(x-0.415205)*(x-0.415403)*(x-0.416)*(x-0.41712)*(x-0.420469)*(x-0.421169)*(x-0.42224)*(x-0.423462)*(x-0.426043)*(x-0.427794)*(x-0.428185)*(x-0.429996)*(x-0.430268)*(x-0.430786)*(x-0.431536)*(x-0.433651)*(x-0.434245)*(x-0.43471)*(x-0.436153)*(x-0.436176)*(x-0.436358)*(x-0.436527)*(x-0.439032)*(x-0.443162)*(x-0.443205)*(x-0.443533)*(x-0.443905)*(x-0.444059)*(x-0.445881)*(x-0.446684)*(x-0.450288)*(x-0.450345)*(x-0.451849)*(x-0.451875)*(x-0.451959)*(x-0.452991)*(x-0.455589)*(x-0.456151)*(x-0.457091)*(x-0.457361)*(x-0.458253)*(x-0.458475)*(x-0.459442)*(x-0.460324)*(x-0.461207)*(x-0.466815)*(x-0.469303)*(x-0.470601)*(x-0.472131)*(x-0.472887)*(x-0.477371)*(x-0.47761)*(x-0.477757)*(x-0.478232)*(x-0.478729)*(x-0.479573)*(x-0.480766)*(x-0.483238)*(x-0.48486)*(x-0.485817)*(x-0.486073)*(x-0.486938)*(x-0.487082)*(x-0.488621)*(x-0.490493)*(x-0.492801)*(x-0.492814)*(x-0.493253)*(x-0.494826)*(x-0.49514)*(x-0.495942)*(x-0.49662)*(x-0.498071)*(x-0.498618)*(x-0.499517)*(x-0.499547)*(x-0.499985)*(x-0.504292)*(x-0.504987)*(x-0.506222)*(x-0.506409)*(x-0.506584)*(x-0.509383)*(x-0.510155)*(x-0.510867)*(x-0.511165)*(x-0.512831)*(x-0.512877)*(x-0.512933)*(x-0.513095)*(x-0.514571)*(x-0.514579)*(x-0.515207)*(x-0.516645)*(x-0.518148)*(x-0.518721)*(x-0.524082)*(x-0.525003)*(x-0.525239)*(x-0.52559)*(x-0.526553)*(x-0.527262)*(x-0.527955)*(x-0.529041)*(x-0.530035)*(x-0.530315)*(x-0.530661)*(x-0.530832)*(x-0.531268)*(x-0.53168)*(x-0.532473)*(x-0.532642)*(x-0.534181)*(x-0.535716)*(x-0.535879)*(x-0.536299)*(x-0.536505)*(x-0.537225)*(x-0.539456)*(x-0.542824)*(x-0.544581)*(x-0.544769)*(x-0.545155)*(x-0.545432)*(x-0.545958)*(x-0.546096)*(x-0.547281)*(x-0.549878)*(x-0.550351)*(x-0.551489)*(x-0.551713)*(x-0.552406)*(x-0.553807)*(x-0.555047)*(x-0.556742)*(x-0.556963)*(x-0.558673)*(x-0.559449)*(x-0.560239)*(x-0.562021)*(x-0.563021)*(x-0.564093)*(x-0.565117)*(x-0.565119)*(x-0.567122)*(x-0.5682)*(x-0.568781)*(x-0.569105)*(x-0.570728)*(x-0.572554)*(x-0.573472)*(x-0.57426)*(x-0.574277)*(x-0.575505)*(x-0.576442)*(x-0.577891)*(x-0.578305)*(x-0.578794)*(x-0.579045)*(x-0.580142)*(x-0.58199)*(x-0.583649)*(x-0.585094)*(x-0.585369)*(x-0.586209)*(x-0.587351)*(x-0.589527)*(x-0.589838)*(x-0.590304)*(x-0.591407)*(x-0.592333)*(x-0.597668)*(x-0.598071)*(x-0.599386)*(x-0.599562)*(x-0.603002)*(x-0.603089)*(x-0.6032)*(x-0.604287)*(x-0.604521)*(x-0.604632)*(x-0.606142)*(x-0.606966)*(x-0.607747)*(x-0.607787)*(x-0.607986)*(x-0.608247)*(x-0.61171)*(x-0.613477)*(x-0.614544)*(x-0.614772)*(x-0.615033)*(x-0.615933)*(x-0.616288)*(x-0.617058)*(x-0.618904)*(x-0.620063)*(x-0.620068)*(x-0.621196)*(x-0.624503)*(x-0.625585)*(x-0.626458)*(x-0.62664)*(x-0.627909)*(x-0.628181)*(x-0.630044)*(x-0.631486)*(x-0.631731)*(x-0.631889)*(x-0.632039)*(x-0.63365)*(x-0.634402)*(x-0.634759)*(x-0.636109)*(x-0.636557)*(x-0.637457)*(x-0.638011)*(x-0.639655)*(x-0.639885)*(x-0.640568)*(x-0.642429)*(x-0.643382)*(x-0.644105)*(x-0.644769)*(x-0.645907)*(x-0.646773)*(x-0.64701)*(x-0.64758)*(x-0.648044)*(x-0.648082)*(x-0.648084)*(x-0.648248)*(x-0.652729)*(x-0.655258)*(x-0.656222)*(x-0.656402)*(x-0.660815)*(x-0.661047)*(x-0.661051)*(x-0.662778)*(x-0.663545)*(x-0.664755)*(x-0.666146)*(x-0.666309)*(x-0.667194)*(x-0.667241)*(x-0.668047)*(x-0.669792)*(x-0.67177)*(x-0.673201)*(x-0.673462)*(x-0.67359)*(x-0.673856)*(x-0.673865)*(x-0.673869)*(x-0.674064)*(x-0.674135)*(x-0.674558)*(x-0.674721)*(x-0.675066)*(x-0.675352)*(x-0.676531)*(x-0.677404)*(x-0.677592)*(x-0.67876)*(x-0.679875)*(x-0.681565)*(x-0.683512)*(x-0.683976)*(x-0.684538)*(x-0.686925)*(x-0.688837)*(x-0.690481)*(x-0.693561)*(x-0.694685)*(x-0.695969)*(x-0.696089)*(x-0.696433)*(x-0.697142)*(x-0.697728)*(x-0.700351)*(x-0.700361)*(x-0.704517)*(x-0.704656)*(x-0.704873)*(x-0.705834)*(x-0.709195)*(x-0.712504)*(x-0.713523)*(x-0.714656)*(x-0.715909)*(x-0.716616)*(x-0.717083)*(x-0.717439)*(x-0.718862)*(x-0.719924)*(x-0.720264)*(x-0.721253)*(x-0.72221)*(x-0.724086)*(x-0.72506)*(x-0.726156)*(x-0.726692)*(x-0.727149)*(x-0.727159)*(x-0.728526)*(x-0.730202)*(x-0.731028)*(x-0.731093)*(x-0.731546)*(x-0.731677)*(x-0.732514)*(x-0.734349)*(x-0.73449)*(x-0.734686)*(x-0.735315)*(x-0.736015)*(x-0.736836)*(x-0.738876)*(x-0.739243)*(x-0.741556)*(x-0.741615)*(x-0.742166)*(x-0.742685)*(x-0.742915)*(x-0.745043)*(x-0.748922)*(x-0.74947)*(x-0.749642)*(x-0.749971)*(x-0.752191)*(x-0.752587)*(x-0.753806)*(x-0.755459)*(x-0.758886)*(x-0.761235)*(x-0.762275)*(x-0.762652)*(x-0.763644)*(x-0.763751)*(x-0.764762)*(x-0.765765)*(x-0.766326)*(x-0.76733)*(x-0.767466)*(x-0.767549)*(x-0.769912)*(x-0.771984)*(x-0.772618)*(x-0.773275)*(x-0.774816)*(x-0.776892)*(x-0.777317)*(x-0.777639)*(x-0.778653)*(x-0.779509)*(x-0.779555)*(x-0.78223)*(x-0.783443)*(x-0.785611)*(x-0.787839)*(x-0.788259)*(x-0.789112)*(x-0.789238)*(x-0.789243)*(x-0.789573)*(x-0.790769)*(x-0.791207)*(x-0.79175)*(x-0.792115)*(x-0.792119)*(x-0.793165)*(x-0.794318)*(x-0.794745)*(x-0.795059)*(x-0.795875)*(x-0.797861)*(x-0.799143)*(x-0.801205)*(x-0.801645)*(x-0.801715)*(x-0.802323)*(x-0.803327)*(x-0.804451)*(x-0.805949)*(x-0.80824)*(x-0.809266)*(x-0.809947)*(x-0.811008)*(x-0.811097)*(x-0.811761)*(x-0.81228)*(x-0.812401)*(x-0.812661)*(x-0.812887)*(x-0.81316)*(x-0.816125)*(x-0.817646)*(x-0.817854)*(x-0.821218)*(x-0.822714)*(x-0.824262)*(x-0.825029)*(x-0.826521)*(x-0.831785)*(x-0.83215)*(x-0.83259)*(x-0.832837)*(x-0.836698)*(x-0.838034)*(x-0.83981)*(x-0.840224)*(x-0.841532)*(x-0.843564)*(x-0.845864)*(x-0.847884)*(x-0.848836)*(x-0.84892)*(x-0.849303)*(x-0.849669)*(x-0.85298)*(x-0.853044)*(x-0.854169)*(x-0.854303)*(x-0.854906)*(x-0.856426)*(x-0.856489)*(x-0.856497)*(x-0.857771)*(x-0.859321)*(x-0.860567)*(x-0.860694)*(x-0.861271)*(x-0.863048)*(x-0.865538)*(x-0.865867)*(x-0.867489)*(x-0.868588)*(x-0.869151)*(x-0.874973)*(x-0.8797)*(x-0.881195)*(x-0.881855)*(x-0.882347)*(x-0.882429)*(x-0.883613)*(x-0.88497)*(x-0.88531)*(x-0.88565)*(x-0.887192)*(x-0.88724)*(x-0.887655)*(x-0.888974)*(x-0.889148)*(x-0.890057)*(x-0.890199)*(x-0.890448)*(x-0.892255)*(x-0.893697)*(x-0.893703)*(x-0.894026)*(x-0.89558)*(x-0.89659)*(x-0.897687)*(x-0.899619)*(x-0.902304)*(x-0.903347)*(x-0.903674)*(x-0.907162)*(x-0.908012)*(x-0.908104)*(x-0.90831)*(x-0.909295)*(x-0.909564)*(x-0.911603)*(x-0.913398)*(x-0.913459)*(x-0.913823)*(x-0.913859)*(x-0.91399)*(x-0.914847)*(x-0.915059)*(x-0.915321)*(x-0.915392)*(x-0.917253)*(x-0.917742)*(x-0.918114)*(x-0.918464)*(x-0.919986)*(x-0.920787)*(x-0.921178)*(x-0.921748)*(x-0.922314)*(x-0.922463)*(x-0.922909)*(x-0.924575)*(x-0.924903)*(x-0.925343)*(x-0.928393)*(x-0.928561)*(x-0.929567)*(x-0.931135)*(x-0.931209)*(x-0.932049)*(x-0.93297)*(x-0.935288)*(x-0.935619)*(x-0.936327)*(x-0.936639)*(x-0.937192)*(x-0.940637)*(x-0.940852)*(x-0.94152)*(x-0.94171)*(x-0.941974)*(x-0.942239)*(x-0.942686)*(x-0.94374)*(x-0.944043)*(x-0.944185)*(x-0.945365)*(x-0.947249)*(x-0.947264)*(x-0.94731)*(x-0.948671)*(x-0.949184)*(x-0.950258)*(x-0.950696)*(x-0.955136)*(x-0.955147)*(x-0.959076)*(x-0.960061)*(x-0.961005)*(x-0.961887)*(x-0.96372)*(x-0.964379)*(x-0.964726)*(x-0.967569)*(x-0.971034)*(x-0.97169)*(x-0.972592)*(x-0.974952)*(x-0.976375)*(x-0.976643)*(x-0.977838)*(x-0.982556)*(x-0.98531)*(x-0.986497)*(x-0.98895)*(x-0.992625)*(x-0.994032)*(x-0.996076)*(x-0.99636)*(x-0.996499)*(x-0.997759)*(x-0.998811);
}


long double fun(long double x) {
//Formula of Function
//f(x)cos(x^2)*sin(x^2)-0.2 ;
	return cos( pow(x, 2) ) * sin( pow(x, 2) ) - 0.2 ;
}