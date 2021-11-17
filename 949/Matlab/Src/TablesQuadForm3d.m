%Last modified: May 18, 2012
%First four columns=barycentric coordinates
%5th columns = weights, sum(weights)=1
%tetraN: formula integrates exactly polynomials of three variables of
%degree less than or equal to N

%integration degree=1,  N_nodes=1:
tetra1=[1/4 1/4 1/4 1/4 1];


%integration degree=3  N_nodes=4:
tetra3=[0.138196601125011   0.585410196624968   0.138196601125011   0.138196601125011   0.250000000000000
        0.138196601125011   0.138196601125011   0.585410196624968   0.138196601125011   0.250000000000000
        0.138196601125011   0.138196601125011   0.138196601125011   0.585410196624968   0.250000000000000
        0.585410196624969   0.138196601125011   0.138196601125011   0.138196601125011   0.250000000000000];


%integration degree=5  N_nodes=14:
tetra5=[0.092735250310891   0.721794249067326   0.092735250310891   0.092735250310891   0.073493043116362
        0.092735250310891   0.092735250310891   0.721794249067326   0.092735250310891   0.073493043116362
        0.092735250310891   0.092735250310891   0.092735250310891   0.721794249067326   0.073493043116362
        0.721794249067326   0.092735250310891   0.092735250310891   0.092735250310891   0.073493043116362
        0.310885919263301   0.067342242210098   0.310885919263301   0.310885919263301   0.112687925718016
        0.310885919263301   0.310885919263301   0.067342242210098   0.310885919263301   0.112687925718016
        0.310885919263301   0.310885919263301   0.310885919263301   0.067342242210098   0.112687925718016
        0.067342242210098   0.310885919263301   0.310885919263301   0.310885919263301   0.112687925718016
        0.454496295874350   0.045503704125650   0.045503704125650   0.454496295874350   0.042546020777081
        0.454496295874350   0.045503704125650   0.454496295874350   0.045503704125650   0.042546020777081
        0.045503704125650   0.045503704125650   0.454496295874350   0.454496295874350   0.042546020777081
        0.454496295874350   0.454496295874350   0.045503704125650   0.045503704125650   0.042546020777081
        0.045503704125650   0.454496295874350   0.045503704125650   0.454496295874350   0.042546020777081
        0.045503704125650   0.454496295874350   0.454496295874350   0.045503704125650   0.042546020777081];


%integration degree=7  N_nodes=25:
tetra7=[0.214602871259915   0.356191386220254   0.214602871259915   0.214602871259915   0.039922750257870
        0.214602871259915   0.214602871259915   0.356191386220254   0.214602871259915   0.039922750257870
        0.214602871259915   0.214602871259915   0.214602871259915   0.356191386220254   0.039922750257870
        0.356191386220254   0.214602871259915   0.214602871259915   0.214602871259915   0.039922750257870
        0.040673958534611   0.877978124396166   0.040673958534611   0.040673958534611   0.010077211055346
        0.040673958534611   0.040673958534611   0.877978124396166   0.040673958534611   0.010077211055346
        0.040673958534611   0.040673958534611   0.040673958534611   0.877978124396166   0.010077211055346
        0.877978124396166   0.040673958534611   0.040673958534611   0.040673958534611   0.010077211055346
        0.322337890142276   0.032986329573173   0.322337890142275   0.322337890142275   0.055357181543927
        0.322337890142276   0.322337890142275   0.032986329573173   0.322337890142275   0.055357181543927
        0.322337890142276   0.322337890142275   0.322337890142275   0.032986329573173   0.055357181543927
        0.032986329573174   0.322337890142275   0.322337890142275   0.322337890142275   0.055357181543927
        0.063661001875017   0.603005664791649   0.269672331458316   0.063661001875018   0.048214285714286
        0.063661001875017   0.603005664791649   0.063661001875018   0.269672331458316   0.048214285714286
        0.269672331458316   0.603005664791649   0.063661001875018   0.063661001875018   0.048214285714286
        0.063661001875017   0.063661001875018   0.603005664791649   0.269672331458316   0.048214285714286
        0.269672331458316   0.063661001875018   0.603005664791649   0.063661001875018   0.048214285714286
        0.269672331458316   0.063661001875018   0.063661001875018   0.603005664791649   0.048214285714286
        0.063661001875017   0.269672331458316   0.603005664791649   0.063661001875018   0.048214285714286
        0.063661001875017   0.269672331458316   0.063661001875018   0.603005664791649   0.048214285714286
        0.603005664791649   0.269672331458316   0.063661001875018   0.063661001875018   0.048214285714286
        0.063661001875017   0.063661001875018   0.269672331458316   0.603005664791649   0.048214285714286
        0.603005664791649   0.063661001875018   0.269672331458316   0.063661001875018   0.048214285714286
        0.603005664791649   0.063661001875018   0.063661001875018   0.269672331458316   0.048214285714286];


%integration degree=9  N_nodes=45:
tetra9=[0.250000000000000   0.250000000000000   0.250000000000000   0.250000000000000  -0.235962039847756
   0.127470936566639   0.617587190300083   0.127470936566639   0.127470936566639   0.024487896356056
   0.617587190300083   0.127470936566639   0.127470936566639   0.127470936566639   0.024487896356056
   0.127470936566639   0.127470936566639   0.127470936566639   0.617587190300083   0.024487896356056
   0.127470936566639   0.127470936566639   0.617587190300083   0.127470936566639   0.024487896356056
   0.032078830392632   0.903763508822103   0.032078830392632   0.032078830392632   0.003948520639826
   0.903763508822103   0.032078830392632   0.032078830392632   0.032078830392632   0.003948520639826
   0.032078830392632   0.032078830392632   0.032078830392632   0.903763508822103   0.003948520639826
   0.032078830392632   0.032078830392632   0.903763508822103   0.032078830392632   0.003948520639826
   0.450222904356719   0.450222904356719   0.049777095643281   0.049777095643281   0.026305552950737
   0.450222904356719   0.049777095643281   0.450222904356719   0.049777095643281   0.026305552950737
   0.450222904356719   0.049777095643281   0.049777095643281   0.450222904356719   0.026305552950737
   0.049777095643281   0.049777095643281   0.450222904356719   0.450222904356719   0.026305552950737
   0.049777095643281   0.450222904356719   0.049777095643281   0.450222904356719   0.026305552950737
   0.049777095643281   0.450222904356719   0.450222904356719   0.049777095643281   0.026305552950737
   0.316269552601450   0.316269552601450   0.183730447398550   0.183730447398550   0.082980383055059
   0.316269552601450   0.183730447398550   0.316269552601450   0.183730447398550   0.082980383055059
   0.316269552601450   0.183730447398550   0.183730447398550   0.316269552601450   0.082980383055059
   0.183730447398550   0.183730447398550   0.316269552601450   0.316269552601450   0.082980383055059
   0.183730447398550   0.316269552601450   0.183730447398550   0.316269552601450   0.082980383055059
   0.183730447398550   0.316269552601450   0.316269552601450   0.183730447398550   0.082980383055059
   0.513280033360881   0.022917787844817   0.231901089397151   0.231901089397151   0.025442624548102
   0.513280033360881   0.231901089397151   0.022917787844817   0.231901089397151   0.025442624548102
   0.513280033360881   0.231901089397151   0.231901089397151   0.022917787844817   0.025442624548102
   0.022917787844817   0.513280033360881   0.231901089397151   0.231901089397151   0.025442624548102
   0.022917787844817   0.231901089397151   0.513280033360881   0.231901089397151   0.025442624548102
   0.022917787844817   0.231901089397151   0.231901089397151   0.513280033360881   0.025442624548102
   0.231901089397151   0.231901089397151   0.022917787844817   0.513280033360881   0.025442624548102
   0.231901089397151   0.022917787844817   0.513280033360881   0.231901089397151   0.025442624548102
   0.231901089397151   0.513280033360881   0.231901089397151   0.022917787844817   0.025442624548102
   0.231901089397151   0.231901089397151   0.513280033360881   0.022917787844817   0.025442624548102
   0.231901089397151   0.022917787844817   0.231901089397151   0.513280033360881   0.025442624548102
   0.231901089397151   0.513280033360881   0.022917787844817   0.231901089397151   0.025442624548102
   0.193746475248804   0.730313427807538   0.037970048471829   0.037970048471829   0.013432438437685
   0.193746475248804   0.037970048471829   0.730313427807538   0.037970048471829   0.013432438437685
   0.193746475248804   0.037970048471829   0.037970048471829   0.730313427807538   0.013432438437685
   0.730313427807538   0.193746475248804   0.037970048471829   0.037970048471829   0.013432438437685
   0.730313427807538   0.037970048471829   0.193746475248804   0.037970048471829   0.013432438437685
   0.730313427807538   0.037970048471829   0.037970048471829   0.193746475248804   0.013432438437685
   0.037970048471829   0.037970048471829   0.730313427807538   0.193746475248804   0.013432438437685
   0.037970048471829   0.730313427807538   0.193746475248804   0.037970048471829   0.013432438437685
   0.037970048471829   0.193746475248804   0.037970048471829   0.730313427807538   0.013432438437685
   0.037970048471829   0.037970048471829   0.193746475248804   0.730313427807538   0.013432438437685
   0.037970048471829   0.730313427807538   0.037970048471829   0.193746475248804   0.013432438437685
   0.037970048471829   0.193746475248804   0.730313427807538   0.037970048471829   0.013432438437685];
   
return
  