make 
#g++ -shared -fPIC -DPIC *.o -o libgslmatrix.so -L /home/wschrep/ArithmosRelease/Libraries/Arithmos/lib/ -lArithmos -lIMpIeee -lMpIeee -lgmp -lContinuedFraction -lSpecHardware -lRational -lSpecialValue -lEasyval
g++ *.o ../err/*.o ../block/*.o -o libgslmatrix.a  -L /home/wschrep/ArithmosRelease/Libraries/Arithmos/lib/ -lArithmos -lIMpIeee -lMpIeee -lgmp -lContinuedFraction -lSpecHardware -lRational -lSpecialValue -lEasyval 
