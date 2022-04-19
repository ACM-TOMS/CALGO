/**
 * 
 */
package jphase.distributions;

import java.util.ArrayList;

import umontreal.iro.lecuyer.rng.F2NL607;
import umontreal.iro.lecuyer.rng.GenF2w32;
import umontreal.iro.lecuyer.rng.LFSR113;
import umontreal.iro.lecuyer.rng.LFSR258;
import umontreal.iro.lecuyer.rng.MRG31k3p;
import umontreal.iro.lecuyer.rng.MRG32k3aL;
import umontreal.iro.lecuyer.rng.RandRijndael;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.WELL1024;
import umontreal.iro.lecuyer.rng.WELL512;
import umontreal.iro.lecuyer.rng.WELL607;

/**
 * This class models the set of available uniform
 * random variate generators
 * @author Andrés Sarmiento Romero. Universidad de los Andes. (C) 2013
 * 
 */

public class UniformManager {

    /**
     * Name of the random number generators
     */
	public ArrayList<String> streams;
	
	/**
	 * Information of the streamers
	 */
	public ArrayList<String> infos;

	/**
	 * Initializes the manager
	 */
	public UniformManager( ) {
		streams = new ArrayList<String>();

		streams.add("RandRijndael(Default)");
		streams.add("MRG32k3aL");
		streams.add("MRG31k3p");
		streams.add("LFSR113");
		streams.add("LFSR258");
		streams.add("WELL512");
		streams.add("WELL607");
		streams.add("WELL1024");
		streams.add("GenF2w32");
		streams.add("F2NL607");
	}
	
	/**
	 * Returns the name of the streamers
	 * @return streamers 
	 */
	public ArrayList<String> darDistribuciones() {
		return streams;
	}
	
	/**
	 * Returns the random number generator given a name
	 * @param nombre Name of the generator
	 * @return the generator
	 */
	public RandomStream darDistribucion(String nombre) {
		if(nombre.equals("MRG32k3aL")) return new MRG32k3aL();
		else if(nombre.equals("MRG31k3p")) return new MRG31k3p();
		else if(nombre.equals("LFSR113")) return new LFSR113();
		else if(nombre.equals("LFSR258")) return new LFSR258();
		else if(nombre.equals("WELL512")) return new WELL512();
		else if(nombre.equals("WELL607")) return new WELL607();
		else if(nombre.equals("WELL1024")) return new WELL1024();
		else if(nombre.equals("GenF2w32")) return new GenF2w32();
		else if(nombre.equals("F2NL607")) return new F2NL607();
		else return new RandRijndael();
	}

	/**
	 * Returns the info of a generator
	 * @param nombre Name of the generator
	 * @return The info of the generator
	 */
	public String darInfoUnif(String nombre) {
		iniInfo();
		if(nombre.equals("MRG32k3aL")) return infos.get(0);
		else if(nombre.equals("MRG31k3p")) return infos.get(1);
		else if(nombre.equals("LFSR113")) return infos.get(2);
		else if(nombre.equals("LFSR258")) return infos.get(3);
		else if(nombre.equals("WELL512")) return infos.get(4);
		else if(nombre.equals("WELL607")) return infos.get(5);
		else if(nombre.equals("WELL1024")) return infos.get(6);
		else if(nombre.equals("GenF2w32")) return infos.get(7);
		else if(nombre.equals("F2NL607")) return infos.get(8);
		else return  infos.get(9);
	}
	
	/**
	 * Initializes the info of the streamers
	 */
	public void iniInfo(){
		if(infos == null){
			infos = new ArrayList<String>();
			
			infos.add("This backbone generator has a period length of p=2^191.\n" +
				"The values of V , W, and Z are 2^51, 2^76, and 2^127, respectively.\n" +
				"(See RandomStream for their definition.) The seed of the RNG,\n" +
				"and the state of a stream at any given step, are six-dimensional\n" +
				"vectors of 32-bit integers, stored in double. The default initial\n" +
				"seed of the RNG is (12345; 12345; 12345; 12345; 12345; 12345).");
			
			infos.add("The diference between the RNG of class MRG32k3a and this one\n" +
					"is that this one has all its coeficients of the form a = +-2^q +-2^r\n" +
					"This permits a faster implementation than for arbitrary coeficients.");

			infos.add("It has four 32-bit components combined by a bitwise xor. Its period\n" +
					"length is p=2^113. The values of V , W and Z are 2^35, 2^55 and 2^90\n" +
					"respectively (see RandomStream for their denition). The seed of the RNG,\n" +
					"and the state of a stream at any given step, are four-dimensional vectors\n" +
					"of 32-bit integers. The default initial seed of the RNG is (12345; 12345;\n" +
					"12345; 12345). The	nextValue method returns numbers with 32 bits of precision.");

			infos.add("It has five components combined by a bitwise xor. Its period length is p=2^258.\n" +
					"The values of V , W and Z are 2^100, 2^100, 2^200 respectively (see RandomStream\n" +
					"for their definition). The seed of the RNG, and the state of a stream at any\n" +
					"given step, are five-dimensional vectors of 32-bit integers. The default initial\n" +
					"seed of the RNG is (1234567890; 1234567890; 1234567890; 1234567890; 1234567890).\n" +
					"The nextValue method returns numbers with 53 bits of precision. This generator is\n" +
					"fast for 64-bit machines.");

			infos.add("The backbone generator is a Well Equidistributed Long period Linear Random Number\n" +
					"Generator (WELL), proposed by F. Panneton, and which has a state size\n" +
					"of 512 bits and a period length of p=2^512. The values of V , W and Z are 2^150\n" +
					", 2^200 and 2^350 respectively (see RandomStream for their denition).\n" +
					"The seed of the RNG, and the state of a stream at any given step, is a 16-dimensional\n" +
					" vector of 32-bit integers.");

			infos.add("The implemented generator is the WELL607, which has a state size of 607 bits and\n" +
					"a period length of p=2^607. The values of V , W and Z are 2^150, 2^250 and 2^400\n" +
					"respectively (see RandomStream for their denition). The seed of the RNG, and the\n" +
					"state of a stream at any given step, is a 19-dimensional vector of 32-bit integers.\n" +
					" The output of nextValue has 32 bits of precision.");

			infos.add("Has a state size of 1024 bits and a period length of p=2^1024. The values of\n" +
					"V , W and Z are 2^300, 2^400 and 2^700 respectively (see RandomStream for their\n" +
					"definition). The seed of the RNG, and the state of a stream at any given step,\n" +
					"is a 16-dimensional vector of 32-bit integers. The output of nextValue has\n" +
					"32 bits of precision.");

			infos.add("The implemented generator is the GenF2w2_32 proposed by Panneton.\n" +
					"Its state is 25 32-bit words and it has a period length of 2^800-1. The values\n" +
					"of V , W and Z are 2^200, 2^300 and 2^500 respectively (see RandomStream for\n" +
					"their denition). The seed of the RNG, and the state of a stream at any given\n" +
					"step, is a 25-dimensional vector of 32-bits integers. Its nextValue method returns\n" +
					"numbers with 32 bits of precision.");

			infos.add("This nonlinear generator is made up of a small number of components (say n) combined\n" +
					"via addition modulo 1. Each component contains an array already filled with a \"random\"\n" +
					"permutation of {0...s-1} where s is the size of the array. These numbers and the lengths\n" +
					"of the components can be changed by the user. Each call to the generator uses the next\n" +
					"number in each array (or the rst one if we are at the end of the array). By default,\n" +
					"thereare 3 components of lengths 1019, 1021, and 1031, respectively. The non-linear\n" +
					"generator is combined with the WELL using a bitwise XOR operation.");

			infos.add("A block of 128 bits is encrypted by the Rijndael algorithm to generate 128\n" +
					"pseudo-random bits. Those bits are split into four words of 32 bits which are\n" +
					"returned successively by the method nextValue. The unencrypted block is the state\n" +
					"of the generator. It is incremented by 1 at every four calls to nextValue.\n" +
					"Thus, the period is 2^130 and jumping ahead is easy. The values of V , W and Z are\n" +
					"2^40, 2^42 and 2^82, respectively (see RandomStream for their definition).\n" +
					"Seeds/states must be given as 16-dimensional vectors of bytes (8-bit integers).\n" +
					"The default initial seed is a vector lled with zeros.");
		}
	}
}
