package jphase.values;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import cern.colt.list.DoubleArrayList;

/**
 * This class generates instances for writing an reading data
 * to/from system files.
 * @author Andrés Sarmiento Romero
 * 
 */
public class ReaderWriter {

	/**
	 * File writer
	 */
	public static PrintWriter out;
	
	/**
	 * File reader
	 */
	public static BufferedReader in;
	
	/**
	 * Given a route that contains, this method reads a file
	 * located a ´ruta´ and stores it in Double List  
	 */
	public static ArrayList<Double> cargarDatosList(String ruta) {

		ArrayList<Double> x = new ArrayList<Double>( );
		try {
			inicializarLector( ruta );
			
			String linea = in.readLine( );
			while( linea != null && linea != "" ){
				double num = Double.parseDouble(linea);
				if(num != 0){
					x.add(num);	
				}
				linea = in.readLine();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
		return x;
	}
	
	/**
	 * Given a route that contains, this method reads a file
	 * located a ´ruta´ and stores it in Double Vector  
	 */
	public static double[] saveFileToVector(String ruta) {

		ArrayList<Double> x = cargarDatosList(ruta);
		double[] d = new double[x.size()];
		for(int i = 0; i != x.size();i++){
			d[i] = x.get(i);
		}
		return d;
	}

	/**
	 * Given a system file location, initialices class 
	 * reader pointing to the file.
	 */
	public static void inicializarLector ( String ruta ) throws IOException{
		in = new BufferedReader( new FileReader( new File( ruta ) ) );
	}
	
	/**
	 * Given a file writer, uses it to write a data array.
	 * If the printwriter is null writes in console
	 */
	public static void escritor (PrintWriter pw, DoubleArrayList lista){
		if(pw == null){
			for(int i = 0; i != lista.size(); i++){
				System.out.println(lista.get(i));
			}
		}
		else{
			for(int i = 0; i != lista.size(); i++){
				pw.println(lista.get(i));
			}
			pw.close();
		}
	}
}
