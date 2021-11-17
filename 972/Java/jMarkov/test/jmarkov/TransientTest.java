/**
 * TransientTest.java
 * Created: Jul 6, 2005
 */
package jmarkov;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.net.URISyntaxException;
import java.net.URL;

import jmarkov.basic.State;
import jmarkov.solvers.JamaTransientSolver;
import junit.framework.TestCase;
import examples.jmarkov.DriveThru;

/**
 * @author Germán Riaño. Universidad de los Andes. (C) 2005
 */
public class TransientTest extends TestCase {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//junit.swingui.TestRunner.run(TransientTest.class);
		junit.textui.TestRunner.run(TransientTest.class);
	}

	/**
	 * @throws IOException
	 */
	public static void testTransientProbs() throws IOException {
        //TODO: check aginst MarkovExcel!!
		// as in handout:
		DriveThru theDT = new DriveThru(80.0, 120.0, 30.0, 4, 2, 1);
		JamaTransientSolver Transient = new JamaTransientSolver(theDT);
		// DriveThru theDT = new DriveThru(80.0, 120.0, 30.0, 4, 2, 2);
		theDT.setMaxStates(Integer.MAX_VALUE);
		theDT.setDebugLevel(0);
		// int [] entrada = {0};
		double[] tiempos = { .05, .10, .20, .30 };
		State sts [] = theDT.getStates().toStateArray();
		State s0 = sts[0];
		double[][] trans = Transient.getTransientProbs(tiempos, s0);
		URL url = TransientTest.class.getResource("TransientTest.class");
		try {
			File classFile = new File(url.toURI());
			String fileName = classFile.getParent()
					+ System.getProperty("file.separator") + "DTResults.txt";
			File outputFile = new File(fileName);
			PrintWriter outF = new PrintWriter(new BufferedWriter(
					new FileWriter(outputFile)));

			int filas = trans.length;
			int columnas = trans[0].length;
			System.out.println("filas" + filas);
			System.out.println("Columnas" + columnas);

			for (int i = 0; i < trans.length; i++) {
				for (int j = 0; j < trans[0].length; j++) {
					outF.print(trans[i][j] + "  ");
					System.out.print(trans[i][j] + "  ");
				}
				outF.println("  ");
				System.out.println("  ");
			}

			outF.close();
		} catch (IOException e) {
			return;
		} catch (URISyntaxException e) {
			return;
		}
		// theDT.showGUI();
		theDT.printAll();
		theDT.printMOPs();
	}

}
