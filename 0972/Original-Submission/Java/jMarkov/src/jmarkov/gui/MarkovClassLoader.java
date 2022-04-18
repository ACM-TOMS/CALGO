/**
 * 
 */
package jmarkov.gui;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * This class helps to load a SimpleMarkovProcess, State or Event class.
 * This is used by the GUI when the class has been recompiled.
 * @author Germán Riaño. Universidad de los Andes.
 */
public class MarkovClassLoader extends ClassLoader {

	/**
	 * Default constructor
	 */
	public MarkovClassLoader() {
	}

	/**
     * Loads a Class Markov Process/State/Event from a file.
	 * @param file
	 * @param instantiate
	 * @return The Class loaded
	 * @throws Exception
	 */
	public Object loadFromFile(File file, boolean instantiate) throws Exception {
		Object ob = null;
		Class<? extends Object> cl = null;
		try {
			FileInputStream is = new FileInputStream(file);
			int numBytes = is.available();
			byte[] bts = new byte[numBytes];
			is.read(bts);
			try {
				cl = defineClass(null, bts, 0, numBytes);
			} catch (LinkageError e) {
				// already loaded
				String className = file.getName();
				// remove ".class"
				className = className.substring(0, className.length() - 6);
																			
				String path = "";
				cl = this.findLoadedClass(className);
				while (cl == null && path == "") {
					path = file.getParent();
					file = new File(path);
					className = file.getName() + "." + className;
					cl = this.findLoadedClass(className);
				}
			}
			// resolveClass(cl);
			if (instantiate)
				ob = cl.newInstance();
			else
				ob = cl;
		} catch (InstantiationException e) {
			throw new Exception("In order to use this feature the class MUST "
					+ "provide a public constructor with no arguments.", e);
		} catch (IOException e) {
			throw new Exception("Unable to open file " + file, e);
		} catch (IllegalAccessError e) {
			throw new Exception(
					"The SimpleMarkovProcess, State and Event classes don't seem to belong to the same model",
					e);
		} catch (Throwable e) {
			throw new Exception(e);
		}

		return ob;
	}

}
