/*
 * Created on 24/12/2005
 */
package jmarkov.jmdp.solvers;

import java.io.File;

import jmarkov.basic.Action;
import jmarkov.basic.State;

/**
 * This interface define the minimium elements for creating a MPS file.
 * 
 * @author Diego Bello, Germán Riaño - Universidad de Los Andes (C) 2005
 * @param <S>
 *            States Class.
 * @param <A>
 *            Action Class.
 * 
 * 
 */
public interface MpsLpSolver<S extends State, A extends Action> extends LPSolver<S, A> {

    /**
     * Returns the MPS file name.
     * 
     * @return Returns the fileName.
     */
    public abstract String getMpsFileName();

    /**
     * Returns the MPS file.
     * 
     * @return Returns the MPS generated file.
     */
    public abstract File getMpsFile();

    /**
     * Returns the working directory (where the MPS file is located).
     * 
     * @return Returns the MPS File folder.
     */
    public File getWorkingDir();
    

}