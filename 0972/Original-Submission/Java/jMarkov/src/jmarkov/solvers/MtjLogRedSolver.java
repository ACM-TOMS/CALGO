/**
 * MtjLogRedSolver.java
 * Created: 2/09/2006
 */
package jmarkov.solvers;

import jmarkov.GeomProcess;
import jmarkov.basic.exceptions.NotUnichainException;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.VectorEntry;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.DenseLU;

/**
 * This class implements the Logarithmic reduction algorithm to 
 * find the steady state probabilities of a QBD process
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public class MtjLogRedSolver extends GeometricSolver {

    private Matrix R = null;
    private GeomProcess mp = null;

    /**
     * @param mp The given QBD
     */
    public MtjLogRedSolver(GeomProcess mp) {
        super(mp);
        this.mp = mp;
    }

    @Override
    public double[][] getRmatrix() throws NotUnichainException {

    	long startTimeR = System.currentTimeMillis();
    	
        double epsilon = 1e-8;
        Matrix A[] = mp.getAMatrices();
        Matrix MA0 = A[0], MA1 = A[1], MA2 = A[2];
        
        int dimen = MA1.numRows();

        Matrix I = new DenseMatrix(Matrices.identity(dimen));

        DenseMatrix mA1 = (DenseMatrix) MA1.copy();
        DenseMatrix H = (DenseMatrix) MA0.copy();
        DenseMatrix L = (DenseMatrix) MA2.copy();
        
        long startTimeS0 = System.currentTimeMillis();
        mA1.scale(-1);
        long startTimeSolve = System.currentTimeMillis();
        DenseLU mA1LU = DenseLU.factorize(mA1);
        mA1LU.solve(H);
        mA1LU.solve(L);
        
        long stopTimeS0 = System.currentTimeMillis();
        long elapsedTimeS0 = stopTimeS0 - startTimeS0;
        long elapsedTimeSolve = stopTimeS0 - startTimeSolve;
        //System.out.println("Time init mult and solve: "+elapsedTimeS0+" ms");
        
        
        DenseMatrix G = new DenseMatrix(dimen, dimen);
        DenseMatrix T = new DenseMatrix(dimen, dimen);
        DenseMatrix U = new DenseMatrix(dimen, dimen);
        DenseMatrix TA = new DenseMatrix(dimen, dimen);
        
        DenseVector one = new DenseVector(G.numRows());
        
        long startTimeE = System.currentTimeMillis();
        for(VectorEntry e: one){
        	e.set(1.0);
        }
        long stopTimeE = System.currentTimeMillis();
        long elapsedTimeE = stopTimeE - startTimeE;
        //System.out.println("Time creating unit vector e: "+elapsedTimeE+" ms");
        
        DenseVector check = new DenseVector(G.numRows());
        
        R = new DenseMatrix(dimen, dimen);
        
        G.set(L);
        T.set(H);
        double compare = 1;
        long stopTimeSolve = 0;
        long stopTimeInit = System.currentTimeMillis();
        long elapsedTimeInit = stopTimeInit - startTimeR;
        //System.out.println("Time comp R initial section: "+elapsedTimeInit+" ms");
        long elapsedTimeMult = 0;
        long elapsedTimeAdd = 0;
        long elapsedTime2 = 0;
        elapsedTimeSolve = 0;
        while (compare > epsilon) {

        	long startTimeMult = System.currentTimeMillis();
	    	
            H.mult(L, U);// U=HL
            L.multAdd(H, U);// U=HL+LH
            H = (DenseMatrix) H.mult(H, H.copy());// M=(H)^2
            
            long stopTimeMult = System.currentTimeMillis();
            elapsedTimeMult += stopTimeMult - startTimeMult;
            
	        
	        long startTimeAdd = System.currentTimeMillis();
	        U.scale(-1);
	        U.add(I);
	     
	        long stopTimeAdd = System.currentTimeMillis();
	        elapsedTimeAdd += stopTimeAdd - startTimeAdd;
	        
	        startTimeSolve = System.currentTimeMillis();
	        
	        DenseLU ULU = DenseLU.factorize(U);
	        ULU.solve(H);
	        
	        stopTimeSolve = System.currentTimeMillis();
	        elapsedTimeSolve += stopTimeSolve - startTimeSolve;
	        
            L = (DenseMatrix) L.mult(L, L.copy());// M = (L)^2
            
            startTimeSolve = System.currentTimeMillis();
	        ULU.solve(L);
	        stopTimeSolve = System.currentTimeMillis();
	        elapsedTimeSolve += stopTimeSolve - startTimeSolve;
	        
	        startTimeMult = System.currentTimeMillis();
	    	
	        T.multAdd(L, G);// G=G+TL
            TA.set(T);
            //T = new DenseMatrix(dimen, dimen);
            T.zero();
            TA.mult(H, T);// T=TH
	        
            stopTimeMult = System.currentTimeMillis();
            elapsedTimeMult += stopTimeMult - startTimeMult;
	        
	        
            // Procedimiento para el calculo de la norma infinito
            // ||1-G1||_{\infty}
	        long startTime2 = System.currentTimeMillis();
	        
	        check.zero();
	        G.mult(one, check);
            //check.add(one);
	        
            compare = 0.0;
            for (VectorEntry e: check) {
                compare = Math.max(compare, Math.abs(1 - e.get()));
            }
            
            long stopTime2 = System.currentTimeMillis();
	        elapsedTime2 += stopTime2 - startTime2;
	        
	        /*
	        System.out.println("Time iter mult: "+elapsedTimeMult+" ms");
	        System.out.println("Time iter add: "+elapsedTimeAdd+" ms");
	        System.out.println("Time iter solve: "+elapsedTimeSolve+" ms");
	        */
	        

        }
        /*
        System.out.println("Time matrix multplications: "+elapsedTimeMult+" ms");
        System.out.println("Time norm computation: "+elapsedTime2+" ms");
        System.out.println("Time solve: "+elapsedTimeSolve+" ms");
        System.out.println("Time iteration: "+(elapsedTimeSolve+elapsedTimeMult+elapsedTime2)+" ms");
        */

        //long startTimeEnd = System.currentTimeMillis();
        
        U.zero();
        U.set(MA1); // U = A0 * G + A1
        MA0.multAdd(G, U);
        U.scale(-1);
        U.transSolve(MA0.transpose(), R);
        R.transpose(); 
        
        /*
        long stopTimeEnd = System.currentTimeMillis();
        long elapsedTimeEnd = stopTimeEnd - startTimeEnd;
        System.out.println("Time comp R final section: "+elapsedTimeEnd+" ms");
        long stopTimeR = System.currentTimeMillis();
        long elapsedTimeR = stopTimeR - startTimeR;
        System.out.println("Time comp R inside solver: "+elapsedTimeR+" ms");*/

        return Matrices.getArray(R);
    }

    /**
     * @see jmarkov.solvers.Solver#label()
     */
    @Override
    public String label() {
        return "MTJ Logarithmic reduction solver";
    }

    public String description() {
        return "MTJ Logarithmic reduction solver. This solver uses the MTJ p"
                + "ackage to handle matrices.";
    }

}
