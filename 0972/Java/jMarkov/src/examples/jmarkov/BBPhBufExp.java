/**
 * BBPhBufExp.java
 * Created: 28/10/2006
 */
package examples.jmarkov;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import jmarkov.basic.exceptions.NotUnichainException;
import jphase.ContPhaseVar;

import static java.lang.Math.sqrt;

import jphase.DenseContPhaseVar;

/**
 * This class model a Bucket Brigades production line with 
 * phase type service times.
 * @author German Riano. Universidad de los Andes. (C) 2006
 */
public class BBPhBufExp extends BBPhBuf {

    /**
     * @param workers
     * @param bufferCaps
     * @param stations
     * @param velocities
     * @param vars
     * @param description
     * @param printHeads
     * @param output
     * @param debugLevel
     * @throws NotUnichainException
     * @throws IOException
     */
    public static void modelPrinting(int workers, int[] bufferCaps,
            int stations, double[][] velocities, ContPhaseVar[] vars,
            String description, boolean printHeads, BufferedWriter output,
            int debugLevel) throws NotUnichainException, IOException {

        BBPhBuf model1 = new BBPhBuf(workers, stations, vars, velocities,
                bufferCaps);

        model1.setDebugLevel(debugLevel);

        model1.setMaxStates(20000);

        model1.generate();

        if (printHeads) {
            output.write(description + "\n");

            String[] MopsNames = model1.getMOPNames();
            for (String s : MopsNames) {
                output.write(s + "\t");

            }
            output.write("Throughput Rate\n");

        }

        double[] mposAvg1 = model1.getMOPsAvg();

        for (double d : mposAvg1)
            output.write(d + "\t");
        output.write(model1.getResetRate() + "\n");

    }

    /**
     * @param workers
     * @param stations
     * @param velocities
     * @param vars
     * @param buffer
     * @param maxBuff
     * @param file
     * @param description
     * @param debugLevel
     * @throws IOException
     * @throws NotUnichainException
     */
    public static void bufferVariation(int workers, int stations,
            double[][] velocities, ContPhaseVar[] vars, boolean[] buffer,
            int maxBuff, FileWriter file, String description, int debugLevel)
            throws IOException, NotUnichainException {

        BufferedWriter output = new BufferedWriter(file);

        int[] bufferCaps = new int[stations - 1];
        boolean printHeads = true;
        for (int b = 0; b < maxBuff + 1; b++) {
            int n = 0;
            for (boolean i : buffer) {
                if (i)
                    bufferCaps[n] = b;
                n++;
            }
            modelPrinting(workers, bufferCaps, stations, velocities, vars,
                    description, printHeads, output, debugLevel);
            if (printHeads)
                printHeads = false;
        }
        output.close();

    }

    /**
     * 
     */
    public static void bufferVariation() {

        ContPhaseVar[] varsArray = new ContPhaseVar[6];
        String[] descriptions = new String[varsArray.length];

        varsArray[5] = DenseContPhaseVar.Erlang(7.0, 7);
        descriptions[5] = "Erlang distributions, 7 phases, mean 1";
        varsArray[0] = DenseContPhaseVar.Erlang(2.0, 2);
        descriptions[0] = "Erlang distributions, 2 phases, mean 1";
        varsArray[1] = DenseContPhaseVar.expo(1.0);
        descriptions[1] = "Exponential distribution, mean 1";

        double[] c2Array = { 1.5, 2.0, 2.5 };
        double[] lambda = new double[2];

        for (int i = 0; i < c2Array.length; i++) {
            lambda[0] = 1 / (1 + sqrt((c2Array[i] - 1) / 2));
            lambda[1] = 1 / (1 - sqrt((c2Array[i] - 1) / 2));
            varsArray[i + 2] = DenseContPhaseVar.HyperExpo(lambda,
                    new double[] { 0.5, 0.5 });
            descriptions[i + 2] = "HyperExponential distribution, mean 1, 2 phases, c2 "
                    + c2Array[i];
        }
        int stations = 5;
        ContPhaseVar[] vars = new ContPhaseVar[stations];
        int workers = 4;
        boolean[] buffer = new boolean[stations - 1];
        for (int i = 0; i < buffer.length; i++)
            buffer[i] = false;
        int maxBuff = 0;

        double velocities[][] = new double[workers][stations];

        for (int i = 0; i < workers; i++) {
            for (int l = 0; l < stations; l++)
                velocities[i][l] = 1;
        }
        int modelCounter = 1;
        for (int l = 0; l < buffer.length + 1; l++) {

            if (l == 0) {
                buffer[0] = true;
            } else if (l > 0 && l < buffer.length) {
                buffer[l - 1] = false;
                buffer[l] = true;
            } else {
                for (int i = 0; i < buffer.length; i++)
                    buffer[i] = true;
            }

            for (int i = 0; i < varsArray.length; i++) {

                for (int j = 0; j < vars.length; j++)
                    vars[j] = varsArray[i].copy();

                try {
                	System.out.println("./BucketFiles/model" + modelCounter
                            + ".txt");
                    FileWriter file = new FileWriter("./BucketFiles/model"+modelCounter+".txt");
                    String description;
                    if (l < buffer.length) {
                        int pos = l + 2;
                        description = descriptions[i] + "buffer position "
                                + pos;
                    } else {
                        description = descriptions[i]
                                + "all stations with buffer";
                    }

                    if (i == varsArray.length - 1) {
                        int maxBuffE = 4;
                        bufferVariation(workers, stations, velocities, vars,
                                buffer, maxBuffE, file, description, 3);
                    }
                    else{
                    bufferVariation(workers, stations, velocities, vars,
                            buffer, maxBuff, file, description, 3);
                    }

                } catch (Exception e) {
                    System.out.println("The destination file was not found.");
                }
                modelCounter++;
            }
        }

    }

    /**
     * 
     */
    public static void BottelNeckBefore() {

        ContPhaseVar[] varsArray = new ContPhaseVar[6];
        ContPhaseVar[] neckVars = new ContPhaseVar[6];
        String[] descriptions = new String[varsArray.length];

        varsArray[5] = DenseContPhaseVar.Erlang(7.0, 7);
        neckVars[5] = DenseContPhaseVar.Erlang(7.0 / 2.0, 7);
        descriptions[5] = "Erlang distributions, 7 phases, mean 1";
        varsArray[0] = DenseContPhaseVar.Erlang(2.0, 2);
        neckVars[0] = DenseContPhaseVar.Erlang(2.0 / 2.0, 2);
        descriptions[0] = "Erlang distributions, 2 phases, mean 1";
        varsArray[1] = DenseContPhaseVar.expo(1.0);
        neckVars[1] = DenseContPhaseVar.expo(1.0 / 2.0);
        descriptions[1] = "Exponential distribution, mean 1";

        double[] c2Array = { 1.5, 2.0, 2.5 };
        double[] lambda = new double[2];
        double[] lambdaNeck = new double[2];

        for (int i = 0; i < c2Array.length; i++) {
            lambda[0] = 1 / (1 + sqrt((c2Array[i] - 1) / 2));
            lambda[1] = 1 / (1 - sqrt((c2Array[i] - 1) / 2));
            lambdaNeck[0] = 1 / (2.0 * (1 + sqrt((c2Array[i] - 1) / 2)));
            lambdaNeck[1] = 1 / (2.0 * (1 - sqrt((c2Array[i] - 1) / 2)));
            varsArray[i + 2] = DenseContPhaseVar.HyperExpo(lambda,
                    new double[] { 0.5, 0.5 });
            neckVars[i + 2] = DenseContPhaseVar.HyperExpo(lambdaNeck,
                    new double[] { 0.5, 0.5 });
            descriptions[i + 2] = "HyperExponential distribution, mean 1, 2 phases, c2 "
                    + c2Array[i];
        }

        int stations = 5;
        int bottelNeck = stations - 3;
        ContPhaseVar[] vars = new ContPhaseVar[stations];
        int workers = 4;
        boolean[] buffer = new boolean[stations - 1];
        for (int i = 0; i < buffer.length; i++)
            buffer[i] = false;
        buffer[bottelNeck] = true;
        int maxBuff = 10;

        double velocities[][] = new double[workers][stations];

        for (int i = 0; i < workers; i++) {
            for (int l = 0; l < stations; l++)
                velocities[i][l] = 1;
        }
        int modelCounter = 1;

        for (int i = 0; i < varsArray.length; i++) {

            for (int j = 0; j < vars.length; j++) {
                if (j == bottelNeck)
                    vars[j] = neckVars[i].copy();
                else
                    vars[j] = varsArray[i].copy();
            }

            try {
                FileWriter file = new FileWriter(
                        "./BucketFiles/bottelNeckmodelBefore"
                                + modelCounter + ".txt");
                String description = descriptions[i] + "buffer position 3";

                bufferVariation(workers, stations, velocities, vars, buffer,
                        maxBuff, file, description, 3);

            } catch (Exception e) {
                System.out.println("The destination file was not found.");
            }
            modelCounter++;
        }
    }

    /**
     * 
     */
    public static void BottelNeckAfter() {

        ContPhaseVar[] varsArray = new ContPhaseVar[6];
        ContPhaseVar[] neckVars = new ContPhaseVar[6];
        String[] descriptions = new String[varsArray.length];

        varsArray[5] = DenseContPhaseVar.Erlang(7.0, 7);
        neckVars[5] = DenseContPhaseVar.Erlang(7.0 / 2.0, 7);
        descriptions[5] = "Erlang distributions, 7 phases, mean 1";
        varsArray[0] = DenseContPhaseVar.Erlang(2.0, 2);
        neckVars[0] = DenseContPhaseVar.Erlang(2.0 / 2.0, 2);
        descriptions[0] = "Erlang distributions, 2 phases, mean 1";
        varsArray[1] = DenseContPhaseVar.expo(1.0);
        neckVars[1] = DenseContPhaseVar.expo(1.0 / 2.0);
        descriptions[1] = "Exponential distribution, mean 1";

        double[] c2Array = { 1.5, 2.0, 2.5 };
        double[] lambda = new double[2];
        double[] lambdaNeck = new double[2];

        for (int i = 0; i < c2Array.length; i++) {
            lambda[0] = 1 / (1 + sqrt((c2Array[i] - 1) / 2));
            lambda[1] = 1 / (1 - sqrt((c2Array[i] - 1) / 2));
            lambdaNeck[0] = 1 / (2.0 * (1 + sqrt((c2Array[i] - 1) / 2)));
            lambdaNeck[1] = 1 / (2.0 * (1 - sqrt((c2Array[i] - 1) / 2)));
            varsArray[i + 2] = DenseContPhaseVar.HyperExpo(lambda,
                    new double[] { 0.5, 0.5 });
            neckVars[i + 2] = DenseContPhaseVar.HyperExpo(lambdaNeck,
                    new double[] { 0.5, 0.5 });
            descriptions[i + 2] = "HyperExponential distribution, mean 1, 2 phases, c2 "
                    + c2Array[i];
        }

        int stations = 5;
        int bottelNeck = stations - 3;
        ContPhaseVar[] vars = new ContPhaseVar[stations];
        int workers = 4;
        boolean[] buffer = new boolean[stations - 1];
        for (int i = 0; i < buffer.length; i++)
            buffer[i] = false;
        buffer[bottelNeck + 1] = true;
        buffer[bottelNeck] = true;
        int maxBuff = 10;

        double velocities[][] = new double[workers][stations];

        for (int i = 0; i < workers; i++) {
            for (int l = 0; l < stations; l++)
                velocities[i][l] = 1;
        }
        int modelCounter = 1;

        for (int i = 0; i < varsArray.length; i++) {

            for (int j = 0; j < vars.length; j++) {
                if (j == bottelNeck)
                    vars[j] = neckVars[i].copy();
                else
                    vars[j] = varsArray[i].copy();
            }

            try {
                FileWriter file = new FileWriter(
                        "./BucketFiles/bottelNeckmodelAfter"
                                + modelCounter + ".txt");
                String description = descriptions[i] + "buffer position 4";

                bufferVariation(workers, stations, velocities, vars, buffer,
                        maxBuff, file, description, 3);

            } catch (Exception e) {
                System.out.println("The destination file was not found.");
            }
            modelCounter++;
        }
    }

    public static void main(String[] s) {
        //BottelNeckBefore();
        //BottelNeckAfter();
        bufferVariation();
    }
}
