package jmarkov.jmdp.solvers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.StringTokenizer;

import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.DecisionRule;
import jmarkov.basic.Policy;
import jmarkov.basic.Solution;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.DTMDP;

/**
 * This solver solves a discounted infinite horizon MDP by building
 * and solving a linear problem using as interface Xpress-Optimizer.
 * It requires the professional version of XpressMP.
 * @author Diego Bello, Germán Riaño - Universidad de Los Andes (C)
 *         2005
 * @param <S>
 * @param <A>
 */
public class MPSXpressDiscounted<S extends State, A extends Action> extends
        MpsLpDiscountedSolver<S, A> {

    private static final String INVCOMMAS = "\"\"";
    private static final String COMMA = ",";
    private static String fs = System.getProperty("file.separator");
    private String fileName = "";
    private Scanner input;
    private String outFileName = "";

    private boolean showXpressOutput = false;
   

    /**
     * The constructor method receives a problem of the type infinite
     * DTMDP, an interest rate, and the working directory and file
     * name that the user wants as path for MPS Files and solutions
     * from Xpress-Optimizer.
     * @param problem the structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     * @param workingDir The path where will be saved and read the MPS
     *        file.
     * @param fileName The MPS file name to create.
     */
    public MPSXpressDiscounted(DTMDP<S, A> problem, double interestRate,
            String workingDir, String fileName) {
        super(problem, interestRate, workingDir, fileName);

    }

    /**
     * This is the default constructor for MPSXpressDiscounted class,
     * and defines the path for the working directory equals to the
     * MPS folder tied to the project. The constructor method
     * exclusively receives a problem of the type infinite DTMDP , an
     * interest rate that is modified for being used as discount
     * factor.
     * @param problem the structure of the problem of type infinite
     *        DTMDP.
     * @param interestRate A rate which is paid for the use of a
     *        resource.
     */
    public MPSXpressDiscounted(DTMDP<S, A> problem, double interestRate) {
        super(problem, interestRate);
    }

    @Override
    public long getIterations() {
        return 0;
    }

//    /**
//     * @see jmarkov.jmdp.solvers.MpsLpDiscountedSolver#buildSolution()
//     */
//    @Override
//    public Solution<S, A> buildSolution() throws SolverException {
//        buildDecisionRuleAndValueFunction();
//        return null;
//    }

    /**
     * This is where the actual solving takes place. The builder is
     * used to manipulate the MPS file. A reader is returned thah can
     * be use to build the solution.
     * @throws SolverException
     */
    @Override
    public void solveLP() throws SolverException {
        File mpsFile = getMpsFile();
        File commandFile = new File(getWorkingDir(), "commands.bat");
        String runCommand = "";
        String mpsFileShortName = "";
        // long initialTime = System.currentTimeMillis();
        try {
            if (!mpsFile.exists()) {
                throw new IOException("MPS File " + mpsFile.getAbsolutePath()
                        + " does not exist.");
            }

            commandFile = createCommandsFile();

            // if (!commandFile.exists()) {
            // commandFile = createCommandsFile();
            // }
            //            
            // cdWorking();
            mpsFileShortName = cutExt(mpsFile.toString());
            String comFileName = commandFile.getPath();
            // String comFileName = "commands.bat";
            runCommand = "optimizer " + mpsFileShortName + " @" + comFileName;
            Runtime runTime = Runtime.getRuntime();
            Process process = runTime.exec(runCommand);

            InputStream inputStream = process.getInputStream();
            InputStreamReader inputStreamReader = new InputStreamReader(
                    inputStream);
            BufferedReader bufferedReader = new BufferedReader(
                    inputStreamReader);
            String line = null;
            while ((line = bufferedReader.readLine()) != null) {
                if (showXpressOutput) {
                    System.out.println(line);
                }
            }
            process.waitFor();
            int exitValue = process.exitValue();
            if (exitValue != 64) {
                switch (exitValue) {
                case (65):
                    throw new SolverException(
                            "Xpress found the problem to be infeasible.");
                case (66):
                    throw new SolverException(
                            "Xpress found the problem to be unbounded.");
                case (128):
                    throw new SolverException(
                            "Xpress licensing error 4: The maximum number of simultaneous users has been reached."
                                    + " Check that no other copies of Xpress-MP are running, or contact "
                                    + "license@dashoptimization.com if you would like to upgrade to a "
                                    + "multi-user license.");
                default:
                    throw new SolverException("ERROR: Xpress exit value ="
                            + exitValue);
                }
            }
        } catch (IOException e) {
            throw new SolverException(
                    "Error invoking XPresss. This solver requires Xpress professional. The run command was "
                            + runCommand + ". Exception produced: " + e);
        } catch (InterruptedException e) {
            throw new SolverException("XPresss Execution was suspended! : "
                    + ". The run command was " + runCommand
                    + ". Exception produced: " + e);
            // e.printStackTrace();
        }
        // solLPTime = (System.currentTimeMillis() - initialTime);
        outFileName = mpsFileShortName + ".asc";
    }

    /**
     * creates a File wit the XPress commands:
     * 
     * <pre>
     * readprob
     * minim
     * writesol
     * quit
     * </pre>
     * 
     * @return a File with the commands to be read by XPress
     * @throws IOException
     */
    public File createCommandsFile() throws IOException {
        // File workingDirFile = getWorkingDir();
        File commandFile = new File("Commands.bat");
        // if (workingDirFile.exists())
        // commandFile = File.createTempFile("Commands", ".bat",
        // workingDirFile);
        // else
        // commandFile = File.createTempFile("Commands", ".bat");
        // if (!commandFile.canWrite()) {
        // commandFile = File.createTempFile("Commands", ".bat");
        // }
        // commandFile.deleteOnExit();
        Writer bw = new OutputStreamWriter(new FileOutputStream(commandFile));
        bw = new BufferedWriter(bw);
        String mpsFileName = getMpsFile().getCanonicalPath();
        bw.write("readprob " + mpsFileName + "\nminim\nwritesol\nstop\nquit");
        bw.flush();
        return commandFile;
    }

    /**
     * Returns a String path without the extension.
     * @param fil A file
     * @return A string file name with no extension.
     */
    String cutExt(String fileName) {
        int pos = fileName.lastIndexOf('.');
        if (pos == -1)
            pos = fileName.length();
        return fileName.substring(0, pos);
    }

    /** CD to working dir */
    void cdWorking() {
        File workingFolder = getWorkingDir();
        if (workingFolder.isDirectory()) {
            System.setProperty("user.dir", workingFolder.getAbsolutePath());
            String newDir = System.getProperty("user.dir");
            problem.debug(1, "Dir set to: " + newDir);
        }
    }

    @Override
    public String label() {
        return "MPS Xpress Solver (Disc)";
    }

    /**
     * Returns true if the solver the captures the Xpress o
     * @return Returns true if the solver the captures the Xpress
     *         output.
     */
    public boolean showXpressOutput() {
        return showXpressOutput;
    }

    /**
     * St this to true if you want to capture Xpress output.
     * @param showXpressOutput true if Xpress output is to be
     *        captured.
     */
    public void setShowXpressOutput(boolean showXpressOutput) {
        this.showXpressOutput = showXpressOutput;
    }

    /**
     * @throws SolverException
     * @see jmarkov.jmdp.solvers.MpsLpDiscountedSolver
     */
    public void openFile() throws SolverException {
        try {
            input = new Scanner(new File(outFileName));
        } // end try
        catch (IOException fileNotFoundException) {
            throw new SolverException("Error opening file: " + outFileName);
        } // end catch

    } // end method openFile

    /**
     * Close the MPS File
     */
    public void closeFile() {
        if (input != null)
            input.close(); // close file
    } // end method closeFile

    /**
     * A solution provided for a MPS File by Xpress-Optmizer has the
     * next framework. The columns are separated by commas, the names
     * are enclosed in inverted commas. For the Rows the most
     * important information is the name that is located in field 2
     * and the dual value, field 9. For the columns the most important
     * informations is given in field 2, Name and in field 5, the
     * value.
     * @throws SolverException
     */
    @Override
    public Solution<S, A> buildSolution() throws SolverException {
        // TODO: documentar este codigo . Es decir los diferentes
        // pasos.
        DecisionRule<S, A> decisionRule = new DecisionRule<S, A>();
        policy = new Policy<S, A>(decisionRule);
        double gain = 0;
        int counter = 0;
        int numActions = 0;
        boolean obj = false;
        States<S> states = getDiscreteProblem().getAllStates();
        Iterator<S> stateIteratorVF = states.iterator();
        Iterator<S> stateIteratorDR = states.iterator();
        S stateDR = null;
        boolean newState = true;
        Actions<A> actions = null;
        Iterator<A> actionsIterator = null;
        A action = null;
        openFile();
        
        try { // read records from file using Scanner object
            while (input.hasNextLine()) {// display record contents
                String lineAux = input.nextLine();
                StringTokenizer value = new StringTokenizer(lineAux, COMMA);
                if (discountFactor == 1.0 && obj == false) {
                    int pos = 0;
                    while (pos < 4) {
                        value.nextToken();
                        pos++;
                    }
                    gain = new Double(value.nextToken());
                    obj = true;

                } else {
                    if (obj == false) {
                        obj = true;
                    } else {
                        if (value.countTokens() == 10) {
                            S state = stateIteratorVF.next();
                            counter = 0;
                            Double answer = null;
                            while (value.hasMoreTokens()) {
                                value.nextToken();
                                counter++;
                                if (counter == 8) {
                                    String constraintValue = value.nextToken();
                                    answer = new Double(constraintValue);
                                    if (discountFactor == 1.0) {
                                        valueFunction.set(state, gain);
                                    } else {
                                        valueFunction.set(state, answer);
                                    }
                                    break;
                                }

                            }
                        }

                        if (value.countTokens() == 9) {
                            Double answer = null;
                            boolean notTransient = true;
                            if (newState) {
                                stateDR = stateIteratorDR.next();
                                actions = getProblem().feasibleActions(stateDR);
                                actionsIterator = actions.iterator();
                                numActions = actions.size();
                                newState = false;
                            }
                            counter = 0;
                            while (value.hasMoreTokens()) {
                                value.nextToken();
                                counter++;
                                if (counter == 4) {
                                    String variableValue = value.nextToken();
                                    answer = new Double(variableValue);
                                    action = actionsIterator.next();
                                    if (answer > 0) {
                                        decisionRule.set(stateDR, action);
                                        notTransient = true;
                                    }
                                    numActions--;
                                    if (numActions == 0) {
                                        if (!notTransient) {
                                            actionsIterator = actions
                                                    .iterator();
                                            decisionRule.set(stateDR,
                                                    actionsIterator.next());
                                        }
                                        newState = true;
                                    }
                                    break;
                                }

                            }

                        }

                    }
                }
            }
        } catch (NoSuchElementException elementException) {
            input.close();
            throw new SolverException("File improperly formated.",
                    elementException.getCause());
        } // end catch
        catch (IllegalStateException stateException) {
            throw new SolverException("Error reading solution from file: "
                    + fileName, stateException.getCause());
        }

        closeFile();
        policy.setDecisionRule(decisionRule);
        return new Solution<S, A>(valueFunction, policy);
    }

}
