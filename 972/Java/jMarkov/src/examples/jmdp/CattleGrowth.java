package examples.jmdp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import jmarkov.MarkovProcess;
import jmarkov.basic.Action;
import jmarkov.basic.Actions;
import jmarkov.basic.ActionsSet;
import jmarkov.basic.State;
import jmarkov.basic.States;
import jmarkov.basic.StatesSet;
import jmarkov.basic.ValueFunction;
import jmarkov.basic.exceptions.SolverException;
import jmarkov.jmdp.FiniteMDP;
import jmarkov.jmdp.solvers.FiniteSolver;
import Jama.Matrix;

/**
 * The Cattle Growth problem represent a male cattle, the states
 * represent the weight at different ages. This problem decide whether
 * keep, sell or give forage to the animal in order to maximize the
 * profit of his future sell, by minimizing the cost of feeding it.
 * The weight of each animal in the next state depends on the
 * transition probabilities matrix (one for each action taken at each
 * stage).
 * @author Ana Karina Maroso, Germán Riaño
 */
public class CattleGrowth extends FiniteMDP<AnimalWeight, AnimalActions> {

    private double keep[][][] = null;
    // TODO: Sugerencia
    // private double keep[][][] = null;// el primer parametro es la
    // etapa t

    private double forage[][][] = null;
    private double preciomx[] = null;
    private double forageCost = 0;

    private ActionsSet<AnimalActions> actions = new ActionsSet<AnimalActions>();

    /**
     * @param initial
     * @param horizon
     * @param keep
     * @param sell
     * @param forage
     * @param preciomx
     * @param forageCost
     */
    private CattleGrowth(States<AnimalWeight> initial, int horizon,
            double[][][] keep, double[][][] forage, double preciomx[],
            double forageCost) {
        super(initial, horizon);
        // super(initial, horizon);
        this.keep = keep;
        this.forage = forage;
        this.preciomx = preciomx;
        this.forageCost = forageCost;
        // initStates();
        actions.add(AnimalActions.SELL);
        actions.add(AnimalActions.KEEP);
        actions.add(AnimalActions.FORAGE);
    }

    /**
     * Builds a cattle Weight problem with the given initial weight.
     * @param initialWeight
     * @param horizon
     * @param keep
     * @param forage
     * @param preciomx
     * @param forageCost
     */
    public CattleGrowth(int initialWeight[], int horizon, double[][][] keep,
            double[][][] forage, double preciomx[], double forageCost) {
        this(initStates(initialWeight), horizon, keep, forage, preciomx,
                forageCost);

    }

    /**
     * @param initialWeight
     * @param horizon
     * @param keepFile
     * @param forageFile
     * @param preciomx
     * @param forageCost
     * @throws IOException 
     * @throws URISyntaxException 
     * @throws FileNotFoundException 
     */
    public CattleGrowth(int initialWeight[], int horizon, String keepFile,
            String forageFile, double preciomx[], double forageCost) throws FileNotFoundException, IOException {
        this(initStates(initialWeight), horizon, initMatricesFromFile(keepFile,
                horizon), initMatricesFromFile(forageFile, horizon), preciomx,
                forageCost);
    }

    private static States<AnimalWeight> initStates(int weights[]) {
        StatesSet<AnimalWeight> stts = new StatesSet<AnimalWeight>();
        for (int weight : weights) {
            stts.add(new AnimalWeight(weight));
        }
        return stts;
    }

    private static double[][][] initMatricesFromFile(String filePattern,
            int horizon) throws IOException{
        double[][][] result = new double[horizon][][];
        for (int t = 0; t < horizon; t++) {
            Matrix mat = loadJamaMatrix(filePattern + t + ".txt");
            result[t] = mat.getArray();
        }
        return result;
    }

    /**
     * @see jmarkov.jmdp.FiniteMDP#finalCost(jmarkov.basic.State)
     */
    @Override
    public double finalCost(AnimalWeight i) {
        // TODO Auto-generated method stub
        return 0.0;
    }

    /**
     * @param initWeight
     * @return The optimal value for this initial weight
     * @throws SolverException
     */
    public double getValue(int initWeight) throws SolverException {
        ValueFunction<AnimalWeight> vf = getOptimalValueFunction();
        AnimalWeight state = new AnimalWeight(initWeight);
        return vf.get(state);
    }

    @Override
    public double immediateCost(AnimalWeight i, AnimalActions a, int t) {
        // (costoForraje * 65) + (t * costoForraje * 51)
        // TODO Poner costo = -precio de venta cuando la accion es
        // SELL
        if (a == AnimalActions.FORAGE)
            return (forageCost);
        else if (a == AnimalActions.SELL)
            return (-800000);
        return 0.0;
    }

    @Override
    public Actions<AnimalActions> feasibleActions(AnimalWeight st, int t) {
        if ((st.isSold()) || (t == 0)) {
            ActionsSet<AnimalActions> act = new ActionsSet<AnimalActions>();
            act.add(AnimalActions.KEEP);
            return act;
        } else if (st.weight >= 600) {
            ActionsSet<AnimalActions> act = new ActionsSet<AnimalActions>();
            act.add(AnimalActions.SELL);
            return act;
        } else if (t == 4) {
            ActionsSet<AnimalActions> act = new ActionsSet<AnimalActions>();
            act.add(AnimalActions.SELL);
            return act;
        } else
            return actions;
    }

    @Override
    public double prob(AnimalWeight pesoactual, AnimalWeight pesofuturo,
            AnimalActions a, int t) {
        if (a == AnimalActions.KEEP)
            return keep[t][(pesoactual.weight - 160) / 20][(pesofuturo.weight - 160) / 20];
        else if (a == AnimalActions.SELL)
            // return sell[(pesoactual.weight - 40) /
            // 20][(pesofuturo.weight - 40) / 20];
            return (pesofuturo.isSold()) ? 1.0 : 0.0;
        else if (a == AnimalActions.FORAGE)
            return forage[t][(pesoactual.weight - 160) / 20][(pesofuturo.weight - 160) / 20];
        return 0;
    }

    @Override
    public States<AnimalWeight> reachable(AnimalWeight arg0,
            AnimalActions arg1, int t) {
        StatesSet<AnimalWeight> statesSet = new StatesSet<AnimalWeight>();
        // Available inventory upon order receival:
        for (int i = 0; i < 22; i++) {
            statesSet.add(new AnimalWeight(160 + (i * 20)));
        }
        return statesSet;
    }

    /**
     * @param fileName Name of the text file with the data
     * @return A jama Matrix with the data
     * @throws IOException 
     * @throws URISyntaxException 
     * @throws IOException 
     */

    public static Matrix loadJamaMatrix(String fileName) throws  IOException {

        // for(int i = 0; i < 10; i++) {no se bien como implementar
        // este
        // for el me dijo q lo hiciera para leer las diferentes
        // matrices
        // no se si vendria aca.
        Matrix Mat = null;
        File fil = null;
        fil = new File(fileName);
        Mat = Matrix.read(new BufferedReader(new FileReader(fil)));

        return Mat;

    }

    /**
     * Small test program
     * @param a Not used
     * @throws SolverException
     * @throws IOException 
     * @throws URISyntaxException 
     * @throws FileNotFoundException 
     */
    // no se bien como cambiaria esta parte, pero las matrices ahora
    // se deben
    // leer con el modulo de Jama matrix(arriba). Ahora tengo 1 matriz
    // para cada
    // etapa y como es una para keep y otra para forage son 10
    // matrices
    // donde dependen de tres parametros (por eso el cambio en los
    // corchetes
    // el primero es en funcion de t para indicar la etapa de la
    // matriz
    // q se necesite.
    public static void main(String a[]) throws SolverException, FileNotFoundException, URISyntaxException, IOException {
        double forageCost = 200;
        // este error aparece porque cambie la entrada del peso
        // inicial de 40
        // a
        // el arreglo de los posibles estados [160, 180, 200, 220,
        // 240, 260,
        // 280,
        // 300, 320, 340]. Pero yo pienso q no se puede cambiar todos
        // los
        // parametros
        // weight a arreglo porque despues no se podria manejar como
        // un entero
        // al
        // hacer las transiciones en las matrices o al calcular los
        // costos de
        // acuerdo al peso.

        double preciomx[] = null;
        int[] initialWeights = new int[] { 160, 180, 200, 220, 240, 260, 280,
                300, 320, 340 };

        CattleGrowth prob = new CattleGrowth(initialWeights, 5, "FORAGE",
                "KEEP", preciomx, forageCost);
        FiniteSolver<AnimalWeight, AnimalActions> theSolver = new FiniteSolver<AnimalWeight, AnimalActions>(
                prob);
        theSolver.solve();
        prob.printSolution();
    }
}

// no se porque aparece este error si esto esta tal cual lo hicimos
// ese dia nosotros
class AnimalWeight extends State {
    int weight; // weight of Animal

    AnimalWeight(int wght) {
        weight = wght;
    }

    @Override
    public int compareTo(State s) {
        AnimalWeight i = (AnimalWeight) s;
        return weight - i.weight;
    }

    @Override
    public String label() {
        return weight + " Kg";
    }

    /**
     * @return True if the animal was sold.
     */
    public boolean isSold() {
        return weight == 0;
    }

    /**
     * @see jmarkov.basic.State#isConsistent()
     */
    @Override
    public boolean isConsistent() {
        // TODO Auto-generated method stub
        return true;
    }

    /**
     * @see jmarkov.basic.State#computeMOPs(jmarkov.MarkovProcess)
     */
    @Override
    public void computeMOPs(MarkovProcess<?, ?> model) {
        // TODO Auto-generated method stub

    }

}

// tampoco entiedo este error

class AnimalActions extends Action {
    /** Sell the Animal */
    public static final AnimalActions SELL = new AnimalActions("Sell");
    /** Keep all the Animal */
    public static final AnimalActions KEEP = new AnimalActions("Keep");
    /** Keep all the Animal */
    public static final AnimalActions FORAGE = new AnimalActions("Forage");

    String name;

    private AnimalActions(String str) {
        name = str;
    }

    /**
     * @see jmarkov.basic.Action#label()
     */
    @Override
    public String label() {
        return name;
    }

    /**
     * @see java.lang.Comparable#compareTo(Object)
     */
    public int compareTo(Action o) {
        return name.compareTo(o.label());
    }

}
