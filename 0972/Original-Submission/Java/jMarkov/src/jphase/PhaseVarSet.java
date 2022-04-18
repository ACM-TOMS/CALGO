package jphase;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Set of Phase-type variables
 */
public class PhaseVarSet implements Serializable{

	/**
     * 
     */
    private static final long serialVersionUID = -8599741746745927478L;
    
    /**
     * Variables in the set 
     */
    private List<PhaseVar> vars = new ArrayList<PhaseVar>();
    
    /**
     * Names of the variables in the set 
     */
    private List<String> names = new ArrayList<String>();
	
    /**
     * 
     */
    public String name = "Untitled";
	
    /**
     * 
     */
    public String fileName = "";
	
    /**
     * 
     */
    public boolean isDirty = false;
    
    /**
     * Default constructor
     *
     */
	public PhaseVarSet(){}
	
    /**
     * Constructs a new set with specified name 
     * @param nam Set name
     */
	public PhaseVarSet(String nam){
        name = nam;
    }
	
    /**
     * Constructs a new set with specified variables
     * @param vars Set variables
     */
	public PhaseVarSet(PhaseVar vars[]){
		for (int i=0;i < vars.length;i++){
			add(vars[i]);
		}
	}
    /**
     * Constructs a new set with specified name and variables
     * @param nam Set name
     * @param vars Set variables
     */
	public PhaseVarSet(String nam,PhaseVar vars[]){
		this(vars);
		name = nam;
		
	}
	
    /**
     * Builds a unique name for a variable from a proposed name
     * @param proposedName proposed name
     * @return Unique name for a variable from a proposed name
     */
	public String newUniqueVarName (String proposedName){
		if (indexOfName(proposedName) == -1) return proposedName;
		int k=1;
		while (indexOfName(proposedName+k)!=-1) {k++;}
		return proposedName+k;
	}

	/**
     * 
     * @param var
	 */
	public void add(PhaseVar var){
		names.add(newUniqueVarName(var.label()));
		vars.add(var);
		isDirty = true;
	}
	
    /**
     * Returns the variable at index i 
     * @param i Index of the required variables
     * @return Variable at index i
     */
	public PhaseVar varAt(int i){
		return  vars.get(i);
	}
	
    /**
     * Returns the index in the set of the variables with the specified name
     * @param s Name to be evaluated
     * @return Index in the set of the variables with the specified name
     */
	public int indexOfName(String s){
		return names.indexOf(s);
	}
	
    /**
     * Returns the number of variables in the set 
     * @return Number of variables in the set
     */
	public int numVars(){
		return vars.size();
	}
	
    /**
     * Removes the variable with specified name 
     * @param varName name of the variable to remove
     * @return index of the removed variable
     */
	public int remove(String varName){
		int idx = indexOfName(varName);
		if (idx!=-1){
			vars.get(idx);
			isDirty = true;
		}
		return idx;
	}
	
    /**
     * Remove the specified variable
     * @param var variable to remove
     */
	public void remove(PhaseVar var){
		vars.remove(var);
		isDirty = true;
	}
	
	/**
	 * Returns a vector with the means of all elements
	 * @return Vector with the means of all elements
     */
	public double[] getMeans(){
		int tot = this.vars.size();
		double vec[] = new double[tot];
		for (int i=0;i<tot;i++){
			vec[i] = varAt(i).expectedValue();
		}
		return vec;
	}
	
    /**
     * Reads a .sed file with the information of a set 
     * @param fileName File with the variable set
     * @return Set of variables in the file
     * @throws Exception IOException
     */
	public static PhaseVarSet open(String fileName)
		throws Exception
	{
		FileInputStream in = new FileInputStream(fileName);
		ObjectInputStream s = new ObjectInputStream(in);
		PhaseVarSet set = (PhaseVarSet)s.readObject();
		s.close();
		set.fileName = fileName;
		String name = (new File(fileName)).getName();
		if (name.endsWith(".sed")) name = name.substring(0,name.length()-4);
		set.name = name;
		set.isDirty = false;
		return set;
	}
	
    /**
     * Reads a .txt file with the information of a set
     * @param fileName File with the variable set
     * @return Set of variables in the file
     * @throws Exception IOException
     */
	public static PhaseVarSet openTxt(String fileName)
		throws Exception
	{
		//StringTokenizer vtk, ttk;
		PhaseVarSet set = new PhaseVarSet();
		BufferedReader in
			= new BufferedReader(new FileReader(fileName));
		set.name = in.readLine();
		while (in.ready()){
			//String line = in.readLine();
            //TODO: revise parser to ix this line
            //PhaseVar var = new DenseContPhaseVar(line);
            PhaseVar var = new DenseContPhaseVar();
			set.add(var);
		}
		in.close ();
		set.fileName = fileName;
		String name = (new File(fileName)).getName();
		if (name.endsWith(".sed")) name = name.substring(0,name.length()-4);
		set.name = name;
		set.isDirty = false;
		return set;
	}
	
    /**
     * 
     * @throws IOException
     */
	public void save() throws IOException{
		if (fileName=="") 
			fileName= name+".sed";
		save(fileName);
	}
    
    /**
     * 
     * @param fileName
     * @throws IOException
     */	
	public void save(String fileName) throws IOException{
			this.fileName = fileName;
			File archie = new File(fileName);
			String theName = archie.getName();
			int le = theName.length();
			this.name = theName.substring(0,le-4);
			FileOutputStream out = new FileOutputStream(archie);
			ObjectOutputStream s = new ObjectOutputStream(out);
			s.writeObject(this);
			s.flush();
			isDirty = false;
	}
    
    /**
     * Saves the set information in a file
     * @return True if the file could be saved, false elsewhere
     * @throws IOException
     */
	public boolean saveTxt() throws IOException{
		if (fileName=="") 
			fileName= name+".sed";
		return saveTxt(fileName);
	}
	
	/**
     * Saves the set information in a file
     * @param fileName File name
     * @return True if the file could be saved, false elsewhere
     * @throws IOException
	 */
	public boolean saveTxt(String fileName) throws IOException{
		try{
			this.fileName = fileName;
			String theName = (new File(fileName)).getName();
			int le = theName.length();
			this.name = theName.substring(0,le-4);
			BufferedWriter out
				= new BufferedWriter(new FileWriter(fileName));
			out.write (this.toString());
			out.close ();
			isDirty = false;
			return true;
		}
		catch (Exception e){
			return false;
		}
	}

	
	@Override
	public String toString(){
		String stg = name +"\n";
		int n = vars.size();
		for (int i=0;i<n;i++){
			stg += varAt(i).toString();
			stg += "\n";
		}
		return stg;
	}
	
}
