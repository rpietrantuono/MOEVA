package util;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import causalOptimization.Main;
import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import solutionType.CauseEffectSolutionType;

public class SolutionUtils {

public static ArrayList<String> getNumericRepresentation(Solution sol) throws Exception {
	
		ArrayList<String> sol_string = new ArrayList<String>();
		if (sol.getType() instanceof CauseEffectSolutionType) {
			for(int j=0; j< ((Sources)(sol.getDecisionVariables()[0]))._sourcesList.size();j++)
				sol_string.add(((Sources)(sol.getDecisionVariables()[0]))._sourcesList.get(j));
			for(int j=0; j< ((Targets)(sol.getDecisionVariables()[1]))._targetsList.size();j++)
				sol_string.add(((Targets)(sol.getDecisionVariables()[1]))._targetsList.get(j));
		} 
		else
			throw new Exception("Solutions are not instances of CauseEffectSolutionType type"); 
		return  sol_string;
		
	}



public static double getDoubleRepresentation(String s) throws Exception {
	int numericRepresentation = Integer.parseInt(s);
	double maxSolutionSize = Main.possibleSources.size() + Main.possibleTargets.size();
	if(numericRepresentation > maxSolutionSize || numericRepresentation<0)
		throw new Exception("Conversion Error: the provided value is outside the limits of the numerical representation"); 
	return Integer.parseInt(s)/maxSolutionSize ; 	
	
}

public static String getStringRepresentation(Double d) throws Exception {
	
	int maxSolutionSize = Main.possibleSources.size() + Main.possibleTargets.size();
	long numericRepresentation = Math.round(d*maxSolutionSize);
	if(numericRepresentation > maxSolutionSize || numericRepresentation<0)
		throw new Exception("Conversion Error: the provided value is outside the limits of the numerical representation: "+d); 
	return String.valueOf(numericRepresentation);
	
}


public static ArrayList<String> getSources(Solution sol){
	
	ArrayList<String> sol_string = new ArrayList<String>();
	for(int j=0; j< ((Sources)(sol.getDecisionVariables()[0]))._sourcesList.size();j++)
		sol_string.add(((Sources)(sol.getDecisionVariables()[0]))._sourcesList.get(j));
	
	return sol_string;
	}


public static ArrayList<String> getTargets(Solution sol){
	
	ArrayList<String> sol_string = new ArrayList<String>();
	for(int j=0; j< ((Targets)(sol.getDecisionVariables()[1]))._targetsList.size();j++)
		sol_string.add(((Targets)(sol.getDecisionVariables()[1]))._targetsList.get(j));
	
	return sol_string;
	}



public static Solution createSolutionFromString(ArrayList<String> solution_string) {

	Variable[] variables = new Variable[2];
	
	ArrayList<String> sources = new ArrayList<String>();
	ArrayList<String> targets = new ArrayList<String>();
	for (int i=0; i<solution_string.size(); i++ ) { 
		if (Integer.parseInt(solution_string.get(i)) <= Main.targetLimitNumber)
			sources.add(solution_string.get(i));
		else
			targets.add(solution_string.get(i));
	}
	variables[0] = new Sources(sources);
	variables[1] = new Targets(targets);
	
	return new Solution(Main.myProblemString,variables); 
}


public static long binomial(int n, int k) {
    if ((n == k) || (k == 0))
        return 1;
    else
        return binomial(n - 1, k) + binomial(n - 1, k - 1);
}

public static void writeMultipleSolutionVariablesToFile(String filename, SolutionSet solutions) {
	for (int i=0; i<solutions.size();i++)
		writeSolutionVariablesToFile(filename, solutions.get(i)); 
}

public static void writeSolutionVariablesToFile(String filename, Solution indiv) {
	FileWriter output=null; 
	PrintWriter pw=null;
	try {
		output  = new FileWriter(filename, true);
		pw = new PrintWriter(new BufferedWriter(output));	
	} catch (IOException e) {e.printStackTrace();}
	
	String line="";
	String line_names="";
	ArrayList<String> sources = new ArrayList<String>();
	ArrayList<String> targets= new ArrayList<String>();
	sources = ((Sources)(indiv.getDecisionVariables()[0]))._sourcesList;
	targets = ((Targets)(indiv.getDecisionVariables()[1]))._targetsList;
	sources.sort(Comparator.comparing(Integer::parseInt));
	targets.sort(Comparator.comparing(Integer::parseInt));
	
	for(int j=0; j< sources.size();j++){
		line = line + sources.get(j) + " ";
		line_names = line_names + Main.getNameFromKey(sources.get(j)) + " ";
	}//end for sources
		
	for(int j=0; j< targets.size();j++){
		line = line + targets.get(j) + " "; //targets
		line_names = line_names + Main.getNameFromKey(targets.get(j)) + " ";
		}
	
	pw.println("\n");
	pw.println(line_names);
	pw.println(line);
	pw.flush();
	pw.close();
		
}

public static void writeSolutionObjectivesToFile(String filename, Solution indiv) {
	FileWriter output=null; 
	PrintWriter pw=null;
	try {
		output  = new FileWriter(filename, true);
		pw = new PrintWriter(new BufferedWriter(output));	
	} catch (IOException e) {e.printStackTrace();}
	
	String line = ""+indiv.getObjective(0);//+" "+indiv.getObjective(1);
		
	pw.println(line);
	pw.flush();
	pw.close();		
}


public static void writePopulationStatisticsToFile(String filename, SolutionSet population) {
	
	
	FileWriter output=null; 
	PrintWriter pw=null;
	try {
		output  = new FileWriter(filename, true);
		pw = new PrintWriter(new BufferedWriter(output));	
	} catch (IOException e) {e.printStackTrace();}
	
	String line="";
	
	ArrayList<Double> sol_set_plaus = new ArrayList<Double>();
	//ArrayList<Double> sol_set_nov = new ArrayList<Double>();
	for (int i=0 ; i< population.size(); i++) { 
		sol_set_plaus.add(population.get(i).getObjective(0));
		//sol_set_nov.add(population.get(i).getObjective(1));
	}
	Collections.sort(sol_set_plaus);
	//Collections.sort(sol_set_nov);
	
	//double averagePlausibility =sol_set_plaus.stream().mapToDouble(val -> val).average().orElse(0.0);
	double averagePlausibility = mean(sol_set_plaus); 
	//double averageNovelty = mean(sol_set_nov); 
	
	double medianPlausibility; 
	if (sol_set_plaus.size()>1)
		medianPlausibility = median(sol_set_plaus);
	else
		medianPlausibility =  sol_set_plaus.get(0);
	
//	double medianPlausibility = median(sol_set_plaus);
	//double medianNovelty = median(sol_set_nov);

	double worstPlausibility = Collections.max(sol_set_plaus); 
	//double worstNovelty = Collections.max(sol_set_nov); 
	
	double bestPlausibility = Collections.min(sol_set_plaus);
	//double bestNovelty = Collections.min(sol_set_nov);
	
	//THIS SHOULD BE PRINTED OFFLINE, at file initialization line = line + "Average Plausibility, Average Novelty, Median Plausibility, Median Novelty, Best Plausibility, Best Novelty, Worst Plausibility, Worst Novelty"; 
	line = averagePlausibility+", " +medianPlausibility+", "+bestPlausibility+", "+worstPlausibility;
	pw.println(line);
	pw.flush();
	pw.close();		
	
	
	
}



public static void writeStatisticsToFile(String filename, ArrayList<Double> values) {
	
	
	
	// FOR A BI-OBJECTIVE solution
	FileWriter output=null; 
	PrintWriter pw=null;
	try {
		output  = new FileWriter(filename, true);
		pw = new PrintWriter(new BufferedWriter(output));	
	} catch (IOException e) {e.printStackTrace();}
	
	String line="";
	
	
	if(values.get(0)==-1) {
		line = " , , , , ";
		pw.println(line);
		pw.flush();
		pw.close();	
		return;
	}
	Collections.sort(values);
	
	double average = mean(values); 
	
	double median; 
	if (values.size()>1)
		median = median(values);
	else
		median =  values.get(0);
	
	double worst = Collections.max(values); 
	
	double best = Collections.min(values);
	
	double std=0;
	try {
		std = sd(values);
	} catch (Exception e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	//line = line + "Average, Median , Best, Worst, Std";
	//pw.println(line);
	line = average+", " +median+", "+best+", "+worst+", "+std;
	pw.println(line);
	pw.flush();
	pw.close();		
}


private static double sum (List<Double> a){
    if (a.size() > 0) {
        double sum = 0;

        for (Double i : a) {
            sum += i;
        }
        return sum;
    }
    return 0;
}
private static double mean (List<Double> a){
    double sum = sum(a);
    double mean = 0;
    mean = sum / (a.size() * 1.0);
    return mean;
}
private static double median (List<Double> a){
    int middle = a.size()/2;
    if (a.size() % 2 == 1) {
        return a.get(middle);
    } else {
       return (a.get(middle-1) + a.get(middle)) / 2.0;
    }
}


private static double sd (List<Double> a) throws Exception{
    int sum = 0;
    double mean = mean(a);
    if (a.size()==1)
    	return  a.get(0); 
    if (a.size()==0)
    	throw new Exception("Standard deviation computation: the input vector is empty");
    for (Double i : a)
        sum += Math.pow((i - mean), 2);
    return Math.sqrt( sum / ( a.size() - 1 ) ); // sample
}

public static double[] Quartiles(double[] val) {
    double ans[] = new double[3];

    for (int quartileType = 1; quartileType < 4; quartileType++) {
        float length = val.length + 1;
        double quartile;
        float newArraySize = (length * ((float) (quartileType) * 25 / 100)) - 1;
        Arrays.sort(val);
        if (newArraySize % 1 == 0) {
            quartile = val[(int) (newArraySize)];
            } else {
            int newArraySize1 = (int) (newArraySize);
            quartile = (val[newArraySize1] + val[newArraySize1 + 1]) / 2;
             }
        ans[quartileType - 1] =  quartile;
    }
    return ans;
}

} 
