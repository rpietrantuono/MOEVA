package causalOptimization;

import util.SolutionUtils;
import jmetal.core.*;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList;
import jmetal.util.Ranking;
import jmetal.util.archive.CrowdingArchive;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.CrowdingDistanceComparator;
import jmetal.util.comparators.DominanceComparator;
import jmetal.util.comparators.EpsilonDominanceComparator;
import jmetal.operators.creativeAbduction.*; 
import jmetal.operators.factualAbduction.*;
import jmetal.operators.analogicalAbduction.*;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

//import org.apache.commons.lang3.RandomUtils;
import java.util.Random;

import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;

public class EVA extends Algorithm{
  
	public EVA(Problem problem){
		super (problem);
   }
	
	
	

  /** Execute the algorithm 
   * @throws JMException */
  public SolutionSet execute() throws JMException, ClassNotFoundException {

	  double starting_time=System.nanoTime();
	  //Init the parameters
	  
    int populationSize_factual, populationSize_analogical, populationSize_creative, maxEvaluations, evaluations, maxEvaluationPerCycle, externalArchiveSize;
    Operator factualAbductionOperator, analogicalAbductionOperator, creativeAbductionOperator, selectionOperator_factual, selectionOperator_analogical, selectionOperator_creative;
    SolutionSet population_factual, population_analogical, population_creative, external_archive_analogical ;
    SolutionSet offspringPopulation_factual, offspringPopulation_analogical, offspringPopulation_creative;
    
    
    
    // Read the parameters
    populationSize_factual    = ((Integer)getInputParameter("populationSize_factual")).intValue();
    populationSize_analogical    = ((Integer)getInputParameter("populationSize_analogical")).intValue();
    populationSize_creative    = ((Integer)getInputParameter("populationSize_creative")).intValue();
    maxEvaluations    = ((Integer)getInputParameter("maxEvaluations")).intValue();
    maxEvaluationPerCycle  = ((Integer)getInputParameter("maxEvaluationPerCycle")).intValue();
    externalArchiveSize = ((Integer)getInputParameter("externalArchiveSize")).intValue(); 
    /*maxSources = ((Integer)getInputParameter("maxSources")).intValue();
    maxTargets = ((Integer)getInputParameter("maxTargets")).intValue();
    minSources = ((Integer)getInputParameter("minSources")).intValue();
    minTargets = ((Integer)getInputParameter("minTargets")).intValue();
     */
    
//  Comparator comparator = new ObjectiveComparator(0) ; // Single objective comparator
    
  // Create the dominator for equadless and dominance: csMOPSO
  Distance distance           = new Distance();
  int archiveSize_ = populationSize_factual+populationSize_analogical+populationSize_creative;
  double eta_ = 0.0075;
  NonDominatedSolutionList eArchive_      = new NonDominatedSolutionList(new EpsilonDominanceComparator(eta_));
  SolutionSet union_factual, union_analogical, union_creative;
    
    // ALl this should be moved in the main and passed as parameter
    try {
		Files.deleteIfExists(Paths.get(Main.baseDir + "/data/sol_objectives_factual"+Main.nameToAppendToFiles+".txt"));
		Files.deleteIfExists(Paths.get(Main.baseDir + "/data/sol_variables_factual"+Main.nameToAppendToFiles+".txt"));
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/statistics_factual"+Main.nameToAppendToFiles+".csv"));
	    
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/sol_objectives_analogical"+Main.nameToAppendToFiles+".txt")); 
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/sol_variables_analogical"+Main.nameToAppendToFiles+".txt"));
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/statistics_analogical"+Main.nameToAppendToFiles+".csv"));
	    
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/sol_objectives_creative"+Main.nameToAppendToFiles+".txt")); 
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/sol_variables_creative"+Main.nameToAppendToFiles+".txt"));
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/statistics_creative"+Main.nameToAppendToFiles+".csv"));
	    
	    Files.deleteIfExists(Paths.get(Main.baseDir + "/data/statistics_by_generation"+Main.nameToAppendToFiles+".csv"));
	    
    } catch (IOException e1) {
		// TODO Auto-generated catch block
		e1.printStackTrace();
	}
    String filename_obj_factual =  Main.baseDir + "/data/sol_objectives_factual"+Main.nameToAppendToFiles+".txt"; 
    String filename_var_factual =  Main.baseDir + "/data/sol_variables_factual"+Main.nameToAppendToFiles+".txt";
    String filename_statistics_factual =  Main.baseDir + "/data/statistics_factual"+Main.nameToAppendToFiles+".csv";
    
    String filename_obj_analogical =  Main.baseDir + "/data/sol_objectives_analogical"+Main.nameToAppendToFiles+".txt"; 
    String filename_var_analogical =  Main.baseDir + "/data/sol_variables_analogical"+Main.nameToAppendToFiles+".txt";
    String filename_statistics_analogical =  Main.baseDir + "/data/statistics_analogical"+Main.nameToAppendToFiles+".csv";
    
    String filename_obj_creative =  Main.baseDir + "/data/sol_objectives_creative"+Main.nameToAppendToFiles+".txt"; 
    String filename_var_creative =  Main.baseDir + "/data/sol_variables_creative"+Main.nameToAppendToFiles+".txt";
    String filename_statistics_creative =  Main.baseDir + "/data/statistics_creative"+Main.nameToAppendToFiles+".csv";
    
    String filename_statistics_generation =  Main.baseDir + "/data/statistics_by_generation"+Main.nameToAppendToFiles+".csv";
    
    // Read the operators
    factualAbductionOperator = operators_.get("factualAbduction");
    analogicalAbductionOperator  = operators_.get("analogicalAbduction");
    creativeAbductionOperator = operators_.get("creativeAbduction");        
    selectionOperator_factual = operators_.get("selection_factual");
    selectionOperator_analogical = operators_.get("selection_analogical");
    selectionOperator_creative = operators_.get("selection_creative");
    
    // Initialize the variables
    population_factual 			= new SolutionSet(populationSize_factual);
    population_analogical 			= new SolutionSet(populationSize_analogical);
    population_creative 			= new SolutionSet(populationSize_creative);
    offspringPopulation_factual 		= new SolutionSet(populationSize_factual);
    offspringPopulation_analogical 		= new SolutionSet(populationSize_analogical);
    offspringPopulation_creative 		= new SolutionSet(populationSize_creative);
    
    external_archive_analogical = new SolutionSet(externalArchiveSize); 
    
    evaluations        = 0;                        
	
	
    // Create the initial population  
    // System.out.println("****DEBUG ****IN ABLA.execute()");
    
    ArrayList<String> allSourcesCurrentPop_factual =  new ArrayList<String>();
    ArrayList<String> allTargetsCurrentPop_factual =  new ArrayList<String>();
    ArrayList<String> allSourcesCurrentPop_analogical =  new ArrayList<String>();
    ArrayList<String> allTargetsCurrentPop_analogical =  new ArrayList<String>();
    ArrayList<String> allSourcesCurrentPop_creative =  new ArrayList<String>();
    ArrayList<String> allTargetsCurrentPop_creative =  new ArrayList<String>();

    System.out.println("****DEBUG ****Initialize solutions from external archive ");

    //intialize external archive
   ArrayList<String> partialSolutionString, finalSolutionString; 
    for (int i = 0; i < externalArchiveSize; i++) {
    	Solution individual = Main.listOfExternalSolutions.get(i);
    	if (individual.isInternal == true) {
    		System.out.println("Error in the creation of external domain solutions: the variable internal_extenral domain must be set to false");
    		System.exit(-1); 
    	}
    	/*
    	//manage the case of multiple target: select the best one
      	partialSolutionString = new ArrayList<String>(); 
    	//finalSolutionString = new ArrayList<String>();
    	partialSolutionString.addAll(SolutionUtils.getSources(individual));
    	//double bestPlaus = -1; //it could be best fitness, but for the external archive we use the most plausible ones in the analogical operator
    	int maxIndex = 0;
    	int targetSize= SolutionUtils.getTargets(individual).size();
    	if (targetSize>1) {
    		partialSolutionString.add(SolutionUtils.getTargets(individual).get(Main.ran.nextInt(targetSize)));
    		/*for (int j = 0; j<SolutionUtils.getTargets(individual).size();j++) {
	    		partialSolutionString.add(SolutionUtils.getTargets(individual).get(j));
	    		Solution partialSolution = SolutionUtils.createSolutionFromString(partialSolutionString);
	    		double plaus=-1;
				try {plaus = Main.plaus.plausibilityEvaluation(partialSolution, false);} catch (Exception e) {e.printStackTrace();}
	    		if (plaus >= bestPlaus) {
	    			bestPlaus = plaus;
	    			maxIndex = j;
	    		}
	    		partialSolutionString.remove(partialSolutionString.size()-1);
	    	}
	    	// *  /  
    	}
    	
    	Solution partialSolution = SolutionUtils.createSolutionFromString(partialSolutionString);
    	//finalSolutionString.addAll(SolutionUtils.getSources(individual));
    //	finalSolutionString.add(SolutionUtils.getTargets(individual).get(maxIndex));
    	
    /*	((CauseEffectSolutionType)individual.getType()).setOperator("EXTERNAL");
    	problem_.evaluate(individual);   
        problem_.evaluateConstraints(individual); 
    	((CauseEffectSolutionType)individual.getType()).setInternal_domain(false);
    	//individual.setObjective(0, 1);
    	//individual.setObjective(1, 1);
        external_archive_analogical.add(individual);
        //individual.setLocation(i);
         * 
         */
    //	printSolution(partialSolution);
        //((CauseEffectSolutionType)partialSolution.getType()).setOperator("EXTERNAL");
    	//((CauseEffectSolutionType)partialSolution.getType()).setInternal_domain(false);
    	//partialSolution.isInternal = false;
    	//problem_.evaluate(partialSolution);   
        //problem_.evaluateConstraints(partialSolution); 
    	//individual.setObjective(0, 1);
    	//individual.setObjective(1, 1);
    	problem_.evaluate(individual);
        external_archive_analogical.add(individual);
    }
    System.out.println("\nEVA \n ");

    int generations = 1; 
    
    
    
    System.out.println("\nCreating initial factual population\n ");
    Main.solutionTypeFlag = "factual"; 
    for (int i = 0; i < populationSize_factual; i++){
    	AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
		double start = System.nanoTime(); 
		
    	Solution individual = new Solution(problem_); //RANDOM.
   //   ((CauseEffectSolutionType)individual.getType()).setOperator("INITIAL");
      individual.isInternal = true;     
      problem_.evaluate(individual);
      problem_.evaluateConstraints(individual);
      individual.setOperator("factual");

      times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
		
      //printSolution(individual);
      SolutionUtils.writeSolutionObjectivesToFile(filename_obj_factual, individual);
      SolutionUtils.writeSolutionVariablesToFile(filename_var_factual, individual);
    
     // System.out.println("Factual solution - Plausibility: "+individual.getObjective(0));
     // System.out.println("Factual solution - Novelty: "+individual.getObjective(1)+"\n");
      evaluations++;

      ArrayList<String> sources =  new ArrayList<String>();
      sources.addAll(SolutionUtils.getSources(individual));

      ArrayList<String> targets = new ArrayList<String>();
      targets.addAll(SolutionUtils.getTargets(individual));
      for(int j=0; j<sources.size();j++) if (!allSourcesCurrentPop_factual.contains(sources.get(j))) allSourcesCurrentPop_factual.add(sources.get(j)); 
      for(int j=0; j<targets.size();j++) if (!allTargetsCurrentPop_factual.contains(targets.get(j))) allTargetsCurrentPop_factual.add(targets.get(j));

      
      population_factual.add(individual);      
      individual.setLocation(i);
    }   
    //FileWriter output=null; PrintWriter pw=null;
//	try { output  = new FileWriter(filename_statistics_factual, true); pw = new PrintWriter(new BufferedWriter(output));} catch (IOException e) {e.printStackTrace();}
//	pw.println("Average Plausibility, Average Novelty, Median Plausibility, Median Novelty, Best Plausibility, Best Novelty, Worst Plausibility, Worst Novelty");
//	pw.flush();
//	pw.close();
	
    //SolutionUtils.writePopulationStatisticsToFile(filename_statistics_factual, population_factual);
    
    
    System.out.println("\nCreating initial analogical population\n ");
    Main.solutionTypeFlag = "analogical"; 
    for (int i = 0; i < populationSize_analogical; i++){
    	AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
		double start = System.nanoTime(); 
		
		Solution individual = new Solution(problem_); //RANDOM.
    //  ((CauseEffectSolutionType)individual.getType()).setOperator("INITIAL");
      individual.isInternal = true; 
      individual.setOperator("analogical");

      problem_.evaluate(individual);
      problem_.evaluateConstraints(individual); 
      times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
		
   //   printSolution(individual);
      SolutionUtils.writeSolutionObjectivesToFile(filename_obj_analogical, individual);
      SolutionUtils.writeSolutionVariablesToFile(filename_var_analogical, individual);
    //  System.out.println("Analogical solution - Plausibility: "+individual.getObjective(0));
     // System.out.println("Analogical solution - Novelty: "+individual.getObjective(1)+"\n");
      evaluations++;
      
      ArrayList<String> sources =  new ArrayList<String>();
      sources.addAll(SolutionUtils.getSources(individual));
      ArrayList<String> targets = new ArrayList<String>();
      targets.addAll(SolutionUtils.getTargets(individual));
      for(int j=0; j<sources.size();j++) if (!allSourcesCurrentPop_analogical.contains(sources.get(j))) allSourcesCurrentPop_analogical.add(sources.get(j)); 
      for(int j=0; j<targets.size();j++) if (!allTargetsCurrentPop_analogical.contains(targets.get(j))) allTargetsCurrentPop_analogical.add(targets.get(j));
      
      population_analogical.add(individual);      
      individual.setLocation(i);
    }   
  //  output=null;  pw=null;
//	try { output  = new FileWriter(filename_statistics_analogical, true); pw = new PrintWriter(new BufferedWriter(output));} catch (IOException e) {e.printStackTrace();}
//	pw.println("Average Plausibility, Average Novelty, Median Plausibility, Median Novelty, Best Plausibility, Best Novelty, Worst Plausibility, Worst Novelty");
//	pw.flush();
//	pw.close();
	
    //SolutionUtils.writePopulationStatisticsToFile(filename_statistics_analogical, population_analogical);
    
    
    
    System.out.println("\nCreating initial creative population\n ");
    Main.solutionTypeFlag = "creative"; 
    for (int i = 0; i < populationSize_creative; i++){
    	AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
		double start = System.nanoTime(); 
		
		Solution individual = new Solution(problem_); //RANDOM.
    //  ((CauseEffectSolutionType)individual.getType()).setOperator("INITIAL");
		individual.isInternal = true; 
	      individual.setOperator("creative");

		problem_.evaluate(individual);
      problem_.evaluateConstraints(individual);
      times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
 
    //  printSolution(individual);
      SolutionUtils.writeSolutionObjectivesToFile(filename_obj_creative, individual);
      SolutionUtils.writeSolutionVariablesToFile(filename_var_creative, individual);
    //  System.out.println("Creative solution - Plausibility: "+individual.getObjective(0));
     // System.out.println("Creative solution - Novelty: "+individual.getObjective(1)+"\n");
      evaluations++;
      
      ArrayList<String> sources =  new ArrayList<String>();
      sources.addAll(SolutionUtils.getSources(individual));
      ArrayList<String> targets = new ArrayList<String>();
      targets.addAll(SolutionUtils.getTargets(individual));
      for(int j=0; j<sources.size();j++) if (!allSourcesCurrentPop_creative.contains(sources.get(j))) allSourcesCurrentPop_creative.add(sources.get(j)); 
      for(int j=0; j<targets.size();j++) if (!allTargetsCurrentPop_creative.contains(targets.get(j))) allTargetsCurrentPop_creative.add(targets.get(j));
      
      population_creative.add(individual);      
      individual.setLocation(i);
    }   
  //  output=null; pw=null;
//	try { output  = new FileWriter(filename_statistics_creative, true); pw = new PrintWriter(new BufferedWriter(output));} catch (IOException e) {e.printStackTrace();}
//	pw.println("Average Plausibility, Average Novelty, Median Plausibility, Median Novelty, Best Plausibility, Best Novelty, Worst Plausibility, Worst Novelty");
//	pw.flush();
//	pw.close();
	
   // SolutionUtils.writePopulationStatisticsToFile(filename_statistics_creative, population_creative);
    
    
    
    /**** START CYCLE ***/
    offspringPopulation_factual = offspringPopulation_analogical = offspringPopulation_creative  = null;
    //population.sort(comparator) ;
    /*Solution[] firstPopulation = new Solution[populationSize_factual];
    Solution[] secondPopulation =  new Solution[populationSize_analogical]; 
    Solution[] thirdPopulation = new Solution[populationSize_creative];
    */ 
    
    Solution chosenSolution = new Solution(); 
	ArrayList<String> offspring_allSourcesCurrentPop_factual =  new ArrayList<String>();
    ArrayList<String> offspring_allTargetsCurrentPop_factual =  new ArrayList<String>();
    ArrayList<String> offspring_allSourcesCurrentPop_analogical =  new ArrayList<String>();
    ArrayList<String> offspring_allTargetsCurrentPop_analogical =  new ArrayList<String>();
    ArrayList<String> offspring_allSourcesCurrentPop_creative =  new ArrayList<String>();
    ArrayList<String> offspring_allTargetsCurrentPop_creative =  new ArrayList<String>();
    
    
    //System.out.println("Initial population created: Press Any Key To Continue...");
    //new java.util.Scanner(System.in).nextLine();

    
    SolutionSet union_temp_generations = ((SolutionSet) population_factual).union(population_analogical); 
    SolutionSet population_generation = ((SolutionSet) union_temp_generations).union(population_creative);
    try {
    SolutionUtils.writePopulationStatisticsToFile(filename_statistics_generation, population_generation);
    Main.sim.computeReferenceSetDistance_by_generation(population_generation, Main.referenceSet,"distance_RefSet_generations_"+Main.technique, Main.run, generations, true);      
	  } catch (Exception e) {e.printStackTrace();}
    generations++;
    
    
    
    /******START GENERATIONS beyond the first one *****/
    
    
    
    
    while (evaluations < maxEvaluations){    
		offspringPopulation_factual = new SolutionSet(populationSize_factual);
		offspringPopulation_analogical= new SolutionSet(populationSize_analogical);
		offspringPopulation_creative = new SolutionSet(populationSize_creative );
	
		System.out.println("\nEvaluations: "+evaluations);
		
		System.out.println("\nUsing factual abduction operator (Run "+Main.run+")...\n");
		for (int i = 0; i < (populationSize_factual); i++) { 
			AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
			double start = System.nanoTime(); 
			System.out.println("\nSolution - Factual: "+(i+1)+" out of "+populationSize_factual);


			chosenSolution = (Solution) selectionOperator_factual.execute(population_factual);
    	//	System.out.println("Chosen solution "+SolutionUtils.getNumericRepresentation(chosenSolution));
    	//	System.out.println("Chosen solution obj function values: "+chosenSolution);
    		Solution individual= (Solution)((BasicFactualAbduction)factualAbductionOperator).execute(chosenSolution, allSourcesCurrentPop_factual, allTargetsCurrentPop_factual);
    		individual.isInternal = true; 
    	      individual.setOperator("factual");

    		problem_.evaluate(individual);
    		problem_.evaluateConstraints(individual);
    		times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
    		
    		offspringPopulation_factual.add(individual);
    		evaluations++;
    
    		//printSolution(individual);
    		SolutionUtils.writeSolutionObjectivesToFile(filename_obj_factual, individual);
    	    SolutionUtils.writeSolutionVariablesToFile(filename_var_factual, individual);
    	    System.out.println("Plausibility: "+individual.getObjective(0));
    	//    System.out.println("Novelty: "+individual.getObjective(1));
    	      
    		ArrayList<String> sources =  new ArrayList<String>();
    		sources.addAll(SolutionUtils.getSources(individual));
    		ArrayList<String> targets = new ArrayList<String>();
    		targets.addAll(SolutionUtils.getTargets(individual));
    		for(int j=0; j<sources.size();j++) if (!offspring_allSourcesCurrentPop_factual.contains(sources.get(j))) offspring_allSourcesCurrentPop_factual.add(sources.get(j)); 
    		for(int j=0; j<targets.size();j++) if (!offspring_allTargetsCurrentPop_factual.contains(targets.get(j))) offspring_allTargetsCurrentPop_factual.add(targets.get(j));
    	}
		
	
		if(populationSize_factual>0) {
	    /***UPDATE BY NSGA_II ***/
	    // Create the solutionSet union of solutionSet and offSpring
	      union_factual = ((SolutionSet) population_factual).union(offspringPopulation_factual);

	      // Ranking the union
	      Ranking ranking_factual = new Ranking(union_factual);

	      int remain = populationSize_factual;
	      int index = 0;
	      SolutionSet front_factual = null;
	      population_factual.clear();

	      // Obtain the next front
	      front_factual = ranking_factual.getSubfront(index);

	      while ((remain > 0) && (remain >= front_factual.size())) {
	        //Assign crowding distance to individuals
	        distance.crowdingDistanceAssignment(front_factual, problem_.getNumberOfObjectives());
	        //Add the individuals of this front
	        for (int k = 0; k < front_factual.size(); k++) {
	          population_factual.add(front_factual.get(k));
	        } // for

	        //Decrement remain
	        remain = remain - front_factual.size();

	        //Obtain the next front
	        index++;
	        if (remain > 0) {
	        	front_factual = ranking_factual.getSubfront(index);
	        } // if        
	      } // while

	      // Remain is less than front(index).size, insert only the best one
	      if (remain > 0) {  // front contains individuals to insert                        
	        distance.crowdingDistanceAssignment(front_factual, problem_.getNumberOfObjectives());
	        front_factual.sort(new CrowdingComparator());
	        for (int k = 0; k < remain; k++) {
	          population_factual.add(front_factual.get(k));
	        } // for

	        remain = 0;
	      } // if                               
	//	SolutionUtils.writePopulationStatisticsToFile(filename_statistics_factual, offspringPopulation_factual);
		}
		System.out.println("\nUsing analogical abduction operator (Run "+Main.run+")...\n");
    	for (int i = 0; i < (populationSize_analogical); i++) { 
    		AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
			double start = System.nanoTime(); 
			
			System.out.println("\nSolution - Analogical: "+(i+1)+" out of "+populationSize_analogical);
    		chosenSolution = (Solution) selectionOperator_analogical.execute(external_archive_analogical);
    		Solution individual = (Solution)((BasicAnalogicalAbduction)analogicalAbductionOperator).execute(chosenSolution, population_analogical, allSourcesCurrentPop_analogical, allTargetsCurrentPop_analogical);
    		
    		individual.isInternal = true;
    	      individual.setOperator("analogical");

    		problem_.evaluate(individual);
    		problem_.evaluateConstraints(individual); 
    		times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
    		
    		offspringPopulation_analogical.add(individual);
    		evaluations++;
    		
    		
    		//printSolution(individual);
    		SolutionUtils.writeSolutionObjectivesToFile(filename_obj_analogical, individual);
      		SolutionUtils.writeSolutionVariablesToFile(filename_var_analogical, individual);
    		System.out.println("Plausibility: "+individual.getObjective(0));
    	  //  System.out.println("Novelty: "+individual.getObjective(1));
    	      
    		ArrayList<String> sources =  new ArrayList<String>();
    		sources.addAll(SolutionUtils.getSources(individual));
    		ArrayList<String> targets = new ArrayList<String>();
    		targets.addAll(SolutionUtils.getTargets(individual));
    		for(int j=0; j<sources.size();j++) if (!offspring_allSourcesCurrentPop_analogical.contains(sources.get(j))) offspring_allSourcesCurrentPop_analogical.add(sources.get(j)); 
    		for(int j=0; j<targets.size();j++) if (!offspring_allTargetsCurrentPop_analogical.contains(targets.get(j))) offspring_allTargetsCurrentPop_analogical.add(targets.get(j));
	    }
	  //  SolutionUtils.writePopulationStatisticsToFile(filename_statistics_analogical, offspringPopulation_analogical);
	 
    	
    	 if(populationSize_analogical>0) { 
	    /***UPDATE BY NSGA_II ***/
	    // Create the solutionSet union of solutionSet and offSpring
	      union_analogical = ((SolutionSet) population_analogical).union(offspringPopulation_analogical);

	      // Ranking the union
	      Ranking ranking_analogical = new Ranking(union_analogical);

	      int remain = populationSize_analogical;
	      int index = 0;
	      SolutionSet front_analogical= null;
	      population_analogical.clear();

	      // Obtain the next front
	      front_analogical = ranking_analogical.getSubfront(index);

	      while ((remain > 0) && (remain >= front_analogical.size())) {
	        //Assign crowding distance to individuals
	        distance.crowdingDistanceAssignment(front_analogical, problem_.getNumberOfObjectives());
	        //Add the individuals of this front
	        for (int k = 0; k < front_analogical.size(); k++) {
	          population_analogical.add(front_analogical.get(k));
	        } // for

	        //Decrement remain
	        remain = remain - front_analogical.size();

	        //Obtain the next front
	        index++;
	        if (remain > 0) {
	        	front_analogical = ranking_analogical.getSubfront(index);
	        } // if        
	      } // while

	      // Remain is less than front(index).size, insert only the best one
	      if (remain > 0) {  // front contains individuals to insert                        
	        distance.crowdingDistanceAssignment(front_analogical, problem_.getNumberOfObjectives());
	        front_analogical.sort(new CrowdingComparator());
	        for (int k = 0; k < remain; k++) {
	          population_analogical.add(front_analogical.get(k));
	        } // for

	        remain = 0;
	      } // if      
	      
    	 }
    	 
    	System.out.println("\nUsing hypothetical-cause abduction operator (Run "+Main.run+")...\n");
    	for (int i = 0; i < (populationSize_creative); i++) {
    		AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
			double start = System.nanoTime(); 
			
			System.out.println("\nSolution - Creative: "+(i+1)+" out of "+populationSize_creative);
    		chosenSolution = (Solution)selectionOperator_creative.execute(population_creative);
    //		System.out.println("Chosen solution creative "+SolutionUtils.getNumericRepresentation(chosenSolution));
    // 		System.out.println("Chosen solution obj function values: "+chosenSolution);
    		Solution individual = (Solution)((BasicCreativeAbduction)creativeAbductionOperator).execute(chosenSolution, allSourcesCurrentPop_creative, allTargetsCurrentPop_creative);
    		
    		individual.isInternal = true; 
    	      individual.setOperator("creative");

    		offspringPopulation_creative.add(individual);
    		problem_.evaluate(individual);
    		problem_.evaluateConstraints(individual); 
    		times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
    		evaluations++;
    		
    	//	printSolution(individual);
    		SolutionUtils.writeSolutionObjectivesToFile(filename_obj_creative, individual);
    	    SolutionUtils.writeSolutionVariablesToFile(filename_var_creative, individual);
    		System.out.println("Plausibility: "+individual.getObjective(0));
    	 //   System.out.println("Novelty: "+individual.getObjective(1));
    	      
    		ArrayList<String> sources =  new ArrayList<String>();
    		sources.addAll(SolutionUtils.getSources(individual));
    		ArrayList<String> targets = new ArrayList<String>();
    		targets.addAll(SolutionUtils.getTargets(individual));
    		for(int j=0; j<sources.size();j++) if (!offspring_allSourcesCurrentPop_creative.contains(sources.get(j))) offspring_allSourcesCurrentPop_creative.add(sources.get(j)); 
    		for(int j=0; j<targets.size();j++) if (!offspring_allTargetsCurrentPop_creative.contains(targets.get(j))) offspring_allTargetsCurrentPop_creative.add(targets.get(j));
    	}
    //	SolutionUtils.writePopulationStatisticsToFile(filename_statistics_creative, offspringPopulation_creative);
    if(populationSize_creative>0) {
	    
	    /***UPDATE BY NSGA_II ***/
	    // Create the solutionSet union of solutionSet and offSpring
	      union_creative = ((SolutionSet) population_creative).union(offspringPopulation_creative);

	      // Ranking the union
	      Ranking ranking_creative = new Ranking(union_creative);

	      int remain = populationSize_creative;
	      int index = 0;
	      SolutionSet front_creative = null;
	      population_creative.clear();

	      // Obtain the next front
	      front_creative = ranking_creative.getSubfront(index);

	      while ((remain > 0) && (remain >= front_creative.size())) {
	        //Assign crowding distance to individuals
	        distance.crowdingDistanceAssignment(front_creative, problem_.getNumberOfObjectives());
	        //Add the individuals of this front
	        for (int k = 0; k < front_creative.size(); k++) {
	          population_creative.add(front_creative.get(k));
	        } // for

	        //Decrement remain
	        remain = remain - front_creative.size();

	        //Obtain the next front
	        index++;
	        if (remain > 0) {
	        	front_creative = ranking_creative.getSubfront(index);
	        } // if        
	      } // while

	      // Remain is less than front(index).size, insert only the best one
	      if (remain > 0) {  // front contains individuals to insert                        
	        distance.crowdingDistanceAssignment(front_creative, problem_.getNumberOfObjectives());
	        front_creative.sort(new CrowdingComparator());
	        for (int k = 0; k < remain; k++) {
	          population_creative.add(front_creative.get(k));
	        } // for

	        remain = 0;
	      } // if      
    	
    }
  	      allSourcesCurrentPop_factual.clear();
	      allTargetsCurrentPop_factual.clear();
	      allSourcesCurrentPop_analogical.clear();
	      allTargetsCurrentPop_analogical.clear();
	      allSourcesCurrentPop_creative.clear();
	      allTargetsCurrentPop_creative.clear();
	     	  
	      for (int j = 0; j < population_factual.size(); j++) {
	    	  ArrayList<String> sources =  new ArrayList<String>();
	    	  sources.addAll(SolutionUtils.getSources(population_factual.get(j)));
	    	  ArrayList<String> targets = new ArrayList<String>();
	    	  targets.addAll(SolutionUtils.getTargets(population_factual.get(j)));
	    	  for(int i=0; i<sources.size();i++) if (!allSourcesCurrentPop_factual.contains(sources.get(i))) allSourcesCurrentPop_factual.add(sources.get(i)); 
	    	  for(int i=0; i<targets.size();i++) if (!allTargetsCurrentPop_factual.contains(targets.get(i))) allTargetsCurrentPop_factual.add(targets.get(i));
	      }
	      for (int j = 0; j < population_analogical.size(); j++) {
	    	  ArrayList<String> sources =  new ArrayList<String>();
	    	  sources.addAll(SolutionUtils.getSources(population_analogical.get(j)));
	    	  ArrayList<String> targets = new ArrayList<String>();
	    	  targets.addAll(SolutionUtils.getTargets(population_analogical.get(j)));
	    	  for(int i=0; i<sources.size();i++) if (!allSourcesCurrentPop_analogical.contains(sources.get(i))) allSourcesCurrentPop_analogical.add(sources.get(i)); 
	    	  for(int i=0; i<targets.size();i++) if (!allTargetsCurrentPop_analogical.contains(targets.get(i))) allTargetsCurrentPop_analogical.add(targets.get(i));
	      }
	      
	      for (int j = 0; j < population_creative.size(); j++) {
	    	  ArrayList<String> sources =  new ArrayList<String>();
	    	  sources.addAll(SolutionUtils.getSources(population_creative.get(j)));
	    	  ArrayList<String> targets = new ArrayList<String>();
	    	  targets.addAll(SolutionUtils.getTargets(population_creative.get(j)));
	    	  for(int i=0; i<sources.size();i++) if (!allSourcesCurrentPop_creative.contains(sources.get(i))) allSourcesCurrentPop_creative.add(sources.get(i)); 
	    	  for(int i=0; i<targets.size();i++) if (!allTargetsCurrentPop_creative.contains(targets.get(i))) allTargetsCurrentPop_creative.add(targets.get(i));
	      }
	      
			/*
			 * for (int i = 0; i < offspring_allSourcesCurrentPop_factual.size(); i++)
			 * allSourcesCurrentPop_factual.add(offspring_allSourcesCurrentPop_factual.get(i
			 * )) ;
			 * 
			 * for (int i = 0; i < offspring_allSourcesCurrentPop_analogical.size(); i++)
			 * allSourcesCurrentPop_analogical.add(offspring_allSourcesCurrentPop_analogical
			 * .get(i)) ;
			 * 
			 * for (int i = 0; i < offspring_allSourcesCurrentPop_creative.size(); i++)
			 * allSourcesCurrentPop_creative.add(offspring_allSourcesCurrentPop_creative.get
			 * (i)) ;
			 * 
			 * for (int i = 0; i < offspring_allTargetsCurrentPop_factual.size(); i++)
			 * allTargetsCurrentPop_factual.add(offspring_allTargetsCurrentPop_factual.get(i
			 * )) ;
			 * 
			 * for (int i = 0; i < offspring_allTargetsCurrentPop_analogical.size(); i++)
			 * allTargetsCurrentPop_analogical.add(offspring_allTargetsCurrentPop_analogical
			 * .get(i)) ;
			 * 
			 * for (int i = 0; i < offspring_allTargetsCurrentPop_creative.size(); i++)
			 * allTargetsCurrentPop_creative.add(offspring_allTargetsCurrentPop_creative.get
			 * (i)) ;
			 */
	      
	      offspring_allSourcesCurrentPop_factual.clear();
	      offspring_allSourcesCurrentPop_analogical.clear();
	      offspring_allSourcesCurrentPop_creative.clear();
	      offspring_allTargetsCurrentPop_factual.clear();
		  offspring_allTargetsCurrentPop_analogical.clear(); 
	      offspring_allTargetsCurrentPop_creative.clear();
	     
	      // Print Solutions for this generation.  
	      union_temp_generations = null; 
	      union_temp_generations = ((SolutionSet) population_factual).union(population_analogical); 
	      population_generation = null;
	      population_generation =((SolutionSet) union_temp_generations).union(population_creative);
	      try {
	      SolutionUtils.writePopulationStatisticsToFile(filename_statistics_generation, population_generation);
	      Main.sim.computeReferenceSetDistance_by_generation(population_generation, Main.referenceSet,"distance_RefSet_generations_"+Main.technique, Main.run, generations, true);      
		  } catch (Exception e) {e.printStackTrace();}
	      generations++;  
      //population.sort(comparator) ;
    } // END While
     
    SolutionSet union_temp = ((SolutionSet) population_factual).union(population_analogical); 
    SolutionSet population = ((SolutionSet) union_temp).union(population_creative);
   
    Ranking ranking = new Ranking(population);//Useless. Ranking has already been done before
    //return ranking.getSubfront(0);  // just the first one
    return population;
  }
 

private void printSolution(Solution indiv) {
	
	String line="";
	String line_names="";
	ArrayList<String> sources = new ArrayList<String>();
	ArrayList<String> targets= new ArrayList<String>();
	sources = ((Sources)(indiv.getDecisionVariables()[0]))._sourcesList;
	targets = ((Targets)(indiv.getDecisionVariables()[1]))._targetsList;
	for(int j=0; j< sources.size();j++){
		line = line + sources.get(j) + " ";
		line_names = line_names + Main.getNameFromKey(sources.get(j)) + " ";
	}//end for sources
		
	for(int j=0; j< targets.size();j++){
		line = line + targets.get(j) + " "; //targets
		line_names = line_names + Main.getNameFromKey(targets.get(j)) + " ";
		}
	System.out.println("Solution, symbolic representation "+line_names);
	System.out.println("Solution, numerical representation "+line);
}


}
