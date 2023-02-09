package causalOptimization;

import util.SolutionUtils;
import jmetal.core.*;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.NonDominatedSolutionList;
import jmetal.util.Ranking;
import jmetal.util.comparators.EpsilonDominanceComparator;
import jmetal.util.comparators.ObjectiveComparator;
import java.util.Comparator;

import jmetal.operators.creativeAbduction.*; 
import jmetal.operators.factualAbduction.*;
import jmetal.operators.analogicalAbduction.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;

//import org.apache.commons.lang3.RandomUtils;
import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;

public class EVAWeighted extends Algorithm{
  
	public EVAWeighted(Problem problem){
		super (problem);
   }
	
	
	

  /** Execute the algorithm 
   * @throws JMException */
  public SolutionSet execute() throws JMException, ClassNotFoundException {

	  //Init the parameters
	  
	  double learningRate=Main.plaus.quartiles[1]/Main.MAX_sources;//0.5;
	  if(learningRate>1) //it can happen when multiple values per variable are admitted (the maximum (and the upper quartile) can be greater than MaxSources)
		  learningRate=1;
	  
	  
    int populationSize_factual, populationSize_analogical, populationSize_creative, maxEvaluations, evaluations, maxEvaluationPerCycle, externalArchiveSize;
    Operator factualAbductionOperator, analogicalAbductionOperator, creativeAbductionOperator, selectionOperator_factual, selectionOperator_analogical, selectionOperator_creative;
    SolutionSet population_factual, population_analogical, population_creative, external_archive_analogical ;
    SolutionSet offspringPopulation;

    
    
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
    int populationSize = populationSize_factual+populationSize_analogical+populationSize_creative;

  Comparator comparator = new ObjectiveComparator(0) ; // Single objective comparator
    
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
    
    external_archive_analogical = new SolutionSet(externalArchiveSize); 
    
    evaluations        = 0;                        
	
	
    // Create the initial population  
    // System.out.println("****DEBUG ****IN ABLA.execute()");
    
    ArrayList<String> allSourcesCurrentPop =  new ArrayList<String>();
    ArrayList<String> allTargetsCurrentPop =  new ArrayList<String>();
    
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
      times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
		
      individual.setOperator("factual");
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
      for(int j=0; j<sources.size();j++) if (!allSourcesCurrentPop.contains(sources.get(j))) allSourcesCurrentPop.add(sources.get(j)); 
      for(int j=0; j<targets.size();j++) if (!allTargetsCurrentPop.contains(targets.get(j))) allTargetsCurrentPop.add(targets.get(j));
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
      problem_.evaluate(individual);
      problem_.evaluateConstraints(individual); 
      times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
		
      individual.setOperator("analogical");
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
      for(int j=0; j<sources.size();j++) if (!allSourcesCurrentPop.contains(sources.get(j))) allSourcesCurrentPop.add(sources.get(j)); 
      for(int j=0; j<targets.size();j++) if (!allTargetsCurrentPop.contains(targets.get(j))) allTargetsCurrentPop.add(targets.get(j));
      
      population_analogical.add(individual);      
      individual.setLocation(i);
      }   
  //  output=null;  pw=null;
//	try { output  = new FileWriter(filename_statistics_analogical, true); pw = new PrintWriter(new BufferedWriter(output));} catch (IOException e) {e.printStackTrace();}
//	pw.println("Average Plausibility, Average Novelty, Median Plausibility, Median Novelty, Best Plausibility, Best Novelty, Worst Plausibility, Worst Novelty");
//	pw.flush();
//	pw.close();
	
    //SolutionUtils.writePopulationStatisticsToFile(filename_statistics_analogical, population_analogical);
    
    
    
    System.out.println("\nCreating initial hypothetical-cause population\n ");
    Main.solutionTypeFlag = "creative"; 
    for (int i = 0; i < populationSize_creative; i++){
    	AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
		double start = System.nanoTime(); 
		
		Solution individual = new Solution(problem_); //RANDOM.
    //  ((CauseEffectSolutionType)individual.getType()).setOperator("INITIAL");
		individual.isInternal = true; 
		problem_.evaluate(individual);
      problem_.evaluateConstraints(individual);
      times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
 
      individual.setOperator("creative");
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
      for(int j=0; j<sources.size();j++) if (!allSourcesCurrentPop.contains(sources.get(j))) allSourcesCurrentPop.add(sources.get(j)); 
      for(int j=0; j<targets.size();j++) if (!allTargetsCurrentPop.contains(targets.get(j))) allTargetsCurrentPop.add(targets.get(j));
      
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
	
    //System.out.println("Initial population created: Press Any Key To Continue...");
    //new java.util.Scanner(System.in).nextLine();

    SolutionSet union_temp = ((SolutionSet) population_factual).union(population_analogical);
    SolutionSet population = ((SolutionSet) union_temp).union(population_creative);

    SolutionSet population_generation = ((SolutionSet) union_temp).union(population_creative);
    try {
    SolutionUtils.writePopulationStatisticsToFile(filename_statistics_generation, population_generation);
    Main.sim.computeReferenceSetDistance_by_generation(population_generation, Main.referenceSet,"distance_RefSet_generations_"+Main.technique, Main.run, generations, true);      
	  } catch (Exception e) {e.printStackTrace();}
    generations++;
    
    
    offspringPopulation =  new SolutionSet(populationSize_factual+populationSize_analogical+populationSize_creative);
    Solution chosenSolution = new Solution(); 
	ArrayList<String> offspring_allSourcesCurrentPop=  new ArrayList<String>();
    ArrayList<String> offspring_allTargetsCurrentPop=  new ArrayList<String>();

    /******START GENERATIONS beyond the first one *****/
    
    
    double weightFactual, weightAnalogical, weightCreative; 
    weightFactual    =((double)populationSize_factual)/populationSize;
    weightAnalogical =((double)populationSize_analogical)/populationSize;
    weightCreative = ((double)populationSize_creative)/populationSize;
    
    System.out.println("weightFactual "+weightFactual);
    System.out.println("weightAnalogical "+weightAnalogical);
    System.out.println("weightCreative "+weightCreative);

    
    while (evaluations < maxEvaluations){    
		
		System.out.println("\nEvaluations: "+evaluations);
		System.out.println("\n populationSize: "+populationSize);
		for (int i = 0; i < (populationSize ); i++) {

			double weight = Main.ran.nextDouble(); 

			if (weight <= weightFactual) {

			System.out.println("\nUsing factual abduction operator (Run "+Main.run+")...\n");
			AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
			double start = System.nanoTime(); 
			System.out.println("\nSolution - Factual: "+(i+1));
			chosenSolution = (Solution) selectionOperator_factual.execute(population);
    		Solution individual= (Solution)((BasicFactualAbduction)factualAbductionOperator).execute(chosenSolution, allSourcesCurrentPop, allTargetsCurrentPop);
    		individual.isInternal = true; 
    		individual.setOperator("factual");
    		problem_.evaluate(individual);
    		problem_.evaluateConstraints(individual);
    		times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 		
    		offspringPopulation.add(individual);
    		evaluations++;
    		System.out.println("\nEvaluations: "+evaluations);
    		ArrayList<String> sources =  new ArrayList<String>();
    		sources.addAll(SolutionUtils.getSources(individual));
    		ArrayList<String> targets = new ArrayList<String>();
    		targets.addAll(SolutionUtils.getTargets(individual));
    		for(int j=0; j<sources.size();j++) if (!offspring_allSourcesCurrentPop.contains(sources.get(j))) offspring_allSourcesCurrentPop.add(sources.get(j)); 
    		for(int j=0; j<targets.size();j++) if (!offspring_allTargetsCurrentPop.contains(targets.get(j))) offspring_allTargetsCurrentPop.add(targets.get(j));
    	}
		
		if (weight > weightFactual && weight <= (weightFactual+weightAnalogical)) {
			System.out.println("\nUsing analogical abduction operator (Run "+Main.run+")...\n");
			AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
			double start = System.nanoTime(); 
			System.out.println("\nSolution - Analogical: "+(i+1));
    		chosenSolution = (Solution) selectionOperator_analogical.execute(external_archive_analogical);
    		Solution individual = (Solution)((BasicAnalogicalAbduction)analogicalAbductionOperator).execute(chosenSolution, population, allSourcesCurrentPop, allTargetsCurrentPop);
    		
    		individual.isInternal = true; 
    		individual.setOperator("analogical");
    		problem_.evaluate(individual);
    		problem_.evaluateConstraints(individual); 
    		times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
    		
    		offspringPopulation.add(individual);
    		evaluations++;
    		System.out.println("\nEvaluations: "+evaluations);
    		ArrayList<String> sources =  new ArrayList<String>();
    		sources.addAll(SolutionUtils.getSources(individual));
    		ArrayList<String> targets = new ArrayList<String>();
    		targets.addAll(SolutionUtils.getTargets(individual));
    		for(int j=0; j<sources.size();j++) if (!offspring_allSourcesCurrentPop.contains(sources.get(j))) offspring_allSourcesCurrentPop.add(sources.get(j)); 
    		for(int j=0; j<targets.size();j++) if (!offspring_allTargetsCurrentPop.contains(targets.get(j))) offspring_allTargetsCurrentPop.add(targets.get(j));
	    }
	  //  SolutionUtils.writePopulationStatisticsToFile(filename_statistics_analogical, offspringPopulation_analogical);
	 
		if (weight > (weightFactual+weightAnalogical) && weight <= (weightFactual+weightAnalogical+weightCreative)) { //weightCreative is derived as 1 minus the other two weights. Could be removed 
			System.out.println("\nUsing creative abduction operator ...\n");
			AdditionalInfo ai = new AdditionalInfo(); ArrayList<Double> times = new ArrayList<Double>();
			double start = System.nanoTime(); 
			System.out.println("\nSolution - Creative: "+(i+1));
    		chosenSolution = (Solution)selectionOperator_creative.execute(population);
    		Solution individual = (Solution)((BasicCreativeAbduction)creativeAbductionOperator).execute(chosenSolution, allSourcesCurrentPop, allTargetsCurrentPop);
    		individual.isInternal = true; 
    		individual.setOperator("creative");
    		offspringPopulation.add(individual);
    		problem_.evaluate(individual);
    		problem_.evaluateConstraints(individual); 
    		times.add(System.nanoTime()-start); ai.setTimeVector(times); individual.setAdditionalInfo(ai); 
    		evaluations++;
			System.out.println("\nEvaluations: "+evaluations);

    		ArrayList<String> sources =  new ArrayList<String>();
    		sources.addAll(SolutionUtils.getSources(individual));
    		ArrayList<String> targets = new ArrayList<String>();
    		targets.addAll(SolutionUtils.getTargets(individual));
    		for(int j=0; j<sources.size();j++) if (!offspring_allSourcesCurrentPop.contains(sources.get(j))) offspring_allSourcesCurrentPop.add(sources.get(j)); 
    		for(int j=0; j<targets.size();j++) if (!offspring_allTargetsCurrentPop.contains(targets.get(j))) offspring_allTargetsCurrentPop.add(targets.get(j));
    	}
    //	SolutionUtils.writePopulationStatisticsToFile(filename_statistics_creative, offspringPopulation_creative);
   
		}
	    /***UPDATE BY NSGA_II ***/
	    // Create the solutionSet union of solutionSet and offSpring
    	SolutionSet union = ((SolutionSet) population).union(offspringPopulation);
    	
    	/*for (int i=0;i<union.size();i++) {
    		System.out.println(" union "+union.get(i).getOperator());
    	}*/
	      // Ranking the union
	      Ranking ranking = new Ranking(union);

	      union.sort(comparator);
	      
	      
	      for (int ind = 0; ind <union.size(); ind++) {
	    	  System.out.println("solutions ordered "+union.get(ind).getObjective(0)); 
	      }
	      
	      population.clear();
	      for (int k = 0; k < populationSize; k++) {
   			population.add(union.get(k));
   			}
	      offspringPopulation.clear();
	
		    int countFactual=0;
		    int countAnalogical=0;
		    int countCreative=0;
	      //update the weights
	      for (int k = 0; k < population.size(); k++) {
	    	 // System.out.println("operator "+population.get(k).getOperator());
	    	  if (population.get(k).getOperator().equals("factual")) {
	    		  countFactual++;	    		  
	    		  
	    	  }
	    	  if (population.get(k).getOperator().equals("analogical")) {
	    		  countAnalogical++;
	    		  
	    		  
	    	  }
	    	  if (population.get(k).getOperator().equals("creative")) {
	    		  countCreative++; //for double checking, can be removed
	    		  
	    	  }
	      }
	      
	      System.out.println("\n population size "+population.size());

	      System.out.println("\ncountFactual "+countFactual);
	      System.out.println("countAnalogical "+countAnalogical);
	      System.out.println("countCreative "+countCreative);
	      
	      //guards:
	      int normalization = population.size();
	      if (countFactual<3){
	    	  normalization  = normalization  + (3-countFactual);
	    	  countFactual=3;
	      }
	      if (countAnalogical<3) {
	    	  normalization  = normalization  + (3-countAnalogical);
	    	  countAnalogical=3; 
	      }
	      if (countCreative<3) {
	    	  normalization  = normalization  + (3-countCreative);
	    	  countCreative=3;
	      }
	      
	      if (normalization != (countFactual+countAnalogical+countCreative)) {
	    	  System.out.println("Exit, normalization condition broken ");
	    	  System.exit(0);
	      }
	      
	      System.out.println("\ncountFactual "+countFactual);
	      System.out.println("countAnalogical "+countAnalogical);
	      System.out.println("countCreative "+countCreative);
	      
	      System.out.println("normalization "+normalization);
	     
	      weightFactual= (((double)countFactual)/normalization)*learningRate + (1-learningRate)*weightFactual; 
	      weightAnalogical =((double)countAnalogical)/normalization *learningRate + (1-learningRate)*weightAnalogical;
	      //creative is 1 minus the previous two
	      weightCreative  =((double)countCreative)/normalization *learningRate + (1-learningRate)*weightCreative  ;
	      
	      
	      
	      System.out.println("\n weightFactual "+weightFactual);
	      System.out.println("weightAnalogical "+weightAnalogical);
	      System.out.println("weightCreative "+weightCreative);

	      
  	      allSourcesCurrentPop.clear();
	      allTargetsCurrentPop.clear();
	      	      
	      for (int j = 0; j < population.size(); j++) {
	    	  ArrayList<String> sources =  new ArrayList<String>();
	    	  sources.addAll(SolutionUtils.getSources(population.get(j)));
	    	  ArrayList<String> targets = new ArrayList<String>();
	    	  targets.addAll(SolutionUtils.getTargets(population.get(j)));
	    	  for(int i=0; i<sources.size();i++) if (!allSourcesCurrentPop.contains(sources.get(i))) allSourcesCurrentPop.add(sources.get(i)); 
	    	  for(int i=0; i<targets.size();i++) if (!allTargetsCurrentPop.contains(targets.get(i))) allTargetsCurrentPop.add(targets.get(i));
	      }
			/*
			 * for (int i = 0; i < offspring_allSourcesCurrentPop.size(); i++)
			 * allSourcesCurrentPop.add(offspring_allSourcesCurrentPop.get(i)) ;
			 * 
			 * for (int i = 0; i < offspring_allTargetsCurrentPop.size(); i++)
			 * allTargetsCurrentPop.add(offspring_allTargetsCurrentPop.get(i)) ;
			 */
	        
	      offspring_allSourcesCurrentPop.clear();
	      offspring_allTargetsCurrentPop.clear();
		  
	     
	      // Print Solutions for this generation.  
	      population_generation.clear(); 
	      population_generation  = (SolutionSet)population_generation.union(population); 
	      try {
	      SolutionUtils.writePopulationStatisticsToFile(filename_statistics_generation, population_generation);
	      Main.sim.computeReferenceSetDistance_by_generation(population_generation, Main.referenceSet,"distance_RefSet_generations_"+Main.technique, Main.run, generations, true);      
		  } catch (Exception e) {e.printStackTrace();}
	      generations++;  
      //population.sort(comparator) ;
    } // END While
     
    
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
