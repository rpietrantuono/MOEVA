package jmetal.operators.analogicalAbduction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import org.apache.commons.lang3.RandomUtils;

import causalOptimization.Main;
import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import solutionType.CauseEffectSolutionType;
import util.SolutionUtils;

public class BasicAnalogicalAbduction extends AnalogicalAbduction  {

	private double analogicalNoveltyIndex_ ;
	private int numberOfMaxAttempts_;

	public BasicAnalogicalAbduction(HashMap parameters)  {
		super(parameters);
		analogicalNoveltyIndex_ = 0.5; 
		if (parameters.get("analogicalNoveltyIndex") != null)
			analogicalNoveltyIndex_ = (double)parameters.get("analogicalNoveltyIndex") ;
		numberOfMaxAttempts_ = 5;  //This should be made customisable 
		
	}

	@Override
	public Object execute(Object object) throws JMException {
		return null;
	}

	/**
	 * @param chosenSol: The chosen solution by the selection operator in the analogy starting domain (SOURCE)
	 * @param population: current population
	 * @param allSources: The list of possible sources taken from the current population
	 * @param allTargets:  The list of possible targets taken from the current population
	 * @return generatedSolution: a new Solution
	 * @throws JMException
	 */
	public Object execute(Object chosenSol, Object population, ArrayList<String> allSources, ArrayList<String> allTargets) throws JMException {

		
		long start_time = System.currentTimeMillis();
		
		Solution generatedSolution = new Solution();
	
	
		//System.out.println("Analogical nov index: "+analogicalNoveltyIndex_);
		
		//ANALOGY: from SORUCE domain (from another (internal/external) archive) to TARGET domain
		// In this case, we don't work on an existing solution in the current pop; but we build it from scratch, 
		//considering the solution in the SOURCE domain
		
		Solution chosenSolution_SourceDomain  = (Solution)chosenSol;
		//String chosenTarget_SourceDomain = SolutionUtils.getTargets(chosenSolution_SourceDomain).get(Main.ran.nextInt((SolutionUtils.getTargets(chosenSolution_SourceDomain)).size()));
		ArrayList<String> sources = SolutionUtils.getSources(chosenSolution_SourceDomain);
		//String chosenTarget = selectRandomFromCurrentPop(allTargets);
		String chosenTarget =null; 
//		if (RandomUtils.nextDouble(0,1) < analogicalNoveltyIndex_)
			chosenTarget = selectElementBinaryTournament(allTargets);
			//chosenTarget = selectBest(allTargets);
//		else 
//			chosenTarget = selectElementBinaryTournament(causeEffect.Main.possibleTargets);

		//	System.out.println("chosenTarget "+chosenTarget);
		
		
		//build relations of chosen solution : characterize the SOURCE solution. 
		//1. Tell how many solution per group. 
		ArrayList<Integer> numberOfSourcesPerGroup_SourceDomain= new ArrayList<Integer>(); // for the sources of the chosen solution (in the SOURCE domain)
		ArrayList<ArrayList<String>> sourcesPerGroup_SourceDomain = new ArrayList<ArrayList<String>>();
		ArrayList<Integer> maxKDegreePerGroup_SourceDomain= new ArrayList<Integer>();
		//2. vectors to map sources in the target domain
		ArrayList<ArrayList<String>> sourcesPerGroup_TargetDomain = new ArrayList<ArrayList<String>>();
		ArrayList<Integer> maxKDegreePerGroup_TargetDomain= new ArrayList<Integer>();
		ArrayList<ArrayList<String>> sourcesPerGroupCurrentPop_TargetDomain = new ArrayList<ArrayList<String>>();
		ArrayList<Integer> numberOfSourcesPerGroupCurrentPop_TargetDomain = new ArrayList<Integer>();//for all the sources of the population (in the TARGET domain)
		
		
		// initialize
		for (int j=0; j< Main.groupLimitNumber.size(); j++) { 
			numberOfSourcesPerGroup_SourceDomain.add(0);  
			numberOfSourcesPerGroupCurrentPop_TargetDomain.add(0);
			maxKDegreePerGroup_SourceDomain.add(0);
			maxKDegreePerGroup_TargetDomain.add(0);
			sourcesPerGroup_SourceDomain.add(new ArrayList<String>());
			sourcesPerGroupCurrentPop_TargetDomain.add(new ArrayList<String>());
			sourcesPerGroup_TargetDomain.add(new ArrayList<String>());
		}
		int previousLimit; 
		
		// SOURCE DOMAIN 

		//LIMIT THE NUMBER TO  MIN_SOURCES, MAX_SOURCES
		//if (sources.size()>Main.MAX_sources) {
			//for (int i=0 ;i<(sources.size()- Main.MAX_sources); i++) 
			//	sources.remove(RandomUtils.nextInt(0, sources.size() ));
		while (sources.size()>Main.MAX_sources) 
			sources.remove(Main.ran.nextInt(sources.size()));
			//}
		//if (sources.size()<Main.MIN_sources) {
			//for (int i=0 ;i<(Main.MIN_sources - sources.size()); i++) 
			//	sources.add(allSources.get(RandomUtils.nextInt(0, allSources.size())));
		while (sources.size()<Main.MIN_sources) 
			sources.add(allSources.get(Main.ran.nextInt(allSources.size())));
		//}
		
		
		Collections.sort(sources, Comparator.comparing(Integer::valueOf));
		
	//	System.out.println("\n[***DEBUG****] Size of allSources of solution chosen in the SOURCE domain: "+sources.size());
	//	System.out.println("\n[***DEBUG****] allSources of solution chosen in the SOURCE domain sorted: "+sources);
		for (int i=0; i<sources.size();i++) { 
		  previousLimit = 0; 
		  for (int j=0; j<Main.groupLimitNumber.size(); j++) { 
			  if (Integer.parseInt(sources.get(i)) <= Main.groupLimitNumber.get(j) && Integer.parseInt(sources.get(i)) > previousLimit) { 
				  numberOfSourcesPerGroup_SourceDomain.set(j, numberOfSourcesPerGroup_SourceDomain.get(j) + 1); 
				  sourcesPerGroup_SourceDomain.get(j).add(sources.get(i));
				  break; 
				  } 
			  previousLimit = Main.groupLimitNumber.get(j); 
			  } 
		  }
		
			
		//print
/*		for (int i=0; i< Main.groupLimitNumber.size(); i++) {
			System.out.println("numberOfSourcesPerGroup_SourceDomain "+numberOfSourcesPerGroup_SourceDomain.get(i));
			System.out.println("sourcesPerGroup_SourceDomain "+sourcesPerGroup_SourceDomain.get(i));
		}
*/		
	  
	  // TARGET DOMAIN
		Collections.sort(allSources, Comparator.comparing(Integer::valueOf));
	//	System.out.println("\n[***DEBUG****] allSources in the TARGET domain (current population) sorted: "+allSources);
		for (int i=0; i<allSources.size();i++) {
			previousLimit = 0; 
			for (int j=0; j< Main.groupLimitNumber.size(); j++) {
				if (Integer.parseInt(allSources.get(i)) <= Main.groupLimitNumber.get(j) && Integer.parseInt(allSources.get(i)) > previousLimit) {
					numberOfSourcesPerGroupCurrentPop_TargetDomain.set(j, numberOfSourcesPerGroupCurrentPop_TargetDomain.get(j) + 1);
					sourcesPerGroupCurrentPop_TargetDomain.get(j).add(allSources.get(i)); 
					break;
				}
				previousLimit = Main.groupLimitNumber.get(j); 
				}
			}
		
		//print
	/*	for (int i=0; i< Main.groupLimitNumber.size(); i++) {
			System.out.println("numberOfSourcesPerGroupCurrentPop_TargetDomain "+numberOfSourcesPerGroupCurrentPop_TargetDomain.get(i));
			System.out.println("sourcesPerGroupCurrentPop_TargetDomain "+sourcesPerGroupCurrentPop_TargetDomain.get(i));
		}
	*/	
		
		
	
		
		/* Build solution */
		ArrayList<String> solution_string = new ArrayList<String>();
		solution_string.add(null); // artificial, then will be removed
		previousLimit = 0;
		boolean alreadyIn=false;
		for (int j= 0; j< numberOfSourcesPerGroup_SourceDomain.size(); j++){  // for each group j 
	//		System.out.println("[***DEBUG****] numberOfSourcesPerGroup_SourceDomain "+numberOfSourcesPerGroup_SourceDomain.get(j));
			for (int k = 0; k < numberOfSourcesPerGroup_SourceDomain.get(j); k++) {  // for each element k in the group
				//select the same number of sources as numberOfSourcesPerGroup: with prob analogicalNoveltyIndex_ from allSources current pop, 1-analogicalNoveltyIndex_ from ontology (only of that type). Just the non-already selected
				String sourceToAdd = null; 
				alreadyIn=false;
				if (Main.ran.nextDouble() < analogicalNoveltyIndex_) { 
	//				System.out.println("[***DEBUG****] Selection from current population");
					if (solution_string.size() < allSources.size()){ 
						ArrayList<String> temporarySources = new ArrayList<String>(sourcesPerGroupCurrentPop_TargetDomain.get(j));
						temporarySources.removeAll(solution_string); 
						
						
						long start_time_while =System.currentTimeMillis();
						while(solution_string.contains(sourceToAdd)&&(alreadyIn==false)) { 
							//sourceToAdd = selectFromCurrentPop(allSources, previousLimitPop, numberOfSourcesPerGroupCurrentPop_TargetDomain.get(j));
							//sourceToAdd = selectRandomFromCurrentPop(sourcesPerGroupCurrentPop_TargetDomain.get(j));
							//sourceToAdd = selectElementBinaryTournament(sourcesPerGroupCurrentPop_TargetDomain.get(j));
							if (!temporarySources.isEmpty())
								sourceToAdd = selectElementBinaryTournament(temporarySources);
							else
								alreadyIn=true;
						}
					}
				}
				else {
			//		System.out.println("[***DEBUG****] Selection from ontology");
					if (solution_string.size()!=Main.numericSources.size()){  // it can happen, for small population, that solution_string = all sources. In this case, do nothing
						
						long start_time_while =System.currentTimeMillis();
						int iteratations=1; //to remove 
						while(solution_string.contains(sourceToAdd)) {
							//sourceToAdd = selectRandomFromOntology(previousLimit, (int)Main.groupLimitNumber.get(j)); // Main.numericSources.get(RandomUtils.nextInt(previousLimit,Main.groupLimitNumber.get(j))); //These correspond to elements of that group type (they are stored in incremental order)
							  sourceToAdd = selectBinaryTournamentFromOntology(previousLimit, (int)Main.groupLimitNumber.get(j));
							  iteratations++;
						}

	//					System.out.println("[***DEBUG****] source to add: "+sourceToAdd);
					}
				}
	//			System.out.println("[***DEBUG****] Seelcted source: "+sourceToAdd);
				if (alreadyIn==false) {
					sourcesPerGroup_TargetDomain.get(j).add(sourceToAdd);
					solution_string.add(sourceToAdd);
				}
			}
			previousLimit = Main.groupLimitNumber.get(j);
			//previousLimitPop = numberOfSourcesPerGroupCurrentPop_TargetDomain.get(j);
			//maxKDegreePerGroup_TargetDomain.set(j, Main.plaus.getMaxKDegreeWithAutoExclusion(sourcesPerGroup_TargetDomain.get(j),true, 0)); //internal (target) domain
		}
		
	/*	
		// ADJUST THE GENERATED SOLUTION TO SATISFY ORDINAL CONSTRAINT 
		for (int index = 0; index < numberOfMaxAttempts_; index++) {
			// return -1 if the relation holds; return the index of the group to improve, if the relation does not hold
			int groupToChange =checkOrdinalConstraint(sourcesPerGroup_SourceDomain, sourcesPerGroup_TargetDomain, maxKDegreePerGroup_SourceDomain, maxKDegreePerGroup_TargetDomain);  
			if(groupToChange !=-1 && sourcesPerGroupCurrentPop_TargetDomain.get(groupToChange).size()!=0) {
				System.out.println("Adjust the solution to match ordinal constraints... ");
				System.out.println("sourcesPerGroupCurrentPop_TargetDomain size: "+sourcesPerGroupCurrentPop_TargetDomain.size() );
				for (int indextemp =0; indextemp<sourcesPerGroupCurrentPop_TargetDomain.size();indextemp++ )System.out.println("sourcesPerGroupCurrentPop_TargetDomain: "+sourcesPerGroupCurrentPop_TargetDomain.get(indextemp) );
				System.out.println("Group to change: "+groupToChange);
				// if -1, do nothing. 
				//If no, try to select sources to satisfy the constraint for at most N times
				String sourceToAdd = null; int indexTemp=0; boolean select=true;
				while(select==true && indexTemp < numberOfMaxAttempts_)  {
						sourceToAdd = selectBestFromCurrentPop(sourcesPerGroupCurrentPop_TargetDomain.get(groupToChange));
						if (!solution_string.contains(sourceToAdd))
							select=false;
						indexTemp++;
				}
				if(select==false) {
					System.out.println("Chosen source: "+sourceToAdd);
					// repalce an element
					int indexOfElementToReplace = -1; int minSupport = Integer.MAX_VALUE; //assume a support cannot be gerater than this 
					for (int indexElement=0; indexElement<sourcesPerGroup_TargetDomain.get(groupToChange).size();indexElement++) {
						int [] itemToMatch = new int[1];
						itemToMatch[0]  = Integer.parseInt(sourcesPerGroup_TargetDomain.get(groupToChange).get(indexElement)); 
	 					int support  = Main.itemsetTree.getSupportOfItemset(itemToMatch);
	 					if (support < minSupport) {
	 						minSupport = support;
	 						indexOfElementToReplace = indexElement;
	 					}
					}
					// REPLACE THE ELEMENT WITH MINIMAL SUPPORT WITH THE BEST OF THE GROUP
					sourcesPerGroup_TargetDomain.get(groupToChange).set(indexOfElementToReplace, sourceToAdd);
				}//else do nothing
			}
		}
		
		*/
		
			// Build the final solution
		solution_string.clear();
		for (int ind=0; ind<sourcesPerGroup_TargetDomain.size(); ind++)
			if(!sourcesPerGroup_TargetDomain.get(ind).contains(null))
				solution_string.addAll(sourcesPerGroup_TargetDomain.get(ind));
		
		solution_string.add(chosenTarget);  // add the target
		//System.out.println(solution_string);
		generatedSolution = SolutionUtils.createSolutionFromString(solution_string);
		
		
		((CauseEffectSolutionType)generatedSolution.getType()).setOperator("ANALOGICAL");
		
		return generatedSolution ;
	}

	private String selectBinaryTournamentFromOntology(int minLimit, int maxLimit) {
		//ArrayList<Integer> n_values_index  = new ArrayList<Integer>();
		int index=0;
		for (int i=0; i < Main.possibleSourcePerType.size(); i++) {
			if (getMaximumFromNumericString(Main.possibleSourcePerType.get(i)) <= maxLimit && getMinimumFromNumericString(Main.possibleSourcePerType.get(i))  >  minLimit) { 
		//		n_values_index.add(i);
				index = i;
				break; 
			}
		}
		//int selectedIndex = Main.ran.nextInt(n_values_index.size());  
		//return selectElementBinaryTournament(Main.possibleSourcePerType.get(n_values_index.get(selectedIndex)));
		return selectElementBinaryTournament(Main.possibleSourcePerType.get(index));
	}
	
	private String selectRandomFromOntology(int minLimit, int maxLimit) {
		ArrayList<Integer> n_values_index  = new ArrayList<Integer>();
		for (int i=0; i < Main.possibleSourcePerType.size(); i++) { 
			if (getMaximumFromNumericString(Main.possibleSourcePerType.get(i)) <= maxLimit && getMinimumFromNumericString(Main.possibleSourcePerType.get(i))  >  minLimit) 
				n_values_index.add(i); 
		}
		int selectedIndex = Main.ran.nextInt(n_values_index.size());  
		return Main.possibleSourcePerType.get(n_values_index.get(selectedIndex)).get(Main.ran.nextInt(Main.possibleSourcePerType.get(n_values_index.get(selectedIndex)).size()));
	}
	
	private int getMaximumFromNumericString(ArrayList<String> array) {
		ArrayList<Integer> numbers = new ArrayList<Integer>();
		for(int i = 0; i < array.size(); i++) 
		   numbers.add(Integer.parseInt(array.get(i)));   
		return Collections.max(numbers);
	}

	private int getMinimumFromNumericString(ArrayList<String> array) {
		ArrayList<Integer> numbers = new ArrayList<Integer>();
		for(int i = 0; i < array.size(); i++) 
		   numbers.add(Integer.parseInt(array.get(i)));   
		return Collections.min(numbers);
	}
	
	private int checkOrdinalConstraint(ArrayList<ArrayList<String>> sourcesPerGroup_SourceDomain,
			ArrayList<ArrayList<String>> sourcesPerGroup_TargetDomain,
			ArrayList<Integer> maxKDegreePerGroup_SourceDomain, ArrayList<Integer> maxKDegreePerGroup_TargetDomain) {
		
			// Select a pair randomly
			int firstGroupIndex = Main.ran.nextInt(maxKDegreePerGroup_SourceDomain.size());
			int secondGroupIndex = (firstGroupIndex + RandomUtils.nextInt(1, maxKDegreePerGroup_SourceDomain.size()))%maxKDegreePerGroup_SourceDomain.size();
			
			//Check if the ordinality relation holds; if not, return the index of the group with minimum MaxKDegree
			if ( (maxKDegreePerGroup_SourceDomain.get(firstGroupIndex)> maxKDegreePerGroup_SourceDomain.get(secondGroupIndex) && maxKDegreePerGroup_TargetDomain.get(firstGroupIndex) <= maxKDegreePerGroup_TargetDomain.get(secondGroupIndex)) )// || (maxKDegreePerGroup_SourceDomain.get(firstGroupIndex) <= maxKDegreePerGroup_SourceDomain.get(secondGroupIndex) && maxKDegreePerGroup_TargetDomain.get(firstGroupIndex) <= maxKDegreePerGroup_TargetDomain.get(secondGroupIndex)) )
				return firstGroupIndex; 		
			if ( (maxKDegreePerGroup_SourceDomain.get(firstGroupIndex)< maxKDegreePerGroup_SourceDomain.get(secondGroupIndex) && maxKDegreePerGroup_TargetDomain.get(firstGroupIndex) >= maxKDegreePerGroup_TargetDomain.get(secondGroupIndex)) ) //|| (maxKDegreePerGroup_SourceDomain.get(firstGroupIndex) <= maxKDegreePerGroup_SourceDomain.get(secondGroupIndex) && maxKDegreePerGroup_TargetDomain.get(firstGroupIndex) <= maxKDegreePerGroup_TargetDomain.get(secondGroupIndex)) )
				return secondGroupIndex;
			return -1;
		}

	private String selectRandomFromCurrentPop(ArrayList<String> elements) {
		 String source = new String(); 
		 source = elements.get(Main.ran.nextInt(elements.size()));
		return source;
	}

	private String selectBest(ArrayList<String> elements) {
		 
	//	System.out.println("Select the best from current population to adjust the solution ");
		 String element = new String(); 
		 int maxSupport = 0;
		 int indexOfChosenElement  = 0; //by default is the first element 
		 
		 for (int indexElement=0; indexElement<elements.size();indexElement++) {
				int [] itemToMatch = new int[1];
				itemToMatch[0]  = Integer.parseInt(elements.get(indexElement)); 
				int support  = Main.itemsetTree.getSupportOfItemset(itemToMatch);
				if (support >= maxSupport) {
					maxSupport = support;
					indexOfChosenElement = indexElement;
				}
			}
		 element = elements.get(indexOfChosenElement); 		 
		 return element; 
	}
	
	
	private String selectElementBinaryTournament(ArrayList<String> elementsList) {
		String element1 = new String();
		String element2 = new String();
	    element1 = elementsList.get(Main.ran.nextInt(elementsList.size()));
	    element2 = elementsList.get(Main.ran.nextInt(elementsList.size()));

	//     for(int h=0; h<elementsList.size();h++) System.out.println("eleemntlist h: "+elementsList.get(h)); 
	    
	    if (elementsList.size() ==1)
	    	return elementsList.get(0);
	    
	    if (elementsList.size() ==2) {
	    	if(elementsList.get(0).equals(elementsList.get(1))){
	    		return elementsList.get(0);
	    	}
	    	else {
	     	    element1 = elementsList.get(0);
	    	    element2 = elementsList.get(1);
	    	}
	    }
	    if (elementsList.size() > 2)
	    	while (element1.equals(element2))
	    		element2 = elementsList.get(Main.ran.nextInt(elementsList.size()));
	    
	    int[] item = new int[1]; 
	    item[0] = Integer.parseInt(element1); 
	    int support1 = causalOptimization.Main.itemsetTree.getSupportOfItemset(item);
	    
	    /*
	    ArrayList<Double> similarities=null;
		try {
			similarities = causalOptimization.Main.sim.computeElementKBSimilarity(element1);
		} catch (Exception e) {e.printStackTrace();}
	    double total_sim =0; 
		  for (int i =0; i<similarities.size();i++) total_sim  += similarities.get(i); 
		  double average_sim = total_sim/similarities.size();  // in 0 e 1
		double novelty1 = 1 - average_sim;
		*/
	    /*System.out.println("support1: "+support1);
	    System.out.println("element1: "+element1);
	    System.out.println("novelty1: "+novelty1);
	    */
	    item = new int[1];
	    item[0] = Integer.parseInt(element2); 
	    int support2 = causalOptimization.Main.itemsetTree.getSupportOfItemset(item);
	    
	   /* similarities.clear();
		try {
			similarities = causalOptimization.Main.sim.computeElementKBSimilarity(element2);
		} catch (Exception e) {e.printStackTrace();}
	    total_sim =0; 
		for (int i =0; i<similarities.size();i++) total_sim  += similarities.get(i); 
		average_sim = total_sim/similarities.size();  // in 0 e 1
		double novelty2 = 1 - average_sim;
		*/
		/*
	    System.out.println("support2: "+support2);
	    System.out.println("element2: "+element2);
	    System.out.println("novelty2: "+novelty2);
	    */
		
		  if (support1 > support2) 
			  return element1; 
		  else if (support1 < support2) 
			  return
		  element2; 
		  else if (Main.ran.nextDouble()<0.5) 
			  return element1;
		  else 
			  return element2;
		   
	    
	    /*
	    int flag = compare(element1, support1, novelty1, element2, support2, novelty2);
	    if (flag == -1) //It means that support and novelty of element 1 are smaller than those of element 2
	      return element2;
	    else if (flag == 1)
	      return element1;
	    else
	      if (Main.ran.nextDouble()<0.5)
	        return element1;
	      else
	        return element2;
	        */
	}

	
	int compare(String s1, double obj1_s1, double obj2_s1, String s2, double obj1_s2, double obj2_s2) {
		
		if (s1==null)
		      return 1;
		    else if (s2 == null)
		      return -1;
		    
		if (obj1_s1 <= obj1_s2 && obj2_s1 <= obj2_s2) {
			if (obj1_s1 == obj1_s2 && obj2_s1 == obj2_s2) {
				return 0;
			} else {
				return -1;
				}
		} else if (obj1_s1 >= obj1_s2 && obj2_s1 >= obj2_s2) {
			return 1;
		}
		else { // they are discordant
			return 0;
		}
	    
	}
	
}
