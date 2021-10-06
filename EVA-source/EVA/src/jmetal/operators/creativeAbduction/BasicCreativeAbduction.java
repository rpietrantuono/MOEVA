package jmetal.operators.creativeAbduction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import causalOptimization.Main;
import jmetal.core.Solution;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import solutionType.CauseEffectSolutionType;
import util.SolutionUtils;

public class BasicCreativeAbduction extends CreativeAbduction{
	
	private int creativeChangeIndex_ ;
	private double creativeNoveltyIndex_;
	public BasicCreativeAbduction(HashMap<String, Object> parameters) {
		super(parameters);
		creativeChangeIndex_ = 1; 
		System.out.println(parameters.size());
		
		if (parameters.get("creativeChangeIndex") != null)
			creativeChangeIndex_ = (int)parameters.get("creativeChangeIndex") ;
	
		if (parameters.get("creativeNoveltyIndex") != null)
			creativeNoveltyIndex_ = (double)parameters.get("creativeNoveltyIndex") ;
		
	}


	@Override
	public Object execute(Object object) throws JMException {
		return null;
	}
	
	/**
	 * @param object: The chosen solution by the selection operator 
	 * @param allSources: The list of possible sources taken from the current population
	 * @param allTargets:  The list of possible targets taken from the current population
	 * @return generatedSolution: a new Solution
	 * @throws JMException
	 */
	public Object execute(Object object, ArrayList<String> allSources, ArrayList<String> allTargets) throws JMException {

	
	//	System.out.println("creativeChangeIndex_ "+creativeChangeIndex_);
//		System.out.println("creativeNoveltyIndex "+creativeNoveltyIndex_);
		Solution chosenSolution  = (Solution)object;
		
		//select one target in the solution.  
		//String currentTarget = SolutionUtils.getTargets(chosenSolution).get(RandomUtils.nextInt(0, (SolutionUtils.getTargets(chosenSolution)).size()));
		String chosenTarget = null; 
	//	if (Main.ran.nextDouble() < creativeNoveltyIndex_)
			chosenTarget = selectElementBinaryTournament(allTargets);
	//	else 
	//		chosenTarget = selectElementBinaryTournament(Main.possibleTargets);

			//	String chosenTarget= SolutionUtils.getTargets(chosenSolution).get(Main.ran.nextInt((SolutionUtils.getTargets(chosenSolution)).size()));
		
		ArrayList<String> sources = SolutionUtils.getSources(chosenSolution);
		
		// the starting population is generated randomly, for the creative case contains elements in Ontoology
		//number of changes to do 
		int numberOfChanges = Main.ran.nextInt(creativeChangeIndex_)+1; //at least one change
	//	System.out.println("Number of changes: "+numberOfChanges);
	//	System.out.println("[***DEBUG***] sources: "+sources);
	//	System.out.println("[***DEBUG***] All sources: "+allSources);
		
		int minChoice=0, maxChoice= 2;
		boolean alreadyIn=false;
		for (int i=0; i< numberOfChanges; i++) {
			int action ;
			if (sources.size()>=causalOptimization.Main.MAX_sources)
				minChoice=1; //,Max 
			else
				minChoice=0;
			if (sources.size()<=causalOptimization.Main.MIN_sources)
				maxChoice=1;
			else
				maxChoice=2;
			
			alreadyIn=false;
			//action = RandomUtils.nextInt(minChoice,maxChoice+1); //Add = 0, Modify = 1, Delete = 2
			action = Main.ran.nextInt(maxChoice+1-minChoice) + minChoice; //Add = 0, Modify = 1, Delete = 2
			if (action == 0 && sources.size()<causalOptimization.Main.MAX_sources) { //ADD source
				String sourceToAdd = sources.get(0);
				if (Main.ran.nextDouble() < creativeNoveltyIndex_ && sources.size() < allSources.size()) { 
				//	if (sources.size() < allSources.size()){ if this does not happen, add from ontology
						ArrayList<String> temporarySources = new ArrayList<String>(allSources);
						temporarySources.removeAll(sources); 
						while(sources.contains(sourceToAdd)&&(alreadyIn==false)) { 
							//sourceToAdd = allSources.get(RandomUtils.nextInt(0,allSources.size()));
							//sourceToAdd = selectElementBinaryTournament(allSources);
							if (!temporarySources.isEmpty())
								sourceToAdd = selectElementBinaryTournament(temporarySources);
								//sourceToAdd = selectBestFromCurrentPop(temporarySources);
							else
								alreadyIn=true;
						}
	//					System.out.println("[***DEBUG***] source to add, from current population: "+sourceToAdd);
					//}
				}
				else {
					//from ontology
					if (sources.size()<causalOptimization.Main.numericSources.size()){ // it can happen, for small population, that solution_string = all sources. In this case, do nothing
						ArrayList<String> temporarySources = new ArrayList<String>(causalOptimization.Main.numericSources);
						temporarySources.removeAll(sources);
						while(sources.contains(sourceToAdd)&&(alreadyIn==false)) { 
							//sourceToAdd = selectRandomFromOntology(0,Main.numericSources.size());
							if (!temporarySources.isEmpty()) {
							//	sourceToAdd = selectBinaryTournamentFromOntology(0,Main.numericSources.size());
								sourceToAdd = selectElementBinaryTournament(Main.possibleSourcePerType.get(Main.ran.nextInt(Main.possibleSourcePerType.size())));
							}
								//sourceToAdd = selectBestFromCurrentPop(temporarySources);
							else
								alreadyIn=true;
						}
	//				System.out.println("[***DEBUG***] source to add, from ontology: "+sourceToAdd);					
					}
				}
				if (alreadyIn==false)
					sources.add(sourceToAdd);
			
			}
			if (action == 1) {
				String sourceToAdd = sources.get(0);
				String sourceToRemove = null;
				if (Main.ran.nextDouble() < creativeNoveltyIndex_) { 
					ArrayList<String> temporarySources = new ArrayList<String>(allSources);
					temporarySources.removeAll(sources); 
					while(sources.contains(sourceToAdd)&&(alreadyIn==false)) {
						//sourceToAdd = allSources.get(RandomUtils.nextInt(0,allSources.size())); //we cna change this random selection in something better, like in analogical
						//sourceToAdd = selectElementBinaryTournament(allSources);
						if (!temporarySources.isEmpty())
							sourceToAdd = selectElementBinaryTournament(temporarySources);
							//sourceToAdd = selectBestFromCurrentPop(temporarySources);
						else {
							alreadyIn=true;}
	//					System.out.println("source to add, from current population: "+sourceToAdd);
					}
				//	}
				}
				else {
					//from ontology
					if (sources.size()< causalOptimization.Main.numericSources.size()){ // it can happen, for small population, that solution_string = all sources. In this case, do nothing
						ArrayList<String> temporarySources = new ArrayList<String>(causalOptimization.Main.numericSources);
						temporarySources.removeAll(sources);
						while(sources.contains(sourceToAdd)&&(alreadyIn==false)) {
							//sourceToAdd = selectRandomFromOntology(0,Main.numericSources.size()); //Main.numericSources.get(RandomUtils.nextInt(0,Main.numericSources.size()));
							if (!temporarySources.isEmpty()) {
								//sourceToAdd = selectElementBinaryTournament(temporarySources);
								//sourceToAdd = selectBestFromCurrentPop(temporarySources);
								sourceToAdd = selectElementBinaryTournament(Main.possibleSourcePerType.get(Main.ran.nextInt(Main.possibleSourcePerType.size())));
							//	sourceToAdd = selectBinaryTournamentFromOntology(0,Main.numericSources.size());
								
								}
							else {
								alreadyIn=true;}
		//					System.out.println("selected from ontology "+sourceToAdd);
						}
		//				System.out.println("source to add, selected from ontology: "+sourceToAdd);
					}
				}
				if (alreadyIn==false) {
					sources.add(sourceToAdd);
				// REMOVE
				//		do{
					//sourceToRemove = sources.get(RandomUtils.nextInt(0,sources.size()));
		//			System.out.println("REMOVE ");
		//			System.out.println("sources "+sources);
					sourceToRemove  = selectElementToRemoveBinaryTournament(sources);
					//		}while(sourceToRemove.equals(sourceToAdd));  //with 2 solutions, it may hang
					//		System.out.println("[***DEBUG***] source to remove: "+sourceToRemove);
				sources.remove(sourceToRemove);
				}
			}
			
			if (action == 2 && sources.size()>causalOptimization.Main.MIN_sources) {
				if (sources.size()>1) {
		//			System.out.println("sources "+sources);
					sources.remove(selectElementToRemoveBinaryTournament(sources)); 
					//sources.remove(RandomUtils.nextInt(0,sources.size()));
				}
			}
		}
		
		ArrayList<String> solution_string = new ArrayList<String>();
		solution_string.addAll(sources);
		solution_string.add(chosenTarget);
		//create the new solution 
		Solution generatedSolution = SolutionUtils.createSolutionFromString(solution_string);
	//	System.out.println("ORIGINAL SOLUTION");
	//	System.out.println(SolutionUtils.getNumericRepresentation(chosenSolution));
	//	System.out.println("\nSOLUTION AFTER CREATIVE ABDUCTION");
	//	System.out.println(solution_string);
		
		((CauseEffectSolutionType)generatedSolution.getType()).setOperator("CREATIVE");
		
		return generatedSolution ;
				
				
	}
	

	private String selectBestFromCurrentPop(ArrayList<String> elements) {
		 
		//	System.out.println("Select the best from current population to adjust the solution ");
			 String source = new String(); 
			 int maxSupport = 0;
			 int indexOfChosenElement  = 0; //by default is the first element 
			 
			 for (int indexElement=0; indexElement<elements.size();indexElement++) {
					int [] itemToMatch = new int[1];
					itemToMatch[0]  = Integer.parseInt(elements.get(indexElement)); 
					int support  = causalOptimization.Main.itemsetTree.getSupportOfItemset(itemToMatch);
					if (support >= maxSupport) {
						maxSupport = support;
						indexOfChosenElement = indexElement;
					}
				}
			 source = elements.get(indexOfChosenElement); 		 
			 return source; 
		}
	
	private String selectBinaryTournamentFromOntology(int minLimit, int maxLimit) {
		ArrayList<Integer> n_values_index  = new ArrayList<Integer>();
		for (int i=0; i < causalOptimization.Main.possibleSourcePerType.size(); i++) { 
			if (getMaximumFromNumericString(causalOptimization.Main.possibleSourcePerType.get(i)) <= maxLimit && getMinimumFromNumericString(Main.possibleSourcePerType.get(i))  >  minLimit) 
				n_values_index.add(i); 
		}
	//	System.out.println("size of n value index creative  "+n_values_index.size());
		int selectedIndex = Main.ran.nextInt(n_values_index.size());  
		return selectElementBinaryTournament(causalOptimization.Main.possibleSourcePerType.get(n_values_index.get(selectedIndex)));
	}
	
	private String selectRandomFromOntology(int minLimit, int maxLimit) {
		ArrayList<Integer> n_values_index  = new ArrayList<Integer>();
		for (int i=0; i < causalOptimization.Main.possibleSourcePerType.size(); i++) { 
			if (getMaximumFromNumericString(causalOptimization.Main.possibleSourcePerType.get(i)) <= maxLimit && getMinimumFromNumericString(Main.possibleSourcePerType.get(i))  >  minLimit) 
				n_values_index.add(i); 
		}
		int selectedIndex = Main.ran.nextInt(n_values_index.size());  
		return causalOptimization.Main.possibleSourcePerType.get(n_values_index.get(selectedIndex)).get(Main.ran.nextInt(Main.possibleSourcePerType.get(n_values_index.get(selectedIndex)).size()));
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


	private String selectElementToRemoveBinaryTournament(ArrayList<String> elementsList) {
		//THE SAME AS THE PREVIOUS, BUT WITH INVERTED OUTCOME (IT IS RETURNED THE WORST ELEMENT AMONG TWO 
		String element1 = new String();
		String element2 = new String();
	    element1 = elementsList.get(Main.ran.nextInt(elementsList.size()));
	    element2 = elementsList.get(Main.ran.nextInt(elementsList.size()));

	//    System.out.println("Element To Remove: elementlist size: "+elementsList.size());
	    //for(int h=0; h<elementsList.size();h++) System.out.println("eleemntlist h: "+elementsList.get(h)); 
	    
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
	    	while (element1.equals(element2)) {
	    		int index_ = Main.ran.nextInt(elementsList.size());
	    		element2 = elementsList.get(index_);
	    		//System.out.println("select element cycle - element to remove - factual "+index_+"; el list "+elementsList);	    	
	    		}
	    
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
		/*
	    System.out.println("support1: "+support1);
	    System.out.println("element1: "+element1);
	    System.out.println("novelty1: "+novelty1);
	    */
	    item = new int[1];
	    item[0] = Integer.parseInt(element2); 
	    int support2 = causalOptimization.Main.itemsetTree.getSupportOfItemset(item);
	    /*
	    similarities.clear();
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
	    	return element2; 
	    else if (support1 < support2) 
	    	return element1; 
	    else if (Main.ran.nextDouble()<0.5) 
	    	return element1; 
	    else return
	  		  element2;
	  		    
	    
		/*
		 int flag = compare(element1, support1, novelty1, element2, support2, novelty2); // IT IS THE OPPOSITE OF THE PREVIOUS FUNCTION 
		 if (flag == -1)
			 return element1; 
		 else if (flag == 1) 
			 	return element2; 
		 	  else if
		 				(Main.ran.nextDouble()<0.5) 
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
	/*
	 * private String selectElementBinaryTournament(ArrayList<String> allTargets) {
	 * String target1, target2; target1 =
	 * allTargets.get(PseudoRandom.randInt(0,allTargets.size()-1)); target2 =
	 * allTargets.get(PseudoRandom.randInt(0,allTargets.size()-1));
	 * 
	 * if (allTargets.size() >= 2) while (target1.equals(target2)) target2 =
	 * allTargets.get(PseudoRandom.randInt(0,allTargets.size()-1));
	 * 
	 * int[] item = new int[1]; item[0] = Integer.parseInt(target1); int support1 =
	 * Main.itemsetTree.getSupportOfItemset(item); item[0] =
	 * Integer.parseInt(target2); int support2 =
	 * Main.itemsetTree.getSupportOfItemset(item);
	 * 
	 * if (support1 > support2) return target1; else if (support1 < support2) return
	 * target2; else if (PseudoRandom.randDouble()<0.5) return target1; else return
	 * target2; }
	 */

}
