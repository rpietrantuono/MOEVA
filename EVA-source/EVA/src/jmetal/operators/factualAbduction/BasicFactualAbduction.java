package jmetal.operators.factualAbduction;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import causalOptimization.Main;
import jmetal.core.Solution;
import jmetal.util.JMException;
//import jmetal.util.PseudoRandom;
import solutionType.CauseEffectSolutionType;
import util.SolutionUtils;

public class BasicFactualAbduction extends FactualAbduction{

	private int factualChangeIndex_ ;
	private double factualNoveltyIndex_;
	
	public BasicFactualAbduction(HashMap<String, Object> parameters) {
		super(parameters);
		factualChangeIndex_ = 1; 
		System.out.println(parameters.size());
		
		if (parameters.get("factualChangeIndex") != null)
			factualChangeIndex_ = (int)parameters.get("factualChangeIndex") ;
		
		if (parameters.get("factualNoveltyIndex") != null)
			factualNoveltyIndex_ = (double)parameters.get("factualNoveltyIndex") ;
	}


	@Override
	
	public Object execute(Object object) throws JMException {
		return null;
	}
	
	/*
	 * @object: chosen solution by the selection operator.
	 * @allSources: all the sources in the KB
	 */
	/**
	 * @param object: The chosen solution by the selection operator 
	 * @param allSources: The list of possible sources taken from the current population
	 * @param allTargets:  The list of possible targets taken from the current population
	 * @return generatedSolution: a new Solution
	 * @throws JMException
	 */
	public Object execute(Object object, ArrayList<String> allSources, ArrayList<String> allTargets) throws JMException {

		Solution generatedSolution = new Solution();
	
//	System.out.println("Factual change index"+factualChangeIndex_);
//	System.out.println("Factual nov index"+factualNoveltyIndex_);

	
		//Consider the current population (selection is be done by a selection operator: best, worst, torunament)
		Solution chosenSolution  = (Solution)object;
		
		//Select one target in the solution.  
		//String currentTarget = SolutionUtils.getTargets(chosenSolution).get(RandomUtils.nextInt(0, (SolutionUtils.getTargets(chosenSolution)).size()));
	//	System.out.println("factualNoveltyIndex_ "+factualNoveltyIndex_);
		String chosenTarget = null; 
		if (Main.ran.nextDouble() < factualNoveltyIndex_)
			chosenTarget = selectElementBinaryTournament(allTargets);
		else 
			chosenTarget = selectElementBinaryTournament(Main.possibleKBtargets);		
	//	String chosenTarget= SolutionUtils.getTargets(chosenSolution).get(Main.ran.nextInt((SolutionUtils.getTargets(chosenSolution)).size()));
		
	//	System.out.println("Chosen target: "+chosenTarget);
		ArrayList<String> sources = SolutionUtils.getSources(chosenSolution);
		
		// the starting population is generated randomly, but for the factual abduction case it contains only elements in KB
		// for the others, can contain elements in Ontology
		int numberOfChanges = Main.ran.nextInt(factualChangeIndex_) +1; //at least one change
	//	System.out.println("Number of changes: "+numberOfChanges);
	//	System.out.println("[***DEBUG***] sources: "+sources);
	//	System.out.println("[***DEBUG***] All sources: "+allSources);
		
		int minChoice=0, maxChoice= 2;
		boolean alreadyIn = false;
		for (int i=0; i< numberOfChanges; i++) {
			int action ;
			if (sources.size()>=causalOptimization.Main.MAX_sources)
				minChoice=1;
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
				if (Main.ran.nextDouble() < factualNoveltyIndex_ && sources.size() < allSources.size()) { // it happens, for small population, that sources = all sources. In this case, add from KB
				//	if (sources.size() < allSources.size()){ // it happens, for small population, that sources = all sources. In this case, do nothing
						ArrayList<String> temporarySources = new ArrayList<String>(allSources);
						temporarySources.removeAll(sources);
					//	System.out.println("temporarySources allsources "+temporarySources);
					//	System.out.println("temporarySources sources "+sources);
						while(sources.contains(sourceToAdd)&&(alreadyIn==false)) { // IT IS NOT NECESSARY
							//sourceToAdd = allSources.get(RandomUtils.nextInt(0,allSources.size())); //we cna change this random selection in something better, like in analogical
							if (!temporarySources.isEmpty())
								sourceToAdd = selectElementBinaryTournament(temporarySources);
								//sourceToAdd = selectBestFromCurrentPop(temporarySources);
							else
								alreadyIn=true;
		//				System.out.println("[***DEBUG***]source to add, Selected from current population: "+sourceToAdd);
						}
					//}
				}
				else {
					//from KB
					if (sources.size()<causalOptimization.Main.numericSources.size()){ // it can happen, for small population, that solution_string = all sources. In this case, do nothing
						ArrayList<String> temporarySources = new ArrayList<String>(causalOptimization.Main.possibleKBsources);
						temporarySources.removeAll(sources);
						while(sources.contains(sourceToAdd)&&(alreadyIn==false)) {
							//sourceToAdd = selectElementBinaryTournament(causeEffect.Main.possibleKBsources);// COULD BE DONE FOR THE TARGET TOO
							//sourceToAdd = selectElementBinaryTournament(temporarySources);// COULD BE DONE FOR THE TARGET TOO
							if (!temporarySources.isEmpty()) {
								//sourceToAdd = selectElementBinaryTournament(temporarySources);
								sourceToAdd = selectElementBinaryTournament(Main.possibleKBsourcesPerType.get(Main.ran.nextInt(Main.possibleKBsourcesPerType.size())));
								//sourceToAdd = selectBinaryTournamentFromKB(0,Main.numericSources.size());
					//			System.out.println("selected from kb "+sourceToAdd);
							}
								//sourceToAdd = selectBestFromCurrentPop(temporarySources);
								//causeEffect.Main.possibleKBsources.get(RandomUtils.nextInt(0,causeEffect.Main.possibleKBsources.size()));
							else
								alreadyIn=true;
						}
					}
				}
				if (alreadyIn==false)
					sources.add(sourceToAdd);
			}
			if (action == 1) {
				String sourceToAdd = sources.get(0);
				String sourceToRemove = null; 
				if (Main.ran.nextDouble() < factualNoveltyIndex_) {
					ArrayList<String> temporarySources = new ArrayList<String>(allSources);
					temporarySources.removeAll(sources); 
					//System.out.println("temporarySources allsources "+temporarySources);
				//	System.out.println("temporarySources sources "+sources);
					while(sources.contains(sourceToAdd)&&(alreadyIn==false)) {
						//sourceToAdd = allSources.get(RandomUtils.nextInt(0,allSources.size())); //we cna change this random selection in something better, like in analogical
						//sourceToAdd = selectElementBinaryTournament(allSources); 
						if (!temporarySources.isEmpty())
							sourceToAdd = selectElementBinaryTournament(temporarySources);
							//sourceToAdd = selectBestFromCurrentPop(temporarySources);
						else
							alreadyIn=true;
					}
	//				System.out.println("[***DEBUG***] source to add,Selected from current population : "+sourceToAdd);
				}
				else {
					//from KB
					if (sources.size()<causalOptimization.Main.numericSources.size()){ // it can happen, for small population, that solution_string = all sources. In this case, do nothing
						ArrayList<String> temporarySources = new ArrayList<String>(causalOptimization.Main.possibleKBsources);
						temporarySources.removeAll(sources);
				///		System.out.println("temporarySources allsources "+temporarySources);
				//		System.out.println("temporarySources sources "+sources);
							while(sources.contains(sourceToAdd)&&(alreadyIn==false)) { 
							//sourceToAdd = selectElementBinaryTournament(causeEffect.Main.possibleKBsources);
							if (!temporarySources.isEmpty()) {
								sourceToAdd = selectElementBinaryTournament(Main.possibleKBsourcesPerType.get(Main.ran.nextInt(Main.possibleKBsourcesPerType.size())));
								//sourceToAdd = selectBinaryTournamentFromKB(0,Main.numericSources.size());
							//	System.out.println("selected from KB"+sourceToAdd);
							//	sourceToAdd = selectElementBinaryTournament(temporarySources);//causeEffect.Main.possibleKBsources.get(RandomUtils.nextInt(0,causeEffect.Main.possibleKBsources.size()));
							//	sourceToAdd = selectBestFromCurrentPop(temporarySources);
							}
							else
								alreadyIn=true;
						}
					//System.out.println("[***DEBUG***] source to add, Selected from  KB: "+sourceToAdd);					
					}
				}
				if (alreadyIn==false) {
					sources.add(sourceToAdd);
				//REMOVE
			//	do {
					sourceToRemove = selectElementToRemoveBinaryTournament(sources);//;sources.get(RandomUtils.nextInt(0,sources.size()));
					//	}while(sourceToRemove.equals(sourceToAdd)) ; //with 2 solutions, it may hang
					//		System.out.println("[***DEBUG***] source to remove: "+sourceToRemove);
					sources.remove(sourceToRemove);
				}
			}
			if (action == 2 && sources.size()>causalOptimization.Main.MIN_sources) {
				String sourceToRemove = null;
				if (sources.size()>1) {
			//		System.out.println("sources "+sources);
					sourceToRemove = selectElementToRemoveBinaryTournament(sources); 
		//			System.out.println("[***DEBUG***] source to remove: "+sourceToRemove);
					sources.remove(sourceToRemove);
					//sources.remove(RandomUtils.nextInt(0,sources.size()));
				}
			}
		
		}
		
		ArrayList<String> solution_string = new ArrayList<String>();
		solution_string.addAll(sources);
		solution_string.add(chosenTarget);
		//create the new solution 
		generatedSolution = SolutionUtils.createSolutionFromString(solution_string);
	//	System.out.println("ORIGINAL SOLUTION");
	//	System.out.println(SolutionUtils.getNumericRepresentation(chosenSolution));
	//	System.out.println("\nSOLUTION AFTER FACTUAL ABDUCTION");
	//	System.out.println(solution_string);
		
		((CauseEffectSolutionType)generatedSolution.getType()).setOperator("FACTUAL");
		
		return generatedSolution ;
	}

	
	private String selectBinaryTournamentFromKB(int minLimit, int maxLimit) {
		ArrayList<Integer> n_values_index  = new ArrayList<Integer>();
		for (int i=0; i < causalOptimization.Main.possibleKBsourcesPerType.size(); i++) { 
			if (getMaximumFromNumericString(causalOptimization.Main.possibleKBsourcesPerType.get(i)) <= maxLimit && getMinimumFromNumericString(Main.possibleKBsourcesPerType.get(i))  >  minLimit && causalOptimization.Main.possibleKBsourcesPerType.get(i).size()>0) 
				n_values_index.add(i); 
		}
		int selectedIndex = Main.ran.nextInt(n_values_index.size());  
		return selectElementBinaryTournament(causalOptimization.Main.possibleKBsourcesPerType.get(n_values_index.get(selectedIndex)));
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
}
