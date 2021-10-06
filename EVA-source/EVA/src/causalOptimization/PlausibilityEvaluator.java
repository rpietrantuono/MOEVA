package causalOptimization;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import javax.print.attribute.standard.Media;
import javax.xml.transform.Source;

import com.sun.jdi.Value;
import com.sun.source.tree.Tree;

import jmetal.core.Solution;
import jmetal.core.Variable;
import myFrequentpatterns.negFIN.AlgoNegFIN;
import myItemsetTree.AssociationRuleIT;
import myItemsetTree.HashTableIT;
import util.SolutionUtils;
import org.apache.commons.math3.random.EmpiricalDistribution;

import ca.pfv.spmf.algorithms.associationrules.TopKRules_and_TNR.Database;
import ca.pfv.spmf.algorithms.associationrules.TopKRules_and_TNR.Transaction;
import ca.pfv.spmf.patterns.itemset_array_integers_with_count.Itemset;
import ca.pfv.spmf.patterns.itemset_array_integers_with_count.Itemsets;
import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;


public class PlausibilityEvaluator {
	
	double average_KB_sol_size;
	double std_KB_sol_size;
	double KB_sample_size;
	double UCL,LCL; // CL at 95% of confidence
	double[] quartiles; 
	//TDistribution tDist =null; 
	//EmpiricalDistribution eDist = null;
	public enum Evaluator {RAW, RAW_RESCALED, ROOT, LOG, RESCALED};
	Evaluator evaluator; 
	int maxDegree; 
	public static int count = 0; 
	
	public void setKB(double[] sizes){
		
	//	if (sizes.length<10){  
			//tDist = new TDistribution(KB_sample_size-1);
		//	UCL = average_KB_sol_size + tDist.inverseCumulativeProbability(1-significance/2)*Math.sqrt(std_KB_sol_size/KB_sample_size); 
		// 	LCL = average_KB_sol_size - tDist.inverseCumulativeProbability(1-significance/2)*Math.sqrt(std_KB_sol_size/KB_sample_size);
	//		System.out.println("A KB under 10 is not enough to use the 'Raw rescaled evaluator'. Will use the default evaluator");
//			evaluator = evaluator.RAW;
//		}
//		else {
			evaluator = evaluator.RAW_RESCALED;
			//eDist = new EmpiricalDistribution();
			//eDist.load(sizes);
		
			average_KB_sol_size = 0;
			std_KB_sol_size = 0; 
			KB_sample_size = sizes.length;
			for (int i = 0; i < sizes.length;i++ ) 
				average_KB_sol_size  = average_KB_sol_size + sizes[i];
			
			average_KB_sol_size = average_KB_sol_size/KB_sample_size;
			
			for (int i = 0; i < sizes.length;i++ )
				std_KB_sol_size = std_KB_sol_size + Math.pow((sizes[i] - average_KB_sol_size),2); 
					
			std_KB_sol_size = Math.sqrt(std_KB_sol_size/(KB_sample_size-1)); 
			
			quartiles = SolutionUtils.Quartiles(sizes);   
			
			System.out.println("average_KB_sol_size "+average_KB_sol_size);
			System.out.println("std_KB_sol_size "+std_KB_sol_size);
			System.out.println("KB_sample_size "+KB_sample_size);
			System.out.println("Quartiles "+quartiles[0]);
			System.out.println("Quartiles "+quartiles[1]);
			System.out.println("Quartiles "+quartiles[2]);
			
			//System.out.println("prob at average "+eDist.cumulativeProbability(Math.round(0))); 
			//.probability(average_KB_sol_size));
//		}
 		
	}
	
	public double evaluate(Solution solution, boolean internal) throws Exception {
		double x = 0;
		switch(evaluator) {
			case RAW: x = rawPlausibilityEvaluation(solution, internal);
			break; 
		
			case RAW_RESCALED: x = rawRescaledPlausibilityEvaluation(solution, internal);
			break;
			
			case ROOT: x = rootPlausibilityEvaluation(solution, internal);
			break;
			
			case LOG: x = logarithmicPlausibilityEvaluation(solution, internal);
			break;
			
			case RESCALED: x = rescaledPlausibilityEvaluation(solution, internal);
			break;
			
			default:  x = plausibilityEvaluation(solution, internal);
			break;
		
		}
		return x;
		
	}
	
	public double logarithmicPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  
		  ArrayList<String> targets = new ArrayList<String>();
		  targets.addAll(SolutionUtils.getTargets(solution));
		  ArrayList<String> sources =  new ArrayList<String>();
		  sources.addAll(SolutionUtils.getSources(solution));
		
		  //System.out.println("Size of the source "+sources.size());
		  int N = sources.size();
		  double x = plausibilityEvaluation(solution, internal); 
		  
		  
		  // ADJUST BASED ON SOURCE SIZE - TO DO
		  //LOGARITHMIC 
		 // System.out.println("result "+Math.log(1+Math.pow(x, 1.0/N))/Math.log(2));
		  
		  return  Math.log(1+Math.pow(x, 1.0/N))/Math.log(2);
		  
		  
	  }
	  
	  public double rootPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  
		  ArrayList<String> targets = new ArrayList<String>();
		  targets.addAll(SolutionUtils.getTargets(solution));
		  ArrayList<String> sources =  new ArrayList<String>();
		  sources.addAll(SolutionUtils.getSources(solution));
		
		  int N = sources.size();
		  double x = plausibilityEvaluation(solution, internal); 
		  
		  
		 return  Math.pow(x, 1.0/N);
		  
		  
	  }

	  
	  public double rescaledPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  
		  ArrayList<String> targets = new ArrayList<String>();
		  targets.addAll(SolutionUtils.getTargets(solution));
		  ArrayList<String> sources =  new ArrayList<String>();
		  sources.addAll(SolutionUtils.getSources(solution));
		
		  int N = sources.size();
		  double x = plausibilityEvaluation(solution, internal); 
		  
		  
		  double rootPlaus =  Math.pow(x, 1.0/N);
		  return rootPlaus * N/Main.MAX_sources;  
		  
	  }

	  
	  public double associationRulesEvaluation(Solution solution, boolean internal) throws Exception {
			
	  	
//			ArrayList<String> sol_string = SolutionUtils.getNumericRepresentation(solution);  
			ArrayList<String> targets = new ArrayList<String>();
			targets.addAll(SolutionUtils.getTargets(solution));
			ArrayList<String> sources =  new ArrayList<String>();
			sources.addAll(SolutionUtils.getSources(solution));
			  
			System.out.println("solution to evaluate: "+sources+" with target: "+targets);
			//String internalKnowledgeBaseFile = Main.baseDir + "/internalKBFile_test.txt";
			
			 
			//Integer[] sInt = Arrays.stream(s).boxed().toArray(Integer[]::new);
			Integer[] sInt  = new Integer[sources.size()]; //+targets.size()]; 
		
			for (int ind = 0; ind < sources.size(); ind++) 
				sInt[ind] = Integer.parseInt(sources.get(ind)); 
			
			Database database_targets_copy = new Database();
		
			for(int k = 0; k< ((Database)Main.database_map.get(targets.get(0))).size();k++) {
				List<Integer> list = (List<Integer>) ((Database)Main.database_map.get(targets.get(0))).getTransactions().get(k).getItems();
			//	ArrayList<Integer> arrayList = new ArrayList<Integer>(list);
				List<String> newList = new ArrayList<>(list.size());
				for (Integer myInt : list) { newList.add(String.valueOf(myInt)); }
				String[] itemsCopy = new String[newList.size()] ;
				newList.toArray(itemsCopy );
				//((Database)Main.database_map.get(targets.get(0))).getTransactions().get(k).getItems().toArray(itemsCopy);
				//itemsCopy = 
				database_targets_copy.addTransaction(itemsCopy);
			}
			
			
			for(int i = 0; i<database_targets_copy.getTransactions().size(); i++) {
				//System.out.println("transactions "+database_targets_copy.getTransactions().get(i).getItems());
				database_targets_copy.getTransactions().get(i).getItems().retainAll(Arrays.asList(sInt));
			}
			 
		//	long start  = System.currentTimeMillis();
			//double minsup = 0.4; // means a minsup of 2 transaction (we used a relative support)
			double minsup = 1.0/database_targets_copy.getTransactions().size();
			
			// Applying the algorithm
			AlgoNegFIN algorithm = new AlgoNegFIN();
			
			if (database_targets_copy.size()==0 ) {
				algorithm.outputCount = 0;
			}
			else if (sInt.length==1) {
				algorithm.outputCount = 1;
			}	
			else {
			
			//algorithm.runAlgorithm("output_db.txt", minsup, "output_db_algo.txt"); // PER CONFRONTO
			algorithm.runAlgorithm_memory(database_targets_copy, minsup);
			algorithm.printStats();
		
			System.out.println("count: "+algorithm.outputCount);
			System.out.println("run: "+Main.run);
			}
			int N = sources.size();
			double scale_factor; 
			if (N>average_KB_sol_size)
					scale_factor = average_KB_sol_size/N;
			 else
					scale_factor = N/average_KB_sol_size;
			 	
		
			return algorithm.outputCount*scale_factor;
	  
	  }
	  
	  
	  public double rawPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  
		  myItemsetTree.ItemsetTree tree = null;  
			  tree = Main.itemsetTree;
		  
	    	int support;
	    	Main.supportMap.clear();
	    	
//			ArrayList<String> sol_string = SolutionUtils.getNumericRepresentation(solution);  
			ArrayList<String> targets = new ArrayList<String>();
			targets.addAll(SolutionUtils.getTargets(solution));
			ArrayList<String> sources =  new ArrayList<String>();
			sources.addAll(SolutionUtils.getSources(solution));
			  
			System.out.println("\nPlausibility evaluation...");
			System.out.println("sources: " + sources);
			System.out.println("targets: " + targets);
				  
		
			  int N = sources.size();//SolutionUtils.getSolutionSize(example); 
			  
			  
			  ArrayList<ArrayList<String>> combinationsEvaluated = new ArrayList<ArrayList<String>>(); 
			  ArrayList<ArrayList<String>> combinationsToEvaluate = new ArrayList<ArrayList<String>>();
			  ArrayList<ArrayList<String>> combinationsToSkip = new ArrayList<ArrayList<String>>();
			  ArrayList<ArrayList<String>>[] allCombinations = (ArrayList<ArrayList<String>>[])new ArrayList[N];
			  int count=0;
			  maxDegree = 0;
			  //TEMPORARY CHECK
		//	  if (targets.size()>1) throw new Exception("Solutions are built with one target at a time"); // MOVE THIS
			  int targetUnderEvaluation = Integer.parseInt(targets.get(0)); 
			  //for (int i=0; i<targets.size();i++){
			  //		  int targetUnderEvaluation = Integer.parseInt(targets.get(i));
			  //INITIALIZE
			  ArrayList<String> temp = new ArrayList<String>(1);
			  for (int init=0; init<N; init++) {
					  temp.add(sources.get(init));
					  combinationsToEvaluate.add(new ArrayList<String>(temp));
					  temp.clear();
				  }

			  
			  for (int k=1; k<=N;k++){  //trascuro K=0
				//  System.out.println("K: "+k);
				  //System.out.println("binomial: "+binomial(N, k));
				  int numberOfEvaluations = combinationsToEvaluate.size();
				  				  
				  for (int j=0; j<numberOfEvaluations ;j++) { 
					  /*** WITH TARGET ***/
					  int [] itemsToMatch = new int[k+1]; // because of appending the target					  
					  for (int h=0; h<k ;h++) itemsToMatch[h] = Integer.parseInt(combinationsToEvaluate.get(j).get(h));
					  itemsToMatch[k] = targetUnderEvaluation;   //append the target
					  /*** WITHOUT TARGET *** 
					  int [] itemsToMatch = new int[k]; 					  
					  for (int h=0; h<k ;h++) itemsToMatch[h] = Integer.parseInt(combinationsToEvaluate.get(j).get(h));
					  */
					  //for (int h=0; h<=k ;h++) System.out.println(" itemsToMatch[h] "+h+" : "+ itemsToMatch[h]);
					  support = tree.getSupportOfItemset(itemsToMatch);
				//	  System.out.println(+j+"th combination Under Evaluation, k = "+k+" : "+combinationsToEvaluate.get(j));
				//	  System.out.println("support: "+support); 
					  if(support!=0){
						  count++;
						  maxDegree=k;
						  if(internal)
							  Main.supportMap.put(combinationsToEvaluate.get(j), support); // TO RETRIVE, USE comparison (or a three-entry map with "Entry")
						 
						  addUniqueParentsToEvaluate(combinationsToEvaluate.get(j), combinationsToEvaluate, combinationsToSkip,N, sources); 
						  	// add N-k parents (append k parents to the solution); CHECK FOR DUPLICATES
						  
			//			  System.out.println("support  "+support);
			//			  System.out.println("count  "+count);
			//			  System.out.println("combination "+combinationsToEvaluate.get(j));
						  
						  combinationsEvaluated.add(combinationsToEvaluate.get(j));
					  }
					  else{
			//			  System.out.println("\nsupport  "+support);
			//			  System.out.println("count  "+count);
			//			  System.out.println("combination "+combinationsToEvaluate.get(j));
						
						  addElementsToSkip(combinationsToEvaluate.get(j), combinationsToSkip, N, sources);
					  }
			  		}
				  for (int init=0; init<combinationsToEvaluate.size(); init++) {
					  if(combinationsToEvaluate.get(init).size()<=k) {
						  combinationsToEvaluate.remove(init);
						  init--;
					  }
				  }
			//	  System.out.println("combinations To Evaluate in the next step (of k-degree) "+combinationsToEvaluate.size());
			//	  for (int init=0; init<combinationsToEvaluate.size(); init++) 
			//		  System.out.println("combinations "+combinationsToEvaluate.get(init));
				  
				mergeCombToEvaluateToSkip(combinationsToEvaluate, combinationsToSkip); 
			//	System.out.println("combinations To Evaluate in the next step after merging (of k-degree) "+combinationsToEvaluate.size());
				  //combinationsToEvaluate.clear();
				//combinationsToSkip.clear();
			  }
			  System.out.println("\n ***+ Plausibility Evaluation: COUNT: "+count);
	//		  System.out.println("****Plausibility Evaluation: combinations Evaluated (combinations with Support): "+combinationsEvaluated.size()+"\n"); 
		//	  for (int a=0; a< combinationsEvaluated.size();a++ ) System.out.println("combinationsEvaluated  "+combinationsEvaluated.get(a));
		/*	  Iterator it = Main.supportMap.entrySet().iterator();
			  System.out.println("\n Map size: "+Main.supportMap.size());
			  while (it.hasNext()) {
				  Map.Entry pair = (Map.Entry)it.next();
				  System.out.println("Key=value: "+pair.getKey()+" = " + pair.getValue());
			  }
			  */
			  //ArrayList<String> listToSearch = new ArrayList<String>();
			  //listToSearch .add("2");
			  //listToSearch .add("13");
			  //System.out.println("\n return support of 2, 13: "+supportMap.get(listToSearch ));
			  //}
			  return ((double)count); 
	  }

	  
	  
 public double fast_rawPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  
		  myItemsetTree.ItemsetTree tree = null;  
		  /*if(internal)*/
			  tree = Main.itemsetTree;
		  /*else
			  tree = COPConfigurator.itemsetTree_external; 
			*/  
		  
	    	int support;
	    	Main.supportMap.clear();
	    	
//			ArrayList<String> sol_string = SolutionUtils.getNumericRepresentation(solution);  
			ArrayList<String> targets = new ArrayList<String>();
			targets.addAll(SolutionUtils.getTargets(solution));
			ArrayList<String> sources =  new ArrayList<String>();
			sources.addAll(SolutionUtils.getSources(solution));
			
			
			  
			System.out.println("\nPlausibility evaluation...");
			System.out.println("sources: " + sources);
			System.out.println("targets: " + targets);
			
			int N = sources.size();//SolutionUtils.getSolutionSize(example); 
			  
			ArrayList<ArrayList<String>>[] combinationsToSkip = (ArrayList<ArrayList<String>>[])new ArrayList[N];
			  //ArrayList<ArrayList<String>>[] allCombinations = (ArrayList<ArrayList<String>>[])new ArrayList[N];
			ArrayList<ArrayList<String>>[] combinationsToSkipSupportNull = (ArrayList<ArrayList<String>>[])new ArrayList[N];

			count = 0;
			maxDegree = 0;

			//intialize combinationsToSkip
			for (int init = 0 ; init < N; init ++) {
				combinationsToSkip[init] = new ArrayList<ArrayList<String>>();
				combinationsToSkip[init].add(new ArrayList<String>());
				combinationsToSkipSupportNull[init] = new ArrayList<ArrayList<String>>();
				combinationsToSkipSupportNull[init].add(new ArrayList<String>());
			}
			
			int targetUnderEvaluation = Integer.parseInt(targets.get(0)); 
				  
			boolean all_parents_have_null_support = true; 
			boolean already_counted = false;
			int with_supports_null=0;
			for (int i=0; i<sources.size(); i++) {
				ArrayList<String> combination = new ArrayList<String>(); 
				combination.add(sources.get(i)); 
				System.out.println("base root source under assessment "+sources.get(i));
				already_counted = false;
				for (int skip_index=0; skip_index<combinationsToSkip[0].size(); skip_index++ ) {   //Some of the recursive calls could update this, so it's not sure that it is empty 
						if(combinationsToSkip[0].get(skip_index).containsAll(combination)) {
							already_counted = true; 
							break;
						}
					}
					
				for (int skip_index=0; skip_index<combinationsToSkipSupportNull[0].size(); skip_index++ ) {    //( DIVIDERE IL CombToSkip VETTOER PER LIVELLI (k) per evitare di controllare sempre
					if(combinationsToSkipSupportNull[0].get(skip_index).containsAll(combination)  && combination.containsAll(combinationsToSkipSupportNull[0].get(skip_index)) ) {
							already_counted = true; 
							with_supports_null++;
							System.out.println("to skip support null : "+combination);
							break;
						}
				
				}
				 
				 
				  if(!already_counted) {
					  if(verifyParentSupport(combination, sources ,targetUnderEvaluation, tree, combinationsToSkip, combinationsToSkipSupportNull)!=0) { //Recursive function
						  all_parents_have_null_support = false; 
						  //break;
					  }
				  }
			  }
			  if (all_parents_have_null_support || with_supports_null==N) {
				  count = (int) (count + Math.pow(2,(0))-1); //CHECK. K=0; count = 0;
			  }
		
			  //mergeRemoveDuplicates();
			  for(int ind = combinationsToSkip.length-1 ; ind >= 0; ind--) {
				  System.out.println("\n\n Merging phase. Comb to skips: "+combinationsToSkip[ind]);
				  System.out.println("\n\n Merging phase. skips null: "+combinationsToSkipSupportNull[ind]);
				  // SDP  rule 
				 // sum of sets =  A1+ (not(A1)AND A2) + (not(A1, A2)AND A3) +
				//F1+F2 ..
					//	  where Fr  = CrAr, Cr = 1 all'inizio, e uguale a Cr-1 * Ar-1  al passo R
			  }
			  
			  ArrayList<ArrayList<String>> combinationsExamined = new ArrayList<ArrayList<String>>(); 
			  //ArrayList<Double> finalCount = new ArrayList<Double>(); 
			  double finalCount =0;
			  boolean first=true; 
			  for (int outer_index = combinationsToSkip.length-1;outer_index>=0; outer_index-- ) {
				  for (int inner_index = 0;inner_index <combinationsToSkip[outer_index].size(); inner_index++ ) {
					  if (combinationsToSkip[outer_index].get(inner_index).size()!=0) {
						  ArrayList<String> comb = combinationsToSkip[outer_index].get(inner_index);
						  if (first) {
							  //combinationsCleaned.add(comb); 
							  //finalCount.add(Math.pow(2, comb.size())-1); // POTREBBE SERVIRE SOLO una varibile
							  combinationsExamined.add(comb); 
							  finalCount = Math.pow(2, comb.size())-1; // POTREBBE SERVIRE SOLO una varibile
							  System.out.println("comb  "+comb);
							  System.out.println("first "+finalCount); 
							  
							  first = false; 
						  }
						  else {
							  finalCount = finalCount + (Math.pow(2, comb.size()) - 1);
							   for (int ind_examined =0; ind_examined<combinationsExamined.size();ind_examined++) {
								   //deterine the intersection
								   int n_intersections = intersection(comb, combinationsExamined.get(ind_examined)).size();
								   System.out.println("\n comb  "+comb);
									  System.out.println("comb examined "+combinationsExamined.get(ind_examined)); 
									  System.out.println("intersections "+n_intersections);  
									  
								   finalCount = finalCount - (Math.pow(2,n_intersections));
								   System.out.println("updated final count  "+finalCount);
							   }
							  combinationsExamined.add(comb); 
						  } 
				  }
			  }
			  }
			  
			  System.out.println("Count "+count);
			  System.out.println("Count updated: "+finalCount);
			  			  
			  System.exit(0);
			  return ((double)count); 
	  }
 
			  public ArrayList<String> intersection(ArrayList<String> list1, ArrayList<String> list2) {

					ArrayList<String> list = new ArrayList<String>();

				     for (String s : list1) {
				         if(list2.contains(s)) {
				             list.add(s);
				         }
				     }

				     return list;
				}

			  
	  private int verifyParentSupport(ArrayList<String> combination, ArrayList<String> sources, int targetUnderEvaluation, myItemsetTree.ItemsetTree tree, ArrayList<ArrayList<String>>[] combinationsToSkip, ArrayList<ArrayList<String>>[] combinationsToSkipSupportNull) {

		  		Collections.sort(combination);
				int N = sources.size();
				int k = combination.size();
				int [] itemsToMatch = new int[combination.size()+1]; // because of appending the target					  
				for (int h=0; h<combination.size() ;h++) itemsToMatch[h] = Integer.parseInt(combination.get(h));
				itemsToMatch[combination.size()] = targetUnderEvaluation;   //append the target
				int support = tree.getSupportOfItemset(itemsToMatch);
		//		if(parentCombination.size()==sources.size()) {
				//}
				System.out.println("\n\n Verify support: "+combination);
					
				System.out.println(" support  "+support);
				ArrayList<String> sourceMinusComb = new ArrayList<String>(sources);
				sourceMinusComb.removeAll(combination);
				 
				if(support!=0) {
					 boolean all_parents_have_null_support = true; 
					 boolean already_counted = false;
					 int with_supports_null = 0;
					 int parents_already_conuted = 0;
					for(int parents_index = 0; parents_index < (N - k);parents_index++) {
						
						ArrayList<String> parentCombination = new ArrayList<String>(combination);
						  // build parent
						parentCombination.add(sourceMinusComb.get(parents_index));
						Collections.sort(parentCombination);
						System.out.println("\nVerify parent suppot - comb parent: "+parentCombination);
						System.out.println("parent index: "+parents_index);
						System.out.println("k: "+k);
						already_counted=false;
						if(k+1<N) {
							outer: //verify
							for (int k1=k+1; k1<N; k1++) {
								for (int skip_index=0; skip_index<combinationsToSkip[k1].size(); skip_index++ ) {   
									if(combinationsToSkip[k1].get(skip_index).containsAll(parentCombination)) {
											already_counted = true; 
											parents_already_conuted++;
											System.out.println("already counted: "+parentCombination);
											//addElementsToSkipFast(parentCombination,combinationsToSkip, k, sources);
											break outer;
										}
								}
							}
						}
						
						
						for (int skip_index=0; skip_index<combinationsToSkipSupportNull[k].size(); skip_index++ ) {    //( DIVIDERE IL CombToSkip VETTOER PER LIVELLI (k) per evitare di controllare sempre
							if(combinationsToSkipSupportNull[k].get(skip_index).containsAll(parentCombination)  && parentCombination.containsAll(combinationsToSkipSupportNull[k].get(skip_index)) ) {
									already_counted = true; 
									with_supports_null++;
									System.out.println("to skip support null : "+parentCombination);
									break;
								}
						
						}
						
						
						if(!already_counted) { // continue the exploration of the tree
							  if(verifyParentSupport(parentCombination, sources ,targetUnderEvaluation, tree, combinationsToSkip, combinationsToSkipSupportNull)!=0) { //Recursive function
								  all_parents_have_null_support = false; 
								  //break;
							  }
						}
					  }
						
					  if ((all_parents_have_null_support || with_supports_null==N-k) && parents_already_conuted<(N-k)) {
						  System.out.println("\n CONTO: aggiungo "+(Math.pow(2,(k))-1));
						  System.out.println("combinzioe che port al conteggio "+combination);
						  if (with_supports_null==N-k)
							  System.out.println("contato perché tutti i parent erano in comb. to skip null ");
						  count = (int) (count + Math.pow(2,(k))-1); //CHECK
						  addElementsToSkipFast(combination,combinationsToSkip, k-1, sources);
						  
					  }
				}
				else {
					addElementsToSkipFast(combination, combinationsToSkipSupportNull, k-1, sources); 
				}
				System.out.println("combinationsToSkip: "+combinationsToSkip[k-1]);
				System.out.println("combinationsToSkipSupportNull: "+combinationsToSkipSupportNull[k-1]);
				System.out.println("at k : "+k);

				return support;
	  }

	public double rawRescaledPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  
		  int N = (SolutionUtils.getSources(solution).size());
		  //double x = fast_rawPlausibilityEvaluation(solution, internal); // 
		  double x =  rawPlausibilityEvaluation(solution, internal);
		  
		  
		 // System.out.println("x "+x);
		  //System.exit(0);
		  double scale_factor; 
		  if (N>average_KB_sol_size)
				scale_factor = average_KB_sol_size/N;
		  else
				scale_factor = N/average_KB_sol_size;
		  return x*scale_factor;
		  
		  
	
	}


	  public double complementedPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  // Divide plausibility by the number of "represented" solutions: (size/Maximum Degree) (or 2^(size/MD))
		  // this discourages the aggregation effect
		  int N = (SolutionUtils.getSources(solution).size());
		  double x = rawPlausibilityEvaluation(solution, internal);
		  System.out.println("N "+N);
		  System.out.println("rawPlausibilityEvaluation(solution, internal) "+x);
		  System.out.println("Max Degree "+maxDegree);
		  if(maxDegree!=0) {
			System.out.println("((double)N)/maxDegree) "+(double)N/maxDegree);
			System.out.println("Math.pow((((double)N)/maxDegree),2) "+Math.pow(2,(double)N/maxDegree));
			x = x/(Math.pow(2,(((double)N)/maxDegree)));
		  }
		  return x;
	
	
	  }
	  
	  
	  public double rescaledPoweredPlausibilityEvaluation(Solution solution, boolean internal) throws Exception {
		  int N = (SolutionUtils.getSources(solution).size());
		  //double x = fast_rawPlausibilityEvaluation(solution, internal); // 
		  double x =  rawPlausibilityEvaluation(solution, internal);
		  double scale_factor; 
		  if (N>average_KB_sol_size)
				scale_factor = average_KB_sol_size/N;
		  else
				scale_factor = N/average_KB_sol_size;
		  
		  scale_factor = Math.pow(scale_factor,Main.MAX_sources/average_KB_sol_size);
		  return x*scale_factor;
		  
		  
		  
	
	  }
	  
	  
	public double plausibilityEvaluation(Solution solution, boolean internal) throws Exception {
			double x = rawPlausibilityEvaluation(solution, internal);
			int N = (SolutionUtils.getSources(solution).size());
			return x/(Math.pow(2, N)-1);
	    }
	
		private static void addElementsToSkip(ArrayList<String> childCombination, 
				ArrayList<ArrayList<String>> combinationsToSkip, int N, ArrayList<String> sources) {
		
			ArrayList<String> newCombination = new ArrayList<String>(childCombination);
			//ArrayList<ArrayList<String>> parents = new ArrayList<ArrayList<String>>();
			boolean duplicate = false;
			for (int k1=0; k1< combinationsToSkip.size(); k1++ ) {
				if(combinationsComparison(newCombination, combinationsToSkip.get(k1))) {
					duplicate  = true; 
					break;
				}
			}
			if (!duplicate) combinationsToSkip.add(new ArrayList<String>(newCombination));
		}
		
		private static boolean addElementsToSkipFast(ArrayList<String> childCombination, 
				ArrayList<ArrayList<String>>[] combinationsToSkip, int k, ArrayList<String> sources) {
		
			ArrayList<String> newCombination = new ArrayList<String>(childCombination);
			boolean duplicate = false;
			for (int k1=0; k1< combinationsToSkip[k].size(); k1++ ) {
				if(combinationsComparison(newCombination, combinationsToSkip[k].get(k1))) {
					duplicate  = true; 
					System.out.println("Duplicate in add elements to skip. ");
					return false;
					//break;
				}
			}
			if (!duplicate) {
				combinationsToSkip[k].add(new ArrayList<String>(newCombination));
				return true;
			}
			return false;
		}
		
	private static void addUniqueParentsToEvaluate(ArrayList<String> childCombination, ArrayList<ArrayList<String>> combinationsToEvaluate,
				ArrayList<ArrayList<String>> combinationsToSkip, int N, ArrayList<String> sources) {
			
			
			ArrayList<String> newCombination = new ArrayList<String>(childCombination);
			ArrayList<ArrayList<String>> parents = new ArrayList<ArrayList<String>>();
			
			// build parents
			boolean duplicate = false;
			boolean toSkip = false; 
			for (int i=0; i<N; i++) {
				if (!newCombination.contains(sources.get(i))) {  
					duplicate = toSkip =  false; 
					newCombination.add(sources.get(i));
					
					for (int k1=0; k1< combinationsToEvaluate.size(); k1++ ) {
						if(combinationsComparison(newCombination, combinationsToEvaluate.get(k1))) { 
							duplicate = true;
							break;
			//				System.out.println("duplicate "+k1);
			//				for (int index=0; index<combinationsToEvaluate.get(k1).size(); index++)System.out.println(" "+combinationsToEvaluate.get(k1).get(index)); 
						}
					}
					for (int k2 = 0; k2 < combinationsToSkip.size(); k2++) {
						if (combinationsToSkip.get(k2).contains(sources.get(i))) {
							toSkip=true;
							break;
		//					System.out.println("TO SKIP "+newCombination);
							}
						//for (int ind2 = 0; ind2< combinationsToSkip.get(ind1).size(); ind2++) {
					}
					
					if (duplicate == false && toSkip == false) {
						parents.add(new ArrayList<String>(newCombination));
						combinationsToEvaluate.add(new ArrayList<String>(newCombination));
					}
				}
				newCombination = new ArrayList<String>(childCombination);
			}
		/*	for (int i=0; i<parents.size(); i++) {
				System.out.println("PARENT "+i);
				System.out.println(parents.get(i));
			//	for (int j=0; j<parents.get(i).size(); j++) {
			//		System.out.println(parents.get(i));
			//	}
			}
*/			
			
		}
	private static void mergeCombToEvaluateToSkip(ArrayList<ArrayList<String>> combinationsToEvaluate,
			ArrayList<ArrayList<String>> combinationsToSkip) {
	//	System.out.println("Combinations to evaluate "+combinationsToEvaluate);	
//		System.out.println("Combinations to evaluate size "+combinationsToEvaluate.size());
	//	System.out.println("Combinations to skip "+combinationsToSkip);
//		System.out.println("Combinations to skip size "+combinationsToSkip.size());
		for (int j=0; j<combinationsToSkip.size(); j++ ) {
			int i = 0;
			if (!combinationsToEvaluate.isEmpty()) {
				while (i<combinationsToEvaluate.size()) { 
					if(combinationsToEvaluate.get(i).containsAll(combinationsToSkip.get(j))) {
						combinationsToEvaluate.remove(i);
		//				System.out.println("Combinations to evaluate IN "+combinationsToEvaluate);
	//					System.out.println("Combinations to evaluate SIZE IN "+combinationsToEvaluate.size());
					}
					else {
						i++;
					}
					if (combinationsToEvaluate.isEmpty()) {
						break;
					}
				}
			}
		}
		
		/*
		 * int i=0;
		 * 
		 * while (i<combinationsToEvaluate.size()) { for (int j=0;
		 * j<combinationsToSkip.size(); j++ ) {
		 * if(combinationsToEvaluate.get(i).containsAll(combinationsToSkip.get(j)))
		 * {//.get(k1))) { combinationsToEvaluate.remove(i);
		 * System.out.println("Combinations to evaluate IN "+combinationsToEvaluate);
		 * System.out.println("Combinations to evaluate SIZE IN "+combinationsToEvaluate
		 * .size()); } } i++; if(combinationsToEvaluate.isEmpty()) break; }
		 */
		//System.out.println("Combinations to evaluate AFTER "+combinationsToEvaluate);	
		
		/*for (int i=0; i<combinationsToEvaluate.size(); i++ ) {
			for (int j=0; j<combinationsToSkip.size(); j++ ) {
					if(combinationsToEvaluate.get(j).containsAll(combinationsToSkip.get(j))) {//.get(k1))) {
						combinationsToEvaluate.remove(j);
						i--;
						System.out.println("Merged, index "+i);
						System.out.println("Size of comb to ev. "+combinationsToEvaluate.size());
					}
			}
		}*/
			
	}
	private static boolean combinationsComparison(ArrayList<String> combination1, ArrayList<String> combination2) {
		// return true if one contains the same elements (in whatever order) than the other. Assume no duplicates in a list
		int count = 0;
//		int min;
		if (combination1.size() != combination2.size()) {
			//System.out.println("Not comparable, they are of different sizes "); 
			return false;
		}
//			return false; min = combination1.size();
//		else
//			min  = combination2.size();
		
		for (int i=0; i<combination1.size();i++) { 
			for (int j=0; j<combination2.size();j++) {
				if(combination1.get(i).equals(combination2.get(j))) 
					count++;
			}
		}
		if (count==combination1.size())
			return true;
		else
			return false; 
	}
/*
	public Integer getMaxKDegreeWithAutoExclusion(ArrayList<String> elements, boolean internal, int minSupport) {
		
		// the elements are excluded from the counting: in other words, a support > 1 is considered instesad than >0
		
		 myItemsetTree.ItemsetTree tree = null;
		 if(internal)
			 tree = COPConfigurator.itemsetTree;
		 else
			 tree = COPConfigurator.itemsetTree_external;
			
		 //tree.printTree();
		 
		  int N = elements.size(); 
//		  System.out.println("MaxKDegree computation: size of the group of elements: "+N);
		  
		  ArrayList<ArrayList<String>> combinationsEvaluated = new ArrayList<ArrayList<String>>(); 
		  ArrayList<ArrayList<String>> combinationsToEvaluate = new ArrayList<ArrayList<String>>();
		  ArrayList<ArrayList<String>> combinationsToSkip = new ArrayList<ArrayList<String>>();
		  ArrayList<ArrayList<String>>[] allCombinations = (ArrayList<ArrayList<String>>[])new ArrayList[N];
		  int count=0;
	
		  ArrayList<String> temp = new ArrayList<String>(1);
		  for (int init=0; init<N; init++) {
				  temp.add(elements.get(init));
				  combinationsToEvaluate.add(new ArrayList<String>(temp));
				  temp.clear();
			  }
		  
		  int maxk_degree  = N; 
		  for (int k=1; k < N; k++){  //trascuro K=0
			  count = 0; 
			  int numberOfEvaluations = combinationsToEvaluate.size();
			  for (int j=0; j<numberOfEvaluations ;j++) { 
				  int [] itemsToMatch = new int[k];
				  for (int h=0; h<k ;h++) itemsToMatch[h] = Integer.parseInt(combinationsToEvaluate.get(j).get(h)); 
				  //for (int h=0; h<k ;h++) System.out.println(" itemsToMatch[h] "+ itemsToMatch[h]);				  
				  int support = tree.getSupportOfItemset(itemsToMatch);
				  if(support>minSupport){ // so as to exclude the "elements" themselves from the counting
					  count++;
					  addUniqueParentsToEvaluate(combinationsToEvaluate.get(j), combinationsToEvaluate, combinationsToSkip,N, elements); 
					  combinationsEvaluated.add(combinationsToEvaluate.get(j));
				  }
				  else{
					  addElementsToSkip(combinationsToEvaluate.get(j), combinationsToSkip, N, elements);
				  }
		  		}
			  for (int init=0; init<combinationsToEvaluate.size(); init++) {
				  if(combinationsToEvaluate.get(init).size()<=k) {
					  combinationsToEvaluate.remove(init);
					  init--;
				  }
			  }
			mergeCombToEvaluateToSkip(combinationsToEvaluate, combinationsToSkip); 
			//System.out.println("combinations To Evaluate in the next step after merging (of k-degree)"+combinationsToEvaluate.size());
			  //combinationsToEvaluate.clear();
			//combinationsToSkip.clear();
			if (count == 0) {
				maxk_degree = k; 
				return maxk_degree;   
			}
		  }	
		  
		return maxk_degree;
	}
	*/
	
	
}
