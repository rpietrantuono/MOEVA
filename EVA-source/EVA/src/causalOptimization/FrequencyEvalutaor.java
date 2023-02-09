package causalOptimization;

import java.io.IOException;
import java.util.PriorityQueue;

import myItemsetTree.ItemsetTree;
import util.MyAlgoTopKClassRules;
import ca.pfv.spmf.algorithms.associationrules.TopKRules_and_TNR.ClassRuleG;
import ca.pfv.spmf.algorithms.associationrules.TopKRules_and_TNR.Database;
import jmetal.core.Solution;


public class FrequencyEvalutaor {

	public ItemsetTree buildItemSetTree(String filepath){
		
		// Applying the algorithm to build the itemset tree
		ItemsetTree itemsetTree = new ItemsetTree();
		// method to construct the tree from a set of transactions in a file
		try {
			itemsetTree.buildTree(filepath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// print the statistics about the tree construction time and print the tree in the console
		itemsetTree.printStatistics();
		//itemsetTree.printTree();
		return itemsetTree;
	}
	
	
	
	public PriorityQueue<ClassRuleG> minFixedConsequentRules(Solution sol, int k, Database database, double minConfidence, int[] itemToBeUsedAsConsequent, Integer... maxAntecedentSize) throws IOException{
		
		// TOP K ASSOCIATION RULES WITH FIXED CONSEQUENT 
	    Integer value = maxAntecedentSize.length > 0 ? maxAntecedentSize[0] : -1;
	    int maxAntSize = value.intValue();
	    
		MyAlgoTopKClassRules algo = new MyAlgoTopKClassRules();
		
		if(maxAntSize!=-1)
			algo.setMaxAntecedentSize(maxAntSize);  // optional
		
		algo.runAlgorithm(k, minConfidence, database, itemToBeUsedAsConsequent);
		algo.printStats();
		algo.writeResultTofile(Main.baseDir+"/output.txt");   // to save results to file
		return algo.getKRules(); 
	
				
	}
	
	
	
	
}
