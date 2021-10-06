package solutionType;

import java.util.ArrayList;

import org.apache.commons.lang3.RandomUtils;

import causalOptimization.Main;
import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;
import jmetal.core.Problem;
import jmetal.core.SolutionType;
import jmetal.core.Variable;

public class CauseEffectSolutionType extends SolutionType {

	public ArrayList<String> _sources; 
	public ArrayList<String> _targets;
	
	public int numberOfVariables;    
	private int numberOfSources; 
	private int numberOfTargets;
	
	private boolean internal_domain; //true:internal; false:external 
	
	private String OPERATOR; //operator that generated the solution
	
	public CauseEffectSolutionType(Problem problem){
		super(problem);
		setInternal_domain(true);
		OPERATOR = "";
	}
	
	public CauseEffectSolutionType(Problem problem, int n_sources, int n_targets){
		super(problem); 
		numberOfSources = n_sources;
		numberOfTargets= n_targets;
		setInternal_domain(true);
	}

	public String getOperator(){
		return OPERATOR;
	}
	public void setOperator(String _operator){
		OPERATOR = _operator;
	}
	
	@Override
	public Variable[] createVariables() throws ClassNotFoundException { 
		Variable [] variables = new Variable[2];

		//ArrayList<String> sources = new ArrayList<String>();
		int sourceSize = Main.MIN_sources + Main.ran.nextInt(Main.MAX_sources-Main.MIN_sources+1);
		
/*		ArrayList<Integer> n_values_index = new ArrayList<Integer>(causeEffect.Main.possibleSources.size());
		for (int ind = 0; ind < causeEffect.Main.possibleSources.size(); ind++) 
			n_values_index.add(ind);
		Collections.shuffle(n_values_index);	
		for (int k = 0; k< sourceSize; k++)
			sources.add(causeEffect.Main.possibleSources.get(n_values_index.get(k)));
		variables[0] = new Sources(sources); // COPY CONSTRUCTOR
*/		
		try {
			variables[0] = new Sources(sourceSize);
		} catch (Exception e) {e.printStackTrace();} 
			
		//	ArrayList<String> targets = new ArrayList<String>();
		int targetSize = Main.MIN_targets + Main.ran.nextInt(Main.MAX_targets - Main.MIN_targets +1 );
		/*
		ArrayList<Integer> n_values_index_targets = new ArrayList<Integer>(causeEffect.Main.possibleTargets.size());
		for (int ind = 0; ind < causeEffect.Main.possibleTargets.size(); ind++) 
			n_values_index_targets.add(ind);
		Collections.shuffle(n_values_index_targets);
		for (int k = 0; k< targetSize; k++) {
			targets.add(causeEffect.Main.possibleTargets.get(n_values_index_targets.get(k)));
			
		int existing_sol_index;  
		// SELECT FROM A KNOWN TARGET. IF NOT ABDUCTIVe, THIS COULD BE CHANGED. 	
		for (int k = 0; k< targetSize; k++) {
			existing_sol_index = RandomUtils.nextInt(0,Main.listOfInternalSolutions.size());
			String targetToAdd = ((Targets)(Main.listOfInternalSolutions.get(existing_sol_index).getDecisionVariables()[1]))._targetsList.get(k);
			targets.add(targetToAdd);
		}
		 */
		try {
			variables[1] = new Targets(targetSize);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return variables ;
	}

	public boolean isInternal_domain() {
		return internal_domain;
	}

	public void setInternal_domain(boolean internal_domain) {
		this.internal_domain = internal_domain;
	}
		
}
