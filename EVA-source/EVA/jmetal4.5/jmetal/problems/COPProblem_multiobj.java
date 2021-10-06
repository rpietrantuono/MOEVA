
package jmetal.problems;

import java.util.ArrayList;

import causalOptimization.Main;
import causalOptimization.CopVariable.Sources;
import causalOptimization.CopVariable.Targets;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.Variable;
import jmetal.util.JMException;
import solutionType.CauseEffectSolutionType;


public class COPProblem_multiobj extends Problem {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	public COPProblem_multiobj(String solutionType){
	
		numberOfVariables_   = 2;  
	    numberOfObjectives_  = 2; 
	    numberOfConstraints_ = 0;  // These are possible inconsistencies not managed in the solution construction phase by the operators
	    
	    problemName_ = "COPProblem";
	    solutionType_ = new CauseEffectSolutionType(this) ;
	}
	
	  /** 
	  * Evaluates a solution 
	  * @param solution The solution to evaluate
	   * @throws JMException 
	  */        
	
	
		  public void evaluate(Solution solution) throws JMException {
		 
		  Variable[] decisionVariables = new Variable[this.numberOfVariables_];
		  decisionVariables[0] = (Sources)(solution.getDecisionVariables()[0]);
		  decisionVariables[1] = (Targets)(solution.getDecisionVariables()[1]);
		  
		  double [] f = new double[numberOfObjectives_];
		  // For external solutions (analogical abduction), we assume the novelty does not matter (it is always 1, i.e. diverse) and plausibility is always 1 (as we consdier teh entire soltion) 
		  if (!Main.technique.equals("COP")) {
			  try {
				  if(!(solution.isInternal)){
					  solution.setObjective(0, 1); //Main.plaus.plausibilityEvaluation(solution, false));  
					  solution.setObjective(1, 1);
					  return;
				  }
			  }catch (Exception e) {e.printStackTrace();}
		  }
		  
		  /* PALUSIBILITY EVALUATION */
		  double plausibility =0;
		  double novelty = 0;  
		  ArrayList<Double> similarities = null; 
		  try {
			  //plausibility = Main.plaus.plausibilityEvaluation(solution, true);
			  plausibility = Main.plaus.rootPlausibilityEvaluation(solution, true); 
			  //plausibility = Main.plaus.logarithmicPlausibilityEvaluation(solution, true);
			  similarities = Main.sim.computeKBSimilarity(solution);
		  }catch (Exception e) { e.printStackTrace();	} 
		  double total_sim =0; 
		  for (int i =0; i<similarities.size();i++) total_sim  += similarities.get(i); 
		  double average_sim = total_sim/similarities.size();  // in 0 e 1
		  novelty = 1 - average_sim;
		  solution.setObjective(0, -plausibility);  
		  solution.setObjective(1, -novelty);
	}
	 
		
			public void evaluateConstraints(Solution solution) throws JMException {
				
				
			/*	ArrayList<Double> similarities = null;
				  try {
					  similarities = Main.sim.computeKBSimilarity(solution);
				  }catch (Exception e) { e.printStackTrace();	} 
				  double total_sim =0; 
				  for (int i =0; i<similarities.size();i++) total_sim  += similarities.get(i); 
				  double average_sim = total_sim/similarities.size();  // in 0 e 1
				  double novelty = 1 - average_sim;
				  double desiredLimit = Main.DesiredLimit; 
					
					double [] constraint = new double[this.getNumberOfConstraints()];
					constraint[0] = novelty-desiredLimit;
					
					double total = 0.0;
				    int number = 0;
				    for (int i = 0; i < this.getNumberOfConstraints(); i++)
				      if (constraint[i]<0.0){
				        total+=constraint[i];
				        number++;
				      }
				    if (total!=0) System.out.println("constraint violated "+total);
				    
				    solution.setOverallConstraintViolation(total);    
				    solution.setNumberOfViolatedConstraint(number);
				    */
				  /*
				
				// for our probelm number of constraints is 0, so this method does nothing
		    	Variable[] decisionVariables = new Variable[this.numberOfVariables_];
		    	double [] constraint = new double[this.getNumberOfConstraints()];
				 
				decisionVariables[0] = (Sources)(solution.getDecisionVariables()[0]);			
				decisionVariables[1] = (Targets)(solution.getDecisionVariables()[1]);
				
				ArrayList<String> componentsString = Main.getKeyContainingName("COMPONENT.NAME");  
				ArrayList<String> componentsProblemString = Main.getKeyContainingName("COMPONENT.PROBLEM"); 
				
				double total = 0.0;
				int number = 0;
				//CONSTRAINT 1: for each component name there should be one component problem and vice-versa
				ArrayList<String> sources = SolutionUtils.getSources(solution); 
				int numberOfComponentNames = 0;
				int numberOfComponentProblems = 0;
				for (int i = 0; i < componentsString.size(); i++) 
					if (sources.contains(componentsString.get(i)))
						numberOfComponentNames++;
				for (int i = 0; i < componentsProblemString.size(); i++) 
					if (sources.contains(componentsProblemString.get(i)))
						numberOfComponentProblems++;
				
				int violations = Math.abs(numberOfComponentNames - numberOfComponentProblems);  
				
				System.out.println("sources "+sources);
				System.out.println("number of component names "+numberOfComponentNames);
				System.out.println("number of component problems "+numberOfComponentProblems);
				
				if (violations>0){
					constraint[0] = -violations ;
					total+=constraint[0];		
					number++;
					System.out.println("constraint violation "+constraint[0]);
					System.out.println("total "+total);
					System.out.println("number "+number);
				}
				  solution.setOverallConstraintViolation(total);    
				    solution.setNumberOfViolatedConstraint(number);
				  	
				//OTHERS ...?
			*/	
			    /*for (int i = 0; i < this.getNumberOfConstraints(); i++) {
				    	if (constraint[i]<0.0){
				    		total+=constraint[i];
				    		number++;
				    	}
			    } */
				
			  
			    /* UNCONSTRAINED*/ 
			    solution.setOverallConstraintViolation(0);    
			    solution.setNumberOfViolatedConstraint(0);
		             
		  } // evaluateConstraints

	    

	  
}
