//  NonUniformMutation.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.operators.mutation;

import jmetal.core.Solution;
import jmetal.encodings.solutionType.ArrayRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.wrapper.XReal;
import solutionType.CauseEffectSolutionType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;

import org.apache.commons.lang3.RandomUtils;

import causalOptimization.Main;
import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;

/**
  * This class implements a non-uniform mutation operator.
  */
public class CoP_NonUniformMutation extends Mutation{
  /**
   * Valid solution types to apply this operator 
   */
	private static final List VALID_TYPES = Arrays.asList(CauseEffectSolutionType.class) ;
	/**
   * perturbation_ stores the perturbation value used in the Non Uniform 
   * mutation operator
   */
  private Double perturbation_ = null;
  
  /**
   * maxIterations_ stores the maximun number of iterations. 
   */
  private Integer maxIterations_ = null;    
  
  /**
   * currentIteration_ stores the iteration in which the operator is going to be
   * applied
   */
  private Integer currentIteration_ = null;
           
  private Double mutationProbability_ = null;

  /** 
  * Constructor
  * Creates a new instance of the non uniform mutation
  */
  public CoP_NonUniformMutation(HashMap<String, Object> parameters) {
  	super(parameters) ;
  	if (parameters.get("probability") != null)
  		mutationProbability_ = (Double) parameters.get("probability") ;  		
  	if (parameters.get("perturbation") != null)
  		perturbation_ = (Double) parameters.get("perturbation") ;  		
  	if (parameters.get("maxIterations") != null)
  		maxIterations_ = (Integer) parameters.get("maxIterations") ;  		
  } // NonUniformMutation
            

  /**
  * Constructor
  * Creates a new instance of the non uniform mutation
  */
  //public NonUniformMutation(Properties properties){
  //   this();
  //} // NonUniformMutation


  /**
  * Perform the mutation operation
  * @param probability Mutation probability
  * @param solution The solution to mutate
   * @throws JMException 
  */
  public void doMutation(double probability, Solution solution) throws JMException {                
	  
	  
	  double lowerSourceLimit = 0; 
		 double upperSourceLimit = ((double)Main.possibleSources.size())/(Main.possibleSources.size()+Main.possibleTargets.size());
		 double lowerTargetLimit = ((double)Main.possibleSources.size())/(Main.possibleSources.size()+Main.possibleTargets.size())+ 0.5/(Main.possibleTargets.size()+Main.possibleTargets.size());
		 double upperTargetLimit = 1;
		
	//	 System.out.println("Upper source limit "+upperSourceLimit); 
	//	 System.out.println("low target limit "+lowerTargetLimit);
		
		  for (int j = 0; j < causalOptimization.Main.MAX_sources; j++){ //REPLACE WITH "NUMBER OF GENES"
			  if (PseudoRandom.randDouble() < probability) {
				  double rand = PseudoRandom.randDouble();
				  double tmp ; 
			  	  double current_value = ((Sources)(solution.getDecisionVariables()[0]))._doubleValueList.get(j);
				  if (rand <= 0.5) {
			          tmp = delta(upperSourceLimit - current_value,
			                      perturbation_.doubleValue());
			          tmp += current_value;
			        }
				  else {
			          tmp = delta(lowerSourceLimit - current_value,
			                      perturbation_.doubleValue());
			          tmp += current_value;
				  }
				  
				String stringRepresentation; 
			  	if (tmp < lowerSourceLimit){
		  		  ((Sources)(solution.getDecisionVariables()[0]))._doubleValueList.set(j, lowerSourceLimit);
		  		  ((Sources)(solution.getDecisionVariables()[0]))._sourcesFullList.set(j, "0");
		  		  System.out.println("<0 ");
		  		  continue;
			  	}
			  	if (tmp > upperSourceLimit){
			  		((Sources)(solution.getDecisionVariables()[0]))._doubleValueList.set(j, upperSourceLimit);
			  		((Sources)(solution.getDecisionVariables()[0]))._sourcesFullList.set(j, String.valueOf(Main.possibleSources.size()));
			  		System.out.println(">limit"
			  				+ " ");
		  		  continue;
			  	}
			 
			  	try {
					stringRepresentation = util.SolutionUtils.getStringRepresentation(tmp);
					((Sources)(solution.getDecisionVariables()[0]))._sourcesFullList.set(j, stringRepresentation);
					
				} catch (Exception e) {e.printStackTrace();}
			  }
		  }
		  try {
			((Sources)(solution.getDecisionVariables()[0])).updateSource();
		} catch (Exception e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		  for (int j = 0; j < causalOptimization.Main.MAX_targets; j++){ //REPLACE WITH "NUMBER OF GENES"
			  if (PseudoRandom.randDouble() < probability) {
				  double rand = PseudoRandom.randDouble();
				  double tmp; 
				  double current_value = ((Targets)(solution.getDecisionVariables()[1]))._doubleValueList.get(j);
				  if (rand <= 0.5) {
			          tmp = delta(upperTargetLimit - current_value,
			                      perturbation_.doubleValue());
			          tmp += current_value;
			        }
				  else {
			          tmp = delta(lowerTargetLimit - current_value,
			                      perturbation_.doubleValue());
			          tmp += current_value;
				  }
				  String stringRepresentation; 
			  	  if (tmp < lowerTargetLimit){
			  		  ((Targets)(solution.getDecisionVariables()[1]))._doubleValueList.set(j, lowerTargetLimit);
			  		  ((Targets)(solution.getDecisionVariables()[1]))._targetsFullList.set(j, "0"); 
			  		  continue;
			  	  }
			  	  if (tmp > upperTargetLimit){
			  		  ((Targets)(solution.getDecisionVariables()[1]))._doubleValueList.set(j, upperTargetLimit);
			  		((Targets)(solution.getDecisionVariables()[1]))._targetsFullList.set(j, String.valueOf(Main.possibleSources.size()+Main.possibleTargets.size()));
			  		  continue;
			  	  }
			  	  try {
					stringRepresentation = util.SolutionUtils.getStringRepresentation(tmp);
					((Targets)(solution.getDecisionVariables()[1]))._targetsFullList.set(j, stringRepresentation);
			  	  } catch (Exception e) {e.printStackTrace();}
			  }
		  }
		  try {
			((Targets)(solution.getDecisionVariables()[1])).updateTarget();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
		      
	  // ADJUST IN CASE OF DUPLICATES
			// Remove duplicates if any
		  
		  LinkedHashSet<String> set = new LinkedHashSet<String>(((Sources)(solution.getDecisionVariables()[0]))._sourcesList);
		  ((Sources)(solution.getDecisionVariables()[0]))._sourcesList.clear();
			 ((Sources)(solution.getDecisionVariables()[0]))._sourcesList.addAll(set);

/*		  
			LinkedHashSet<String> hs1 = new LinkedHashSet<String>();
			for (int i=0; i< ((Sources)(solution.getDecisionVariables()[0]))._sourcesList.size() ; i++)
				hs1.add(((Sources)(solution.getDecisionVariables()[0]))._sourcesList.get(i));
		
		 ((Sources)(solution.getDecisionVariables()[0]))._sourcesList.clear();
		  Iterator it1 = hs1.iterator();
		
			while (it1.hasNext()) 
		    	((Sources)(solution.getDecisionVariables()[0]))._sourcesList.add(String.valueOf(it1.next()));
*/			
		    while (((Sources)(solution.getDecisionVariables()[0]))._sourcesList.size()<Main.MIN_sources) {
		    	String sourceToAdd; 
		    	do {
		    		sourceToAdd = Main.possibleSources.get(Main.ran.nextInt(Main.possibleSources.size()));
		    	}
		    	while (((Sources)(solution.getDecisionVariables()[0]))._sourcesList.contains(sourceToAdd));
		    	((Sources)(solution.getDecisionVariables()[0]))._sourcesList.add(sourceToAdd);
		    }	    
		    
		    // FOR THE TARGET IS NOT NEEDED IF MAX_TARGET=1. IF NOT, THIS SHOUlD BE MANAGED
		    
  } // doMutation
    

  /**
   * Calculates the delta value used in NonUniform mutation operator
   */
  private double delta(double y, double bMutationParameter) {
    double rand = PseudoRandom.randDouble();
    int it,maxIt;
    it    = currentIteration_.intValue();
    maxIt = maxIterations_.intValue();
        
    return (y * (1.0 - 
                Math.pow(rand,
                         Math.pow((1.0 - it /(double) maxIt),bMutationParameter)
                         )));
  } // delta

  /**
  * Executes the operation
  * @param object An object containing a solution
  * @return An object containing the mutated solution
   * @throws JMException 
  */
  public Object execute(Object object) throws JMException {
    Solution solution = (Solution)object;
    
		if (!VALID_TYPES.contains(solution.getType().getClass())) {
      Configuration.logger_.severe("NonUniformMutation.execute: the solution " +
      		solution.getType() + "is not of the right type");

      Class cls = java.lang.String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".execute()") ;
    } // if  
    
  	if (getParameter("currentIteration") != null)
  		currentIteration_ = (Integer) getParameter("currentIteration") ;  		

    doMutation(mutationProbability_,solution);
        
    return solution;    
  } // execute
} // NonUniformMutation
