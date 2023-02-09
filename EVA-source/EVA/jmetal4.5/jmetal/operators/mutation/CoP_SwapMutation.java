//  SwapMutation.java
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
import jmetal.encodings.solutionType.PermutationSolutionType;
import jmetal.encodings.variable.Permutation;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.lang3.RandomUtils;

import causalOptimization.Main;
import causalOptimization.variable.Sources;

/**
 * This class implements a swap mutation. Swap Mutation with an element of the Ontology
 */
public class CoP_SwapMutation extends Mutation{
  /**
   * Valid solution types to apply this operator 
   */
  private static final List VALID_TYPES = Arrays.asList(solutionType.CauseEffectSolutionType.class) ;
  
  private Double mutationProbability_ = null ;

  /** 
   * Constructor
   */
  public CoP_SwapMutation(HashMap<String, Object> parameters) {    
  	super(parameters) ;
  	
  	if (parameters.get("probability") != null)
  		mutationProbability_ = (Double) parameters.get("probability") ;  		
  } // Constructor


  /**
   * Constructor
   */
  //public SwapMutation(Properties properties) {
  //  this();
  //} // Constructor

  /**
   * Performs the operation
   * @param probability Mutation probability
   * @param solution The solution to mutate
   * @throws JMException 
   */
  public void doMutation(double probability, Solution solution) throws JMException {   
    String sourceList[];
    int sourceListLength  ;
	    if (solution.getType().getClass() == solutionType.CauseEffectSolutionType.class) {

	      sourceListLength = ((Sources)(solution.getDecisionVariables()[0]))._sourcesList.size();
	    //  permutation = ((Permutation)solution.getDecisionVariables()[0]).vector_ ;
	      sourceList = new String[sourceListLength]; 
	      
	      for (int i=0; i<sourceListLength; i++) {
	    	sourceList[i] = ((Sources)(solution.getDecisionVariables()[0]))._sourcesList.get(i);
		//	System.out.println("source list "+((Sources)(solution.getDecisionVariables()[0]))._sourcesList.get(i));
	      }

	      if (PseudoRandom.randDouble() < probability) {
	        int pos ;
	        //int pos2 ;

	        pos = PseudoRandom.randInt(0,sourceListLength-1) ;
	        //pos2 = PseudoRandom.randInt(0,sourceListLength-1) ;

	       /* while (pos1 == pos2) {
	          if (pos1 == (sourceListLength - 1)) 
	            pos2 = PseudoRandom.randInt(0, sourceListLength- 2);
	          else 
	            pos2 = PseudoRandom.randInt(pos1, sourceListLength- 1);
	        } // while
	        // swap
	         * 
	         */
	        String sourceToAdd; 
		    do {
		    	sourceToAdd = Main.possibleSources.get(Main.ran.nextInt(Main.possibleSources.size()));
		    }
		    while (((Sources)(solution.getDecisionVariables()[0]))._sourcesList.contains(sourceToAdd));
		    
	        sourceList[pos] = sourceToAdd; //COPConfigurator.possibleSources.get(RandomUtils.nextInt(0,COPConfigurator.possibleSources.size())); 
//	        int temp = permutation[pos1];
//	        permutation[pos1] = permutation[pos2];/
	        //permutation[pos2] = temp;    
	      } // if
	    } // if
	    else  {
	      Configuration.logger_.severe("SwapMutation.doMutation: invalid type. " +
	          ""+ solution.getDecisionVariables()[0].getVariableType());

	      Class cls = java.lang.String.class;
	      String name = cls.getName(); 
	      throw new JMException("Exception in " + name + ".doMutation()") ;
	    }
	    
	    for (int i=0; i<sourceListLength; i++) 
	    	((Sources)(solution.getDecisionVariables()[0]))._sourcesList.set(i,sourceList[i]);
	    
	    
	    

  } // doMutation

  /**
   * Executes the operation
   * @param object An object containing the solution to mutate
   * @return an object containing the mutated solution
   * @throws JMException 
   */
  public Object execute(Object object) throws JMException {
    Solution solution = (Solution)object;
    
		if (!VALID_TYPES.contains(solution.getType().getClass())) {
			Configuration.logger_.severe("SwapMutation.execute: the solution " +
					"is not of the right type. The type should be 'Binary', " +
					"'BinaryReal' or 'Int', but " + solution.getType() + " is obtained");

			Class cls = java.lang.String.class;
			String name = cls.getName();
			throw new JMException("Exception in " + name + ".execute()");
		} // if 

    
    this.doMutation(mutationProbability_, solution);
    return solution;
  } // execute  
} // SwapMutation
