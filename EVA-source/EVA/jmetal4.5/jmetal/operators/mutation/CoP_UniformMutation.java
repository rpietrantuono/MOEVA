//  UniformMutation.java
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
import util.SolutionUtils;

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
 * This class implements a uniform mutation operator.
 */
public class CoP_UniformMutation extends Mutation{
  /**
   * Valid solution types to apply this operator 
   */
  private static final List VALID_TYPES = Arrays.asList(CauseEffectSolutionType.class) ;
  /**
   * Stores the value used in a uniform mutation operator
   */
  private Double perturbation_;
  
  private Double mutationProbability_ = null;

  /** 
   * Constructor
   * Creates a new uniform mutation operator instance
   */
  public CoP_UniformMutation(HashMap<String, Object> parameters) {
  	super(parameters) ;
  	
  	if (parameters.get("probability") != null)
  		mutationProbability_ = (Double) parameters.get("probability") ;  		
  	if (parameters.get("perturbation") != null)
  		perturbation_ = (Double) parameters.get("perturbation") ;  		

  } // UniformMutation


  /**
   * Constructor
   * Creates a new uniform mutation operator instance
   */
  //public UniformMutation(Properties properties) {
  //  this();
  //} // UniformMutation


  /**
  * Performs the operation
  * @param probability Mutation probability
  * @param solution The solution to mutate
   * @throws JMException 
  */
  public void doMutation(double probability, Solution solution) throws JMException {  
  	//XReal x = new XReal(solution) ; 

	 
	 double lowerSourceLimit = 0; 
	 double upperSourceLimit = ((double)Main.possibleSources.size())/(Main.possibleSources.size()+Main.possibleTargets.size());
	 double lowerTargetLimit = ((double)Main.possibleSources.size())/(Main.possibleSources.size()+Main.possibleTargets.size())+ 0.5/(Main.possibleTargets.size()+Main.possibleTargets.size());
	 double upperTargetLimit = 1;
	
	 
	  for (int j = 0; j < causalOptimization.Main.MAX_sources; j++){ //REPLACE WITH "NUMBER OF GENES"
		  if (PseudoRandom.randDouble() < probability) {
			  double rand = PseudoRandom.randDouble();
			  double tmp = (rand - 0.5)*perturbation_.doubleValue();
			  String stringRepresentation; 
		  	  double current_value = ((Sources)(solution.getDecisionVariables()[0]))._doubleValueList.get(j);
		  	tmp += current_value;
		  	if (tmp < lowerSourceLimit){
	  		  ((Sources)(solution.getDecisionVariables()[0]))._doubleValueList.set(j, lowerSourceLimit);
	  		((Sources)(solution.getDecisionVariables()[0]))._sourcesFullList.set(j, "0");
	  		continue;
	  	  }
	  	  if (tmp > upperSourceLimit){
	  		  ((Sources)(solution.getDecisionVariables()[0]))._doubleValueList.set(j, upperSourceLimit);
	  		  ((Sources)(solution.getDecisionVariables()[0]))._sourcesFullList.set(j, String.valueOf(Main.possibleSources.size()));
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
			  double tmp = (rand - 0.5)*perturbation_.doubleValue();
			  String stringRepresentation; 
		  	  double current_value = ((Targets)(solution.getDecisionVariables()[1]))._doubleValueList.get(j);
		  	  tmp += current_value;
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

	  /*LinkedHashSet<String> hs1 = new LinkedHashSet<String>();
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
	    
	
  } // doMutation
  
  /**
  * Executes the operation
  * @param object An object containing the solution to mutate
   * @throws JMException 
  */
  public Object execute(Object object) throws JMException {
    Solution solution = (Solution )object;
    
		if (!VALID_TYPES.contains(solution.getType().getClass())) {
      Configuration.logger_.severe("UniformMutation.execute: the solution " +
          "is not of the right type. The type should be 'Real', but " +
          solution.getType() + " is obtained");

      Class cls = java.lang.String.class;
      String name = cls.getName(); 
      throw new JMException("Exception in " + name + ".execute()") ;
    } // if 
    
    doMutation(mutationProbability_,solution);
        
    return solution;
  } // execute                  
} // UniformMutation
