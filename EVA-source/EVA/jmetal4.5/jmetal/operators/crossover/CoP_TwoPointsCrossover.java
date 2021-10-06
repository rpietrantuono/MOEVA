//  TwoPointsCrossover.java
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

package jmetal.operators.crossover;

import jmetal.core.Solution;
import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;

import causalOptimization.Main;
import causalOptimization.variable.Sources;
import causalOptimization.variable.Targets;

/**
 * This class allows to apply a two points crossover operator using two parent
 * solutions. 
 * NOTE: the type of the solutions must be Causal (CauseEffect)
 * The operator is applied only to SOURCE (and not to TARGET) for this version. 
 */
public class CoP_TwoPointsCrossover extends Crossover {

  /**
   * Valid solution types to apply this operator 
   */
  private static final List VALID_TYPES = Arrays.asList(solutionType.CauseEffectSolutionType.class) ;

  private Double crossoverProbability_ = null;

	/**
	 * Constructor
	 * Creates a new intance of the two point crossover operator
	 */
	public CoP_TwoPointsCrossover(HashMap<String, Object> parameters) {
		super(parameters) ;
		
  	if (parameters.get("probability") != null)
  		crossoverProbability_ = (Double) parameters.get("probability") ;  		
	} // TwoPointsCrossover


	/**
	 * Constructor
	 * @param A properties containing the Operator parameters
	 * Creates a new intance of the two point crossover operator
	 */
	//public TwoPointsCrossover(Properties properties) {
	//	this();
	//}


	/**
	 * Perform the crossover operation
	 * @param probability Crossover probability
	 * @param parent1 The first parent
	 * @param parent2 The second parent
	 * @return Two offspring solutions
	 * @throws Exception 
	 */
	public Solution[] doCrossover(double   probability, 
			Solution parent1, 
			Solution parent2) throws Exception {

		Solution [] offspring = new Solution[2];

		offspring[0] = new Solution(parent1);
		offspring[1] = new Solution(parent2);
		
		if (parent1.getType().getClass() == solutionType.CauseEffectSolutionType.class) {
				if (PseudoRandom.randDouble() < probability) {
					int crosspoint1        ;
					int crosspoint2        ;
					int sourceListLength1 , sourceListLength2 ;
					int parent1Vector[]    ;
					int parent2Vector[]    ;
					int offspring1Vector[] ;
					int offspring2Vector[] ;

				
//					if (((Sources)(parent1.getDecisionVariables()[0]))._sourcesList.size() < ((Sources)(parent2.getDecisionVariables()[0]))._sourcesList.size())
					sourceListLength1 = ((Sources)(parent1.getDecisionVariables()[0]))._sourcesList.size(); 
				    //else 
					sourceListLength2 = ((Sources)(parent2.getDecisionVariables()[0]))._sourcesList.size();
				    
					 parent1Vector 	= new int[sourceListLength1]; 
				     parent2Vector 	= new int[sourceListLength2];
				     offspring1Vector = new int[sourceListLength1];
				     offspring2Vector = new int[sourceListLength2];
				    
					
					for (int i=0; i<sourceListLength1; i++) {
						parent1Vector[i] =  Integer.parseInt(((Sources)(parent1.getDecisionVariables()[0]))._sourcesList.get(i));
						offspring1Vector[i] =  Integer.parseInt(((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.get(i));
					}
					for (int i=0; i<sourceListLength2; i++) {
						parent2Vector[i] =  Integer.parseInt(((Sources)(parent2.getDecisionVariables()[0]))._sourcesList.get(i));
				    	offspring2Vector[i] =  Integer.parseInt(((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.get(i));
				    	}

			/*		for (int i=0; i<sourceListLength1; i++) 
						System.out.println("OFFSPRING 1 "+offspring1Vector[i]);
					
					for (int i=0; i<sourceListLength2; i++) 
				    	System.out.println("OFFSPRING 2 "+offspring2Vector[i]);
			*/		
					
					// STEP 1: Get two cutting points WITH REFERENCE TO THE SHORTEST LENGTH
					int min = Math.min(sourceListLength1, sourceListLength2)-2;
					crosspoint1 = PseudoRandom.randInt(0,min) ;
					crosspoint2 = PseudoRandom.randInt(0,min) ;

				      if (sourceListLength1<Main.MIN_sources || sourceListLength2<Main.MIN_sources)
				    	  throw new Exception("COP Crossover ERROR: list < min ");
				/*     System.out.println("crossover point 1 "+crosspoint1);
				     System.out.println("crossover point 2 "+crosspoint2);
				     System.out.println("sourceListLength1 "+sourceListLength1);
				     System.out.println("sourceListLength2 "+sourceListLength2);*/
					while (crosspoint2 == crosspoint1)  {
						crosspoint2 = PseudoRandom.randInt(0,min) ;
					}
					if (crosspoint1 > crosspoint2) {
						int swap ;
						swap        = crosspoint1 ;
						crosspoint1 = crosspoint2 ;
						crosspoint2 = swap          ;
					} // if

				/*	System.out.println("cutting point1: "+ crosspoint1);
				      System.out.println("cutting point2: "+ crosspoint2);

				      
				      System.out.println("source list length: "+ sourceListLength1);
				      System.out.println("source list length: "+ sourceListLength2);
				  */    
					// STEP 2: Obtain the first child
					int m = 0;
					for(int j = 0; j < sourceListLength2; j++) {
						boolean exist = false;
						int temp = parent2Vector[j];
						for(int k = crosspoint1; k <= crosspoint2; k++) {
							if (temp == offspring1Vector[k]) {
								exist = true;
								break;
							} // if
						} // for
						if (!exist) {
							if (m == crosspoint1)
								m = crosspoint2 + 1;
						//	System.out.println("m: "+m);
					//		System.out.println("offspring1Vector[m] "+offspring1Vector[m]);
				//			System.out.println("temp "+temp);
							offspring1Vector[m++] = temp;
						} // if
						if (m==sourceListLength1)
							break;
					} // for

					// STEP 3: Obtain the second child
					m = 0;
					for(int j = 0; j < sourceListLength1; j++) {
						boolean exist = false;
						int temp = parent1Vector[j];
						for(int k = crosspoint1; k <= crosspoint2; k++) {
							if (temp == offspring2Vector[k]) {
								exist = true;
								break;
							} // if
						} // for
						if(!exist) {
							if (m == crosspoint1)
								m = crosspoint2 + 1;
							offspring2Vector[m++] = temp;
						} // if
						if (m==sourceListLength2)
							break;
					} // for
					
					// Remove duplicates if any
					LinkedHashSet<String> hs1 = new LinkedHashSet<String>();
					for (int i=0; i< offspring1Vector.length; i++)
						hs1.add(String.valueOf(offspring1Vector[i]));
					
					LinkedHashSet<String> hs2= new LinkedHashSet<String>();
					for (int i=0; i< offspring2Vector.length; i++)
						hs2.add(String.valueOf(offspring2Vector[i]));
					
				//	if there is a duplicate, the hash list is smaller.   
				//	sourceList  must be cleaned and re-created from scratch
					
					((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.clear();
					((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.clear();
					((Targets)(offspring[0].getDecisionVariables()[1]))._targetsList.clear();
					((Targets)(offspring[1].getDecisionVariables()[1]))._targetsList.clear();
					
					
				//	Iterator it1 = hs1.iterator();
				//	Iterator it2 = hs2.iterator();
					((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.addAll(hs1);
					((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.addAll(hs1);

				/*	while (it1.hasNext()) 
				    	((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.add(String.valueOf(it1.next()));
					//for (int i =0; i< sourceListLength2; i++)
				    while (it2.hasNext())
				    	((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.add(String.valueOf(it2.next()));
				*/	
				    while (((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.size()<Main.MIN_sources) {
				    	String sourceToAdd; 
				    	do {
				    		sourceToAdd = Main.possibleSources.get(Main.ran.nextInt(Main.possibleSources.size()));
				    	}
				    	while (((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.contains(sourceToAdd));
				    	((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.add(sourceToAdd);
				    }
				    while (((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.size()<Main.MIN_sources) {
				    	String sourceToAdd; 
				    	do {
				    		sourceToAdd = Main.possibleSources.get(Main.ran.nextInt(Main.possibleSources.size()));
				    	}
				    	while (((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.contains(sourceToAdd));
				    	((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.add(sourceToAdd);
				    }
				    
				    ((Targets)(offspring[0].getDecisionVariables()[1]))._targetsList.add(((Targets)(parent1.getDecisionVariables()[1]))._targetsList.get(Main.ran.nextInt(((Targets)(parent1.getDecisionVariables()[1]))._targetsList.size()))); 
				    ((Targets)(offspring[1].getDecisionVariables()[1]))._targetsList.add(((Targets)(parent2.getDecisionVariables()[1]))._targetsList.get(Main.ran.nextInt(((Targets)(parent1.getDecisionVariables()[1]))._targetsList.size())));
				    
				   // System.out.println("crossover done");
				    
				    //System.out.println("crossover done");
				    //System.out.println("size "+((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.size());
				    //System.out.println("size "+((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.size());
				/*	for (int i =0; i< hs1.size(); i++) 
						System.out.println("First SOL: " +((Sources)(offspring[0].getDecisionVariables()[0]))._sourcesList.get(i));
					for (int i =0; i< hs2.size(); i++)
						System.out.println("SECOND SOL: " +((Sources)(offspring[1].getDecisionVariables()[0]))._sourcesList.get(i));
				*/
					// if
					
					
			    }
		} // if
			else
			{
				Configuration.logger_.severe("TwoPointsCrossover.doCrossover: invalid " +
						"type" + 
						parent1.getDecisionVariables()[0].getVariableType());
				Class cls = java.lang.String.class;
				String name = cls.getName(); 
				throw new JMException("Exception in " + name + ".doCrossover()") ; 
			}

		
	    
		return offspring;                                                                                      
	} // makeCrossover

	/**
	 * Executes the operation
	 * @param object An object containing an array of two solutions 
	 * @return An object containing an array with the offSprings
	 * @throws JMException 
	 */
	public Object execute(Object object) throws JMException {
		Solution [] parents = (Solution [])object;
		Double crossoverProbability ;

    if (!(VALID_TYPES.contains(parents[0].getType().getClass())  &&
        VALID_TYPES.contains(parents[1].getType().getClass())) ) {

			Configuration.logger_.severe("TwoPointsCrossover.execute: the solutions " +
					"are not of the right type. The type should be 'CauseEffect', but " +
					parents[0].getType() + " and " + 
					parents[1].getType() + " are obtained");
		} // if 

		crossoverProbability = (Double)getParameter("probability");

		if (parents.length < 2)
		{
			Configuration.logger_.severe("TwoPointsCrossover.execute: operator needs two " +
			"parents");
			Class cls = java.lang.String.class;
			String name = cls.getName(); 
			throw new JMException("Exception in " + name + ".execute()") ;      
		}

		Solution[] offspring = null;
		try {
			offspring = doCrossover(crossoverProbability_,
					parents[0],
					parents[1]);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return offspring; 
	} // execute

} // TwoPointsCrossover
