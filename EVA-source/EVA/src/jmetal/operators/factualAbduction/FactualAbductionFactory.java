
package jmetal.operators.factualAbduction;

import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.util.HashMap;

/**
 * Class implementing a factory for abductive operators.
 */
public class FactualAbductionFactory {
    
  /**
   * Gets a crossover operator through its name.
   * @param name Name of the operator
   * @return The operator
   */
	public static FactualAbduction getFactualAbductionOperator(String name, HashMap parameters) throws JMException {
	    if (name.equalsIgnoreCase("BasicFactualAbduction"))
	      return new BasicFactualAbduction(parameters);
	    else {
	      Configuration.logger_.severe("FactualAbductionFactory.getFactualAbductionOperator. " +
	          "Operator '" + name + "' not found ");
	      throw new JMException("Exception in " + name + ".getFactualAbductionOperator()") ;
	    } // else        
	  } // 
	} // 

