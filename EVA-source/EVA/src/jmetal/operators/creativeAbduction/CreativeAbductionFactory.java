package jmetal.operators.creativeAbduction;

import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.util.HashMap;

/**
 * Class implementing a factory for abductive operators.
 */
public class CreativeAbductionFactory {
    
  /**
   * Gets a crossover operator through its name.
   * @param name Name of the operator
   * @return The operator
   */
	public static CreativeAbduction getCreativeAbductionOperator(String name, HashMap parameters) throws JMException {
	    if (name.equalsIgnoreCase("BasicCreativeAbduction"))
	      return new BasicCreativeAbduction(parameters);
	    else {
	      Configuration.logger_.severe("CreativeAbductionFactory.getCreativeAbductionOperator. " +
	          "Operator '" + name + "' not found ");
	      throw new JMException("Exception in " + name + ".getCreativeAbductionOperator()") ;
	    } // else        
	  } // 
	} // 

