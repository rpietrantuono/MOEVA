//  StandardStudy.java
//
//  Authors:
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

package jmetal.experiments.studies;

import jmetal.core.Algorithm;
import jmetal.experiments.Experiment;
import jmetal.experiments.Settings;
import jmetal.experiments.settings.*;
import jmetal.experiments.util.Friedman;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class implementing a typical experimental study. Five algorithms are 
 * compared when solving the ZDT, DTLZ, and WFG benchmarks, and the hypervolume,
 * spread and additive epsilon indicators are used for performance assessment.
 */
public class MyStandardStudy extends Experiment {

  /**
   * Configures the algorithms in each independent run
   * @param problemName The problem to solve
   * @param problemIndex
   * @throws ClassNotFoundException 
   */
  public void algorithmSettings(String problemName, 
  		                          int problemIndex, 
  		                          Algorithm[] algorithm) throws ClassNotFoundException {
    try {
      int numberOfAlgorithms = algorithmNameList_.length;

      HashMap[] parameters = new HashMap[numberOfAlgorithms];

      for (int i = 0; i < numberOfAlgorithms; i++) {
        parameters[i] = new HashMap();
      } // for

      if (!paretoFrontFile_[problemIndex].equals("")) {
        for (int i = 0; i < numberOfAlgorithms; i++)
          parameters[i].put("paretoFrontFile_", paretoFrontFile_[problemIndex]);
        } // if

        algorithm[0] = new NSGAII_Settings(problemName).configure(parameters[0]);
        algorithm[1] = new SPEA2_Settings(problemName).configure(parameters[1]);
        algorithm[2] = new MOCell_Settings(problemName).configure(parameters[2]);
        algorithm[3] = new RandomSearch_Settings(problemName).configure(parameters[3]);
     //   algorithm[3] = new SMPSO_Settings(problemName).configure(parameters[3]);
      //  algorithm[4] = new GDE3_Settings(problemName).configure(parameters[4]);
      } catch (IllegalArgumentException ex) {
      Logger.getLogger(MyStandardStudy.class.getName()).log(Level.SEVERE, null, ex);
    } catch (IllegalAccessException ex) {
      Logger.getLogger(MyStandardStudy.class.getName()).log(Level.SEVERE, null, ex);
    } catch  (JMException ex) {
      Logger.getLogger(MyStandardStudy.class.getName()).log(Level.SEVERE, null, ex);
    }
  } // algorithmSettings

  /**
   * Main method
   * @param args
   * @throws JMException
   * @throws IOException
   */
  public static void main(String[] args) throws JMException, IOException {
    MyStandardStudy exp = new MyStandardStudy();

    exp.experimentName_ = "MyStandardStudy";
    exp.algorithmNameList_ = new String[]{
                                "NSGAII", "IBEA", "MOCell", "RandomSearch"};
  
    exp.problemList_ = new String[]{"ReliabilityProblem"};
    exp.paretoFrontFile_ = new String[]{"ReliabilityProblem.pf"};

    
    
    exp.indicatorList_ = new String[]{"HV", "IGD"};//, "EPSILON"};

    int numberOfAlgorithms = exp.algorithmNameList_.length;

    exp.experimentBaseDirectory_ = "/Users/robertopietrantuono/Documents/RelAlloc/" +
                                   exp.experimentName_;
    exp.paretoFrontDirectory_ = "/Users/robertopietrantuono/Documents/RelAlloc/data/paretoFronts";

    exp.algorithmSettings_ = new Settings[numberOfAlgorithms];

    exp.independentRuns_ = 30;

    exp.initExperiment();

    // Run the experiments
    int numberOfThreads ;
    exp.runExperiment(numberOfThreads = 4) ;

    exp.generateQualityIndicators() ;

    // Generate latex tables
    exp.generateLatexTables() ;

    // Configure the R scripts to be generated
    int rows  ;
    int columns  ;
    String prefix ;
    String [] problems ;
    boolean notch ;

    
    rows = 1 ;
    columns = 1 ;
    prefix = new String("ReliabilityProblem");
    problems = new String[]{"ReliabilityProblem"};
    
    exp.generateRBoxplotScripts(rows, columns, problems, prefix, notch = true, exp) ;
    exp.generateRWilcoxonScripts(problems, prefix, exp) ;

    
    // Applying Friedman test
    Friedman test = new Friedman(exp);
    //test.executeTest("EPSILON");
    test.executeTest("HV");
    test.executeTest("IGD");
  } // main
} // StandardStudy


