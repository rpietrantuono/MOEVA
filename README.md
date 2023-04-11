# The Evolutionary Abduction (EVA) Repository
## Source

This repository contains the material of the following work: 


Roberto Pietrantuono. 2023. An Evolutionary Strategy for Automatic Hypotheses Generation inspired by Abductive Reasoning. In *Genetic and Evolutionary Computation Conference Companion* (GECCO ’23 Companion), July 15–19, 2023, Lisbon, Portugal. ACM, New York, NY, USA, 5 pages. https://doi.org/10.1145/3583133.3590568

The work proposes an evolutionary algorithm to solve a specific class of optimization problems, called Causal Combinatorial Optimization Problems (CCOP), making use of operators that emulate the abductive human reasoning. The algorithm automatically advances hypotheses for potential causes for a given effect, hence allowing to find the most plausible (accoridng to a defined plausibility metric) explanations for a certain event. The application is on 4 datasets (in medial, decision-support system, and safety engineering domains). Results show that the algorithm can formulate hypotheses that are equal or very simialr to really-occurred events. 

## Description
The repository contains the artefacts required to run EVA. The code exploits  jMetal (https://github.com/jMetal), a Java framework to develop and experiment evolutionary algorithms. 

The repository also contains the **Appendix** (Appendix.pdf) with details about setting not reported in the paper and a description of the ASRS dataset. 

The artefacts include: 
- EVA.jar. The executable JAR file to run the algorithm on a benchmark problem given as input. 
Usage: 
If the EVA algorithm has to be run, then type: 
java -jar EVA.jar <configuration.txt> r=<number of runs (optional, default is 1)> technique=EVA

where: 
* "configuration.txt" is a file in which the user can set all the parameters required for the run. Examples are in the BENCHMARKS folder. To use it in a benchmark different from the experimented ones, namely for replicating our study or for comparing EVA with other strategies, just copy and customise one of the available files from the BENCHMARKS folder. In that file, the text after the "#" character are comments describing the meaning of every parameter; 
* r is the number of runs (optional, default is 1); 

- BENCHMARKS. It contains the four benchmarks problems. For each of them, there is: an instance of the configuration.txt file mentioned above; the dataset (e.g., medical.csv for the first problem); the OntologySource.txt and OntologyTargets.txt files that enumerates the variables and their possible values (in this format: 1:1, 1:2, ...2:1, 2:2, ..., meaning: variable1: value1, variable1:value2). It is also reported a "description" folder that has the original dataset and a description of the attributes, which we used as "variables:values".  

Finally, if a MOEA strategy among those used in the paper as baseline need to be run, then refer to the CCOP repository (https://github.com/rpietrantuono/CCOP).
