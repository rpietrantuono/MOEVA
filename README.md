# The Evolutionary Abduction (EVA) Repository
## Source

This repository contains the material of the following work: 

"Automated Hypotheses Generation via Evolutionary Abduction", submitted for review to the  The Tenth International Conference on Learning Representations, 2022.

## Description
The repository contains the artefacts required to run EVA, the algorithm for evolutionary abduction implemented to solve Combinatorial Causal Optimization Problem (CCOP), as well as to reproduce the results reported in the paper (on four benchmark datasets). The code exploits the jMetal framework (https://github.com/jMetal), a Java framework to develop and experiment evolutionary algorithms. 

The artefacts include: 
- EVA.jar. The executable JAR file to run the algorithm on a benchmark problem given as input. 
Usage: 
If the EVA algorithm or the RANDOM baseline strategy have to be run, then type: 
java -jar EVA.jar <configuration.txt> r=<number of runs (optional, default is 1)> technique=EVA (or =RAN)

where: 
* "configuration.txt" is a file in which the user can set all the parameters required for the run. Examples are in the BENCHMARKS folder. To use it in a benchmark different from the experimented ones, namely for replicating our study or for comparing EVA with other strategies, just copy and customise one of the available files from the BENCHMARKS folder. In that file, the text after the "#" character are comments describing the meaning of every parameter; 
* r is the number of runs (optional, default is 1); 

If a GB baseline strategy has to be run, then type: 
java -jar ../../EVA.jar ./BEST_01/configuration.txt r=10 technique=GB w=weights.txt

where: 
* "configuration.txt" is the same as above; 
* r is the same as above; 
* technique = GB is to specify that a GB technique will be run
* w gives the weights.txt file as input, specifying the weights obtained by the casual structure discovery algorithm (in the paper, they were FGES, GFCI, RFCI). The file should contain a list of comma-separated double values, one per each cause variable, representing the probability for the cause variable to be causally related to the effect variable (these are obtained, in our paper, from the mentioned causal structure discovery algorithms). 

Note: in order to reproduce the same result of the paper, there are <run.sh> files in the sub-folders of /RESULTS, as described in the following. Those scripts run all the experiments related to a given dataset, with either EVA or with the baselines strategies.  Alternatively, to run one single experiment, you can use EVA.jar directly (with the commands above): in that case, to have the same results of the paper, leave the "seed" parameter set to 1 in the configuration file, as well as the "KB size" and "external KB size" unchanged; then, run EVA with r=10 repetitions (the seed is incremented by 1 at every repetition). Pre-set configuration files can also be found in the folders related to the experiment to repeat (again under /RESULTS)

- EVA-SOURCE. It contains the source code. The code includes not only the EVA algorithm but also the baselines. 

- BENCHMARKS. It contains the four benchmarks problems. For each of them, there is: an instance of the configuration.txt file mentioned above; the dataset (e.g., medical.csv for the first problem); the OntologySource.txt and OntologyTargets.txt files that enumerates the variables and their possible values (in this format: 1:1, 1:2, ...2:1, 2:2, ..., meaning: variable1: value1, variable1:value2). It is also reported a "description" folder that has the original dataset and a description of the attributes, which we used as "variables:values".  

- RESULTS. It contains the results obtained for the four benchmarks described in the paper. There is one folder for EVA and one for all the baselines called GB. Inside, there is one folder per benchmark problem. Within these, there is a folder FINAL_ICLR, with results of the main paper, obtained in the six configurations (BEST, WORST) x (Novelty constraint = 0.1, 0.4, 0.7). In the EVA folders, you will also find TUNING_ICLR phase (with results of tuning, in which the BEST and WORST configurations have been found (see Appendix of the article)). All the results can be re-obtained (overwriting the existing ones) by simply running the <run.sh> (or tuning.sh for TUNING) scripts within these folders. 
All the results (for EVA as well as for the four GB and RANDOM strategies adopted) are in the form of text files reporting the distances obtained for each algorithm with respect to the test dataset (*Ref*.txt files) as well as with respect to the knowledge base (*KB*.txt files), for every repetition. Other statistics are available too. 
