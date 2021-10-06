package causalOptimization;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.regex.Pattern;

import causalOptimization.CopVariable.Sources;
import causalOptimization.CopVariable.Targets;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.core.Variable;
import jmetal.problems.COPProblem_multiobj;

public class RecomputeDistanceFromFile {

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		Main.myProblemString =new COPProblem_multiobj("CauseEffectSolutionType"); 
		SimilarityEvalutaor sim = new SimilarityEvalutaor();
		FrequencyEvalutaor freq = new FrequencyEvalutaor();
		Main.baseDir = "/Users/robertopietrantuono/Dropbox/Papers/in_preparazione/TEVC-LEARNING/GECCO/MEDICAL";
		String tempDir = "/MyStandardStudyCOP-100-10000-6min"; 
		Main.inputKB_CSVFile = Paths.get(Main.baseDir + "/medical.csv"); 
		Main.listOfInternalSolutions = Main.getInternalSolutionsFromCSV(0, 364); 
		
		SolutionSet solutions= readSolutionsFromFile(Main.baseDir +tempDir + "/solutions.txt"); 
        
		int startingPoint = 365;
		Main.referenceSet = readFromDataSetFile(Main.baseDir + tempDir+"/dataset.txt", startingPoint);
		
		Main.internalKnowledgeBaseFile = Main.baseDir + tempDir + "/internalKBFile.txt";
		Main.itemsetTree_dataset = freq.buildItemSetTree(Main.baseDir + tempDir + "/dataset.txt");
		Main.itemsetTree = freq.buildItemSetTree(Main.internalKnowledgeBaseFile);
		
	    System.out.println("Computing test set distacne...");
	    System.out.println("COPConfigurator.referenceSet  "+Main.referenceSet.size());
	    System.out.println("solutions..."+solutions.size());
	    System.out.println("COPConfigurator.listOfInternalSolutions..."+Main.listOfInternalSolutions.size());
	    
	    sim.computeReferenceSetDistance(solutions, Main.referenceSet ,"closestItem_RefSet_modified",0, true);
        sim.computeReferenceSetDistance(solutions, Main.listOfInternalSolutions,"closestItem_KBSet_modified",0, true);
        
	}

	private static SolutionSet readFromDataSetFile(String filename, int startingPoint) throws IOException {
		ArrayList<String> fileLines = (ArrayList<String>)Files.readAllLines(Paths.get(filename), StandardCharsets.UTF_8);
		COPProblem_multiobj myProblem = new COPProblem_multiobj("CauseEffectSolutionType");
		SolutionSet solSet = new SolutionSet(fileLines.size());
		
		ArrayList<String> sources = new ArrayList<String>();
		ArrayList<String> targets  = new ArrayList<String>();
		StringTokenizer st ;
		for (int i=startingPoint-1; i<fileLines.size();i++) {
			st = new StringTokenizer(fileLines.get(i), " ");   
			while (st.hasMoreTokens()) {
				sources.add(st.nextToken().trim()); 
			}
			//target
			targets.add(sources.get(sources.size()-1));
			sources.remove(sources.size()-1); //remvoce the last, it is a target
			Variable[] variables_string = new Variable[2];
			//ArrayList<String> sol_string = new ArrayList<String>();
			variables_string[0] = new Sources(sources); 
			variables_string[1] =  new Targets(targets);
			Solution ind = new Solution(myProblem, variables_string);
			ind.isInternal=true;
			solSet.add(ind);
			/*** CLEAR ****/
			sources.clear();
			targets.clear();
			}
		
		return solSet;
	}

	private static SolutionSet readSolutionsFromFile(String filename) throws IOException {
		ArrayList<String> fileLines = (ArrayList<String>)Files.readAllLines(Paths.get(filename), StandardCharsets.UTF_8);
		COPProblem_multiobj myProblem = new COPProblem_multiobj("CauseEffectSolutionType");
		SolutionSet solSet = new SolutionSet(fileLines.size());
		ArrayList<String> sources = new ArrayList<String>();
		ArrayList<String> targets  = new ArrayList<String>();
		StringTokenizer st ;
		for (int i=0; i<fileLines.size();i++) {
			String[] cells = fileLines.get(i).split("],");
			String[] cells_1 = cells[0].split(Pattern.quote("["));
			st = new StringTokenizer(cells_1[1], ",");   
			while (st.hasMoreTokens()) {
				sources.add(st.nextToken().trim()); 
			}
			
			//target
			targets.add(sources.get(sources.size()-1));
			sources.remove(sources.size()-1); //remvoce the last, it is a target
			Variable[] variables_string = new Variable[2];
			//ArrayList<String> sol_string = new ArrayList<String>();
			variables_string[0] = new Sources(sources); 
			variables_string[1] =  new Targets(targets);
			Solution ind = new Solution(myProblem, variables_string);
			ind.isInternal=true;
			solSet.add(ind);
			/*** CLEAR ****/
			sources.clear();
			targets.clear();
			
			
		}
		System.out.println("fsize "+fileLines.size());
		
		return solSet;
		
	}

	
}
