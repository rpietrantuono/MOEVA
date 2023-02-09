package causalOptimization.variable;

import java.util.ArrayList;
import java.util.Collections;

import causalOptimization.Main;

//import org.apache.commons.lang3.RandomUtils;
//import java.util.Random;

import jmetal.core.Variable;
import util.SolutionUtils;

public class Sources extends Variable {
	
	public ArrayList<String> _sourcesList; // Both To Be Made "Private"
	public ArrayList<Double> _doubleValueList;  // double representation of the solution. Size = MaxSources, padded with '0' values
	public ArrayList<String> _sourcesFullList;  // Used by PSO 
	
	
	public Sources(ArrayList<String> sourcesList) {
		super();
		this._sourcesList = new ArrayList<String>();
		for (int i=0; i<sourcesList.size();i++) _sourcesList.add(sourcesList.get(i));
	}

	public Sources(Sources sources) throws Exception {
		this._sourcesList = new ArrayList<String>();
		this._doubleValueList = new ArrayList<Double>(Main.MAX_sources); //will put '0' for the Max_source - n_sources empty values
		this._sourcesFullList = new ArrayList<String>(Main.MAX_sources); //will put '0' for the Max_source - n_sources empty values
	
		for (int i=0; i<sources._sourcesList.size();i++) {
			this._sourcesList.add(sources._sourcesList.get(i));
			this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(sources._sourcesList.get(i)));
			this._sourcesFullList.add(sources._sourcesList.get(i)); 
		}
		for(int i=0; i<Main.MAX_sources - this._sourcesList.size();i++) { 
			this._doubleValueList.add(0.0); // Pad with '0'
			this._sourcesFullList.add("0");
		}
		
		}

	
	
	public Sources(int n_sources) throws Exception {
		//n_sources = Math.max(n_sources , Main.MAX_sources);
		this._sourcesList = new ArrayList<String>(n_sources);
		this._doubleValueList = new ArrayList<Double>(Main.MAX_sources); //will put '0' for the Max_source - n_sources empty values
		this._sourcesFullList = new ArrayList<String>(Main.MAX_sources); //will put '0' for the Max_source - n_sources empty values
		
		// Sources taken from knowledge base
		if (Main.solutionTypeFlag.equals("factual")) {
			if(n_sources>Main.possibleKBsourcesPerType.size())
				n_sources = Main.possibleKBsourcesPerType.size();  //limit n_sources in this case 
			
			//ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleKBsources.size());
			ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleKBsourcesPerType.size());
			System.out.println("possible kb source type "+Main.possibleKBsourcesPerType);
		    for (int ind = 0; ind < Main.possibleKBsourcesPerType.size(); ind++) 
		    	n_values_index.add(ind);
		    Collections.shuffle(n_values_index, Main.ran);
		    for (int k = 0; k< n_sources; k++) {
		    	String soruceToAdd = (Main.possibleKBsourcesPerType.get(n_values_index.get(k))).get(Main.ran.nextInt(Main.possibleKBsourcesPerType.get(n_values_index.get(k)).size()));
		    	//		System.out.println("source in constr "+soruceToAdd);
		    	this._sourcesList.add(soruceToAdd);
		    	this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(soruceToAdd));
		    	this._sourcesFullList.add(soruceToAdd);
		    }
			}
		// Sources taken from ontology
	
		// two-level selection is needed: first select the source type and then the specific source (so as the probability of selection is not influenced by how many elements per source are present
		if (Main.solutionTypeFlag.equals("analogical") || Main.solutionTypeFlag.equals("creative")) {
			ArrayList<Integer> n_values_index_external = new ArrayList<Integer>(Main.possibleSourcePerType.size());
			for (int ind = 0; ind < Main.possibleSourcePerType.size(); ind++) 
				n_values_index_external.add(ind);
		    Collections.shuffle(n_values_index_external, Main.ran);
		    for (int k = 0; k< n_sources; k++) {
		    	String sourceToAdd = (Main.possibleSourcePerType.get(n_values_index_external.get(k))).get(Main.ran.nextInt(Main.possibleSourcePerType.get(n_values_index_external.get(k)).size()));
		    	this._sourcesList.add(sourceToAdd);
		    	this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(sourceToAdd));
		    	this._sourcesFullList.add(sourceToAdd);
		    }
		}
			/*
			ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleSources.size());//from ontology
		    for (int ind = 0; ind < Main.possibleSources.size(); ind++) 
		    	n_values_index.add(ind);
		    Collections.shuffle(n_values_index);
		    for (int k = 0; k< n_sources; k++)
		    	this._sourcesList.add(Main.possibleSources.get(n_values_index.get(k)));
			}
		*/
	    
		// two-level selection is needed: first select the source type and then the specific source (so as the probability of selection is not influenced by how many elements per source are present
		if (Main.solutionTypeFlag.equals("graph")) {
			double r;
			int duplicate = 0; 
			int max_duplicate = 10;
			double cum=0;
			int index,element=0;
			for (int k = 0; k< n_sources; k++) {
				r = Main.ran.nextDouble();
				cum=0; 
				index=0;
				for(index=0; index< Main.graph_weights.length; index++){
					cum = cum + Main.graph_weights[index];
					if(r<=cum)
						break;
				}
				
				element=0;
				cum=0;
				r = Main.ran.nextDouble();
				for (element=0; element<Main.possibleSourcePerType.get(index).size(); element++) {
					cum = cum + Main.distributions[index].get(element);
					if(r<=cum)
						break;
				}
				
				String sourceToAdd = (Main.possibleSourcePerType.get(index)).get(element);
			//	System.out.println(" Main.possibleSourcePerType idnex  "+Main.possibleSourcePerType.get(index));
				
			//	System.out.println("source to add: "+sourceToAdd);
				if(_sourcesList.contains(sourceToAdd)) {
					//System.out.println("duplicate");
					duplicate++;
					//System.out.println("solutio " + _sourcesList);
					k--;
				}else {
					this._sourcesList.add(sourceToAdd);
				   	this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(sourceToAdd));
				   	this._sourcesFullList.add(sourceToAdd);
				}
				if (duplicate> max_duplicate)   // TO improve this check
					break;
			    }
		}
		
		if (!Main.solutionTypeFlag.equals("factual") && !Main.solutionTypeFlag.equals("analogical") && !Main.solutionTypeFlag.equals("creative")&& !Main.solutionTypeFlag.equals("graph")) 
			throw new Exception("Creation of random solution: solution type is not correct");
		
		for(int i=0; i<Main.MAX_sources - this._sourcesList.size();i++) { 
			this._doubleValueList.add(0.0); // Pad with '0'
			this._sourcesFullList.add("0");
		}
	}

		
	public ArrayList<String> get_sourcesList() {
		return _sourcesList;
	}

	public void set_sourcesList(ArrayList<String> _sourcesList) {
		this._sourcesList = _sourcesList;
	}

	@Override
	public Variable deepCopy() {
		try {
			return new Sources(this);
		} catch (Exception e) {e.printStackTrace();}
		return null;
	}

	public void updateSource() throws Exception {
		_sourcesList.clear();
		ArrayList<Integer> zeroIndexes = new ArrayList<Integer>();
		//int duplicates=0;
		for (int i= 0; i< this._sourcesFullList.size(); i++) {
			if(!this._sourcesFullList.get(i).equals("0")) {
				if (!this._sourcesList.contains(this._sourcesFullList.get(i))) {
					_sourcesList.add(this._sourcesFullList.get(i));
				}
				else {
					//duplicates++;
					zeroIndexes .add(i);  // Pad with '0'
				}
			}
		else {
			// remember the 0 elements, can be useful later
			zeroIndexes .add(i); 
			}
		}
		
		int minSize = _sourcesList.size() - Main.MIN_sources;  
		if (minSize<0) {
			ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleSources.size());
			for (int ind = 0; ind < Main.possibleSources.size(); ind++) 
		    	n_values_index.add(ind);
		    Collections.shuffle(n_values_index, Main.ran);
			for (int k = 0; k < Math.abs(minSize); k++) {
		    	this._sourcesList.add(Main.possibleSources.get(n_values_index.get(k)));
		    	this._doubleValueList.set(zeroIndexes .get(k), SolutionUtils.getDoubleRepresentation(Main.possibleSources.get(n_values_index.get(k))));
		    	this._sourcesFullList.set(zeroIndexes .get(k), Main.possibleSources.get(n_values_index.get(k)));
		    }
		    	
		}
		
	}
}
