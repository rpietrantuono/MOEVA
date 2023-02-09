package causalOptimization.variable;

import java.util.ArrayList;
import java.util.Collections;

import causalOptimization.Main;

//import org.apache.commons.lang3.RandomUtils;

import jmetal.core.Variable;
import util.SolutionUtils;

public class Targets extends Variable {
	
	public ArrayList<String> _targetsList;
	public ArrayList<Double> _doubleValueList;  // double representation of the solution. Size = MaxSources, padded with '0' values
	public ArrayList<String> _targetsFullList;
	
	public Targets(ArrayList<String> targetsList) {
		super();
		this._targetsList = new ArrayList<String>();
		for (int i=0; i<targetsList.size();i++) _targetsList.add(targetsList.get(i));
	}

	//COSTRUTTORE DI COPIA
	public Targets(Targets targets) throws Exception {
		this._targetsList = new ArrayList<String>();
		this._doubleValueList = new ArrayList<Double>(Main.MAX_targets); //will put '0' for the Max_source - n_sources empty values
		this._targetsFullList = new ArrayList<String>(Main.MAX_targets); //will put '0' for the Max_source - n_sources empty values
	
		for (int i=0; i<targets._targetsList.size();i++) {
			this._targetsList.add(targets._targetsList.get(i));
			this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(targets._targetsList.get(i)));
			this._targetsFullList.add(targets._targetsList.get(i)); 
		}
		for(int i=0; i<Main.MAX_targets - this._targetsList.size();i++) { 
			this._doubleValueList.add(0.0); // Pad with '0'
			this._targetsFullList.add("0");
		}
		
		}

	
	
	public Targets(int n_targets) throws Exception {
		//n_targets= Math.max(n_targets, Main.MAX_targets);
		this._targetsList = new ArrayList<String>(n_targets);
		this._doubleValueList = new ArrayList<Double>(Main.MAX_targets); //will put '0' for the Max_source - n_sources empty values
		this._targetsFullList = new ArrayList<String>(Main.MAX_targets); //will put '0' for the Max_source - n_sources empty values
	
		// Select from targest in the KB
		if (Main.solutionTypeFlag.equals("factual")||Main.solutionTypeFlag.equals("analogical") || Main.solutionTypeFlag.equals("creative")) {
			int existing_sol_index, target_index;  
			for (int k = 0; k< n_targets; k++) {
				existing_sol_index = Main.ran.nextInt(Main.listOfInternalSolutions.size());
//				for (int h=0; h< ((Targets)(Main.listOfInternalSolutions.get(existing_sol_index).getDecisionVariables()[1]))._targetsList.size();h++)
//					System.out.println("Selected solution targets "+((Targets)(Main.listOfInternalSolutions.get(existing_sol_index).getDecisionVariables()[1]))._targetsList.get(h));
				target_index = Main.ran.nextInt(((Targets)(Main.listOfInternalSolutions.get(existing_sol_index).getDecisionVariables()[1]))._targetsList.size());;
				String targetToAdd = ((Targets)(Main.listOfInternalSolutions.get(existing_sol_index).getDecisionVariables()[1]))._targetsList.get(target_index);
				_targetsList.add(targetToAdd);
				this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(targetToAdd));
				this._targetsFullList.add(targetToAdd);
			}
		}
		
		/*
		// Select from targest in the ontology
		if (Main.solutionTypeFlag.equals("analogical") || Main.solutionTypeFlag.equals("creative")) {
			ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleTargets.size());
		    for (int ind = 0; ind < Main.possibleTargets.size(); ind++) 
		    	n_values_index.add(ind);
		    Collections.shuffle(n_values_index, Main.ran);
		    for (int k = 0; k< n_targets; k++) {
		    	this._targetsList.add(Main.possibleTargets.get(n_values_index.get(k)));
		    	this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(Main.possibleTargets.get(n_values_index.get(k))));
		    	this._targetsFullList.add(Main.possibleTargets.get(n_values_index.get(k)));
		    }
		}
		*/
		if (Main.solutionTypeFlag.equals("graph")) {
			// target is unique in our tests. Not needed a random selection of the column, just the row. 
			double r = Main.ran.nextDouble();
			double cum=0;
			int index=0;
			for (index=0; index<Main.possibleTargets.size(); index++) {
				cum = cum + Main.distributions[Main.distributions.length-1].get(index);
				if(r<=cum)
					break;
			}
			String targetToAdd = (Main.possibleTargets.get(index)); //assume one value
			
	    	
		    this._targetsList.add(targetToAdd);
		    this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(targetToAdd));
		    this._targetsFullList.add(targetToAdd);
		}
		
		if (!Main.solutionTypeFlag.equals("factual") && !Main.solutionTypeFlag.equals("analogical") && !Main.solutionTypeFlag.equals("creative")&& !Main.solutionTypeFlag.equals("graph")) 
			throw new Exception("Creation of random solution: solution type is not correct");
		
		
		for(int i=0; i<Main.MAX_targets - this._targetsList.size();i++) { 
			this._doubleValueList.add(0.0); // Pad with '0'
			this._targetsFullList.add("0");
		}
			
	}

	

	public ArrayList<String> get_targetsList() {
		return _targetsList;
	}

	public void set_targetsList(ArrayList<String> _targetsList) {
		this._targetsList = _targetsList;
	}

	@Override
	public Variable deepCopy() {
		try {
			return new Targets(this);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	public void updateTarget() throws Exception {
		_targetsList.clear();
		ArrayList<Integer> zeroIndexes = new ArrayList<Integer>();
		for (int i= 0; i< this._targetsFullList.size(); i++) {
			if(!this._targetsFullList.get(i).equals("0")) {
				if (!this._targetsList.contains(this._targetsFullList.get(i)))
				_targetsList.add(this._targetsFullList.get(i));
			}
			else {
				// remember the 0 elements, can be useful later
				zeroIndexes .add(i); 
			}
		}
		int minSize = _targetsList.size() - Main.MIN_targets;  
		if (minSize<0) {
			ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleTargets.size());
			for (int ind = 0; ind < Main.possibleTargets.size(); ind++) 
		    	n_values_index.add(ind);
		    Collections.shuffle(n_values_index, Main.ran);
			for (int k = 0; k < Math.abs(minSize); k++) {
		    	this._targetsList.add(Main.possibleTargets.get(n_values_index.get(k)));
		    	this._doubleValueList.set(zeroIndexes .get(k), SolutionUtils.getDoubleRepresentation(Main.possibleTargets.get(n_values_index.get(k))));
		    	this._targetsFullList.set(zeroIndexes .get(k), Main.possibleTargets.get(n_values_index.get(k)));
		    }
		    	
		}
	}
}
