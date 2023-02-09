package causalOptimization.CopVariable;

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
	
		ArrayList<Integer> n_values_index = new ArrayList<Integer>(Main.possibleTargets.size());
	    for (int ind = 0; ind < Main.possibleTargets.size(); ind++) 
	    	n_values_index.add(ind);
	    Collections.shuffle(n_values_index, Main.ran);
	    for (int k = 0; k< n_targets; k++) {
	    	this._targetsList.add(Main.possibleTargets.get(n_values_index.get(k)));
	    	this._doubleValueList.add(SolutionUtils.getDoubleRepresentation(Main.possibleTargets.get(n_values_index.get(k))));
	    	this._targetsFullList.add(Main.possibleTargets.get(n_values_index.get(k)));
	    }
		
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
