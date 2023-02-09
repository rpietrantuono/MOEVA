package util;

import java.util.ArrayList;

import ca.pfv.spmf.patterns.cluster.DoubleArray;

public abstract class ElementDistance extends Distance {

	
	
	public abstract int calculateDistance(int element1, ArrayList<Integer> elementVector); 

	@Override
	public abstract String getName() ;

}
