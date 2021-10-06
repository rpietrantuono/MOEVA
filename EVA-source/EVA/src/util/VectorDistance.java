package util;

import java.util.ArrayList;

import ca.pfv.spmf.patterns.cluster.DoubleArray;

public abstract class VectorDistance extends Distance {

	public abstract Double calculateDistance(ArrayList<Integer> vector1, ArrayList<Integer> vector2) ;

	@Override
	public abstract String getName() ;

}
