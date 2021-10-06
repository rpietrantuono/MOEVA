package util;

import java.util.ArrayList;

public abstract class MatrixDistance extends Distance {

	
	public abstract double calculateDistance(ArrayList<Double> vector1, ArrayList<Object> vectorOfVectors) ;
	
	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return null;
	}

}
