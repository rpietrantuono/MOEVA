package util;

import java.util.ArrayList;

public class VectorDistanceJaccardImpl extends VectorDistance {

	// compare two vectors of solutions. Can be used to comptue the distance of a solution from the KB  
	static String NAME = "VECTOR_JACCARD";

	@Override
	public Double calculateDistance(ArrayList<Integer> vector1,
			ArrayList<Integer> vector2) {
		
		double count11 = 0;	  // count of M11
		double count10or01or11 = 0; // count of M01, M10 and M11
		
		// for each position in the vector

		for(int i=0; i< vector1.size();i++){
			// if it is not  two 0s
			if(vector1.get(i) != 0  || vector2.get(i) != 0) {
				// if it is two 1s
				if(vector1.get(i) == 1  && vector2.get(i) == 1) {
					count11++;
				}
				// increase the count of not two 0s
				count10or01or11++;
			}
		}
		
		return (1 - (count11 / count10or01or11));
	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return NAME;
	}

}
