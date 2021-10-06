package util;

import java.util.ArrayList;

public class ElementDistanceBinaryImpl extends ElementDistance{

	/*Distance of a property of an entity with respect to a vector of properties 
	 * (it can be used both to compute the distance from the knowledge base and from the ontology)
	 * */ 
	 
	static String NAME = "ELEMENT_BINARY";
	
	@Override
	public int calculateDistance(int element, ArrayList<Integer> elementVector) {
		
		//binary comparison: exist or not exist
		for (int i = 0; i < elementVector.size();)
			if (elementVector.get(i)==element)
				return 1;
			else return 0;
		return 0;
	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return NAME;
	}

	
}
