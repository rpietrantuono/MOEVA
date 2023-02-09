package util;

import ca.pfv.spmf.algorithms.clustering.distanceFunctions.DistanceCorrelation;
import ca.pfv.spmf.algorithms.clustering.distanceFunctions.DistanceCosine;
import ca.pfv.spmf.algorithms.clustering.distanceFunctions.DistanceEuclidian;
import ca.pfv.spmf.algorithms.clustering.distanceFunctions.DistanceFunction;
import ca.pfv.spmf.algorithms.clustering.distanceFunctions.DistanceJaccard;
import ca.pfv.spmf.algorithms.clustering.distanceFunctions.DistanceManathan;
import ca.pfv.spmf.patterns.cluster.DoubleArray;

public abstract class Distance {
		
	public  abstract String getName();
		
		
		/**
		 * This method returns the distance function having a given name
		 * @param name the name  (euclidian, manathan, cosine, correlation,...)
		 * @return the distance function
		 */
		public static Distance getDistanceFunctionByName(String name){
			if(ElementDistanceBinaryImpl.NAME.equals(name)) {
				return new ElementDistanceBinaryImpl();
			}else if(MatrixDistanceBinaryImpl.NAME.equals(name)) {
				return new MatrixDistanceBinaryImpl();
			}else if(VectorDistanceJaccardImpl.NAME.equals(name)) {
				return new VectorDistanceJaccardImpl();
			}
			return null;
		}

}
