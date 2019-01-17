package jalgs.jseg;

import java.util.List;

public class AgglomCluster{
	//a simple agglomerative clustering algorithm
	//concept, find closest distance between two points or point and cluster
	//then merge those points or clusters
	//need to keep a hierarchy as well as cluster distances
	//cluster distance can be mean, min, or max
	//max number of cluster levels will be n-1 (points have ever increasing distance)
	//min will be 1 (all points are equidistant)
	//identify each point with an id (0 to n-1) and each cluster with an id (in order of discovery)
	//need to return opposite trajectory: start at single cluster and list splits
	public static List<List<float[]>> cluster(List<float[]> points){
		return null;
	}

}
