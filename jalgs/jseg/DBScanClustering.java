package jalgs.jseg;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class DBScanClustering{
	public int minPts;
	public double eps;
	
	public DBScanClustering(double eps,int minPts){
		this.minPts=minPts;
		this.eps=eps;
	}
	
	public List<List<float[]>> cluster(List<float[]> points){
		List<List<float[]>> clusters=new ArrayList<List<float[]>>();
		Map<float[],Integer> visited=new HashMap<float[],Integer>();
		//the integer in the map is the point status: 0 is noise, 1 is part of cluster
		for(int i=0;i<points.size();i++){
			float[] point=points.get(i);
			if(visited.get(point)!=null){
				continue;
			}
			List<float[]> neighbors=getNeighbors(point,points);
			if(neighbors.size()>=minPts){
				List<float[]> cluster=new ArrayList<float[]>();
				clusters.add(expandCluster(cluster,point,neighbors,points,visited));
			} else {
				visited.put(point,new Integer(0));
			}
		}
		return clusters;
	}
	
	public List<float[]> expandCluster(List<float[]> cluster,float[] point,List<float[]> neighbors,List<float[]> points,Map<float[],Integer> visited){
		cluster.add(point);
		visited.put(point,new Integer(1));
		List<float[]> seeds=new ArrayList<float[]>(neighbors);
		int index=0;
		while(index<seeds.size()){
			float[] current=seeds.get(index);
			Integer pStatus=visited.get(current);
			if(pStatus==null){
				List<float[]> currentNeighbors=getNeighbors(current,points);
				if(currentNeighbors.size()>=minPts){
					seeds=merge(seeds,currentNeighbors);
				}
			}
			if(pStatus==null || pStatus.intValue()!=1){
				visited.put(current,new Integer(1));
				cluster.add(current);
			}
			index++;
		}
		return cluster;
	}
	
	public List<float[]> getNeighbors(float[] point,List<float[]> points){
		List<float[]> neighbors=new ArrayList<float[]>();
		for(int i=0;i<points.size();i++){
			float[] neighbor=points.get(i);
			if(point!=neighbor && distance(neighbor,point)<eps){
				neighbors.add(neighbor);
			}
		}
		return neighbors;
	}
	
	public double distance(float[] point1,float[] point2){
		double d2=0.0;
		for(int i=0;i<point1.length;i++) d2+=(point2[i]-point1[i])*(point2[i]-point1[i]);
		return Math.sqrt(d2);
	}
	
	public List<float[]> merge(List<float[]> one,List<float[]> two){
		Set<float[]> oneSet=new HashSet<float[]>(one);
		for(float[] item:two){
			if(!oneSet.contains(item)){
				one.add(item);
			}
		}
		return one;
	}
	
}
