package jguis;

import java.util.ArrayList;
import java.util.List;

import skeleton_analysis.AnalyzeSkeleton_;
import skeleton_analysis.Edge;
import skeleton_analysis.Graph;
import skeleton_analysis.Point;
import skeleton_analysis.SkeletonResult;
import skeleton_analysis.Vertex;
import ij.IJ;
import ij.ImagePlus;

public class SkeletonUtilities{
	
	public AnalyzeSkeleton_ as;
	public SkeletonResult skelResult;
	public byte[][] data;
	public int width,height,slices;
	public double[][] dists;
	public int[][] next;
	
	public SkeletonUtilities(ImagePlus imp){
		as=new AnalyzeSkeleton_();
		AnalyzeSkeleton_.calculateShortestPath=true;
		as.setup("",imp);
		skelResult=as.run(AnalyzeSkeleton_.NONE,false,true,null,true,false);
	}
	
	public List<List<int[]>> getShortestPaths(){
		ArrayList<Point>[] sppoints=as.getShortestPathPoints();
		List<List<int[]>> pts=new ArrayList<List<int[]>>();
		for(int i=0;i<sppoints.length;i++){
			List<int[]> path=new ArrayList<int[]>();
			for(int j=0;j<sppoints[i].size();j++){
				Point temp=sppoints[i].get(j);
				path.add(new int[]{temp.x,temp.y,temp.z});
			}
			pts.add(path);
		}
		return pts;
	}
	
	public double[][] getSpStartPositions(){
		return skelResult.getSpStartPosition();
	}
	
	//here we get the shortest path ourselves
	public List<List<int[]>> getShortestPaths2(){
		Graph[] graphs=as.getGraphs();
		List<List<int[]>> pts=new ArrayList<List<int[]>>();
		for(int i=0;i<graphs.length;i++){
			List<int[]> temp=getShortestPath(graphs[i]);
			if(temp!=null && temp.size()>0) pts.add(temp);
		}
		return pts;
	}
	
	/*******************
	 * this is my attempt at the warshall algorithm, I pull extensively from the Fiji version written by Huub Hovens
	 * @param graph
	 * @return
	 */
	public List<int[]> getShortestPath(Graph graph){
		List<Vertex> vertices=graph.getVertices();
		if(vertices.size()==1) return null;
		List<Edge> edges=graph.getEdges();
		int nvertices=vertices.size();
		dists=new double[nvertices][nvertices];
		next=new int[nvertices][nvertices];
		for(int i=0;i<dists.length;i++){
			for(int j=0;j<dists[i].length;j++){
				dists[i][j]=Double.POSITIVE_INFINITY;
				//if(i==j) dists[i][j]=0.0f;
				next[i][j]=-1;
			}
		}
		for(Edge edge:edges){
			Vertex v1=edge.getV1();
			Vertex v2=edge.getV2();
			int row=vertices.indexOf(v1);
			int col=vertices.indexOf(v2);
			dists[row][row]=0.0f;
			dists[col][col]=0.0f;
			dists[row][col]=edge.getLength();
			dists[col][row]=dists[row][col];
			next[row][row]=-1;
			next[col][col]=-1;
			next[row][col]=row;
			next[col][row]=col;
		}
		for(int k=0;k<nvertices;k++){
			for(int i=0;i<nvertices;i++){
				for(int j=0;j<nvertices;j++){
					if((dists[i][k]+dists[k][j])<dists[i][j]){
						dists[i][j]=dists[i][k]+dists[k][j];
						next[i][j]=next[k][j]; //fiji version
						//next[i][j]=next[i][k]; //wiki version
					}
				}
			}
		}
		//maxPath is the longest shortest (direct) path
		//c is the start vertex and d is the end vertex of that path
		double maxPath=0.0;
		int c=0;
		int d=0;
		for(int i=0;i<nvertices;i++){
			for(int j=0;j<nvertices;j++){
				if(dists[i][j]>maxPath && dists[i][j]!=Double.POSITIVE_INFINITY){
					maxPath=dists[i][j];
					c=i;
					d=j;
				}
			}
		}
		//now walk through and find the best vertex path
		int a=c;
		int b=d;
		List<Vertex> path1=new ArrayList<Vertex>();
		while(b!=a){
			path1.add(vertices.get(b));
			b=next[a][b]; //fiji version
			//a=next[a][b]; //wiki version
		}
		path1.add(vertices.get(a));
		List<int[]> path=new ArrayList<int[]>();
		//now trace the shortest edges between these vertices
		for(int i=1;i<path1.size();i++){
			List<Edge> segedges=new ArrayList<Edge>();
			List<Integer> orientations=new ArrayList<Integer>();
			//find the edges for this segment
			Vertex v1=path1.get(i-1);
			Vertex v2=path1.get(i);
			for(Edge edge:edges){
				if(edge.getV1()==v1 && edge.getV2()==v2){
					segedges.add(edge);
					orientations.add(1);
				} else if(edge.getV1()==v2 && edge.getV2()==v1){
					segedges.add(edge);
					orientations.add(-1);
				}
			}
			if(segedges.size()==0) continue;
			//find the shortest edge for this segment
			int selected=0;
			double shortestlength=segedges.get(0).getLength();
			for(int j=1;j<segedges.size();j++){
				if(segedges.get(j).getLength()<shortestlength){
					shortestlength=segedges.get(j).getLength();
					selected=j;
				}
			}
			//now trace this edge
			List<Point> slabs=segedges.get(selected).getSlabs();
			for (int j=0;j<slabs.size();j++){
				Point p=slabs.get(j);
				if(orientations.get(selected).intValue()<0) p=slabs.get(slabs.size()-j-1);
				path.add(new int[]{p.x,p.y,p.z});
			}
			//for now we will ignore the vertices
		}
		return path;
	}
	
	public String getVertexInfo(int graphindex){
		Graph[] graphs=as.getGraphs();
		Graph graph=graphs[graphindex];
		List<Vertex> vertices=graph.getVertices();
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<vertices.size();i++){
			sb.append(vertices.get(i).pointsToString()+"\n");
		}
		return sb.toString();
	}

}
