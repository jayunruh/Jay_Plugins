package jguis;

import java.awt.Polygon;
import java.util.ArrayList;
import java.util.List;

public class MultiPolygon{
	
	public List<Polygon> contents;
	public Polygon poly;
	
	public MultiPolygon(Polygon[] polys){
		contents=new ArrayList<Polygon>();
		contents.add(polys[0]);
		poly=polys[0];
		for(int i=1;i<polys.length;i++){
			contents.add(polys[i]);
			addPolygon(polys[i]);
		}
	}
	
	public boolean contains(int x,int y){
		for(int i=0;i<contents.size();i++){
			if(contents.get(i).contains(x,y)) return true;
		}
		return false;
	}
	
	public void addPolygon(Polygon p){
		int[] closest=findClosestPair(poly,p);
		poly=combine(poly,p,closest);
		contents.add(p);
	}
	
	public static Polygon combine(Polygon p1,Polygon p2,int[] closest){
		//combines separated polygons with an infinitesimal isthmus between them
		//closest lists the indices of closest approach for the two polygons
		int newnpts=p1.npoints+p2.npoints+2;
		int[] newxpts=new int[newnpts];
		int[] newypts=new int[newnpts];
		//move "clockwise" until closest approach
		for(int i=0;i<=closest[0];i++){
			newxpts[i]=p1.xpoints[i];
			newypts[i]=p1.ypoints[i];
		}
		int pos=closest[0]+1;
		//now move "clockwise" around p2 until its endpoint
		for(int i=closest[1];i<p2.npoints;i++){
			newxpts[pos]=p2.xpoints[i];
			newypts[pos]=p2.ypoints[i];
			pos++;
		}
		//now from the beginning of p2 to its closest point (closest point is duplicated)
		for(int i=0;i<=closest[1];i++){
			newxpts[pos]=p2.xpoints[i];
			newypts[pos]=p2.ypoints[i];
			pos++;
		}
		//and finally from the closest p1 point to its endpoint
		for(int i=closest[0];i<p1.npoints;i++){
			newxpts[pos]=p1.xpoints[i];
			newypts[pos]=p1.ypoints[i];
			pos++;
		}
		return new Polygon(newxpts,newypts,newnpts);
	}
	
	public static int[] findClosestPair(Polygon p1,Polygon p2){
		float mindist2=getDistSqr(p1.xpoints[0],p1.ypoints[0],p2.xpoints[0],p2.ypoints[0]);
		int[] minindices={0,0};
		for(int i=0;i<p1.npoints;i++){
			for(int j=0;j<p2.npoints;j++){
				float temp=getDistSqr(p1.xpoints[i],p1.ypoints[i],p2.xpoints[j],p2.ypoints[j]);
				if(temp<mindist2){
					mindist2=temp;
					minindices[0]=i; minindices[1]=j;
				}
			}
		}
		return minindices;
	}
	
	public static float getDistSqr(float x1,float y1,float x2,float y2){
		return (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);
	}

}
