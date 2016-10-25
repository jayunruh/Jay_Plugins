/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import java.util.ArrayList;
import java.util.List;

public class tracker{
	// these are static utility methods for tracking based on known object
	// positions
	// see the two findblobs classes for methods to find the object positions
	// eventually will add progeny tracking
	public int linkdelay,maxtrajs;
	public float linkrange;
	public float zratio;
	private int nclosed;
	private boolean[] closed;
	private int[] tmleft;

	public tracker(float linkrange,int linkdelay){
		this.linkrange=linkrange;
		this.linkdelay=linkdelay;
		this.zratio=1.0f;
		maxtrajs=100000;
	}

	public List<List<float[]>> track2D(List<List<float[]>> params){
		// the incoming params array must have x and y positions as its first
		// two values
		// each sublist must be frame specific
		// on output the list contains a list of parameter arrays for each
		// trajectory
		List<List<float[]>> trajlist=new ArrayList<List<float[]>>();
		closed=new boolean[maxtrajs];
		tmleft=new int[maxtrajs];
		for(int i=0;i<maxtrajs;i++)
			tmleft[i]=linkdelay;
		nclosed=0;
		for(int i=0;i<params.size();i++){
			add_objects(params.get(i),true,trajlist,i);
		}
		return trajlist;
	}

	public List<List<float[]>> track2D(track_interface ti){
		// the incoming params array must have x and y positions as its first
		// two values
		// this version gets the parameters for each frame through the
		// track_interface
		// on output the list contains a list of parameter arrays for each
		// trajectory
		// the last parameter is the first frame the object is detected
		List<List<float[]>> trajlist=new ArrayList<List<float[]>>();
		closed=new boolean[maxtrajs];
		tmleft=new int[maxtrajs];
		for(int i=0;i<maxtrajs;i++)
			tmleft[i]=linkdelay;
		nclosed=0;
		for(int i=0;i<ti.getNFrames();i++){
			int[] assign=add_objects(ti.getNextFrameParams(),true,trajlist,i);
			ti.put_assignments(assign);
			ti.show_progress(i,ti.getNFrames());
		}
		return trajlist;
	}
	
	public List<List<float[]>> track3D(track_interface ti){
		// the incoming params array must have x and y and z positions as its first three values
		// this version gets the parameters for each frame through the track_interface
		// on output the list contains a list of parameter arrays for each trajectory
		// the last parameter is the first frame the object is detected
		List<List<float[]>> trajlist=new ArrayList<List<float[]>>();
		closed=new boolean[maxtrajs];
		tmleft=new int[maxtrajs];
		for(int i=0;i<maxtrajs;i++)
			tmleft[i]=linkdelay;
		nclosed=0;
		for(int i=0;i<ti.getNFrames();i++){
			int[] assign=add_objects(ti.getNextFrameParams(),false,trajlist,i);
			ti.put_assignments(assign);
			ti.show_progress(i,ti.getNFrames());
		}
		return trajlist;
	}

	private int[] add_objects(List<float[]> params,boolean twoD,List<List<float[]>> trajlist,int currframe){
		int nnewobj=params.size();
		int[] assignments=new int[nnewobj]; // this array contains which trajectory the object was assigned to
		int noldobj=trajlist.size()-nclosed;
		float[][] dist=new float[noldobj][nnewobj];
		int[] indices=new int[noldobj];
		int counter=0;
                //distList added by Chris Wood 10-24-2016
                ArrayList<List<Integer>> distList = new ArrayList<>();
		for(int i=0;i<trajlist.size();i++){
			if(!closed[i]){
				indices[counter]=i;
				List<float[]> traj=trajlist.get(i);
				float[] tp=traj.get(traj.size()-1);
				for(int j=0;j<nnewobj;j++){
					dist[counter][j]=distance(tp,params.get(j),twoD);
					if(dist[counter][j]>linkrange)
						dist[counter][j]=-1.0f;
                                        else {
                                                List<Integer> cj = new ArrayList<>();
                                                cj.add(counter);
                                                cj.add(j);
                                                distList.add(cj);
                                        }
				}
				counter++;
			}
		}
		boolean[] newavail=new boolean[nnewobj];
		for(int i=0;i<nnewobj;i++)
			newavail[i]=true;
		boolean[] oldavail=new boolean[noldobj];
		for(int i=0;i<noldobj;i++)
			oldavail[i]=true;
		// now repeatedly find the global minimum distance, making pairs unavailable as we find them
		// once we can't find any distances below the linkrange we cut our losses and move on
		for(int i=0;i<noldobj;i++){

			int[] pair=find_global_min(dist, distList, oldavail,newavail);
			if(pair==null)
				break;
			float[] oldparams=trajlist.get(indices[pair[0]]).get(0);
			float discindex=oldparams[oldparams.length-1];
			(trajlist.get(indices[pair[0]])).add(add_to_array(params.get(pair[1]),discindex));
			assignments[pair[1]]=indices[pair[0]]+1;
		}
		for(int i=0;i<noldobj;i++){
			if(oldavail[i]){
				if(tmleft[i]<1){
					closed[indices[i]]=true;
					nclosed++;
				}else{
					tmleft[i]--;
				}
			}else{
				tmleft[i]=linkdelay;
			}
		}
		for(int i=0;i<nnewobj;i++){
			if(newavail[i]){
				List<float[]> temp=new ArrayList<float[]>();
				temp.add(add_to_array(params.get(i),currframe));
				trajlist.add(temp);
				assignments[i]=trajlist.size();
			}
		}
		return assignments;
	}

	private float[] add_to_array(float[] arr,float addval){
		float[] temp=new float[arr.length+1];
		System.arraycopy(arr,0,temp,0,arr.length);
		temp[arr.length]=addval;
		return temp;
	}

	private float distance(float[] coords1,float[] coords2,boolean twoD){
		if(twoD){
			return (float)Math.sqrt((coords2[0]-coords1[0])*(coords2[0]-coords1[0])+(coords2[1]-coords1[1])*(coords2[1]-coords1[1]));
		}else{
			return (float)Math.sqrt((coords2[0]-coords1[0])*(coords2[0]-coords1[0])+(coords2[1]-coords1[1])*(coords2[1]-coords1[1])+zratio*zratio*(coords2[2]-coords1[2])*(coords2[2]-coords1[2]));
		}
	}

	private int[] find_global_min(float[][] dist, ArrayList<List<Integer>> distList,
                                    boolean[] oldavail,boolean[] newavail){

		float min=-1.0f;
		int oldmin=0;
		int newmin=0;
		boolean first=true;
		//for(int i=0;i<dist.length;i++){
                //iterate over the distList rather than the whole dist array
                for (List<Integer> xij : distList) {
                        int i = xij.get(0);
                        int j = xij.get(1);
			if(oldavail[i]){
//				for(int j=0;j<dist[i].length;j++){
				if(newavail[j]&&dist[i][j]>=0.0f){
					if(first){
						min=dist[i][j];
						oldmin=i;
						newmin=j;
						first=false;
					}else {
						if(dist[i][j]<min){
							min=dist[i][j];
							oldmin=i;
							newmin=j;
						}
					}
				//}
				}
			}
		}
		if(min>=0.0f){
			oldavail[oldmin]=false;
			newavail[newmin]=false;
			int[] temp={oldmin,newmin};
			return temp;
		}else{
			return null;
		}
	}
}
