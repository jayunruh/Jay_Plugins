/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class sim_common_setups{
	public static final String[] samples={"fluorescein","cy5","egfp","mcherry"};
	public static final float[][] Dvals={{300.0f,300.0f,90.0f,90.0f},{60.0f,60.0f,15.0f,15.0f},{5.0f,5.0f,1.0f,1.0f},{5.0f,5.0f,1.0f,1.0f}};
	public static final float[][] brightvals={{10000.0f,0.0f,3000.0f,0.0f},{1000.0f,10000.0f,300.0f,3000.0f}};
	public static final String[] locations={"solution","cytosol","membranexy","membranexz"};

	public sim_multispecies get_common_setups(int sampleindex,int locationindex,int clustersize){
		int confineindex=1;
		int num=100;
		if(locationindex<2){
			confineindex=0;
			num=200;
		}
		if(locationindex>2){
			confineindex=2;
		}
		sim_setup set=new sim_setup(10.0f,3.2f,64,0,confineindex,0.17f,0.7f);
		float[] tempbright={brightvals[0][sampleindex]*clustersize,brightvals[1][sampleindex]*clustersize};
		sim_species ss=new sim_species(set,Dvals[locationindex][sampleindex],tempbright,new float[1],new int[]{1},false,0);
		return new sim_multispecies(ss,num);
	}

	public sim_multispecies get_common_setups(int sampleindex1,int sampleindex2,int clustersize1,int clustersize2,int locationindex,double finteracting){
		int confineindex=1;
		int num=100;
		if(locationindex<2){
			confineindex=0;
			num=200;
		}
		if(locationindex>2){
			confineindex=2;
		}
		sim_setup set=new sim_setup(10.0f,3.2f,64,0,confineindex,0.17f,0.7f);
		float[] tempbright1={brightvals[0][sampleindex1]*clustersize1,brightvals[1][sampleindex1]*clustersize1};
		float[] tempbright2={brightvals[0][sampleindex2]*clustersize2,brightvals[1][sampleindex2]*clustersize2};
		float[] tempbright3={tempbright1[0]+tempbright2[0],tempbright1[1]+tempbright2[1]};
		sim_species[] sslist=new sim_species[3];
		float dval=Dvals[locationindex][sampleindex1];
		if(Dvals[locationindex][sampleindex2]<dval){
			dval=Dvals[locationindex][sampleindex2];
		}
		sslist[0]=new sim_species(set,Dvals[locationindex][sampleindex1],tempbright1,new float[1],new int[]{1},false,0);
		sslist[1]=new sim_species(set,Dvals[locationindex][sampleindex2],tempbright2,new float[1],new int[]{1},false,0);
		sslist[0]=new sim_species(set,dval,tempbright3,new float[1],new int[]{1},false,0);
		int numinteracting=(int)(finteracting*num);
		int numnotinteracting=num-numinteracting;
		int[] numarray={numnotinteracting,numnotinteracting,numinteracting};
		return new sim_multispecies(sslist,null,numarray);
	}

}
