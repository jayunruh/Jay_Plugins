/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class sim_multispecies{
	// this class is an umbrella class for simulation species
	// it contains the particles (sps) and the species definitions (sslist) as
	// well as species interchange info
	public sim_species[] sslist;
	public float[][] trate,trate2;
	public sim_particle[] sps;
	public int nspecies,nchannels;
	public boolean transitions;
	public rngs random;
	public sim_setup set;
	public boolean bleachbyphoton;

	public sim_multispecies(sim_species ss1,int num1){
		set=ss1.set;
		random=set.random;
		sslist=new sim_species[1];
		sslist[0]=ss1;
		trate=new float[1][1];
		trate[0][0]=0.0f;
		trate2=new float[1][1];
		trate2[0][0]=0.0f;
		transitions=false;
		nspecies=1;
		nchannels=ss1.bright.length;
		sps=new sim_particle[num1];
		sps=ss1.get_molecules(num1);
	}

	public sim_multispecies(sim_species[] ss1,float[][] trate1,int[] num1){
		set=ss1[0].set;
		sslist=ss1;
		nspecies=sslist.length;
		nchannels=sslist[0].bright.length;
		trate2=new float[nspecies][nspecies];
		float tottrate=0.0f;
		if(trate1!=null){
			trate=trate1;
			for(int i=0;i<trate.length;i++){
				for(int j=0;j<trate[0].length;j++){
					tottrate+=trate[i][j];
					trate2[i][j]=trate[i][j]*set.frametime2;
				}
			}
		}
		if(tottrate<=0.0f){
			transitions=false;
		}else{
			transitions=true;
		}
		random=set.random;
		int totnum=0;
		for(int i=0;i<nspecies;i++){
			totnum+=num1[i];
		}
		sps=new sim_particle[totnum];
		int counter=0;
		for(int i=0;i<nspecies;i++){
			for(int j=0;j<num1[i];j++){
				sps[counter]=sslist[i].get_molecule();
				counter++;
			}
		}
	}

	public void add_species(sim_species ss1,float[][] trate1,int num){
		// the trate array must be 2xn in dimension. The first column contains
		// the rates of transition
		// from other species. The second column has the rates of transition to
		// other species
		sim_species[] tempsslist=new sim_species[nspecies+1];
		float[][] temptrate=new float[nspecies+1][nspecies+1];
		for(int i=0;i<nspecies;i++){
			tempsslist[i]=sslist[i];
			System.arraycopy(trate[i],0,temptrate[i],0,nspecies);
			if(trate1!=null){
				temptrate[i][nspecies]=trate1[0][i];
			}
		}
		tempsslist[nspecies]=ss1;
		if(trate1!=null){
			System.arraycopy(trate1[1],0,temptrate[nspecies],0,nspecies);
		}
		sslist=tempsslist;
		trate=temptrate;
		nspecies++;
		for(int i=0;i<nspecies;i++){
			for(int j=0;j<nspecies;j++){
				trate2[i][j]=trate[i][j]*set.frametime2;
			}
		}
		sim_particle[] tempsps=new sim_particle[sps.length+num];
		System.arraycopy(sps,0,tempsps,0,sps.length);
		for(int i=0;i<num;i++){
			tempsps[i+sps.length]=ss1.get_molecule();
		}
		sps=tempsps;
	}

	public void change_setup(sim_setup set1){
		set=set1;
		for(int i=0;i<nspecies;i++){
			for(int j=0;j<nspecies;j++){
				trate2[i][j]=trate[i][j]*set.frametime2;
			}
			sslist[i].change_setup(set);
		}
	}

	public void perform_time_step(){
		for(int i=0;i<sps.length;i++){
			sps[i].perform_time_step(this);
		}
	}

}
