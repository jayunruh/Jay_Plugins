/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class sim_particle{
	public float[] coords;
	public int[] bleach_state;
	public boolean stuck;
	public int speciesID;
	public float[] photonsused;

	public sim_particle(float[] coords1,int[] bleach_state1,boolean stuck1,int speciesID1){
		coords=coords1;
		bleach_state=bleach_state1;
		stuck=stuck1;
		speciesID=speciesID1;
		photonsused=new float[bleach_state.length];
	}

	public void get_from_particle(sim_particle sp){
		coords=sp.coords;
		bleach_state=sp.bleach_state;
		stuck=sp.stuck;
		speciesID=sp.speciesID;
		sp=null;
	}

	public void perform_time_step(sim_multispecies sm){
		// start by looking for interspecies transitions
		transitionlabel: if(sm.transitions){
			float[] trate=sm.trate2[speciesID];
			for(int i=0;i<sm.nspecies;i++){
				if(i!=speciesID&&trate[i]>0.0f){
					if(sm.random.expdev(1.0/trate[i])<1.0){
						speciesID=i;
						break transitionlabel;
					}
				}
			}
		}
		// now look for bleaching events
		for(int i=0;i<bleach_state.length;i++){
			double kbleach=sm.sslist[speciesID].kbleach2[i];
			if(bleach_state[i]>0){
				if(kbleach>0.0){
					boolean bleachbyphoton=sm.bleachbyphoton;
					if(bleachbyphoton){
						// here kbleach is the average inverse number of photons
						// required to bleach
						// calculate the number of photons for this molecule
						if(sm.random.expdev(1.0/kbleach)<photonsused[i]){
							bleach_state[i]--;
							if(sm.sslist[speciesID].bleached_reappear){
								get_from_particle(sm.sslist[speciesID].get_molecule());
							}
						}
					}else{
						// here we bleach irrespective of the photon budget
						if(sm.random.expdev(1.0/kbleach)<1.0){
							bleach_state[i]--;
							if(sm.sslist[speciesID].bleached_reappear){
								get_from_particle(sm.sslist[speciesID].get_molecule());
							}
						}
					}
				}
			}
		}
		// finally move the particle
		sm.sslist[speciesID].mp.do_move_particle(this);
	}

	public float get_brightness(sim_multispecies sm,int channel){
		return sm.sslist[speciesID].bright2[channel]*((float)bleach_state[channel]/(float)sm.sslist[speciesID].bsteps[channel]);
	}
}
