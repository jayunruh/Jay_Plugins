/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class sim_species implements Cloneable{
	// this is a species for a simulation
	// this class defines the behavior for a species of molecules and provides
	// the information needed to
	// perform time steps on the molecules
	public float D,D2;
	public float[] bright,bright2,kbleach,kbleach2;
	public int speciesID;
	public int[] bsteps;
	public boolean bleached_reappear;
	public rngs random;
	public move_particles mp;
	public sim_setup set;

	public sim_species(sim_setup set1,int speciesID1){
		set=set1;
		random=set.random;
		speciesID=speciesID1;
		D=40.0f;
		bright=new float[1];
		bright[0]=10000.0f;
		bsteps=new int[1];
		bsteps[0]=1;
		kbleach=new float[1];
		kbleach[0]=0.0f;
		bright2=new float[1];
		bright2[0]=bright[0]*set.frametime2;
		D2=D*(set.frametime2/(set.pixelsize*set.pixelsize));
		kbleach2=new float[1];
		kbleach2[0]=kbleach[0]*set.frametime2;
		bleached_reappear=false;
		mp=new move_particles(this,set);
	}

	public sim_species(float D1,float bright1,sim_setup set1,int speciesID1){
		set=set1;
		random=set.random;
		D=D1;
		bright=new float[1];
		bright[0]=bright1;
		bsteps=new int[1];
		bsteps[0]=1;
		kbleach=new float[1];
		kbleach[0]=0.0f;
		bright2=new float[1];
		bright2[0]=bright[0]*set.frametime2;
		D2=D*(set.frametime2/(set.pixelsize*set.pixelsize));
		kbleach2=new float[1];
		kbleach2[0]=kbleach[0]*set.frametime2;
		bleached_reappear=false;
		speciesID=speciesID1;
		mp=new move_particles(this,set);
	}

	public sim_species(sim_setup set1,float D1,float[] bright1,float[] kbleach1,int[] bsteps1,boolean bleached_reappear1,int speciesID1){
		set=set1;
		random=set.random;
		D=D1;
		bright=bright1;
		bsteps=bsteps1;
		for(int i=0;i<bsteps.length;i++)
			if(bsteps[i]<1){
				bsteps[i]=1;
			}
		kbleach=kbleach1;
		bright2=new float[bright.length];
		for(int i=0;i<bright.length;i++){
			bright2[i]=bright[i]*set.frametime2;
		}
		D2=D*(set.frametime2/(set.pixelsize*set.pixelsize));
		kbleach2=new float[kbleach.length];
		for(int i=0;i<kbleach.length;i++)
			kbleach2[i]=kbleach[i]*set.frametime2;
		bleached_reappear=bleached_reappear1;
		speciesID=speciesID1;
		mp=new move_particles(this,set);
	}

	public Object Clone() throws CloneNotSupportedException{
		sim_species clone=(sim_species)super.clone();
		clone.bright=bright.clone();
		clone.bright2=bright2.clone();
		clone.mp=(move_particles)mp.Clone();
		return clone;
	}

	public void change_brightnesses(float[] bright1){
		bright=bright1;
		bright2=new float[bright.length];
		for(int i=0;i<bright.length;i++){
			bright2[i]=bright[i]*set.frametime2;
		}
	}

	public void change_D(float D1){
		D=D1;
		D2=D*(set.frametime2/(set.pixelsize*set.pixelsize));
		mp.D2=D2;
	}

	public void change_kbleach(float[] kbleach1){
		kbleach=kbleach1;
		for(int i=0;i<kbleach.length;i++)
			kbleach2[i]=kbleach[i]*set.frametime2;
	}

	public void change_setup(sim_setup set1){
		set=set1;
		bright2=new float[bright.length];
		for(int i=0;i<bright.length;i++){
			bright2[i]=bright[i]*set.frametime2;
		}
		D2=D*(set.frametime2/(set.pixelsize*set.pixelsize));
		for(int i=0;i<kbleach.length;i++)
			kbleach2[i]=kbleach[i]*set.frametime2;
	}

	public sim_particle[] get_molecules(int num){
		sim_particle[] sps=new sim_particle[num];
		for(int i=0;i<num;i++){
			sps[i]=get_molecule();
		}
		return sps;
	}

	public sim_particle get_molecule(){
		if(!set.startcenter){
			float[] coords=new float[3];
			boolean stuck=false;
			int[] bleach_state=bsteps.clone();
			coords[0]=(float)random.unidev(set.fboxpixels,0.0);
			if(set.confineindex<2){
				coords[1]=(float)random.unidev(set.fboxpixels,0.0);
			}else{
				coords[1]=0.5f*set.fboxpixels;
			}
			if(set.confineindex!=1&&set.confineindex!=3){
				coords[2]=(float)random.unidev(set.fboxpixels,0.0);
			}else{
				coords[2]=0.5f*set.fboxpixels;
			}
			if(set.restrict){
				float testx=coords[0]%(2.0f*set.gridsizepixels);
				float testy=coords[1]%(2.0f*set.gridsizepixels);
				if((testx>set.gridsizepixels&&testy>set.gridsizepixels)||(testx<set.gridsizepixels&&testy<set.gridsizepixels)){
					stuck=true;
				}
			}
			return new sim_particle(coords,bleach_state,stuck,speciesID);
		}else{
			float[] coords=new float[3];
			boolean stuck=false;
			int[] bleach_state=bsteps.clone();
			coords[0]=0.5f*set.fboxpixels;
			coords[1]=0.5f*set.fboxpixels;
			coords[2]=0.5f*set.fboxpixels;
			if(set.restrict){
				float testx=coords[0]%(2.0f*set.gridsizepixels);
				float testy=coords[1]%(2.0f*set.gridsizepixels);
				if((testx>set.gridsizepixels&&testy>set.gridsizepixels)||(testx<set.gridsizepixels&&testy<set.gridsizepixels)){
					stuck=true;
				}
			}
			return new sim_particle(coords,bleach_state,stuck,speciesID);
		}
	}

}
