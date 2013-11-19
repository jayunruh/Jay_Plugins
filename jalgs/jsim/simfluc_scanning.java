/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class simfluc_scanning{
	simflucinterface callingclass;
	public boolean tpe;
	public boolean bleachbyphoton;

	public simfluc_scanning(simflucinterface callingclass1){
		callingclass=callingclass1;
		tpe=false;
		bleachbyphoton=false;
	}

	public float[] do_simfluc_point(int nparticles,int npoints){
		// here is the default single point simulation
		sim_setup set=new sim_setup();
		sim_species ss=new sim_species(set,0);
		sim_multispecies sm=new sim_multispecies(ss,nparticles);
		sm.bleachbyphoton=bleachbyphoton;
		sim_calcintensity ci=new sim_calcintensity(set);
		float[] back=new float[1];
		back[0]=0.0f;
		float[] result_traj=new float[npoints];
		for(int i=0;i<npoints;i++){
			sm.perform_time_step();
			result_traj[i]=(ci.calc_point(sm,back))[0];
			show_progress(i,npoints);
		}
		return result_traj;
	}

	public float[][] do_simfluc_point(float boxsize,int boxpixels,float w0,float z0,int npoints,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[] kbleach,int[] bsteps){
		float[][] kbleach2=new float[bright.length][bright[0].length];
		int[][] bsteps2=new int[bright.length][bright[0].length];
		for(int i=0;i<bright.length;i++){
			for(int j=0;j<bright[0].length;j++){
				kbleach2[i][j]=kbleach[i];
				bsteps2[i][j]=bsteps[i];
			}
		}
		return do_simfluc_point(boxsize,boxpixels,w0,z0,npoints,pixeltime,back,yflowrate,confineindex,startcenter,restrict,bleached_reappear,noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,
				offset,readstdev,nparticles,bright,D,trate,kbleach2,bsteps2);
	}

	public float[][] do_simfluc_point(float boxsize,int boxpixels,float w0,float z0,int npoints,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[][] kbleach,int[][] bsteps){
		// here we simulate a point with all options
		sim_setup set=new sim_setup(pixeltime,boxsize,boxpixels,boundaryindex,confineindex,restrict,startcenter,yflowrate,gridsize,dfraction,jumpprob,noiseindex,gain,offset,readstdev,w0,z0);
		int nspecies=nparticles.length;
		int nchannels=bright[0].length;
		boolean[] chuse=new boolean[nchannels];
		int nchused=0;
		for(int i=0;i<nchannels;i++){
			for(int j=0;j<nspecies;j++){
				if(bright[j][i]>0.0f){
					chuse[i]=true;
				}
			}
			if(back[i]>0.0f){
				chuse[i]=true;
			}
			if(chuse[i]){
				nchused++;
			}
		}
		sim_species[] sslist=new sim_species[nspecies];
		for(int i=0;i<nspecies;i++){
			sslist[i]=new sim_species(set,D[i],bright[i],kbleach[i],bsteps[i],bleached_reappear,i);
		}
		sim_multispecies sm=new sim_multispecies(sslist,trate,nparticles);
		sm.bleachbyphoton=bleachbyphoton;
		sim_calcintensity ci=new sim_calcintensity(set);
		float[][] result_traj=new float[nchused][npoints];
		for(int i=0;i<npoints;i++){
			sm.perform_time_step();
			float[] temp;
			if(!tpe){
				temp=(ci.calc_point(sm,back));
			}else{
				temp=(ci.calc_point_tpe(sm,back));
			}
			int counter=0;
			for(int j=0;j<nchannels;j++){
				if(chuse[j]){
					result_traj[counter][i]=temp[j];
					counter++;
				}
			}
			show_progress(i,npoints);
		}
		return result_traj;
	}

	public float[][] do_simfluc_line(float boxsize,int boxpixels,float w0,float z0,int npoints,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[] kbleach,int[] bsteps){
		float[][] kbleach2=new float[bright.length][bright[0].length];
		int[][] bsteps2=new int[bright.length][bright[0].length];
		for(int i=0;i<bright.length;i++){
			for(int j=0;j<bright[0].length;j++){
				kbleach2[i][j]=kbleach[i];
				bsteps2[i][j]=bsteps[i];
			}
		}
		return do_simfluc_line(boxsize,boxpixels,w0,z0,npoints,pixeltime,back,yflowrate,confineindex,startcenter,restrict,bleached_reappear,noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,
				offset,readstdev,nparticles,bright,D,trate,kbleach2,bsteps2);
	}

	public float[][] do_simfluc_line(float boxsize,int boxpixels,float w0,float z0,int nlines,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[][] kbleach,int[][] bsteps){
		// here we simulate a line with all options
		sim_setup set=new sim_setup(pixeltime,boxsize,boxpixels,boundaryindex,confineindex,restrict,startcenter,yflowrate,gridsize,dfraction,jumpprob,noiseindex,gain,offset,readstdev,w0,z0);
		int nspecies=nparticles.length;
		int nchannels=bright[0].length;
		boolean[] chuse=new boolean[nchannels];
		int nchused=0;
		for(int i=0;i<nchannels;i++){
			for(int j=0;j<nspecies;j++){
				if(bright[j][i]>0.0f){
					chuse[i]=true;
				}
			}
			if(back[i]>0.0f){
				chuse[i]=true;
			}
			if(chuse[i]){
				nchused++;
			}
		}
		sim_species[] sslist=new sim_species[nspecies];
		for(int i=0;i<nspecies;i++){
			sslist[i]=new sim_species(set,D[i],bright[i],kbleach[i],bsteps[i],bleached_reappear,i);
		}
		sim_multispecies sm=new sim_multispecies(sslist,trate,nparticles);
		sm.bleachbyphoton=bleachbyphoton;
		sim_calcintensity ci=new sim_calcintensity(set);
		float[][] result_traj=new float[nchused][nlines*boxpixels];
		for(int i=0;i<nlines;i++){
			for(int j=0;j<boxpixels;j++){
				sm.perform_time_step();
				float[] temp;
				if(!tpe){
					temp=(ci.calc_point(sm,back,j));
				}else{
					temp=(ci.calc_point_tpe(sm,back,j));
				}
				int counter=0;
				for(int k=0;k<nchannels;k++){
					if(chuse[k]){
						result_traj[counter][j+i*boxpixels]=temp[k];
						counter++;
					}
				}
			}
			show_progress(i,nlines);
		}
		return result_traj;
	}

	public float[][][] do_simfluc_image(float boxsize,int boxpixels,float w0,float z0,int npoints,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[] kbleach,int[] bsteps){
		float[][] kbleach2=new float[bright.length][bright[0].length];
		int[][] bsteps2=new int[bright.length][bright[0].length];
		for(int i=0;i<bright.length;i++){
			for(int j=0;j<bright[0].length;j++){
				kbleach2[i][j]=kbleach[i];
				bsteps2[i][j]=bsteps[i];
			}
		}
		return do_simfluc_image(boxsize,boxpixels,w0,z0,npoints,pixeltime,back,yflowrate,confineindex,startcenter,restrict,bleached_reappear,noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,
				offset,readstdev,nparticles,bright,D,trate,kbleach2,bsteps2);
	}

	public float[][][] do_simfluc_image(float boxsize,int boxpixels,float w0,float z0,int nframes,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[][] kbleach,int[][] bsteps){
		// here we simulate a line with all options
		sim_setup set=new sim_setup(pixeltime,boxsize,boxpixels,boundaryindex,confineindex,restrict,startcenter,yflowrate,gridsize,dfraction,jumpprob,noiseindex,gain,offset,readstdev,w0,z0);
		int nspecies=nparticles.length;
		int nchannels=bright[0].length;
		boolean[] chuse=new boolean[nchannels];
		int nchused=0;
		for(int i=0;i<nchannels;i++){
			for(int j=0;j<nspecies;j++){
				if(bright[j][i]>0.0f){
					chuse[i]=true;
				}
			}
			if(back[i]>0.0f){
				chuse[i]=true;
			}
			if(chuse[i]){
				nchused++;
			}
		}
		sim_species[] sslist=new sim_species[nspecies];
		for(int i=0;i<nspecies;i++){
			sslist[i]=new sim_species(set,D[i],bright[i],kbleach[i],bsteps[i],bleached_reappear,i);
		}
		sim_multispecies sm=new sim_multispecies(sslist,trate,nparticles);
		sm.bleachbyphoton=bleachbyphoton;
		sim_calcintensity ci=new sim_calcintensity(set);
		float[][][] results=new float[nchused][nframes][boxpixels*boxpixels];
		for(int i=0;i<nframes;i++){
			for(int j=0;j<boxpixels;j++){
				for(int k=0;k<boxpixels;k++){
					sm.perform_time_step();
					float[] temp;
					if(!tpe){
						temp=(ci.calc_point(sm,back,k,j));
					}else{
						temp=(ci.calc_point_tpe(sm,back,k,j));
					}
					int counter=0;
					for(int l=0;l<nchannels;l++){
						if(chuse[l]){
							results[counter][i][k+j*boxpixels]=temp[l];
							counter++;
						}
					}
				}
			}

			show_progress(i,nframes);
		}
		return results;
	}

	public void show_progress(int i,int end){
		callingclass.showprogress(i,end);
	}

}
