/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class simfluc_camera{
	simflucinterface callingclass;
	public boolean bleachbyphoton;

	public simfluc_camera(simflucinterface callingclass1){
		callingclass=callingclass1;
		bleachbyphoton=false;
	}

	public float[][] do_simfluc_camera(int nparticles,int nframes){
		// here is the default camera simulation
		sim_setup set=new sim_setup();
		sim_species ss=new sim_species(set,0);
		sim_multispecies sm=new sim_multispecies(ss,nparticles);
		sm.bleachbyphoton=bleachbyphoton;
		sim_calcintensity ci=new sim_calcintensity(set);
		float[] back=new float[1];
		back[0]=0.0f;
		float[][] result_traj=new float[nframes][set.boxpixels*set.boxpixels];
		for(int i=0;i<nframes;i++){
			sm.perform_time_step();
			result_traj[i]=(ci.calc_image(sm,back))[0];
			show_progress(i,nframes);
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

	public float[][] do_simfluc_line(float boxsize,int boxpixels,float w0,float z0,int nlines,float frametime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[][] kbleach,int[][] bsteps){
		// here we simulate a line with all options
		sim_setup set=new sim_setup(frametime,boxsize,boxpixels,boundaryindex,confineindex,restrict,startcenter,yflowrate,gridsize,dfraction,jumpprob,noiseindex,gain,offset,readstdev,w0,z0);
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
			sm.perform_time_step();
			float[][] temp=ci.calc_line(sm,back);
			int counter=0;
			for(int k=0;k<nchannels;k++){
				if(chuse[k]){
					System.arraycopy(temp[k],0,result_traj[counter],i*boxpixels,boxpixels);
					counter++;
				}
			}
			show_progress(i,nlines);
		}
		return result_traj;
	}

	public float[][][] do_simfluc_lineimage(float boxsize,int boxpixels,float w0,float z0,int npoints,float pixeltime,float[] back,float yflowrate,int confineindex,boolean startcenter,
			boolean restrict,boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,
			float[][] bright,float[] D,float[][] trate,float[] kbleach,int[] bsteps){
		float[][] kbleach2=new float[bright.length][bright[0].length];
		int[][] bsteps2=new int[bright.length][bright[0].length];
		for(int i=0;i<bright.length;i++){
			for(int j=0;j<bright[0].length;j++){
				kbleach2[i][j]=kbleach[i];
				bsteps2[i][j]=bsteps[i];
			}
		}
		return do_simfluc_lineimage(boxsize,boxpixels,w0,z0,npoints,pixeltime,back,yflowrate,confineindex,startcenter,restrict,bleached_reappear,noiseindex,boundaryindex,gridsize,dfraction,jumpprob,
				gain,offset,readstdev,nparticles,bright,D,trate,kbleach2,bsteps2);
	}

	public float[][][] do_simfluc_lineimage(float boxsize,int boxpixels,float w0,float z0,int nframes,float frametime,float[] back,float yflowrate,int confineindex,boolean startcenter,
			boolean restrict,boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,
			float[][] bright,float[] D,float[][] trate,float[][] kbleach,int[][] bsteps){
		// here we simulate a line with all options
		sim_setup set=new sim_setup(frametime,boxsize,boxpixels,boundaryindex,confineindex,restrict,startcenter,yflowrate,gridsize,dfraction,jumpprob,noiseindex,gain,offset,readstdev,w0,z0);
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
		float[][][] result_traj=new float[nchused][nframes][boxpixels*boxpixels];
		for(int i=0;i<nframes;i++){
			for(int j=0;j<boxpixels;j++){
				sm.perform_time_step();
				float[][] temp=ci.calc_line(sm,back,j);
				int counter=0;
				for(int k=0;k<nchannels;k++){
					if(chuse[k]){
						System.arraycopy(temp[k],0,result_traj[counter][i],j*boxpixels,boxpixels);
						counter++;
					}
				}
			}
			show_progress(i,nframes);
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

	public float[][][] do_simfluc_image(float boxsize,int boxpixels,float w0,float z0,int nframes,float frametime,float[] back,float yflowrate,int confineindex,boolean startcenter,boolean restrict,
			boolean bleached_reappear,int noiseindex,int boundaryindex,float gridsize,float dfraction,float jumpprob,float gain,float offset,float readstdev,int[] nparticles,float[][] bright,
			float[] D,float[][] trate,float[][] kbleach,int[][] bsteps){
		// here we simulate a line with all options
		sim_setup set=new sim_setup(frametime,boxsize,boxpixels,boundaryindex,confineindex,restrict,startcenter,yflowrate,gridsize,dfraction,jumpprob,noiseindex,gain,offset,readstdev,w0,z0);
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
		// callingclass.showmessage(""+sm.nspecies+" , "+sm.nchannels);
		// callingclass.showmessage(""+sm.trate2[0][1]+" , "+trate[0][1]);
		sim_calcintensity ci=new sim_calcintensity(set);
		float[][][] result_traj=new float[nchused][nframes][boxpixels*boxpixels];
		for(int i=0;i<nframes;i++){
			sm.perform_time_step();
			float[][] temp=ci.calc_image(sm,back);
			int counter=0;
			for(int k=0;k<nchannels;k++){
				if(chuse[k]){
					result_traj[counter][i]=temp[k];
					counter++;
				}
			}
			show_progress(i,nframes);
		}
		return result_traj;
	}

	public void show_progress(int i,int end){
		callingclass.showprogress(i,end);
	}

}
