/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class sim_setup{
	// this class contains the information about instrument setup
	public float frametime,frametime2,fboxpixels,boxsize,pixelsize;
	public float gridsize,gridsizepixels,dfraction,jumpprob,yflowrate,yflowrate2;
	public float gain,offset,readstdev,w0,w0pixels,z0,z0pixels;
	public int boundaryindex,boxpixels,confineindex,noiseindex;
	public boolean startcenter,restrict;
	public rngs random;

	public static final String[] confine_options={"None","xy plane","xz plane","1D xz"};
	public static final String[] boundary_options={"Periodic","Reflective","New Particle"};
	public static final String[] noise_options={"Poisson","Analog","None"};

	public sim_setup(){
		frametime=10.0f;
		boxpixels=64;
		boxsize=3.2f;
		restrict=false;
		yflowrate=0.0f;
		boundaryindex=0;
		confineindex=0;
		startcenter=false;
		gridsize=0.0f;
		dfraction=0.0f;
		jumpprob=0.0f;
		frametime2=frametime*(float)(1.0e-6);
		fboxpixels=boxpixels;
		pixelsize=boxsize/fboxpixels;
		gridsizepixels=gridsize/pixelsize;
		yflowrate2=yflowrate*(frametime2/pixelsize);
		noiseindex=0;
		gain=0.0f;
		offset=0.0f;
		readstdev=0.0f;
		w0=0.17f;
		z0=0.7f;
		w0pixels=w0/pixelsize;
		z0pixels=z0/pixelsize;
		random=new rngs();
	}

	public sim_setup(float frametime1,float boxsize1,int boxpixels1,int boundaryindex1,int confineindex1,float w01,float z01){
		frametime=frametime1;
		boxpixels=boxpixels1;
		boxsize=boxsize1;
		restrict=false;
		yflowrate=0.0f;
		boundaryindex=boundaryindex1;
		confineindex=confineindex1;
		startcenter=false;
		gridsize=0.0f;
		dfraction=0.0f;
		jumpprob=0.0f;
		frametime2=frametime*(float)(1.0e-6);
		fboxpixels=boxpixels;
		pixelsize=boxsize/fboxpixels;
		gridsizepixels=gridsize/pixelsize;
		yflowrate2=yflowrate*frametime2;
		noiseindex=0;
		gain=0.0f;
		offset=0.0f;
		readstdev=0.0f;
		w0=w01;
		z0=z01;
		w0pixels=w0/pixelsize;
		z0pixels=z0/pixelsize;
		random=new rngs();
	}

	public sim_setup(float frametime1,float boxsize1,int boxpixels1,int boundaryindex1,int confineindex1,boolean restrict1,boolean startcenter1,float yflowrate1,float gridsize1,float dfraction1,
			float jumpprob1,int noiseindex1,float gain1,float offset1,float readstdev1,float w01,float z01){
		frametime=frametime1;
		boxpixels=boxpixels1;
		boxsize=boxsize1;
		restrict=restrict1;
		yflowrate=yflowrate1;
		boundaryindex=boundaryindex1;
		confineindex=confineindex1;
		startcenter=startcenter1;
		gridsize=gridsize1;
		dfraction=dfraction1;
		jumpprob=jumpprob1;
		frametime2=frametime*(float)(1.0e-6);
		fboxpixels=boxpixels;
		pixelsize=boxsize/fboxpixels;
		gridsizepixels=gridsize/pixelsize;
		yflowrate2=yflowrate*frametime2;
		noiseindex=noiseindex1;
		gain=gain1;
		offset=offset1;
		readstdev=readstdev1;
		w0=w01;
		z0=z01;
		w0pixels=w0/pixelsize;
		z0pixels=z0/pixelsize;
		random=new rngs();
	}

}
