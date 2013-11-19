/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class PALMutils{
	float stdev,baseline,gain,maxerr,r;
	gausfunc gf;
	int xsize,ysize;

	public PALMutils(float stdev,float baseline,float gain,int xsize,int ysize){
		this.stdev=stdev;
		this.baseline=baseline;
		this.gain=gain;
		gf=new gausfunc();
		this.xsize=xsize;
		this.ysize=ysize;
		maxerr=2.0f*stdev;
		r=0.625f; // r is the ratio of original pixel size to stdev
	}

	public float[] render_const_size(float[][] molecules,float[] limits,float stdev2,float multiplier){
		// the molecules array has columns with x,y,amp,c2,t,z
		// this method ignores the z and t component
		if(stdev2==0.0f)
			return render_pix_size(molecules,limits,multiplier);
		float xmin=limits[0];
		float xmax=limits[1];
		float ymin=limits[2];
		float ymax=limits[3];
		float scaling=(float)xsize/(xmax-xmin);
		// recalculate ymax to make sure image size is retained appropriately
		ymax=(float)ysize/scaling+ymin;
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float pixelerr=scaling*stdev2;
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*(float)xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*(float)ysize;
			gf.draw_2D_func(image,pixelx,pixely,xsize,ysize,pixelerr,multiplier*molecules[2][i],2.0f);
		}
		return image;
	}

	public float[] render_err_size(float[][] molecules,float[] limits,float multiplier){
		// the molecules array has columns with x,y,amp,c2,t,z
		// this method ignores the z and t component
		float xmin=limits[0];
		float xmax=limits[1];
		float ymin=limits[2];
		float ymax=limits[3];
		float scaling=(float)xsize/(xmax-xmin);
		// recalculate ymax to make sure image size is retained appropriately
		ymax=(float)ysize/scaling+ymin;
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float photons=molecules[2][i]*(stdev*stdev*2.0f*(float)Math.PI);
			photons/=gain;
			float pixelerr=scaling*geterr(photons);
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*(float)xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*(float)ysize;
			float newamp=scaling*scaling/(pixelerr*pixelerr*2.0f*(float)Math.PI);
			gf.draw_2D_func(image,pixelx,pixely,xsize,ysize,pixelerr,multiplier*newamp,2.0f);
		}
		return image;
	}

	public float[] render_pix_size(float[][] molecules,float[] limits,float multiplier){
		// the molecules array has columns with x,y,amp,c2,t,z
		// this method ignores the z and t component
		float xmin=limits[0];
		float xmax=limits[1];
		float ymin=limits[2];
		float ymax=limits[3];
		float scaling=(float)xsize/(xmax-xmin);
		// recalculate ymax to make sure image size is retained appropriately
		ymax=(float)ysize/scaling+ymin;
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float photons=molecules[2][i]*(stdev*stdev*2.0f*(float)Math.PI);
			photons/=gain;
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*(float)xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*(float)ysize;
			int tempx=(int)(pixelx+0.5f);
			int tempy=(int)(pixely+0.5f);
			if(tempx>=0&&tempx<xsize&&tempy>=0&&tempy<ysize){
				image[tempx+tempy*xsize]+=photons;
			}
		}
		return image;
	}

	public float[][][] get_filtered_molecules(float[][] molecules,float[][] filters){
		// here we return two arrays: one filtered and the other antifiltered
		// if the appropriate filtering array is null we don't filter it
		// we take one pass to get the number of filtered items
		boolean[] filtered=new boolean[molecules[0].length];
		int nfiltered=0;
		for(int i=0;i<molecules[0].length;i++){
			for(int j=0;j<filters.length;j++){
				if(filters[j]!=null){
					if(molecules[j][i]<filters[j][0]||molecules[j][i]>filters[j][1]){
						filtered[i]=true; // molecule is outside filter range
						nfiltered++;
						break;
					}
				}
			}
		}
		float[][][] molfilt={new float[molecules.length][molecules[0].length-nfiltered],new float[molecules.length][nfiltered]};
		int counter1=0;
		int counter2=0;
		for(int i=0;i<molecules[0].length;i++){
			if(!filtered[i]){
				for(int j=0;j<molecules.length;j++){
					molfilt[0][j][counter1]=molecules[j][i];
				}
				counter1++;
			}else{
				for(int j=0;j<molecules.length;j++){
					molfilt[1][j][counter2]=molecules[j][i];
				}
				counter2++;
			}
		}
		return molfilt;
	}

	public float[][] get_filtered_molecules2(float[][] molecules,float[][] filters){
		// here we return a single filtered array
		// if the appropriate filtering array is null we don't filter it
		// we take one pass to get the number of filtered items
		boolean[] filtered=new boolean[molecules[0].length];
		int nfiltered=0;
		for(int i=0;i<molecules[0].length;i++){
			for(int j=0;j<filters.length;j++){
				if(filters[j]!=null){
					if(molecules[j][i]<filters[j][0]||molecules[j][i]>filters[j][1]){
						filtered[i]=true; // molecule is outside filter range
						nfiltered++;
						break;
					}
				}
			}
		}
		float[][] molfilt=new float[molecules.length][molecules[0].length-nfiltered];
		int counter1=0;
		for(int i=0;i<molecules[0].length;i++){
			if(!filtered[i]){
				for(int j=0;j<molecules.length;j++){
					molfilt[j][counter1]=molecules[j][i];
				}
				counter1++;
			}
		}
		return molfilt;
	}

	private float geterr(float photons){
		float basephotons=baseline/gain;
		// note that r is the original pixel size divided by stdev
		float var=(stdev*stdev/photons)*(1.0f+r*r/12+4.0f*(float)Math.PI*basephotons*basephotons/(r*photons));
		float err=(float)Math.sqrt(var);
		return Math.min(err,maxerr);
	}

	public static float get_gaus_integral(float amp,float stdev){
		return amp*stdev*stdev*2.0f*(float)Math.PI;
	}

	public static float get_gaus_amp(float integral,float stdev){
		return integral/(stdev*stdev*2.0f*(float)Math.PI);
	}

}
