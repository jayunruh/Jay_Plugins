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
		return render_const_size(molecules,limits,stdev2,multiplier,false);
	}

	public float[] render_const_size(float[][] molecules,float[] limits,float stdev2,float multiplier,boolean ignoreamp){
		// the molecules array has columns with x,y,amp,c2,t,z
		// this method ignores the z and t component
		if(stdev2==0.0f)
			return render_pix_size(molecules,limits,multiplier);
		float xmin=limits[0];
		float xmax=limits[1];
		float ymin=limits[2];
		float ymax=limits[3];
		float scaling=xsize/(xmax-xmin);
		// recalculate ymax to make sure image size is retained appropriately
		ymax=ysize/scaling+ymin;
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float amp=multiplier*molecules[2][i];
			if(ignoreamp) amp=1.0f;
			float pixelerr=scaling*stdev2;
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*ysize;
			gf.draw_2D_func(image,pixelx,pixely,xsize,ysize,pixelerr,amp,2.0f);
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
		float scaling=xsize/(xmax-xmin);
		// recalculate ymax to make sure image size is retained appropriately
		ymax=ysize/scaling+ymin;
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float photons=molecules[2][i]*(stdev*stdev*2.0f*(float)Math.PI);
			photons/=gain;
			float pixelerr=scaling*geterr(photons);
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*ysize;
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
		float scaling=xsize/(xmax-xmin);
		// recalculate ymax to make sure image size is retained appropriately
		ymax=ysize/scaling+ymin;
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float photons=molecules[2][i]*(stdev*stdev*2.0f*(float)Math.PI);
			photons/=gain;
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*ysize;
			int tempx=(int)(pixelx+0.5f);
			int tempy=(int)(pixely+0.5f);
			if(tempx>=0&&tempx<xsize&&tempy>=0&&tempy<ysize){
				image[tempx+tempy*xsize]+=photons;
			}
		}
		return image;
	}
	
	public float[] render_circle_size(float[][] molecules,float[] limits,float multiplier,float radius,boolean ignoreamp){
		// the molecules array has columns with x,y,amp,c2,t,z
		// this method ignores the z and t component
		float xmin=limits[0];
		float xmax=limits[1];
		float ymin=limits[2];
		float ymax=limits[3];
		float scaling=xsize/(xmax-xmin); //in pix/um
		// recalculate ymax to make sure image size is retained appropriately
		ymax=ysize/scaling+ymin;
		float[] circletemp=null;
		float radpix=radius*scaling;
		int halfsize=1+(int)radpix;
		int size=halfsize*2;
		if(radpix>=1.0f){
			circletemp=new float[size*size];
			float rad2=radpix*radpix;
			for(int i=0;i<size;i++){
				int ypos=i-halfsize;
				float ypos2=(float)(ypos*ypos);
				for(int j=0;j<size;j++){
					int xpos=j-halfsize;
					float xpos2=(float)(xpos*xpos);
					if((xpos2+ypos2)<=rad2){
						circletemp[j+i*size]=1.0f;
					}
				}
			}
		}
		float[] image=new float[xsize*ysize];
		for(int i=0;i<molecules[0].length;i++){
			float photons=multiplier;
			if(!ignoreamp){
				photons=multiplier*molecules[2][i]*(stdev*stdev*2.0f*(float)Math.PI);
				photons/=gain;
			}
			float pixelx=((molecules[0][i]-xmin)/(xmax-xmin))*xsize;
			float pixely=((molecules[1][i]-ymin)/(ymax-ymin))*ysize;
			int tempx=(int)(pixelx+0.5f);
			int tempy=(int)(pixely+0.5f);
			if(radpix<1.0f){
				//here we draw a single pixel
    			if(tempx>=0&&tempx<xsize&&tempy>=0&&tempy<ysize){
    				image[tempx+tempy*xsize]+=photons;
    			}
			} else {
				//here we draw the circle template
				for(int k=0;k<size;k++){
					int ypos=tempy-halfsize+k;
					if(ypos<0) continue; if(ypos>=ysize) continue;
					int offset=ypos*xsize;
					int offset2=k*size;
					for(int j=0;j<size;j++){
						int xpos=tempx-halfsize+j;
						if(xpos<0) continue; if(xpos>=xsize) continue;
						image[offset+xpos]+=photons*circletemp[offset2+j];
					}
				}
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
	
	public static float[][] pixels2Molecules(float[] pixels,int width,int height){
		float[][] mols=new float[width][height-1];
		for(int i=0;i<(height-1);i++){
			for(int j=0;j<width;j++){
				mols[j][i]=pixels[j+(i+1)*width];
			}
		}
		return mols;
	}
	
	public static float[] molecules2Pixels(float[][] molecules,float[] metadata){
		int nmols=molecules[0].length;
		int nparams=molecules.length;
		float[] pixels=new float[nparams*(nmols+1)];
		System.arraycopy(metadata,0,pixels,0,metadata.length);
		for(int i=0;i<nmols;i++){
			for(int j=0;j<nparams;j++){
				pixels[j+(i+1)*nparams]=molecules[j][i];
			}
		}
		return pixels;
	}
	
	public static float[][] accumulateMolecules(float[][] molecules,float maxdist,float maxoff){
		// the molecules array has columns with x,y,amp,c2,t,z
		//here we recursively search molecules detected within maxoff time for maxdist
		int nparams=molecules.length;
		int nmols=molecules[0].length;
		float[][] newmols=new float[nparams][nmols];
		int nnewmols=1;
		float maxdist2=maxdist*maxdist;
		for(int i=0;i<nparams;i++) newmols[i][0]=molecules[i][0];
		for(int i=1;i<nmols;i++){
			boolean found=false;
			for(int j=nnewmols-1;j>=0;j--){
				float dt=molecules[4][i]-newmols[4][j];
				if(dt>maxoff){
					break;
				}
				float dist2=getdist2(new float[]{molecules[0][i],molecules[1][i]},new float[]{newmols[0][j],newmols[1][j]});
				if(dist2<=maxdist2){
					newmols[0][j]+=molecules[0][i]; newmols[0][j]*=0.5f;
					newmols[1][j]+=molecules[1][i]; newmols[1][j]*=0.5f;
					newmols[2][j]+=molecules[2][i];
					newmols[3][j]+=molecules[3][i]; newmols[3][j]*=0.5f;
					newmols[4][j]=molecules[4][i];
					newmols[5][j]+=molecules[5][i]; newmols[5][j]*=0.5f;
					found=true;
					break;
				}
			}
			if(!found){
				for(int k=0;k<nparams;k++) newmols[k][nnewmols]=molecules[k][i];
				nnewmols++;
			}
		}
		float[][] newmols2=new float[nparams][nnewmols];
		for(int i=0;i<nparams;i++) System.arraycopy(newmols[i],0,newmols2[i],0,nnewmols);
		return newmols2;
	}
	
	private static float getdist(float[] coords1,float[] coords2){
		float d2=0.0f;
		for(int i=0;i<coords1.length;i++) d2+=(coords2[i]-coords1[i])*(coords2[i]-coords1[i]);
		return (float)Math.sqrt(d2);
	}
	
	private static float getdist2(float[] coords1,float[] coords2){
		float d2=0.0f;
		for(int i=0;i<coords1.length;i++) d2+=(coords2[i]-coords1[i])*(coords2[i]-coords1[i]);
		return d2;
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
