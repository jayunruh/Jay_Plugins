/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import jalgs.jsim.*;
import jguis.*;
import jalgs.jfit.*;

public class nnmf_spectral_unmixing_jru_v1 implements PlugIn {
	float min_val=1.0f;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int type=0;
		if(stack.getPixels(1) instanceof short[]) type=1;
		if(stack.getPixels(1) instanceof float[]) type=2;
		int nchannels=stack.getSize();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int npixels=width*height;
		int nspectra=2;
		float background=100.0f;
		float saturation=65535.0f;
		float lambda=0.0f;
		int preiter=10;
		int iterations=30;
		boolean plot=false;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("#_of_species",nspectra,0);
		gd.addNumericField("background",background,5,15,null);
		gd.addNumericField("saturation level",saturation,5,15,null);
		gd.addNumericField("segregation bias",lambda,5,15,null);
		gd.addNumericField("#_of_iterations",iterations,0);
		gd.addCheckbox("Initialize_From_Plot",plot);
		gd.addCheckbox("Output_Residuals",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		nspectra=(int)gd.getNextNumber();
		background=(float)gd.getNextNumber();
		saturation=(float)gd.getNextNumber();
		lambda=(float)gd.getNextNumber();
		iterations=(int)gd.getNextNumber();
		plot=gd.getNextBoolean();
		boolean outresid=gd.getNextBoolean();
		boolean[] fixed=new boolean[nspectra];
		float sigma=0.1f+((float)nchannels-1.0f)/((float)nspectra+1.0f);
		float centerinc=(float)nchannels/((float)nspectra+1.0f);
		float[] centers=new float[nspectra];
		float[] sigmas=new float[nspectra];
		for(int i=0;i<nspectra;i++){
			centers[i]=centerinc*(float)(i+1);
			sigmas[i]=sigma;
		}
		GenericDialog gd2=new GenericDialog("Fix_Spectra");
		for(int i=0;i<nspectra;i++){
			gd2.addCheckbox("Fix_Spectrum_"+(i+1),fixed[i]);
			gd2.addNumericField("Spectrum_Center_"+(i+1),centers[i],5,15,null);
			gd2.addNumericField("Spectrum_StDev_"+(i+1),sigmas[i],5,15,null);
		}
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		for(int i=0;i<nspectra;i++){
			fixed[i]=gd2.getNextBoolean();
			centers[i]=(float)gd2.getNextNumber();
			sigmas[i]=(float)gd2.getNextNumber();
		}
		boolean[] valid=new boolean[npixels];
		int nvalid=0;
		//at least one channel per pixel must be above background
		//no channels can be saturated
		for(int j=0;j<npixels;j++){
			for(int i=0;i<nchannels;i++){
				if(get_value(stack,i,j,type)>=background){valid[j]=true; nvalid++; break;}
			}
			for(int i=0;i<nchannels;i++){
				if(get_value(stack,i,j,type)>saturation){valid[j]=false; nvalid--; break;}
			}
		}
		float[][] y=new float[nchannels][nvalid]; //this is the data stack
		int counter=0;
		for(int j=0;j<npixels;j++){
			if(valid[j]){
				for(int i=0;i<nchannels;i++){
					y[i][counter]=get_value(stack,i,j,type)-background;
					if(y[i][counter]<min_val) y[i][counter]=min_val;
				}
				counter++;
			}
		}
		float[][] x=new float[nspectra][nvalid]; //x contains the concentrations
		float[][] a=new float[nchannels][nspectra]; //a contains the spectra
		init_x(x);
		if(!init_a(a,plot,centers,sigmas)) return;
		PlotWindow4 pw=new PlotWindow4("Fit Spectra","channel","intensity",transpose(a),null);
		pw.draw();
		for(int i=0;i<preiter;i++){
			IJ.showStatus("Initializing Concentrations");
			update_x(a,x,y,0.0f);
		}
		for(int i=0;i<iterations;i++){
			normalize(a,x);
			update_a(a,x,y,fixed);
			update_x(a,x,y,0.0f);
			float[][] temp=transpose(a);
			for(int j=0;j<nspectra;j++) pw.updateSeries(temp[j],j,true);
			IJ.showProgress(i,iterations);
		}
		ImageStack xstack=new ImageStack(width,height);
		for(int i=0;i<nspectra;i++){
			int counter2=0;
			float[] temp=new float[npixels];
			for(int j=0;j<npixels;j++){
				if(valid[j]){
					temp[j]=x[i][counter2];
					counter2++;
				}
			}
			xstack.addSlice("",temp);
		}
		ImagePlus impout=new ImagePlus("Unmixed",xstack);
		impout.setOpenAsHyperStack(true);
		impout.setDimensions(nspectra,1,1);
		impout=new CompositeImage(impout,CompositeImage.COLOR);
		impout.show();
		if(outresid){
			ImageStack resid=new ImageStack(width,height);
			float[][] spectra=transpose(a);
			//Object[] xdata=jutils.stack2array(xstack);
			for(int i=0;i<nchannels;i++){
				float[] temp2=new float[width*height];
				for(int j=0;j<width*height;j++){
					float fit=background;
					//for(int k=0;k<nspectra;k++) fit+=spectra[k][i]*((float[])xdata[k])[j];
					for(int k=0;k<nspectra;k++) fit+=spectra[k][i]*get_value(xstack,k,j,2);
					temp2[j]=get_value(stack,i,j,type)-fit;
				}
				resid.addSlice("",temp2);
			}
			ImagePlus impr=new ImagePlus("Residuals",resid);
			impr.setOpenAsHyperStack(true);
			impr.setDimensions(nchannels,1,1);
			new CompositeImage(impr,CompositeImage.COLOR).show();
		}
	}

	public float get_value(ImageStack stack,int slice,int pixel,int type){
		if(type==0) return (float)(((byte[])stack.getPixels(slice+1))[pixel]&0xff);
		if(type==1) return (float)(((short[])stack.getPixels(slice+1))[pixel]&0xffff);
		return ((float[])stack.getPixels(slice+1))[pixel];
	}

	public float[][] transpose(float[][] input){
		float[][] output=new float[input[0].length][input.length];
		for(int i=0;i<input.length;i++){
			for(int j=0;j<input[0].length;j++){
				output[j][i]=input[i][j];
			}
		}
		return output;
	}

	public void normalize(float[][] a,float[][] x){
		float[] sum=new float[a[0].length];
		for(int i=0;i<a[0].length;i++){
			for(int j=0;j<a.length;j++) sum[i]+=a[j][i];
		}
		for(int i=0;i<a[0].length;i++){
			for(int j=0;j<a.length;j++) a[j][i]/=sum[i];
			for(int j=0;j<x[0].length;j++) x[i][j]*=sum[i];
		}
	}

	public void update_a(float[][] a,float[][] x,float[][] y,boolean[] fixed){
		int nch=a.length; int nspec=a[0].length; int npix=x[0].length;
		//calculate the estimated signal x*a
		float[][] ax=new float[nch][npix];
		for(int i=0;i<npix;i++){
			for(int j=0;j<nch;j++){
				for(int k=0;k<nspec;k++){
					ax[j][i]+=x[k][i]*a[j][k];
				}
			}
		}
		//calculate the denominator
		float[] sumRowX=new float[nspec];
		for(int i=0;i<nspec;i++){
			for(int j=0;j<npix;j++){
				sumRowX[i]+=x[i][j];
			}
		}
		for(int i=0;i<nspec;i++){
			if(!fixed[i]){
				for(int j=0;j<nch;j++){
					float f=0.0f;
					for(int k=0;k<npix;k++){
						f+=y[j][k]*x[i][k]/ax[j][k];
					}
					a[j][i]*=f/sumRowX[i];
				}
			}
		}
	}

	public void update_x(float[][] a,float[][] x,float[][] y,float lambda){
		int nch=a.length; int nspec=a[0].length; int npix=x[0].length;
		//calculate the estimated signal x*a
		float[][] ax=new float[nch][npix];
		for(int i=0;i<npix;i++){
			for(int j=0;j<nch;j++){
				for(int k=0;k<nspec;k++){
					ax[j][i]+=x[k][i]*a[j][k];
				}
			}
		}
		//calculate the segregation bias
		float[][] segbiasterm=null;
		if(lambda>0.0f){
			segbiasterm=new float[nspec][npix];
			for(int i=0;i<npix;i++){
				double sumx=0.0f; double sumx2=0.0f;
				for(int j=0;j<nspec;j++){
					sumx+=x[j][i];
					sumx2+=x[j][i]*x[j][i];
				}
				sumx/=sumx2;
				sumx2=Math.sqrt(sumx2);
				sumx/=sumx2; //this term is now sumx/(sumx2^1.5)
				for(int j=0;j<nspec;j++){
					segbiasterm[j][i]=lambda*(float)(sumx*x[j][i]-1.0/sumx2);
				}
			}
		}
		//calculate the denominator
		float[] sumColA=new float[nspec];
		for(int i=0;i<nspec;i++){
			for(int j=0;j<nch;j++){
				sumColA[i]+=a[j][i];
			}
		}
		//now do update
		for(int i=0;i<nspec;i++){
			for(int j=0;j<npix;j++){
				float f=0.0f;
				for(int k=0;k<nch;k++){
					f+=y[k][j]*a[k][i]/ax[k][j];
				}
				if(lambda>0.0f) x[i][j]*=(f+segbiasterm[i][j])/sumColA[i];
				else x[i][j]*=f/sumColA[i];
				if(x[i][j]<min_val) x[i][j]=min_val;
			}
		}
	}

	public boolean init_a(float[][] a,boolean plot,float[] centers,float[] sigmas){
		int nch=a.length;
		int nspec=a[0].length;
		int nplotspec=0;
		if(plot){
			ImageWindow[] iw=jutils.selectPlots(false,1);
			if(iw==null) return false;
			float[][] plotspectra=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
			nplotspec=plotspectra.length;
			for(int i=0;i<plotspectra.length;i++){
				float sum=0.0f;
				for(int j=0;j<nch;j++){
					a[j][i]=plotspectra[i][j];
					sum+=a[j][i];
				}
				for(int j=0;j<nch;j++) a[j][i]/=sum;
			}
		}
		for(int i=nplotspec;i<nspec;i++){
			float center=centers[i];
			float sigma=sigmas[i];
			float sum=0.0f;
			for(int j=0;j<nch;j++){
				a[j][i]=(float)Math.exp(-0.5*(j-center)*(j-center)/(sigma*sigma));
				sum+=a[j][i];
			}
			for(int j=0;j<nch;j++){
				a[j][i]/=sum;
			}
		}
		return true;
	}

	public void init_x(float[][] x){
		for(int i=0;i<x.length;i++){
			for(int j=0;j<x[0].length;j++){
				x[i][j]=(float)(0.5*Math.random()+0.5);
			}
		}
	}

}
