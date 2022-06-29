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
import java.io.*;
import jguis.*;
import jalgs.*;
import jalgs.jfit.*;
import java.util.concurrent.*;

public class linear_unmixing_jru_v2 implements PlugIn {
	int species,mainoptionsindex,totspecies,maxlength;
	boolean subback;
	String[] speciesnames;
	float[][] speciesspectra;
	float[] backgroundspectrum;
	int[] speciesindices;
	int nthreads=1;

	public void run(String arg) {
		//this plugin performs linear unmixing and in so doing finds the intensity of multiple components in frequency resolved image
		//the image could be frequency resolved in color or in temporal fourier frequency (lifetime or photoactivation)
		//this version doesn't use the database and allows for selection of the channel region for unmixing
		//include the background as a spectrum if you would like
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Spectral_Image","Spectra_Plot"});
		if(imps==null) return;
		GenericDialog gd2=new GenericDialog("More Options");
		gd2.addCheckbox("Normalize_Spectra",true);
		gd2.addCheckbox("Output_Residuals?",false);
		gd2.addCheckbox("Output_chi^2?",false);
		gd2.addCheckbox("Output_errors (plots only)",false);
		gd2.addCheckbox("Truncate_Negative_Values",true);
		gd2.addNumericField("Start_ch",1,0);
		int tempnch=imps[0].getNChannels();
		if(jutils.isPlot(imps[0].getWindow())){
			tempnch=((int[])jutils.runPW4VoidMethod(imps[0].getWindow(),"getNpts"))[0];
		}
		boolean carpet=false;
		if(tempnch==1){
			tempnch=imps[0].getWidth();
			carpet=true;
		}
		gd2.addNumericField("End_ch",tempnch,0);
		gd2.addNumericField("N_Threads",1,0);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		boolean normspec=gd2.getNextBoolean();
		boolean outres=gd2.getNextBoolean();
		boolean outc2=gd2.getNextBoolean();
		boolean outerrs=gd2.getNextBoolean();
		boolean truncate=gd2.getNextBoolean();
		int startch=(int)gd2.getNextNumber()-1;
		int endch=(int)gd2.getNextNumber()-1;
		nthreads=(int)gd2.getNextNumber();
		if(jutils.isPlot(imps[0].getWindow())){
			Object[] output=exec(imps[0].getWindow(),imps[1].getWindow(),startch,endch,outres,outc2,truncate,outerrs,normspec);
			//if we only have one series, output the data, fit, and all scaled contributions
			int nser=(Integer)jutils.runPW4VoidMethod(imps[0].getWindow(),"getNSeries");
			if(nser==1){
				float[][] data=(float[][])jutils.runPW4VoidMethod(imps[0].getWindow(),"getYValues");
				float[][] spectra=(float[][])jutils.runPW4VoidMethod(imps[1].getWindow(),"getYValues");
				float[][] datafit=new float[2+spectra.length][spectra[0].length];
				datafit[0]=data[0];
				for(int i=0;i<data[0].length;i++){
					for(int j=0;j<spectra.length;j++){
						float temp=((float[][])output[0])[0][j]*spectra[j][i];
						datafit[j+2][i]=temp;
						datafit[1][i]+=temp;
					}
				}
				(new PlotWindow4("Fit","channel","intensity",datafit,null)).draw();
			}
			//output the contributions to a plot window
			PlotWindow4 pw=new PlotWindow4("Concentrations","dye","concentration",(float[][])output[0],null);
			pw.draw();
			int outcounter=1;
			if(outres){
				(new PlotWindow4("Residuals","channel","residual",(float[][])output[outcounter],null)).draw();
				outcounter++;
			}
			if(outc2){
				//just dump the chi squared values to the log window
				float[] c2vals=(float[])output[outcounter];
				for(int i=0;i<c2vals.length;i++){
					IJ.log("c2_"+(i+1)+" = "+c2vals[i]);
				}
				outcounter++;
			}
			if(outerrs){
				//add the error bars to the concentrations plot
				float[][] errstack=(float[][])output[outcounter];
				pw.addErrors(errstack);
			}
		} else if(carpet) {
			Object[] output=exec(imps[0],imps[1].getWindow(),startch,endch,outres,outc2,truncate,outerrs,normspec);
			PlotWindow4 pw=new PlotWindow4("Concentrations","slice","concentration",transpose((float[][])output[0]),null);
			pw.draw();
			int outcounter=1;
			if(outres){
				(new PlotWindow4("Residuals","channel","residual",transpose((float[][])output[outcounter]),null)).draw();
				outcounter++;
			}
			if(outc2){
				//just dump the chi squared values to the log window
				float[] c2vals=(float[])output[outcounter];
				for(int i=0;i<c2vals.length;i++){
					IJ.log("c2_"+(i+1)+" = "+c2vals[i]);
				}
				outcounter++;
			}
			if(outerrs){
				//add the error bars to the concentrations plot
				float[][] errstack=transpose((float[][])output[outcounter]);
				pw.addErrors(errstack);
			}
		} else {
			Object[] outimages=exec(imps[0],imps[1].getWindow(),startch,endch,outres,outc2,truncate,normspec);
			for(int i=0;i<outimages.length;i++) ((ImagePlus)outimages[i]).show();
		}
	}

	public float[][] transpose(float[][] src){
		float[][] out=new float[src[0].length][src.length];
		for(int i=0;i<src[0].length;i++){
			for(int j=0;j<src.length;j++){
				out[i][j]=src[j][i];
			}
		}
		return out;
	}


	public Object[] exec(ImageWindow input,ImageWindow plot,int startch,int endch,boolean outresid,boolean outc2,boolean truncneg,boolean outerrs,boolean normspec){
		//here we unmix a set of plots
		float[][] spectra1=(float[][])jutils.runPW4VoidMethod(plot,"getYValues");
		float[][] spectra=spectra1;
		if(normspec) spectra=normSpectra(spectra1);
		float[][] data=(float[][])jutils.runPW4VoidMethod(input,"getYValues");
		int slices=data.length;
		int channels=data[0].length;
		if(startch<0) startch=0;
		if(endch>=channels) endch=channels-1;
		int fitch=endch-startch+1;
		//start by creating the linear least squares object
		linleastsquares lls=new linleastsquares(spectra,false,startch,endch);
		//now go through pixel by pixel and do the unmixing
		float[][] outstack=new float[slices][spectra.length];
		float[][] residstack=null; float[] c2stack=null; float[][] errstack=null;
		if(outresid) residstack=new float[slices][channels];
		if(outc2) c2stack=new float[slices];
		if(outerrs) errstack=new float[slices][spectra.length];
		for(int j=0;j<slices;j++){
			double[] contributions=null;
			if(outerrs){
				double[][] temp=lls.getfiterrors(data[j],null);
				contributions=temp[0];
				errstack[j]=algutils.convert_arr_float(temp[1]);
			} else {
				contributions=lls.fitdata(data[j],null);
			}
			for(int l=0;l<contributions.length;l++){
				if(truncneg && contributions[l]<0.0f) contributions[l]=0.0f;
				outstack[j][l]=(float)contributions[l];
			}
			if(outresid){
				float[] residcol=lls.get_fresid(contributions,data[j],null);
				for(int l=0;l<channels;l++){residstack[j][l]=residcol[l];}
			}
			if(outc2){
				c2stack[j]=(float)lls.get_c2(contributions,data[j],null);
			}
			IJ.showProgress(j,slices);
		}
		int nout=1;
		if(outresid) nout++;
		if(outc2) nout++;
		if(outerrs) nout++;
		Object[] output=new Object[nout];
		output[0]=outstack;
		int outcounter=1;
		if(outresid){output[outcounter]=residstack; outcounter++;}
		if(outc2){output[outcounter]=c2stack; outcounter++;}
		if(outerrs){output[outcounter]=errstack; outcounter++;}
		return output;
	}

	public static float[][] normSpectra(float[][] spectra){
		float[][] normed=new float[spectra.length][spectra[0].length];
		float[] sums=new float[spectra.length];
		for(int i=0;i<spectra.length;i++){
			sums[i]=jstatistics.getstatistic("Sum",spectra[i],null);
			for(int j=0;j<spectra[i].length;j++){
				normed[i][j]=spectra[i][j]/sums[i];
			}
		}
		return normed;
	}

	public Object[] exec(ImagePlus input,ImageWindow plot,int startch,int endch,boolean outresid,boolean outc2,boolean truncneg,boolean outerrs,boolean normspec){
		//here we unmix a carpet
		float[][] spectra1=(float[][])jutils.runPW4VoidMethod(plot,"getYValues");
		float[][] spectra=spectra1;
		if(normspec) spectra=normSpectra(spectra1);
		//float[][] data=(float[][])jutils.runPW4VoidMethod(input,"getYValues");
		int slices=input.getHeight();
		int channels=input.getWidth();
		float[][] data=new float[slices][channels];
		float[] pixdata=(float[])input.getProcessor().convertToFloat().getPixels();
		for(int i=0;i<slices;i++){
			System.arraycopy(pixdata,i*channels,data[i],0,channels);
		}
		//int slices=data.length;
		//int channels=data[0].length;
		if(startch<0) startch=0;
		if(endch>=channels) endch=channels-1;
		int fitch=endch-startch+1;
		//start by creating the linear least squares object
		linleastsquares lls=new linleastsquares(spectra,false,startch,endch);
		//now go through pixel by pixel and do the unmixing
		float[][] outstack=new float[slices][spectra.length];
		float[][] residstack=null; float[] c2stack=null; float[][] errstack=null;
		if(outresid) residstack=new float[slices][channels];
		if(outc2) c2stack=new float[slices];
		if(outerrs) errstack=new float[slices][spectra.length];
		for(int j=0;j<slices;j++){
			double[] contributions=null;
			if(outerrs){
				double[][] temp=lls.getfiterrors(data[j],null);
				contributions=temp[0];
				errstack[j]=algutils.convert_arr_float(temp[1]);
			} else {
				contributions=lls.fitdata(data[j],null);
			}
			for(int l=0;l<contributions.length;l++){
				if(truncneg && contributions[l]<0.0f) contributions[l]=0.0f;
				outstack[j][l]=(float)contributions[l];
			}
			if(outresid){
				float[] residcol=lls.get_fresid(contributions,data[j],null);
				for(int l=0;l<channels;l++){residstack[j][l]=residcol[l];}
			}
			if(outc2){
				c2stack[j]=(float)lls.get_c2(contributions,data[j],null);
			}
			IJ.showProgress(j,slices);
		}
		int nout=1;
		if(outresid) nout++;
		if(outc2) nout++;
		if(outerrs) nout++;
		Object[] output=new Object[nout];
		output[0]=outstack;
		int outcounter=1;
		if(outresid){output[outcounter]=residstack; outcounter++;}
		if(outc2){output[outcounter]=c2stack; outcounter++;}
		if(outerrs){output[outcounter]=errstack; outcounter++;}
		return output;
	}

	public Object[] exec(ImagePlus input,ImageWindow plot,int startch,int endch,boolean outresid,boolean outc2,boolean truncneg,boolean normspec){
		//here we unmix a hyperstack
		float[][] spectra1=(float[][])jutils.runPW4VoidMethod(plot,"getYValues");
		float[][] spectra=spectra1;
		if(normspec) spectra=normSpectra(spectra1);
		int frames=input.getNFrames();
		int slices=input.getNSlices();
		int channels=input.getNChannels();
		int width=input.getWidth(); int height=input.getHeight();
		if(startch<0) startch=0;
		if(endch>=channels) endch=channels-1;
		int fitch=endch-startch+1;
		ImageStack stack=input.getStack();
		//start by creating the linear least squares object
		linleastsquares lls=new linleastsquares(spectra,false,startch,endch);
		//now go through pixel by pixel and do the unmixing
		ImageStack outstack=new ImageStack(width,height);
		ImageStack residstack=new ImageStack(width,height);
		ImageStack c2stack=new ImageStack(width,height);
		int counter=0;
		ExecutorService executor=null;
		if(nthreads>1) executor=Executors.newFixedThreadPool(nthreads);
		for(int i=0;i<frames;i++){
			counter=0;
			for(int j=0;j<slices;j++){
				Object[] cseries=jutils.get3DCSeries(stack,j,i,frames,slices,channels);
				float[][] contr=new float[spectra.length][width*height];
				float[][] resid=null;
				if(outresid) resid=new float[channels][width*height];
				float[] c2=null;
				if(outc2) c2=new float[width*height];
				for(int n=0;n<height;n++){
					for(int m=0;m<width;m++){
						int k=m+n*width;
						if(nthreads>1){
							Runnable worker=new UnmixSpectrum(lls,cseries,width,height,contr,k,truncneg,resid,c2);
							//worker.run();
							executor.execute(worker);
						} else {
						float[] col=algutils.convert_arr_float(algutils.get_stack_col(cseries,width,height,k,cseries.length));
						double[] contributions=lls.fitdata(col,null);
						for(int l=0;l<contributions.length;l++){
							if(truncneg && contributions[l]<0.0f) contributions[l]=0.0f;
							contr[l][k]=(float)contributions[l];
						}
						if(outresid){
							float[] residcol=lls.get_fresid(contributions,col,null);
							for(int l=0;l<channels;l++){resid[l][k]=residcol[l];}
						}
						if(outc2){
							c2[k]=(float)lls.get_c2(contributions,col,null);
						}
						}
					}
					/*if(nthreads>1){
						executor.shutdown();
						//wait for threads to complete
						while(!executor.isTerminated()){;}
						//System.gc();
						executor=Executors.newFixedThreadPool(nthreads);
					}*/
				}
				/*if(nthreads>1){
					executor.shutdown();
					//wait for threads to complete
					while(!executor.isTerminated()){;}
					//System.gc();
					executor=Executors.newFixedThreadPool(nthreads);
				}*/
				for(int k=0;k<spectra.length;k++) outstack.addSlice("",contr[k]);
				if(outresid) for(int k=0;k<channels;k++) residstack.addSlice("",resid[k]);
				if(outc2) c2stack.addSlice("",c2);
				counter++;
				IJ.showProgress(counter,slices);
				IJ.log("frame "+counter+" out of "+slices+" initiated");
			}
		}
		if(nthreads>1){
			executor.shutdown();
			//wait for threads to complete
			IJ.log("waiting for unmixing threads to finish");
			while(!executor.isTerminated()){;}
		}
		IJ.log("showing image");
		ImagePlus outimp=new ImagePlus("Unmixed Stack",outstack);
		outimp.copyScale(input);
		outimp.setOpenAsHyperStack(true);
		outimp.setDimensions(spectra.length,slices,frames);
		int nout=1;
		if(outresid) nout++;
		if(outc2) nout++;
		Object[] output=new Object[nout];
		int outcounter=0;
		output[0]=new CompositeImage(outimp,CompositeImage.COLOR); outcounter++;
		if(outresid){
			ImagePlus residimp=new ImagePlus("Residuals Stack",residstack);
			residimp.copyScale(input);
			residimp.setOpenAsHyperStack(true);
			residimp.setDimensions(channels,slices,frames);
			output[outcounter]=new CompositeImage(residimp,CompositeImage.GRAYSCALE);
			outcounter++;
		}
		if(outc2){
			ImagePlus c2imp=new ImagePlus("chi^2",c2stack);
			c2imp.copyScale(input);
			c2imp.setOpenAsHyperStack(true);
			c2imp.setDimensions(1,slices,frames);
			output[outcounter]=c2imp;
		}
		return output;
	}
	

}

class UnmixSpectrum implements Runnable{
	linleastsquares lls;
	Object[] cseries;
	int index,width,height;
	float[][] contr;
	float[][] resid;
	float[] c2;
	boolean truncneg;
	
	public UnmixSpectrum(linleastsquares lls,Object[] cseries,int width,int height,float[][] contr,int index,boolean truncneg,float[][] resid,float[] c2){
		this.lls=lls;
		this.cseries=cseries;
		this.contr=contr;
		this.index=index;
		this.truncneg=truncneg;
		this.resid=resid;
		this.c2=c2;
		this.width=width;
		this.height=height;
	}

	public void run(){
		float[] col=algutils.convert_arr_float(algutils.get_stack_col(cseries,width,height,index,cseries.length));
		double[] contributions=lls.fitdata(col,null);
		for(int l=0;l<contributions.length;l++){
			if(truncneg && contributions[l]<0.0f) contributions[l]=0.0f;
			contr[l][index]=(float)contributions[l];
		}
		if(resid!=null){
			float[] residcol=lls.get_fresid(contributions,col,null);
			for(int l=0;l<residcol.length;l++){resid[l][index]=residcol[l];}
		}
		if(c2!=null){
			c2[index]=(float)lls.get_c2(contributions,col,null);
		}
	}
}
