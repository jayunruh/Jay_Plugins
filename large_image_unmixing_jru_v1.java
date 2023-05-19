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
import ij.io.*;
import jguis.*;
import jalgs.*;
import jalgs.jfit.*;
import java.util.concurrent.*;

public class large_image_unmixing_jru_v1 implements PlugIn {
	int species;
	boolean subback;

	public void run(String arg) {
			String[] plotlabels={"Spectra Plot","Background Plot"};
			ImageWindow[] iw=jutils.selectPlots(true,2,plotlabels);
			float[][] tempspectra2=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
			int maxspecies=tempspectra2.length;
			String[] tempnames=new String[maxspecies];
			int[] tempindices=new int[maxspecies];
			for(int i=0;i<maxspecies;i++){tempindices[i]=i; tempnames[i]="series"+i;}
			species=maxspecies;
			int maxlength=tempspectra2[0].length;
			float[] backgroundspectrum=new float[maxlength];
			if(iw[1]!=null){
				backgroundspectrum=((float[][])jutils.runPW4VoidMethod(iw[1],"getYValues"))[0];
			}

			GenericDialog gd2=new GenericDialog("More Options");
			gd2.addNumericField("How_Many_Species_are_Present",species,0);
			gd2.addCheckbox("Subtract_Background?",subback);
			gd2.addCheckbox("Output_chi^2?",false);
			gd2.addCheckbox("Truncate_Negative_Values",true);
			gd2.addNumericField("Start_ch",1,0);
			gd2.addNumericField("End_ch",maxlength,0);
			gd2.addCheckbox("Frame_Batch",false);
			gd2.addNumericField("N_Threads",1,0);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			species=(int)gd2.getNextNumber();
			subback=gd2.getNextBoolean();
			boolean outc2=gd2.getNextBoolean();
			boolean truncate=gd2.getNextBoolean();
			int startch=(int)gd2.getNextNumber()-1;
			int endch=(int)gd2.getNextNumber()-1;
			boolean framebatch=gd2.getNextBoolean();
			int nthreads=(int)gd2.getNextNumber();
			if(species>maxspecies || species==0){IJ.showMessage("Not enough reference spectra"); return;}

			GenericDialog gd3=new GenericDialog("Pick Reference Spectra");
			for(int i=0;i<species;i++){gd3.addChoice("Species_"+(i+1),tempnames,tempnames[tempindices[i]]);}
			gd3.showDialog(); if(gd3.wasCanceled()){return;}
			for(int i=0;i<species;i++){tempindices[i]=gd3.getNextChoiceIndex();}

			OpenDialog od = new OpenDialog("Open Image...", arg);
        			String directory = od.getDirectory();
			String fname=od.getFileName();
			if(fname==null || fname.length()==0){return;}

			LOCI_random_access_file_reader r=new LOCI_random_access_file_reader(directory,fname,0,false);
			if(r==null) return;
			int nch=r.channels;
			int slices=r.slices;
			int frames=r.frames;
			if(nch>maxlength){nch=maxlength;}
			int width=r.width;
			int height=r.height;
			ImageStack resstack=new ImageStack(width,height);
			ImageStack c2stack=null;
			if(outc2) c2stack=new ImageStack(width,height);
			for(int i=0;i<frames*slices;i++){
				for(int j=0;j<species;j++){
					resstack.addSlice(tempnames[tempindices[j]]+i,new float[width*height]);
				}
				if(outc2){
					c2stack.addSlice("",new float[width*height]);
				}
			}
			float[][] tempspectra=new float[species][];
			for(int i=0;i<species;i++){
				tempspectra[i]=tempspectra2[tempindices[i]];
			}
			linleastsquares lls=new linleastsquares(tempspectra,false,startch,endch);
			int totpixels=width*height*slices*frames;
			int counter=0;
			ExecutorService executor=null;
			if(nthreads>1) executor=Executors.newFixedThreadPool(nthreads);
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					float[][] resseries=new float[species][];
					for(int k=0;k<species;k++){
						resseries[k]=((float[])resstack.getPixels(i*slices*species+j*species+k+1));
					}
					float[] c2img=null;
					if(outc2){
						c2img=((float[])c2stack.getPixels(i*slices+j+1));
					}
					Object[] cseries=null;
					if(framebatch){
						cseries=r.getCSeries(j,i,0,0,width,height);
					}
					for(int m=0;m<height;m++){
						float[][] carpet=null;
						if(framebatch){
							carpet=new float[nch][];
							for(int l=0;l<nch;l++) carpet[l]=algutils.convert_arr_float2(algutils.get_subarray(cseries[l],m*width,width));
						} else {
							carpet=r.getCCarpet(m,j,i);
						}
						for(int l=0;l<width;l++){
							float[] spectrum=new float[nch];
							for(int k=0;k<nch;k++) spectrum[k]=carpet[k][l];
							if(subback && backgroundspectrum!=null){
								for(int k=0;k<nch;k++) spectrum[k]-=backgroundspectrum[k];
							}
							if(nthreads>1){
								Runnable worker=new UnmixSpectrum2(lls,spectrum,width,height,resseries,l+m*width,truncate,c2img);
								executor.execute(worker);
							} else {
								double[] contributions=lls.fitdata(spectrum,null);
								if(truncate){for(int k=0;k<species;k++){if(contributions[k]<0.0) contributions[k]=0.0;}}
								for(int k=0;k<species;k++){
									//((float[])resstack.getPixels(i*slices*species+j*species+k+1))[l+m*width]=(float)contributions[k];
									resseries[k][l+m*width]=(float)contributions[k];
								}
								if(outc2){
									//((float[])c2stack.getPixels(i*slices+j+1))[l+m*width]=(float)lls.get_c2(contributions,spectrum,null);
									c2img[l+m*width]=(float)lls.get_c2(contributions,spectrum,null);
								}
							}
							counter++;
							IJ.showProgress(counter,totpixels);
						}
					}
					//IJ.showProgress(j+i*slices,frames*slices);
					IJ.log("slice "+(j+i*slices)+" of "+(frames*slices)+" complete");
				}
			}
			if(nthreads>1){
				executor.shutdown();
				//wait for threads to complete
				IJ.log("waiting for unmixing threads to finish");
				while(!executor.isTerminated()){;}
			}
			ImagePlus imp5=new ImagePlus("Unmixed Stack",resstack);
			imp5.setDimensions(species,slices,frames);
			imp5.setOpenAsHyperStack(true);
			jutils.set_psize(imp5,r.psize);
			jutils.set_pdepth(imp5,r.zsize);
			jutils.set_pinterval(imp5,r.tsize);
			if(species>1){
				new CompositeImage(imp5,CompositeImage.COLOR).show();
			} else {
				imp5.show();
			}
			if(outc2){
				new ImagePlus("chi^2",c2stack).show();
			}
			r.dispose();
	}

	

}

class UnmixSpectrum2 implements Runnable{
	linleastsquares lls;
	float[] spectrum;
	int index,width,height;
	float[][] resseries;
	float[] c2img;
	ImageStack resstack,c2stack;
	boolean truncneg;
	
	public UnmixSpectrum2(linleastsquares lls,float[] spectrum,int width,int height,float[][] resseries,int index,boolean truncneg,float[] c2img){
		this.lls=lls;
		this.spectrum=spectrum;
		this.resseries=resseries;
		this.index=index;
		this.truncneg=truncneg;
		this.c2img=c2img;
		this.width=width;
		this.height=height;
	}

	public void run(){
		double[] contributions=lls.fitdata(spectrum,null);
		if(truncneg){for(int k=0;k<contributions.length;k++){if(contributions[k]<0.0) contributions[k]=0.0;}}
		for(int k=0;k<contributions.length;k++){
			resseries[k][index]=(float)contributions[k];
		}
		if(c2stack!=null){
			c2img[index]=(float)lls.get_c2(contributions,spectrum,null);
		}
	}
}
