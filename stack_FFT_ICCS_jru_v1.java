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
import jalgs.*;
import jalgs.jfft.*;
import ij.measure.*;

public class stack_FFT_ICCS_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean avg_sections,interpolate,isfloat;
		int i,height,width,start,end,size;
		float avg,variance,overall_avg,g0;
		ImageStack stack1,stack2;
		Calibration cal;
		ImagePlus currimp=WindowManager.getCurrentImage();
		if(currimp.getNChannels()>1){
			cal=currimp.getCalibration().copy();
			width=currimp.getWidth();
			height=currimp.getHeight();
			stack1=new ImageStack(width,height);
			stack2=new ImageStack(width,height);
			ImageStack tempstack=currimp.getStack();
			int tempsize=tempstack.getSize();
			size=tempsize/2;
			isfloat=(currimp.getProcessor() instanceof FloatProcessor);
			for(i=0;i<tempsize/2;i++){
				stack1.addSlice("",tempstack.getPixels(2*i+1));
				stack2.addSlice("",tempstack.getPixels(2*i+2));
			}
		} else {			
			int[] wList = WindowManager.getIDList();
			String[] titles = new String[wList.length];
			for(i=0;i<wList.length;i++){
				ImagePlus imp = WindowManager.getImage(wList[i]);
				if(imp!=null){titles[i]=imp.getTitle();}
				else{titles[i]="";}
			}

			GenericDialog gd2=new GenericDialog("ICCS Options");
			gd2.addChoice("Image 1",titles,titles[0]);
			gd2.addChoice("Image 2",titles,titles[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int index1=gd2.getNextChoiceIndex();
			int index2=gd2.getNextChoiceIndex();
			ImagePlus imp1 = WindowManager.getImage(wList[index1]);
			ImagePlus imp2 = WindowManager.getImage(wList[index2]);
			isfloat=(imp1.getProcessor() instanceof FloatProcessor);
			cal=currimp.getCalibration().copy();
			height = imp1.getHeight();
			width = imp1.getWidth();
			stack1 = imp1.getStack();
			size = stack1.getSize();
			stack2=imp2.getStack();
		}
		boolean avg_quad;
		GenericDialog gd = new GenericDialog("Options");
		gd.addCheckbox("Avg Quadrants?",true);
		gd.addCheckbox("Interpolate G(0)?",true);
		gd.addNumericField("Beginning Slice",1.0,0);
		gd.addNumericField("End Slice",(double)size,0);
		boolean brightcorr=false;
		gd.addCheckbox("Brightcorr?",brightcorr);
		boolean pearson=false;
		gd.addCheckbox("Pearson",pearson);
		gd.addCheckbox("Calc_Autocorr",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		avg_sections = gd.getNextBoolean();
		interpolate = gd.getNextBoolean();
		start = (int)gd.getNextNumber();
		if(start<1){start=1;}
		end = (int)gd.getNextNumber();
		if(end>size){end=size;}
		size = end-start+1;
		brightcorr=gd.getNextBoolean();
		pearson=gd.getNextBoolean();
		boolean calcauto=gd.getNextBoolean();
		crosscorr2D cc2D=new crosscorr2D(width,height);
		autocorr2D ac2D=new autocorr2D(width,height);
		float[] ac=new float[width*height];
		float[] ac1=new float[width*height];
		float[] ac2=new float[width*height];
		if(isfloat){
			for(i=start;i<=end;i++){
				//get the appropriate stack image
				float[] pixels1=(float[])stack1.getPixels(i);
				float[] pixels2=(float[])stack2.getPixels(i);
				float[] tempac=cc2D.docrosscorr2D(pixels1,pixels2,true,true,avg_sections,brightcorr);
				float[] tempac1=null;
				float[] tempac2=null;
				if(calcauto){
					tempac1=ac2D.doautocorr2D(pixels1,true,true,avg_sections,brightcorr);
					tempac2=ac2D.doautocorr2D(pixels2,true,true,avg_sections,brightcorr);
				}
				if(pearson){
					float stdev1=jstatistics.getstatistic("StDev",pixels1,null);
					float stdev2=jstatistics.getstatistic("StDev",pixels2,null);
					float avg1=jstatistics.getstatistic("Avg",pixels1,null);
					float avg2=jstatistics.getstatistic("Avg",pixels2,null);
					for(int j=0;j<width*height;j++) tempac[j]*=(avg1*avg2/(stdev1*stdev2));
				}
				for(int j=0;j<width*height;j++){
					ac[j]+=tempac[j]/(float)size;
					if(calcauto){
						ac1[j]+=tempac1[j]/(float)size;
						ac2[j]+=tempac2[j]/(float)size;
					}
				}
				IJ.showProgress(i-start,size);
			}
		} else {
			for(i=start;i<=end;i++){
				//get the appropriate stack image
				short[] pixels1=(short[])stack1.getPixels(i);
				short[] pixels2=(short[])stack2.getPixels(i);
				float[] tempac=cc2D.docrosscorr2D(pixels1,pixels2,true,true,avg_sections,brightcorr);
				float[] tempac1=null;
				float[] tempac2=null;
				if(calcauto){
					tempac1=ac2D.doautocorr2D(pixels1,true,true,avg_sections,brightcorr);
					tempac2=ac2D.doautocorr2D(pixels2,true,true,avg_sections,brightcorr);
				}
				if(pearson){
					float stdev1=jstatistics.getstatistic("StDev",pixels1,null);
					float stdev2=jstatistics.getstatistic("StDev",pixels2,null);
					float avg1=jstatistics.getstatistic("Avg",pixels1,null);
					float avg2=jstatistics.getstatistic("Avg",pixels2,null);
					for(int j=0;j<width*height;j++) tempac[j]*=(avg1*avg2/(stdev1*stdev2));
				}
				for(int j=0;j<width*height;j++){
					ac[j]+=tempac[j]/(float)size;
					if(calcauto){
						ac1[j]+=tempac1[j]/(float)size;
						ac2[j]+=tempac2[j]/(float)size;
					}
				}
				IJ.showProgress(i-start,size);
			}
		}
		if(interpolate){
			ac[(height/2)*width+width/2]=(ac[(height/2)*width+width/2-1]+ac[(height/2)*width+width/2+1])/2.0f;
			if(calcauto){
				ac1[(height/2)*width+width/2]=(ac1[(height/2)*width+width/2-1]+ac1[(height/2)*width+width/2+1])/2.0f;
				ac2[(height/2)*width+width/2]=(ac2[(height/2)*width+width/2-1]+ac2[(height/2)*width+width/2+1])/2.0f;
			}
		}
		//output the autocorrelated files
		if(calcauto){
			ImageStack corrstack=new ImageStack(width,height);
			corrstack.addSlice("Autocorr1",ac1);
			corrstack.addSlice("Autocorr2",ac2);
			corrstack.addSlice("Crosscorr",ac);
			ImagePlus imp10=new ImagePlus("Cross-correlation",corrstack);
			imp10.setOpenAsHyperStack(true);
			imp10.setDimensions(3,1,1);
			imp10.setCalibration(cal);
			new CompositeImage(imp10,CompositeImage.COMPOSITE).show();
		} else {
			ImagePlus imp10=new ImagePlus("Cross-correlation",new FloatProcessor(width,height,ac,null));
			imp10.setCalibration(cal);
			imp10.show();
		}
	}

}
