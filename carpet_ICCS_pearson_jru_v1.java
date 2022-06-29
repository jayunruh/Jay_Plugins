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
import jguis.*;

public class carpet_ICCS_pearson_jru_v1 implements PlugIn {

	public void run(String arg) {
		int width=0; int height=0; Object carpet1=null; Object carpet2=null; float psize=0.0f;
		ImagePlus currimp=WindowManager.getCurrentImage();
		if(currimp.getNChannels()>1){
			width=currimp.getWidth();
			height=currimp.getHeight();
			carpet1=currimp.getStack().getPixels(1);
			carpet2=currimp.getStack().getPixels(2);
			psize=(float)jutils.get_psize(currimp);
		} else {
			int[] wList = WindowManager.getIDList();
			String[] titles = new String[wList.length];
			for(int i=0;i<wList.length;i++){
				ImagePlus imp = WindowManager.getImage(wList[i]);
				if(imp!=null){titles[i]=imp.getTitle();}
				else{titles[i]="";}
			}
			GenericDialog gd=new GenericDialog("Options");
			gd.addChoice("Image 1",titles,titles[0]);
			gd.addChoice("Image 2",titles,titles[0]);
			gd.showDialog(); if(gd.wasCanceled()){return;}
			int index1=gd.getNextChoiceIndex();
			int index2=gd.getNextChoiceIndex();
			ImagePlus imp1 = WindowManager.getImage(wList[index1]);
			ImagePlus imp2 = WindowManager.getImage(wList[index2]);
			height = imp1.getHeight();
			width = imp1.getWidth();
			carpet1=imp1.getProcessor().getPixels();
			carpet2=imp2.getProcessor().getPixels();
			psize=(float)jutils.get_psize(imp1);
		}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Temporal_Correlation?",false);
		gd.addNumericField("Time_per_frame(s)",1.0,5,15,null);
		gd.addCheckbox("Symmetrize spatial corr",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean temporal=gd.getNextBoolean();
		float T=(float)gd.getNextNumber();
		boolean symmetrize=gd.getNextBoolean();
		double p2length=Math.log((double)width)/Math.log(2.0);
		boolean p2=false;
		if(p2length==Math.floor(p2length)) p2=true;
		crosscorr acclass=new crosscorr(width,p2);
		float[] ics1=new float[width];
		float[] xvals=new float[width];
		for(int i=0;i<height;i++){
			float[] temp1=new float[width];
			float[] temp2=new float[width];
			for(int j=0;j<width;j++){
				if(carpet1 instanceof short[]){
					temp1[j]=(float)(((short[])carpet1)[i*width+j]&0xffff);
					temp2[j]=(float)(((short[])carpet2)[i*width+j]&0xffff);
				} else {
					temp1[j]=((float[])carpet1)[i*width+j];
					temp2[j]=((float[])carpet2)[i*width+j];
				}
			}
			float[][] temp3=null;
			if(p2) temp3=acclass.docrosscorr(temp1,temp2,false);
			else temp3=acclass.docrosscorrnofft(temp1,temp2,false);
			float stdev1=jstatistics.getstatistic("StDev",temp1,null);
			float stdev2=jstatistics.getstatistic("StDev",temp2,null);
			for(int j=0;j<width;j++){
				int position=j+width/2;
				if(position>=width){position-=width;}
				ics1[j]+=temp3[0][position]*temp3[1][0]*temp3[1][1]/(stdev1*stdev2*height);
				xvals[j]=psize*(float)(j-width/2);
			}
			IJ.showProgress(i,height);
		}
		if(symmetrize){
			float[] xvals3=new float[width/2];
			float[] ics3=new float[width/2];
			xvals3[0]=xvals[width/2];
			ics3[0]=ics1[width/2];
			for(int i=1;i<width/2;i++){
				xvals3[i]=xvals[i+width/2];
				ics3[i]=0.5f*(ics1[i+width/2]+ics1[width/2-i]);
			}
			(new PlotWindow4("ICCS","xi","G(xi)",xvals3,ics3)).draw();
		} else {
			(new PlotWindow4("ICCS","xi","G(xi)",xvals,ics1)).draw();
		}
		if(temporal){
			p2length=Math.log((double)height)/Math.log(2.0);
			p2=false;
			if(p2length==Math.floor(p2length)) p2=true;
			acclass=new crosscorr(height,p2);
			float[] ics2=new float[height/2];
			float[] xvals2=new float[height];
			for(int i=0;i<width;i++){
				float[] temp1=algutils.convert_arr_float2(algutils.get_image_col(carpet1,width,height,i));
				float[] temp2=algutils.convert_arr_float2(algutils.get_image_col(carpet2,width,height,i));
				float[][] temp3=null;
				if(p2) temp3=acclass.docrosscorr(temp1,temp2,false);
				else temp3=acclass.docrosscorrnofft(temp1,temp2,false);
				float stdev1=jstatistics.getstatistic("StDev",temp1,null);
				float stdev2=jstatistics.getstatistic("StDev",temp2,null);
				for(int j=0;j<height/2;j++){
					ics2[j]+=temp3[0][j]*temp3[1][0]*temp3[1][1]/(stdev1*stdev2*width);
					xvals2[j]=T*(float)j;
				}
				IJ.showProgress(i,height);
			}
			(new PlotWindow4("ITCS","tau","G(tau)",xvals2,ics2)).draw();
		}
	}

}
