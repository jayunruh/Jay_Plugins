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

public class carpet_ICCS_jru_v1 implements PlugIn {

	public void run(String arg) {
		int width=0; int height=0; Object carpet1=null; Object carpet2=null;
		int sublength=100000;
		ImagePlus currimp=WindowManager.getCurrentImage();
		if(currimp.getNChannels()>1){
			width=currimp.getWidth();
			height=currimp.getHeight();
			sublength=height;
			carpet1=currimp.getStack().getPixels(1);
			carpet2=currimp.getStack().getPixels(2);
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
			sublength=height;
		}
		if(sublength>height){sublength=height;}
		int sections=(int)((float)height/(float)sublength);
		crosscorr acclass=new crosscorr(width);
		float[] ics1=new float[width*sections];
		for(int k=0;k<sections;k++){
			for(int i=0;i<sublength;i++){
				float[] temp1=new float[width];
				float[] temp2=new float[width];
				for(int j=0;j<width;j++){
					if(carpet1 instanceof short[]){
						temp1[j]=(float)(((short[])carpet1)[(i+k*sublength)*width+j]&0xffff);
						temp2[j]=(float)(((short[])carpet2)[(i+k*sublength)*width+j]&0xffff);
					} else {
						temp1[j]=((float[])carpet1)[(i+k*sublength)*width+j];
						temp2[j]=((float[])carpet2)[(i+k*sublength)*width+j];
					}
				}
				temp1=acclass.docrosscorr(temp1,temp2);
				if(!(Float.isNaN(temp1[0]) || Float.isInfinite(temp1[0]))){
					for(int j=0;j<width;j++){
						int position=j+width/2;
						if(position>=width){position-=width;}
						ics1[j+k*width]+=temp1[position]/(float)sublength;
					}
				}
				IJ.showProgress(k*sublength+i,sections*sublength);
			}
			//if(extrapolate){ics1[width/2+k*width]=(ics1[width/2-1+k*width]+ics1[width/2+1+k*width])/2.0f;}
		}
		if(sections==1){(new PlotWindow4("ICCS","xi","G(xi)",ics1)).draw();}
		else{
			new ImagePlus("Carpet ICCS",new FloatProcessor(width,sections,ics1,null)).show();
		}
	}

}
