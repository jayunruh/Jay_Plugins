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
import jguis.*;

public class carpet_ICS_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageProcessor ip=imp.getProcessor();
		int width=imp.getWidth();
		int height=imp.getHeight();

		GenericDialog gd=new GenericDialog("Options");
		boolean extrapolate=true;
		gd.addCheckbox("Extrapolate G(0)?",extrapolate);
		int sublength=height;
		gd.addNumericField("Sublength?",sublength,0,10,null);
		boolean center=true;
		gd.addCheckbox("Center Correlation?",center);
		boolean dobrightcorr=false;
		gd.addCheckbox("BrightCorr?",dobrightcorr);
		boolean avgall=false;
		gd.addCheckbox("Use Global Avg?",avgall);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		extrapolate=gd.getNextBoolean();
		sublength=(int)gd.getNextNumber();
		center=gd.getNextBoolean();
		dobrightcorr=gd.getNextBoolean();
		avgall=gd.getNextBoolean();
		if(sublength>height){sublength=height;}
		int sections=(int)((float)height/(float)sublength);
		autocorr acclass=new autocorr(width);
		float[] ics1=new float[width*sections];
		if(ip instanceof ShortProcessor){
			short[] carpet=(short[])ip.getPixels();
			for(int k=0;k<sections;k++){
				float globalavg=0.0f;
				for(int i=0;i<sublength;i++){
					float[] temp2=new float[width];
					for(int j=0;j<width;j++){
						temp2[j]=(float)(carpet[(i+k*sublength)*width+j]&0xffff);
					}
					float[][] temp1=acclass.doautocorr(temp2,dobrightcorr);
					if(avgall){
						float tempmult=temp1[1][0];
						globalavg+=tempmult/(float)sublength;
						for(int j=0;j<width;j++){
							if(dobrightcorr){temp1[0][j]/=tempmult;}
							temp1[0][j]+=1.0f;
							temp1[0][j]*=tempmult*tempmult;
						}
					}
					if(!(Float.isNaN(temp1[0][0]) || Float.isInfinite(temp1[0][0]))){
						for(int j=0;j<width;j++){
							if(center){
								int position=j+width/2;
								if(position>=width){position-=width;}
								ics1[j+k*width]+=temp1[0][position]/(float)sublength;
							} else {
								ics1[j+k*width]+=temp1[0][j]/(float)sublength;
							}
						}
					}
					IJ.showProgress(k*sublength+i,sections*sublength);
				}
				if(extrapolate){
					if(center){
						ics1[width/2+k*width]=(ics1[width/2-1+k*width]+ics1[width/2+1+k*width])/2.0f;
					} else {
						ics1[k*width]=(ics1[1+k*width]+ics1[(width-1)+k*width])/2.0f;
					}
				}
				if(avgall){
					for(int j=0;j<width;j++){
						ics1[j+k*width]/=(globalavg*globalavg);
						ics1[j+k*width]-=1.0f;
						if(dobrightcorr){ics1[j+k*width]*=globalavg;}
					}
				}
			}
		} else {
			float[] carpet=(float[])ip.getPixels();
			for(int k=0;k<sections;k++){
				float globalavg=0.0f;
				for(int i=0;i<sublength;i++){
					float[] temp2=new float[width];
					for(int j=0;j<width;j++){
						temp2[j]=carpet[(i+k*sublength)*width+j];
					}
					float[][] temp1=acclass.doautocorr(temp2,dobrightcorr);
					if(avgall){
						float tempmult=temp1[1][0];
						globalavg+=tempmult/(float)sublength;
						for(int j=0;j<width;j++){
							if(dobrightcorr){temp1[0][j]/=tempmult;}
							temp1[0][j]+=1.0f;
							temp1[0][j]*=tempmult*tempmult;
						}
					}
					if(!(Float.isNaN(temp1[0][0]) || Float.isInfinite(temp1[0][0]))){
						for(int j=0;j<width;j++){
							if(center){
								int position=j+width/2;
								if(position>=width){position-=width;}
								ics1[j+k*width]+=temp1[0][position]/(float)sublength;
							} else {
								ics1[j+k*width]+=temp1[0][j]/(float)sublength;
							}
						}
					}
					IJ.showProgress(k*sublength+i,sections*sublength);
				}
				if(extrapolate){
					if(center){
						ics1[width/2+k*width]=(ics1[width/2-1+k*width]+ics1[width/2+1+k*width])/2.0f;
					} else {
						ics1[k*width]=(ics1[1+k*width]+ics1[(width-1)+k*width])/2.0f;
					}
				}
				if(avgall){
					for(int j=0;j<width;j++){
						ics1[j+k*width]/=(globalavg*globalavg);
						ics1[j+k*width]-=1.0f;
						if(dobrightcorr){ics1[j+k*width]*=globalavg;}
					}
				}
			}
		}
		float[] xvals=new float[width];
		for(int i=0;i<width;i++){
			xvals[i]=(float)(i-width/2);
		}
		if(sections==1){(new PlotWindow4("ICS","xi","G(xi)",xvals,ics1)).draw();}
		else{
			new ImagePlus("Carpet ICS",new FloatProcessor(width,sections,ics1,null)).show();
		}
	}

}
