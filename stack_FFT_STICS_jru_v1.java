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

public class stack_FFT_STICS_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean avg_quad,use3D;
		ImagePlus imp = WindowManager.getCurrentImage();
		int height = imp.getHeight();
		int width = imp.getWidth();
		ImageStack stack = imp.getStack();
		int slices = stack.getSize();
		GenericDialog gd2=new GenericDialog("STICS Type");
		use3D=false;
		gd2.addCheckbox("Use 3D FFT (uses power of 2 slices)?",use3D);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		use3D=gd2.getNextBoolean();
		ImageStack result_stack=new ImageStack(width,height);
		if(!use3D){
			GenericDialog gd = new GenericDialog("STICS Options");
			avg_quad=false;
			gd.addCheckbox("Average Quadrants?",avg_quad);
			int start=1;
			gd.addNumericField("Beginning Slice",start,0);
			int end=slices;
			gd.addNumericField("End Slice",end,0);
			int shiftstart=0;
			gd.addNumericField("Frame Shift Start",(double)shiftstart,0);
			int shiftend=10;
			gd.addNumericField("Frame Shift End",(double)shiftend,0);
			int shiftspace=1;
			gd.addNumericField("Frame Shift Space",(double)shiftspace,0);
			gd.showDialog();
			if(gd.wasCanceled()){return;}
			avg_quad = gd.getNextBoolean();
			start = (int)gd.getNextNumber();
			end = (int)gd.getNextNumber();
			shiftstart = (int)gd.getNextNumber();
			shiftend = (int)gd.getNextNumber();
			shiftspace = (int)gd.getNextNumber();
			crosscorr2D cc2D=new crosscorr2D(width,height);
			if(imp.getProcessor() instanceof FloatProcessor){
				for(int i=shiftstart;i<(shiftend+1);i+=shiftspace){
					int size = end-start+1-i;
					float[] cc=new float[width*height];
					for(int j=start;j<(end-i+1);j++){
						//get the appropriate stack images
						float[] tempfloat1=(float[])stack.getPixels(j);
						float[] tempfloat2=(float[])stack.getPixels(j+i);
						tempfloat1=cc2D.docrosscorr2D(tempfloat1,tempfloat2,true,true,avg_quad,false);
						for(int k=0;k<width*height;k++){
							cc[k]+=tempfloat1[k]/(float)size;
						}
					}
					//if(avg_quad){avg_quadrants(cc,width,height);}
					result_stack.addSlice("shift "+i,cc);
					IJ.showProgress(i-shiftstart,shiftend-shiftstart+1);
				}
			} else {
				for(int i=shiftstart;i<(shiftend+1);i+=shiftspace){
					int size = end-start+1-i;
					float[] cc=new float[width*height];
					for(int j=start;j<(end-i+1);j++){
						//get the appropriate stack images
						short[] tempfloat1=(short[])stack.getPixels(j);
						short[] tempfloat2=(short[])stack.getPixels(j+i);
						float[] tempfloat=cc2D.docrosscorr2D(tempfloat1,tempfloat2,true,true,avg_quad,false);
						for(int k=0;k<width*height;k++){
							cc[k]+=tempfloat[k]/(float)size;
						}
					}
					//if(avg_quad){avg_quadrants(cc,width,height);}
					result_stack.addSlice("shift "+i,cc);
					IJ.showProgress(i-shiftstart,shiftend-shiftstart+1);
				}
			}
		} else {
			GenericDialog gd=new GenericDialog("STICS options");
			int p2length=(int)(Math.log((double)(slices))/Math.log(2.0));
			String[] p2lengths=new String[p2length-1];
			for(int i=0;i<p2length-1;i++){
				p2lengths[i]=""+(int)(Math.pow(2.0,p2length-i));
			}
			gd.addChoice("Analysis length",p2lengths,p2lengths[0]);
			avg_quad=false;
			gd.addCheckbox("Average Quadrants",avg_quad);
			gd.showDialog(); if(gd.wasCanceled()){return;}
			p2length-=gd.getNextChoiceIndex();
			avg_quad=gd.getNextBoolean();
			int size = (int)Math.pow(2.0,p2length);
			int segments=(int)((float)slices/(float)size);
			autocorr3D ac3D=new autocorr3D(width,height,size);
			for(int i=0;i<size/2;i++){result_stack.addSlice("Shift "+i,new float[width*height]);}
			for(int i=0;i<segments;i++){
				Object[] tempstack=new Object[size];
				for(int j=0;j<size;j++){
					tempstack[j]=stack.getPixels(i*size+j+1);
				}
				Object[] tempstack2=ac3D.doautocorr3D(tempstack,true);
				for(int j=0;j<size/2;j++){
					for(int k=0;k<width*height;k++){
						((float[])result_stack.getPixels(j+1))[k]+=((float[])tempstack2[j])[k]/(float)segments;
					}
				}
				IJ.showProgress(i,segments);
			}
			if(avg_quad){
				for(int j=0;j<size/2;j++){
					avg_quadrants((float[])result_stack.getPixels(j+1),width,height);
				}
			}
		}
		//output the STICS file
		ImagePlus imp4 = new ImagePlus("STICS",result_stack);
		imp4.show();
	}

	private void avg_quadrants(float[] pixels,int width,int height){
		for(int i=1;i<height/2;i++){
			float avg=(pixels[i*width]+pixels[(height-i)*width])/2.0f;
		}
		for(int i=1;i<width/2;i++){
			float avg=(pixels[i]+pixels[width-i])/2.0f;
		}
		for(int i=1;i<height/2;i++){
			for(int j=1;j<width/2;j++){
				float avg=(pixels[i*width+j]+pixels[i*width+width-j]+pixels[(height-i)*width+j]+pixels[(height-i)*width+width-j])/4.0f;
				pixels[i*width+j]=avg; pixels[i*width+width-j]=avg; pixels[(height-i)*width+j]=avg; pixels[(height-i)*width+width-j]=avg;
			}
		}
	}

}
