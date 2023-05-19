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

public class stack_temporal_FFT_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin calculates the temporal FFT at each pixel
		//assume that real and imaginary are interlaced for inverse fft
		GenericDialog gd=new GenericDialog("Options");
		boolean inverse=false;
		gd.addCheckbox("Inverse?",inverse);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		inverse=gd.getNextBoolean();

		//get the current image and its info
		ImagePlus imp = WindowManager.getCurrentImage();
		int height = imp.getHeight();
		int width = imp.getWidth();
		ImageStack stack = imp.getStack();
		int slices = stack.getSize();
		if(inverse){slices/=2;}
		ImageProcessor ip=imp.getProcessor();
		Rectangle r = ip.getRoi();
		int[] index=fftutils.get_best_index(slices,false,19);
		po4realfft fft=new po4realfft(index[1],index[0]);
			
		//now that we have the data, calculate the autocorrelation
		if(!inverse){
			int counter=0;
			float[][] rfftmatrix=new float[index[1]][r.width*r.height];
			float[][] ifftmatrix=new float[index[1]][r.width*r.height];
			Object[] original=jutils.stack2array(stack);
			for(int i=r.y;i<(r.y+r.height);i++){
				for(int j=r.x;j<(r.x+r.width);j++){
					float[] temp=algutils.convert_arr_float(algutils.get_stack_col(original,width,height,j,i,index[1]));
					float[] temp4=new float[slices];
					fft.realfft(temp,temp4,false);
					for(int k=0;k<index[1];k++){
						rfftmatrix[k][counter]=temp[k];
						ifftmatrix[k][counter]=temp4[k];
					}
					counter++;
					IJ.showProgress(counter,r.width*r.height);
				}
			}
			//copy the fft to the result stack
			ImageStack result_stack = new ImageStack(r.width,r.height);
			for(int j=0;j<index[1];j++){
				result_stack.addSlice(null,rfftmatrix[j]);
				result_stack.addSlice(null,ifftmatrix[j]);
			}
			(new ImagePlus("Forward FFT",result_stack)).show();
		}
		else {
			int counter=0;
			float[][] rfftmatrix=new float[index[1]][r.width*r.height];
			Object[] original=jutils.stack2array(stack);
			for(int i=r.y;i<(r.y+r.height);i++){
				for(int j=r.x;j<(r.x+r.width);j++){
					float[] temp=algutils.convert_arr_float(algutils.get_stack_col(original,width,height,j,i,2*index[1]));
					float[] temp2 = new float[slices];
					float[] temp4=new float[slices];
					for(int k=0;k<index[1];k++){
						temp2[k]=temp[2*k]; temp4[k]=temp[2*k+1];
					}
					fft.realfft(temp2,temp4,true);
					for(int k=0;k<index[1];k++){
						rfftmatrix[k][counter]=temp2[k];
					}
					counter++;
					IJ.showProgress(counter,r.width*r.height);
				}
			}
			//copy the fft to the result stack
			ImageStack result_stack = new ImageStack(r.width,r.height);
			for(int j=0;j<index[1];j++){
				result_stack.addSlice(null,rfftmatrix[j]);
			}
			(new ImagePlus("Inverse FFT",result_stack)).show();
		}
	}

}


