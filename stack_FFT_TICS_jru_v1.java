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

public class stack_FFT_TICS_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin calculates the temporal correlation at each pixel using an FFT
		int i,height,width,slices,j,k,p2length,size;

		//get the current image and its info
		ImagePlus imp = WindowManager.getCurrentImage();
		height = imp.getHeight();
		width = imp.getWidth();
		ImageStack stack = imp.getStack();
		slices = stack.getSize();
		ImageProcessor ip=imp.getProcessor();
		Rectangle r = ip.getRoi();

		p2length=(int)(Math.log((double)slices)/Math.log(2.0));
		String[] p2lengths=new String[p2length-1];
		for(i=0;i<p2length-1;i++){
			p2lengths[i]=""+(int)(Math.pow(2.0,p2length-i));
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Analysis length",p2lengths,p2lengths[0]);
		boolean binlog=false;
		gd.addCheckbox("log bin?",binlog);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		p2length-=gd.getNextChoiceIndex();
		binlog=gd.getNextBoolean();
		size = (int)Math.pow(2.0,p2length);
		int segments=(int)((float)slices/(float)size);
		if(size>slices){
			IJ.showMessage("Error","Analysis runs past end of stack");
			return;
		}
			
		//now that we have the data, calculate the autocorrelation
		int counter=0;
		Object[] acmatrix=new Object[size/2];
		for(i=0;i<size/2;i++){
			acmatrix[i]=new float[r.height*r.width];
		}
		autocorr acclass=new autocorr(size);
		if(ip instanceof FloatProcessor){
			for(i=r.y;i<r.y+r.height;i++){
				for(j=r.x;j<r.x+r.width;j++){
					for(int l=0;l<segments;l++){
						float[] temp = new float[size];
						for(k=0;k<size;k++){
							float[] temp2=(float[])stack.getPixels(k+l*size+1);
							temp[k]=temp2[i*width+j];
						}
						temp=acclass.doautocorr(temp);
						for(k=0;k<size/2;k++){
							((float[])acmatrix[k])[counter]+=temp[k]/(float)segments;
						}
					}
					counter++;
					IJ.showProgress(counter,r.width*r.height);
				}
			}
		}
		if(ip instanceof ShortProcessor){
			for(i=r.y;i<r.y+r.height;i++){
				for(j=r.x;j<r.x+r.width;j++){
					for(int l=0;l<segments;l++){
						float[] temp = new float[size];
						for(k=0;k<size;k++){
							short[] temp2=(short[])stack.getPixels(k+l*size+1);
							temp[k]=temp2[i*width+j]&0xffff;
						}
						temp=acclass.doautocorr(temp);
						for(k=0;k<size/2;k++){
							((float[])acmatrix[k])[counter]+=temp[k]/(float)segments;
						}
					}
					counter++;
					IJ.showProgress(counter,r.width*r.height);
				}
			}
		}

		//logarithmically bin the autocorrelation if called for
		int newsize=size/2;
		if(binlog){
			binmultilog bml=new binmultilog();
			float[] xvals=bml.getxvals(size/2);
			newsize=xvals.length;
			for(i=0;i<r.width*r.height;i++){
				float[] temp=new float[size/2];
				for(j=0;j<size/2;j++){
					temp[j]=((float[])acmatrix[j])[i];
				}
				temp=bml.dobinmultilog(temp,size/2);
				for(j=0;j<newsize;j++){
					((float[])acmatrix[j])[i]=temp[j];
				}
			}
		}
		//copy the autocorrelation to the result stack
		ImageStack result_stack = new ImageStack(r.width,r.height);
		for(j=0;j<newsize;j++){
			result_stack.addSlice(null,acmatrix[j]);
		}
		ImagePlus imp3=new ImagePlus("Temporal Correlation",result_stack);
		imp3.show();
	}

}


