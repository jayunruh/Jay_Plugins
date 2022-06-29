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

public class fft2D_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean inverse=false;
		gd.addCheckbox("Inverse?",inverse);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		inverse=gd.getNextBoolean();
		ImagePlus oimp = WindowManager.getCurrentImage();
		int oheight = oimp.getHeight();
		int owidth = oimp.getWidth();
		ImageStack ostack = oimp.getStack();
		int osize = ostack.getSize();
		if(inverse){osize/=2;}

		//do the fft of the convolved image
		int[] windex=fftutils.get_best_index(owidth,false,19);
		int[] hindex=fftutils.get_best_index(oheight,false,19);
		po4realfft2D fft=new po4realfft2D(windex[1],hindex[1],windex[0],hindex[0]);
		Object[] real=new Object[osize];
		Object[] im=new Object[osize];
		if(inverse){
			for(int i=0;i<osize;i++){
				real[i]=algutils.get_region2(ostack.getPixels(2*i+1),0,0,windex[1],hindex[1],owidth,oheight);
				im[i]=algutils.get_region2(ostack.getPixels(2*i+2),0,0,windex[1],hindex[1],owidth,oheight);
			}
		}else{
			for(int i=0;i<osize;i++){
				real[i]=algutils.get_region2(ostack.getPixels(i+1),0,0,windex[1],hindex[1],owidth,oheight);
				im[i]=new float[windex[1]*hindex[1]];
			}
		}

		for(int i=0;i<osize;i++) fft.dorealfft2D((float[])real[i],(float[])im[i],inverse);

		ImageStack outstack=new ImageStack(windex[1],hindex[1]);
		if(!inverse){
			for(int i=0;i<osize;i++){
				outstack.addSlice("",real[i]);
				outstack.addSlice("",im[i]);
			}
			(new ImagePlus("FFT",outstack)).show();
		}
		else{
			for(int i=0;i<osize;i++){
				outstack.addSlice("",real[i]);
			}
			(new ImagePlus("iFFT",outstack)).show();
		}
	}

}
