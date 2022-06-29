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
import jalgs.jfft.*;

public class stack_real_FFT_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean shiftcenter=true;
		gd.addCheckbox("Shift Center",shiftcenter);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		shiftcenter=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		ImageStack cstack=new ImageStack(width,height);
		po4realfft2D fft=new po4realfft2D(width,height);
		for(int j=0;j<size;j++){
			float[] rpixels=(float[])stack.getPixels(j+1);
			float[] real=new float[width*height];
			System.arraycopy(rpixels,0,real,0,width*height);
			float[] im=new float[width*height];
			fft.dorealfft2D(real,im,false);
			if(shiftcenter){
				rotate_to_center(real,width,height);
				rotate_to_center(im,width,height);
			}
			cstack.addSlice("",(Object)real);
			cstack.addSlice("",(Object)im);
		}
		ImagePlus imp2=new ImagePlus("FFT",cstack);
		imp2.show();
	}

	public void rotate_to_center(float[] data,int width, int height){
		//here we rearrange so that zero point is in the center
		float[] tempdata=new float[width*height];
		int dumvaly,dumvalx;
		for(int i=0;i<height;i++){
			dumvaly=i+height/2;
			if(dumvaly>=height){dumvaly-=height;}
			for(int j=0;j<width;j++){
				dumvalx=j+width/2;
				if(dumvalx>=width){dumvalx-=width;}
				tempdata[j+i*width]=data[dumvalx+width*dumvaly];
			}
		}
		for(int i=0;i<width*height;i++){
			data[i]=tempdata[i];
		}
	}

}
