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

public class stack_real_3D_FFT_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		po4realfft3D fft=new po4realfft3D(width,height,slices);
		Object[] rdata=new Object[slices];
		Object[] idata=new Object[slices];
		//start by copying the data for the fft
		for(int j=0;j<slices;j++){
			float[] indata=(float[])stack.getPixels(j+1);
			float[] temp=new float[width*height];
			System.arraycopy(indata,0,temp,0,width*height);
			rdata[j]=temp;
			idata[j]=new float[width*height];
		}
		//now do the fft
		fft.dorealfft3D(rdata,idata,false);
		ImageStack cstack=new ImageStack(width,height);
		for(int j=0;j<slices;j++){
			cstack.addSlice("",rdata[j]);
			cstack.addSlice("",idata[j]);
		}
		ImagePlus imp2=new ImagePlus("FFT",cstack);
		imp2.show();
	}

	/*public void rotate_to_center(float[] data,int width, int height){
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
	}*/

}
