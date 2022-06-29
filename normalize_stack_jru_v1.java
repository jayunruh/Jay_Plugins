/*******************************************************************************
 * Copyright (c) 2019 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import jalgs.*;

public class normalize_stack_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we normalize a stack
		//normalize every 3D volume to its maximum (for each channel and frame)
		//the stack will be duplicated and converted to 32 bit
		GenericDialog gd=new GenericDialog("Options");
		String[] normoptions={"Sum","Max"};
		gd.addChoice("Normalization",normoptions,normoptions[0]);
		double multfactor=1.0;
		gd.addNumericField("Multiplication Factor",multfactor,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int normindex=gd.getNextChoiceIndex();
		multfactor=gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		Object[] arrstack=jutils.stack2array(stack);
		int nslices=imp.getNSlices();
		int nchans=imp.getNChannels();
		int nframes=imp.getNFrames();
		float[][] normstat=new float[nchans][nframes];
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nchans;j++){
				Object[] zstack=jutils.get3DZSeries(stack,j,i,nframes,nslices,nchans);
				if(normindex==1) normstat[j][i]=getStackMax(zstack);
				else normstat[j][i]=getStackSum(zstack);
			}
		}
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
					float[] normdata=algutils.convert_arr_float(arrstack[k+j*nchans+i*nslices*nchans]);
					for(int l=0;l<normdata.length;l++){
						normdata[l]*=(multfactor/normstat[k][i]);
					}
					retstack.addSlice("",normdata);
				}
			}
		}
		jutils.create_hyperstack("Normalized Image",retstack,imp).show();
	}

	public float getStackMax(Object[] stack){
		float max=jstatistics.getstatistic("Max",stack[0],null);
		for(int i=1;i<stack.length;i++){
			float temp=jstatistics.getstatistic("Max",stack[i],null);
			if(temp>max) max=temp;
		}
		return max;
	}

	public float getStackSum(Object[] stack){
		float sum=jstatistics.getstatistic("Sum",stack[0],null);
		for(int i=1;i<stack.length;i++){
			sum+=jstatistics.getstatistic("Sum",stack[i],null);
		}
		return sum;
	}
}
