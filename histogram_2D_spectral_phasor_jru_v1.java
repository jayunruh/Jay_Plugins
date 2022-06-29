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
import java.awt.event.*;
import java.awt.image.*;
import ij.plugin.*;
import ij.measure.*;
import java.io.*;
import jguis.*;
import jalgs.*;

public class histogram_2D_spectral_phasor_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int nchans=imp.getNChannels();
		ImagePlus[] imps=stack2gsi(imp.getStack(),nchans,1);
		//double[][] matrix={{0.14,-0.01,-0.02},{0.32,0.12,0.30},{1.0,1.0,1.0}};
		//float[] xpts={0.14f,-0.01f,-0.02f};
		//float[] ypts={0.32f,0.12f,0.30f};
		//double[][] minv=(new matrixsolve()).gjinv2(matrix,3);
		//IJ.log(table_tools.print_double_array(minv));
		//boolean inside=PU.contains(0.012f,0.22f,xpts,ypts);
		//IJ.log(""+inside);
		final Hist2DWindow cw=new Hist2DWindow();
		cw.init(imps[0],imps[1],null,imp,imps[2],4);
		Hist2DWindow.launch_frame("Interactive 2D Histogram",cw);
	}

	ImagePlus[] stack2gsi(ImageStack stack,int nchans,int harmonic){
		int slices=stack.getSize();
		int width=stack.getWidth();
		int height=stack.getHeight();
		int result_slices=1;
		if(slices!=nchans){
			result_slices=(int)((float)slices/(float)nchans);
		}
		ImageStack gstack=new ImageStack(width,height);
		ImageStack sstack=new ImageStack(width,height);
		ImageStack intstack=new ImageStack(width,height);
		float[] cosvals=new float[nchans];
		float[] sinvals=new float[nchans];
		if(nchans==4){
			cosvals=new float[]{1.0f,-1.0f,-1.0f,1.0f};
			sinvals=new float[]{1.0f,1.0f,-1.0f,-1.0f};
		} else {
			for(int i=0;i<nchans;i++){
				cosvals[i]=(float)Math.cos(2.0*Math.PI*((double)(i*harmonic)/(double)nchans));
				sinvals[i]=(float)Math.sin(2.0*Math.PI*((double)(i*harmonic)/(double)nchans));
			}
		}
		for(int k=0;k<result_slices;k++){
			float[] intensity=new float[width*height];
			float[] G=new float[width*height];
			float[] S=new float[width*height];
			if(stack.getProcessor(1) instanceof FloatProcessor){
				for(int i=0;i<(width*height);i++){
					for(int j=0;j<nchans;j++){
						float temp=((float[])stack.getPixels(j+k*nchans+1))[i];
						intensity[i]+=temp;
						G[i]+=temp*cosvals[j];
						S[i]+=temp*sinvals[j];
					}
					if(intensity[i]>0.0f){
						G[i]/=intensity[i];
						S[i]/=intensity[i];
					}
				}
			}
			if(stack.getProcessor(1) instanceof ShortProcessor){
				for(int i=0;i<(width*height);i++){
					for(int j=0;j<nchans;j++){
						float temp=(float)(((short[])stack.getPixels(j+k*nchans+1))[i]&0xffff);
						intensity[i]+=temp;
						G[i]+=temp*cosvals[j];
						S[i]+=temp*sinvals[j];
					}
					if(intensity[i]>0.0f){
						G[i]/=intensity[i];
						S[i]/=intensity[i];
					}
				}
			}
			if(stack.getProcessor(1) instanceof ByteProcessor){
				for(int i=0;i<(width*height);i++){
					for(int j=0;j<nchans;j++){
						float temp=(float)(((byte[])stack.getPixels(j+k*nchans+1))[i]&0xff);
						intensity[i]+=temp;
						G[i]+=temp*cosvals[j];
						S[i]+=temp*sinvals[j];
					}
					if(intensity[i]>0.0f){
						G[i]/=intensity[i];
						S[i]/=intensity[i];
					}
				}
			}
			gstack.addSlice("",G);
			sstack.addSlice("",S);
			intstack.addSlice("",intensity);
		}
		ImagePlus[] imps=new ImagePlus[3];
		imps[0]=new ImagePlus("G",gstack);
		imps[1]=new ImagePlus("S",sstack);
		imps[2]=new ImagePlus("Intensity",intstack);
		return imps;
	}

}
