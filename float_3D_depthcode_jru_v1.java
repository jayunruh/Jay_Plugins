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
import ij.measure.*;
import jalgs.*;
import jguis.*;

public class float_3D_depthcode_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int slices=imp.getNSlices();
		int channels=imp.getNChannels();
		int frames=imp.getNFrames();
		GenericDialog gd=new GenericDialog("Options");
		String[] filters={"avg depth","max depth","min depth"};
		gd.addChoice("Filter?",filters,filters[0]);
		String[] projections={"c","z","t"};
		gd.addChoice("Projection?",projections,projections[1]);
		int mindepth=0;
		gd.addNumericField("Minimum depth",mindepth,0);
		int maxdepth=slices-1;
		if(maxdepth==0) maxdepth=channels-1;
		if(maxdepth==0) maxdepth=frames-1;
		gd.addNumericField("Maximum depth",maxdepth,0);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int filterindex=gd.getNextChoiceIndex();
		int projindex=gd.getNextChoiceIndex();
		mindepth=(int)gd.getNextNumber();
		maxdepth=(int)gd.getNextNumber();
		ImageStack dstack=new ImageStack(width,height);
		ImageStack codestack=new ImageStack(width,height);
		int outchannels=channels;
		int outslices=slices;
		int outframes=frames;
		if(projindex==0){
			outchannels=1;
			for(int f=0;f<frames;f++){
				for(int z=0;z<slices;z++){
					float[] intresult=new float[height*width];
					float[] dresult=new float[height*width];
					Object[] cseries=jutils.get3DCSeries(stack,z,f,frames,slices,channels);
					for(int i=0;i<width*height;i++){
						float[] col=algutils.convert_arr_float(algutils.get_stack_col(cseries,width,height,i,cseries.length));
						dresult[i]=get_depth(col,filterindex);
						intresult[i]=get_proj(col,filterindex);
					}
					int[] colorimage=depthcode(intresult,dresult,(float)mindepth,(float)maxdepth);
					codestack.addSlice("",colorimage);
					dstack.addSlice("",dresult);
					IJ.showProgress(z+f*slices,frames*slices);
				}
			}
		}
		if(projindex==1){
			outslices=1;
			for(int f=0;f<frames;f++){
				for(int c=0;c<channels;c++){
					float[] intresult=new float[height*width];
					float[] dresult=new float[height*width];
					Object[] zseries=jutils.get3DZSeries(stack,c,f,frames,slices,channels);
					for(int i=0;i<width*height;i++){
						float[] col=algutils.convert_arr_float(algutils.get_stack_col(zseries,width,height,i,zseries.length));
						dresult[i]=get_depth(col,filterindex);
						intresult[i]=get_proj(col,filterindex);
					}
					int[] colorimage=depthcode(intresult,dresult,(float)mindepth,(float)maxdepth);
					codestack.addSlice("",colorimage);
					dstack.addSlice("",dresult);
					IJ.showProgress(c+f*channels,frames*channels);
				}
			}
		}
		if(projindex==2){
			outframes=1;
			for(int z=0;z<slices;z++){
				for(int c=0;c<channels;c++){
					float[] intresult=new float[height*width];
					float[] dresult=new float[height*width];
					Object[] tseries=jutils.get3DTSeries(stack,z,c,frames,slices,channels);
					for(int i=0;i<width*height;i++){
						float[] col=algutils.convert_arr_float(algutils.get_stack_col(tseries,width,height,i,tseries.length));
						dresult[i]=get_depth(col,filterindex);
						intresult[i]=get_proj(col,filterindex);
					}
					int[] colorimage=depthcode(intresult,dresult,(float)mindepth,(float)maxdepth);
					codestack.addSlice("",colorimage);
					dstack.addSlice("",dresult);
					IJ.showProgress(c+z*channels,channels*slices);
				}
			}
		}
		Calibration cal=(imp.getCalibration()).copy();
		ImagePlus dimp=new ImagePlus("Depth Image",dstack);
		dimp.copyScale(imp);
		dimp.setOpenAsHyperStack(true);
		dimp.setDimensions(outchannels,outslices,outframes);
		dimp.show();
		ImagePlus imp2=new ImagePlus(projections[projindex]+" "+filters[filterindex]+" projection",codestack);
		imp2.copyScale(imp);
		imp2.setOpenAsHyperStack(true);
		imp2.setDimensions(outchannels,outslices,outframes);
		imp2.show();
	}

	private int[] depthcode(float[] intensity, float[] depth,float mindepth,float maxdepth){
		float maxint=maxvector(intensity);
		float minint=minvector(intensity);
		int[] color_pixels=new int[intensity.length];
		int r,g,b;
		for(int i=0;i<intensity.length;i++){
			float fdepth=(depth[i]-mindepth)/(maxdepth-mindepth);
			float fint=(intensity[i]-minint)/(maxint-minint);
			if(fdepth<0.5f){
				if(fdepth<0.0f){fdepth=0.0f;}
				r=0;
				g=(int)((2.0f*fdepth)*fint*256.9999f);
				b=(int)((1.0f-2.0f*fdepth)*fint*255.9999f);
			} else {
				if(fdepth>1.0f){fdepth=1.0f;}
				b=0;
				r=(int)((2.0f*(fdepth-0.5f))*fint*256.9999f);
				g=(int)((1.0f-2.0f*(fdepth-0.5f))*fint*255.9999f);
			}
			color_pixels[i]=0xff000000 + (r<<16) + (g<<8) + b;
		}
		return color_pixels;
	}

	float get_depth(float[] values,int filterindex){
		if(filterindex==0) return avgindex(values);
		if(filterindex==1) return maxindex(values);
		return minindex(values);
	}

	float get_proj(float[] values,int filterindex){
		if(filterindex==0) return avgvector(values);
		if(filterindex==1) return maxvector(values);
		return minvector(values);
	}

	float minindex(float[] values){
		int length=values.length;
		float minval=values[0];
		int temp=0;
		for(int i=1;i<length;i++){
			if(values[i]<minval){minval=values[i]; temp=i;}
		}
		return (float)temp;
	}

	float maxindex(float[] values){
		int length=values.length;
		float maxval=values[0];
		int temp=0;
		for(int i=1;i<length;i++){
			if(values[i]>maxval){maxval=values[i]; temp=i;}
		}
		return (float)temp;
	}

	float avgindex(float[] values){
		int length;
		float dumfloat=0.0f;
		float dumfloat2=0.0f;
		length=values.length;
		for(int i=0;i<length;i++){
			dumfloat+=values[i]/(float)length;
			dumfloat2+=(values[i]*(float)i)/(float)length;
		}
		dumfloat2/=dumfloat;
		if(dumfloat2>(float)length) dumfloat2=(float)length;
		if(dumfloat2<0.0f) dumfloat2=0.0f;
		return dumfloat2;
	}

	float minvector(float[] values){
		int length=values.length;
		float minval=values[0];
		for(int i=1;i<length;i++){
			if(values[i]<minval){minval=values[i];}
		}
		return minval;
	}

	float maxvector(float[] values){
		int length=values.length;
		float maxval=values[0];
		for(int i=1;i<length;i++){
			if(values[i]>maxval){maxval=values[i];}
		}
		return maxval;
	}

	float avgvector(float[] values){
		int length;
		float dumfloat=0.0f;
		length=values.length;
		for(int i=0;i<length;i++){
			dumfloat+=values[i]/length;
		}
		return dumfloat;
	}

}
