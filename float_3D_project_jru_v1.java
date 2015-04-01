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

public class float_3D_project_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] filters=jstatistics.stats;
		gd.addChoice("Filter?",filters,filters[0]);
		String[] projections={"x","y","z"};
		gd.addChoice("Projection?",projections,projections[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int filterindex=gd.getNextChoiceIndex();
		String stat=filters[filterindex];
		int projindex=gd.getNextChoiceIndex();
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int frames=imp.getNFrames();
		int slices=imp.getNSlices();
		int channels=imp.getNChannels();
		if(slices==1){frames=slices; frames=1;}
		float[] result=new float[height*width];
		int reswidth=width;
		int resheight=height;
		ImageStack resstack=null;
		if(projindex==0){
			reswidth=height;
			resheight=slices;
			resstack=new ImageStack(reswidth,resheight);
			for(int i=0;i<frames;i++){
				for(int j=0;j<channels;j++){
					Object[] substack=jutils.get3DZSeries(stack,j,i,frames,slices,channels);
					float[] proj=get_projx(substack,width,height,stat);
					resstack.addSlice("",proj);
					IJ.showProgress(j+i*channels,channels*frames);
				}
			}
		}
		if(projindex==1){
			reswidth=width;
			resheight=slices;
			resstack=new ImageStack(reswidth,resheight);
			for(int i=0;i<frames;i++){
				for(int j=0;j<channels;j++){
					Object[] substack=jutils.get3DZSeries(stack,j,i,frames,slices,channels);
					float[] proj=get_projy(substack,width,height,stat);
					resstack.addSlice("",proj);
					IJ.showProgress(j+i*channels,channels*frames);
				}
			}
		}
		if(projindex==2){
			reswidth=width;
			resheight=height;
			resstack=new ImageStack(reswidth,resheight);
			for(int i=0;i<frames;i++){
				for(int j=0;j<channels;j++){
					Object[] substack=jutils.get3DZSeries(stack,j,i,frames,slices,channels);
					float[] proj=get_projz(substack,width,height,stat);
					resstack.addSlice("",proj);
					IJ.showProgress(j+i*channels,channels*frames);
				}
			}
		}
		String newname=projections[projindex]+" "+filters[filterindex]+" projection";
		Calibration cal=(imp.getCalibration()).copy();
		if(reswidth>1 && resheight>1){
			ImagePlus imp2=new ImagePlus(newname,resstack);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(channels,1,frames);
			imp2.setCalibration(cal);
			Calibration cal2=imp2.getCalibration();
			if(projindex==0){
				cal2.pixelHeight=1.0f;
			}
			if(projindex==1){
				cal2.pixelWidth=cal2.pixelHeight;
				cal2.pixelHeight=1.0;
			}
			if(imp2.getNChannels()>1){
				CompositeImage ci=new CompositeImage(imp2);
				if(imp.isComposite()){
					LUT[] lut=((CompositeImage)imp).getLuts();
					ci.setLuts(lut);
					ci.resetDisplayRanges();
					ci.show();
				}
			} else {
				imp2.show();
			}
		}
		else{
			float xmin=(float)(-cal.xOrigin*cal.pixelWidth);
			float ymin=(float)((cal.yOrigin-(double)(height-1))*cal.pixelHeight);
			int length=result.length;
			float[] xvals = new float[length];
			if(projindex==0){for(int i=0;i<length;i++){xvals[i]=xmin+(float)cal.pixelWidth*(float)i;}}
			else{for(int i=0;i<length;i++){xvals[i]=ymin+(float)cal.pixelHeight*(float)(length-i-1);}}
			PlotWindow4 pw=new PlotWindow4(newname,"Pixel",stat,xvals,(float[])resstack.getPixels(1));
			pw.draw();
			for(int i=1;i<channels*frames;i++){
				pw.addPoints(xvals,(float[])resstack.getPixels(i+1),true);
			}
		}
	}

	public float[] get_projx(Object[] image,int width,int height,String stat){
		int slices=image.length;
		float[] proj=new float[height*slices];
		for(int j=0;j<slices;j++){
			for(int i=0;i<height;i++){
				Object row=algutils.get_image_row(image[j],width,height,i);
				proj[i+j*height]=jstatistics.getstatistic(stat,row,null);
			}
		}
		return proj;
	}

	public float[] get_projy(Object[] image,int width,int height,String stat){
		int slices=image.length;
		float[] proj=new float[width*slices];
		for(int j=0;j<slices;j++){
			for(int i=0;i<width;i++){
				Object col=algutils.get_image_col(image[j],width,height,i);
				proj[i+j*width]=jstatistics.getstatistic(stat,col,null);
			}
		}
		return proj;
	}

	public float[] get_projz(Object[] image,int width,int height,String stat){
		int slices=image.length;
		float[] proj=new float[width*height];
		for(int j=0;j<height;j++){
			for(int i=0;i<width;i++){
				Object col=algutils.get_stack_col(image,width,height,i,j,slices);
				proj[i+j*width]=jstatistics.getstatistic(stat,col,null);
			}
		}
		return proj;
	}

}
