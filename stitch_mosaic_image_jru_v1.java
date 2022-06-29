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
import ij.plugin.frame.*;
import jguis.*;
import jalgs.*;
import ij.io.*;

public class stitch_mosaic_image_jru_v1 implements PlugIn, FrameInterface, gui_interface {
	boolean showrois;
	int tempchan,tempslice,tempframe,slices,channels,frames,trueframes;
	ImageStack tempstack;
	Object outimg;

	public void run(String arg) {
		/*ImagePlus imp=WindowManager.getCurrentImage();
		ImageWindow[] iws=jutils.selectPlots(false,1,new String[]{"Position_Plot"});
		if(iws==null) return;*/
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Mosaic_Stack","Position_Plot"});
		if(imps==null) return;
		ImageWindow iw=imps[1].getWindow();
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Ignore_scaling",false);
		gd.addCheckbox("Output_file",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean ignorescale=gd.getNextBoolean();
		boolean outfile=gd.getNextBoolean();
		if(outfile){
			SaveDialog sd=new SaveDialog("Save As","",".tif");
			String dir=sd.getDirectory();
			String fname=sd.getFileName();
			if(fname==null || fname.length()==0) return;
			exec(imps[0],iw,ignorescale,dir+fname);
		} else {
			Object[] retvals=exec(imps[0],iw,ignorescale);
			ImagePlus retimp=(ImagePlus)retvals[0];
			if(imps[0].getNChannels()==1) retimp.show();
			else{
				(new CompositeImage(retimp,CompositeImage.COLOR)).show();
			}
		}
	}

	public static float[] getAvgOverlap(float[][] coords,float width,float height){
		float hover=0.0f; int nhpairs=0;
		float vover=0.0f; int nvpairs=0;
		float halfwidth=0.5f*width;
		float halfheight=0.5f*height;
		for(int i=0;i<coords[0].length;i++){
			for(int j=(i+1);j<coords[0].length;j++){
				float xdist=Math.abs(coords[0][j]-coords[0][i]);
				float ydist=Math.abs(coords[1][j]-coords[1][i]);
				if(xdist>halfwidth && xdist<width && ydist<halfheight){
					hover+=xdist; nhpairs++;
				}
				if(xdist<halfwidth && ydist>halfheight && ydist<height){
					vover+=ydist; nvpairs++;
				}
			}
		}
		hover/=(float)nhpairs;
		vover/=(float)nvpairs;
		return new float[]{width-hover,height-vover};
	}

	public Object[] exec(ImagePlus imp,ImageWindow iw,boolean ignorescale){
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected"); if(sel<0) sel=0;
		float psize=(float)jutils.get_psize(imp);
		if(ignorescale) psize=1.0f;
		float[] overlap=getAvgOverlap(new float[][]{xvals[sel],yvals[sel]},(float)imp.getWidth()*psize,(float)imp.getHeight()*psize);
		IJ.log("avg overlap = "+overlap[0]+" , "+overlap[1]);
		return exec(imp,xvals[sel],yvals[sel],overlap[0],overlap[1],psize);
	}

	public Object[] exec(ImagePlus imp,float[] xvals,float[] yvals,float hoverlap,float voverlap,float psize){
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		channels=imp.getNChannels();
		slices=imp.getNSlices();
		frames=imp.getNFrames();
		if(frames==1){
			frames=slices; slices=1;
		}
		trueframes=1;
		if(frames>xvals.length){
			trueframes=(int)(frames/xvals.length);
			frames=xvals.length;
		}
		int typeindex=algutils.get_array_type(stack.getPixels(1));
		stitching sclass=new stitching(width,height,xvals,yvals,psize,typeindex,this);
		ImageStack stack2=new ImageStack(sclass.newwidth,sclass.newheight);
		tempstack=stack;
		for(int k=0;k<trueframes;k++){
			for(int i=0;i<slices;i++){
				for(int j=0;j<channels;j++){
					tempchan=j;
					tempslice=i;
					tempframe=k*frames;
					Object stitched=sclass.stitch_frame(this,frames,true,hoverlap,voverlap);
					stack2.addSlice("",stitched);
				}
			}
		}
		if(showrois){
			int[][] rois=sclass.getRois();
			RoiManager rman=RoiManager.getInstance();
			if(rman==null) rman=new RoiManager();
			for(int i=0;i<rois.length;i++){
				rman.addRoi(new Roi(rois[i][0],rois[i][1],rois[i][2],rois[i][3]));
			}
		}
		ImagePlus imp5=new ImagePlus("Stitched Image",stack2);
		imp5.copyScale(imp);
		imp5.setOpenAsHyperStack(true);
		imp5.setDimensions(channels,slices,trueframes);
		//imp5.show();
		return new Object[]{imp5};
	}

	public void exec(ImagePlus imp,ImageWindow iw,boolean ignorescale,String outpath){
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected"); if(sel<0) sel=0;
		float psize=(float)jutils.get_psize(imp);
		if(ignorescale) psize=1.0f;
		float[] overlap=getAvgOverlap(new float[][]{xvals[sel],yvals[sel]},(float)imp.getWidth()*psize,(float)imp.getHeight()*psize);
		IJ.log("avg overlap = "+overlap[0]+" , "+overlap[1]);
		exec(imp,xvals[sel],yvals[sel],overlap[0],overlap[1],psize,outpath);
		return;
	}

	public void exec(ImagePlus imp,float[] xvals,float[] yvals,float hoverlap,float voverlap,float psize,String outpath){
		//this version writes the output one frame at a time to a file
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		channels=imp.getNChannels();
		slices=imp.getNSlices();
		frames=imp.getNFrames();
		if(frames==1){
			frames=slices; slices=1;
		}
		trueframes=1;
		if(frames>xvals.length){
			trueframes=(int)(frames/xvals.length);
			frames=xvals.length;
		}
		int typeindex=algutils.get_array_type(stack.getPixels(1));
		stitching sclass=new stitching(width,height,xvals,yvals,psize,typeindex,this);
		//ImageStack stack2=new ImageStack(sclass.newwidth,sclass.newheight);
		tempstack=stack;
		Tiff_Writer tw=new Tiff_Writer(imp,sclass.newwidth,sclass.newheight,trueframes*slices*channels,this);
		tw.saveAsTiffStack2(outpath);
		for(int k=0;k<trueframes;k++){
			for(int i=0;i<slices;i++){
				for(int j=0;j<channels;j++){
					tempchan=j;
					tempslice=i;
					tempframe=k*frames;
					Object stitched=sclass.stitch_frame(this,frames,true,hoverlap,voverlap);
					//stack2.addSlice("",stitched);
					tw.saveFrame(stitched);
				}
			}
		}
		tw.endWrite();
	}

	public Object getNextFrame(){
		Object temp=jutils.get3DSlice(tempstack,tempframe,tempslice,tempchan,frames,slices,channels);
		tempframe++;
		return temp;
	}

	public void showMessage(String message){
		IJ.log(message);
	}
	
	public void showProgress(int currpos,int finalpos){
		IJ.showProgress(currpos,finalpos);
	}

}

