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
import ij.io.*;
import jguis.*;
import jalgs.*;
import java.io.*;

public class cat_hyperstacks_jru_v1 implements PlugIn, FrameInterface {
	String directory,extension;
	String[] flist;
	int frameiterator,nfiles,nframes,nslices,nchannels,stacksize,fileiterator;
	int channeliterator,sliceiterator,width,height;
	double scale;
	ImagePlus currimp;
	boolean maxproj;
	//Opener opener;

	public void run(String arg) {
		OpenDialog od = new OpenDialog("Open Image Sequence...", arg);
        		directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Sort_Names_Numerically",false);
		gd.addCheckbox("Sort_Date_Modified",false);
		gd.addStringField("File_Name_Contains:","",10);
		gd.addNumericField("Scale Factor",1.0,5,15,null);
		gd.addCheckbox("Force Z Stack",true);
		gd.addCheckbox("Max_Project",false);
		gd.addCheckbox("Open_Virtual_After",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean sortnum=gd.getNextBoolean();
		boolean sortdatemod=gd.getNextBoolean();
		String filter=gd.getNextString();
		scale=gd.getNextNumber();
		boolean forcez=gd.getNextBoolean();
		maxproj=gd.getNextBoolean();
		boolean openafter=gd.getNextBoolean();
		int period=fname.lastIndexOf('.');
		extension=fname.substring(period+1);
		String[] mask=null;
		if(filter!="" && filter!=null && filter.length()>0){
			mask=new String[2]; mask[1]=filter; mask[0]=extension;
			//IJ.log("masks: "+mask[0]+" , "+mask[1]);
		} else {
			mask=new String[1]; mask[0]=extension;
		}
		if(sortnum){
			flist=(new jdataio()).get_numeric_sorted_string_list(directory,mask,null);
		} else {
			if(sortdatemod){
				flist=(new jdataio()).get_datemod_sorted_string_list(directory,mask,null);
			} else {
				flist=(new jdataio()).get_sorted_string_list(directory,mask,null);
			}
		}
		for(int i=0;i<flist.length;i++){
			IJ.log(""+flist[i]);
		}
		if(flist==null){return;}
		nfiles=flist.length;
		//IJ.showMessage("test");
		//opener=new Opener(); opener.setSilentMode(true);
		if(!extension.equals("lsm")){
			//opener.openImage(directory,flist[0]);
			//currimp=WindowManager.getCurrentImage();
			currimp=(new LOCI_file_reader()).get_loci_imp(directory,flist[0]);
		} else {
			currimp=(new LSM_file_reader()).open(directory,flist[0]);
		}
		if(currimp==null){return;}
		width=currimp.getWidth(); height=currimp.getHeight();
		nchannels=currimp.getNChannels();
		nslices=currimp.getNSlices();
		nframes=currimp.getNFrames();
		if(forcez && nslices==1){
			nslices=nframes;
			nframes=1;
			currimp.setDimensions(nchannels,nslices,nframes);
		}
		stacksize=nchannels*nslices*nframes;
		if(maxproj) stacksize/=nslices;
		SaveDialog sd=new SaveDialog("Save as Tiff...",currimp.getTitle(),".tif");
		String fileName=sd.getFileName();
		String fileDir=sd.getDirectory();
		if(fileName==null){return;}
		Tiff_Writer tw=new Tiff_Writer(currimp,nfiles*stacksize,this);
		frameiterator=0;
		channeliterator=0;
		sliceiterator=0;
		fileiterator=0;
		tw.saveAsTiffStack(fileDir+fileName);
		currimp.changes=false;
		currimp.close();
		if(openafter) IJ.run("TIFF Virtual Stack...","open="+fileDir+fileName);
	}

	public Object getNextFrame(){
		if(channeliterator>=nchannels){
			channeliterator=0;
			sliceiterator++;
		}
		if(sliceiterator>=nslices){
			sliceiterator=0;
			frameiterator++;
		}
		if(maxproj && sliceiterator>0){
			sliceiterator=0;
			frameiterator++;
		}
		if(frameiterator>=nframes){
			currimp.changes=false;
			currimp.close();
			frameiterator=0;
			fileiterator++;
			IJ.showProgress(fileiterator,nfiles-1);
			if(!extension.equals("lsm")){
				//opener.openImage(directory,flist[fileiterator]);
				//currimp=WindowManager.getCurrentImage();
				currimp=(new LOCI_file_reader()).get_loci_imp(directory,flist[fileiterator]);
			} else {
				currimp=(new LSM_file_reader()).open(directory,flist[fileiterator]);
			}
		}
		ImageProcessor ip=null;
		if(maxproj){
			Object pixels=jutils.get3DProjZStat2(currimp.getStack(),frameiterator,channeliterator,nframes,nslices,nchannels,"Max");
			if(pixels instanceof float[]) ip=new FloatProcessor(currimp.getWidth(),currimp.getHeight(),(float[])pixels,null);
			if(pixels instanceof short[]) ip=new ShortProcessor(currimp.getWidth(),currimp.getHeight(),(short[])pixels,null);
			if(pixels instanceof byte[]) ip=new ByteProcessor(currimp.getWidth(),currimp.getHeight(),(byte[])pixels,null);
		} else {
			ip=currimp.getStack().getProcessor(frameiterator*nslices*nchannels+sliceiterator*nchannels+channeliterator+1);
		}
		channeliterator++;
		if(scale!=1.0){
			ImageProcessor ip2=ip.resize((int)(width*scale),(int)(height*scale));
			return ip2.getPixels();
		} else {
			return ip.getPixels();
		}
	}

}
