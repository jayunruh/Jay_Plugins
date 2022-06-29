/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import ij.text.*;
import java.util.*;
import jguis.*;
import ij.io.*;
import jalgs.jseg.*;
import jalgs.*;

public class create_montages_from_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addNumericField("Number_Of_Columns",2,0);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		int ncols=(int)gd2.getNextNumber();
		GenericDialog gd=new GenericDialog("Options");
		for(int i=0;i<ncols;i++){
			gd.addChoice("Image"+(i+1)+"_Column",col_labels,col_labels[0]);
		}
		gd.addChoice("Slice_Name_Column",col_labels,col_labels[0]);
		gd.addNumericField("Bin_By",2,0);
		gd.addCheckbox("Force_Color",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int[] selcols=new int[ncols];
		for(int i=0;i<ncols;i++){
			selcols[i]=gd.getNextChoiceIndex();
		}
		int scol=gd.getNextChoiceIndex();
		int binby=(int)gd.getNextNumber();
		boolean forcecolor=gd.getNextBoolean();
		DirectoryChooser dc=new DirectoryChooser("Choose Directory");
		String dir=dc.getDirectory();
		if(dir==null) return;
		ImageStack stack=null;
		int nchans=1;
		int nslices=1;
		int width=1;
		int height=1;
		for(int i=0;i<listtable.size();i++){
			ImagePlus[] imps=new ImagePlus[ncols];
			List<String> temprow=listtable.get(i);
			for(int j=0;j<ncols;j++){
				String temp="";
				if(temprow.size()>selcols[j]) temp=temprow.get(selcols[j]);
				if(temp.length()>1) imps[j]=IJ.openImage(dir+listtable.get(i).get(selcols[j]));
				else imps[j]=makeBlankImp(imps[0]);
				if(imps[j]==null) break;
				//IJ.log(""+imps[j].getWidth());
			}
			if(imps[ncols-1]==null) continue;
			String name=listtable.get(i).get(scol);
			ImageStack[] stacks=new ImageStack[ncols];
			for(int j=0;j<ncols;j++) stacks[j]=imps[j].getStack();
			if(i==0){
				width=imps[0].getWidth();
				height=imps[0].getHeight();
				stack=new ImageStack(ncols*width/binby,height/binby);
				nchans=imps[0].getNChannels();
				nslices=imps[0].getNSlices();
			}
			for(int j=0;j<stacks[0].getSize();j++){
				Object[] pix=new Object[ncols];
				for(int k=0;k<ncols;k++){
					if(binby!=1) pix[k]=jsmooth.bin2D(stacks[k].getPixels(j+1),width,height,binby,true);
					else pix[k]=stacks[k].getPixels(j+1);	
				}
				if(forcecolor) pix=forceColor(pix);
				//IJ.log(pix[0].toString());
				Object temp=makeMontage(pix,width/binby,height/binby);
				stack.addSlice(name,temp);
			}
			IJ.showProgress(i,listtable.size());
		}
		jutils.create_hyperstack("Montage_Stack",stack,listtable.size(),nslices,nchans,true,null).show();
	}

	public ImagePlus makeBlankImp(ImagePlus template){
		int width=template.getWidth();
		int height=template.getHeight();
		int slices=template.getStack().getSize();
		int dtype=algutils.get_array_type(template.getStack().getPixels(1));
		Object[] blankarray=new Object[slices];
		for(int i=0;i<slices;i++) blankarray[i]=algutils.create_array(width*height,dtype);
		ImageStack stack=jutils.array2stack(blankarray,width,height);
		return new ImagePlus("blank",stack);
	}

	public Object[] forceColor(Object[] imgs){
		//here we convert all images to color images
		//scale 32 bit and 16 bit between min and max
		Object[] newimgs=new Object[imgs.length];
		for(int i=0;i<imgs.length;i++){
			if(imgs[i] instanceof byte[]){
				byte[] temp1=(byte[])imgs[i];
				int[] temp2=new int[temp1.length];
				for(int j=0;j<temp2.length;j++){
					temp2[j]=jutils.rgb2intval(temp1[j],temp1[j],temp1[j]);
				}
				newimgs[i]=temp2;
			} else if(imgs[i] instanceof short[]){
				short[] temp1=(short[])imgs[i];
				float min=jstatistics.getstatistic("Min",temp1,null);
				float max=jstatistics.getstatistic("Max",temp1,null);
				float range=max-min; if(range==0.0f) range=1.0f;
				int[] temp2=new int[temp1.length];
				for(int j=0;j<temp2.length;j++){
					float grayval=255.0f*((float)(temp1[j]&0xffff)-min)/range;
					temp2[j]=jutils.rgb2intval(grayval,grayval,grayval);
				}
				newimgs[i]=temp2;
			} else if(imgs[i] instanceof float[]){
				float[] temp1=(float[])imgs[i];
				float min=jstatistics.getstatistic("Min",temp1,null);
				float max=jstatistics.getstatistic("Max",temp1,null);
				float range=max-min; if(range==0.0f) range=1.0f;
				int[] temp2=new int[temp1.length];
				for(int j=0;j<temp2.length;j++){
					float grayval=255.0f*(temp1[j]-min)/range;
					temp2[j]=jutils.rgb2intval(grayval,grayval,grayval);
				}
				newimgs[i]=temp2;
			} else {
				newimgs[i]=imgs[i];
			}
		}
		return newimgs;
	}

	public Object makeMontage(Object[] imgs,int width,int height){
		if(imgs[0] instanceof byte[]){
			int newwidth=width*imgs.length;
			byte[] newimg=new byte[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((byte[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		} else if(imgs[0] instanceof short[]){
			int newwidth=width*imgs.length;
			short[] newimg=new short[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((short[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		} else if(imgs[0] instanceof int[]){
			int newwidth=width*imgs.length;
			int[] newimg=new int[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((int[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		} else {
			int newwidth=width*imgs.length;
			float[] newimg=new float[newwidth*height];
			for(int i=0;i<height;i++){
				for(int j=0;j<imgs.length;j++){
					System.arraycopy((float[])imgs[j],i*width,newimg,i*newwidth+j*width,width);
				}
			}
			return newimg;
		}
	}

}
