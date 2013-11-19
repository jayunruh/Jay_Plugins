/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.*;
import ij.process.*;
import jalgs.algutils;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

//import org.slf4j.LoggerFactory;
//import ch.qos.logback.classic.Level;
//import ch.qos.logback.classic.Logger;

import loci.formats.*;
import loci.formats.meta.IMetadata;
import loci.plugins.util.*;

public class LOCI_random_access_file_reader{
	// this plugin simply uses the loci library to open files
	public int nseries,width,height,channels,num,slices,frames;
	public float psize,tsize,zsize;
	public String order,name;
	public ImageProcessorReader r;
	public boolean nometa=false;
	
	public LOCI_random_access_file_reader(String directory,String fname,int series,boolean outmeta){
		IMetadata omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			if(series>=nseries)
				series=0;
			r.setSeries(series);
			num=r.getImageCount();
			width=r.getSizeX();
			height=r.getSizeY();
			channels=r.getSizeC();
			slices=r.getSizeZ();
			frames=r.getSizeT();
			order=r.getDimensionOrder();

			if(outmeta&&!nometa){
				Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
				if(globalMeta!=null)
					dumpMetaData(globalMeta);
				Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
				if(seriesMeta!=null)
					dumpMetaData(seriesMeta);
			}
			name=""+fname;
			if(nseries>1&&!nometa)
				name=omexmlMetadata.getImageName(series);
			else if(nseries>1)
				name=name+series;
			psize=1.0f;
			zsize=1.0f;
			tsize=1.0f;
			if(!nometa){
				if(omexmlMetadata.getPixelsPhysicalSizeX(series)!=null)
					psize=omexmlMetadata.getPixelsPhysicalSizeX(series).getValue().floatValue();
				if(omexmlMetadata.getPixelsPhysicalSizeZ(series)!=null)
					zsize=omexmlMetadata.getPixelsPhysicalSizeZ(series).getValue().floatValue();
				if(omexmlMetadata.getPixelsTimeIncrement(series)!=null)
					tsize=omexmlMetadata.getPixelsTimeIncrement(series).floatValue();
			}
		}catch(FormatException e){
			dispose();
		}catch(IOException e){
			dispose();
		}
	}
	
	public void dispose(){
		try{
			r.close();
		}catch(IOException e){
			r=null;
		}
	}
	
	public Object getSubImage(int x,int y,int rwidth,int rheight,int channel,int slice,int frame){
		int index=get_stack_index(channel,slice,frame,channels,slices,frames,order);
		try{
			if(r!=null){
				ImageProcessor ip=r.openProcessors(index,x,y,rwidth,rheight)[0];
				//Object pixels=r.openPlane(index,x,y,rwidth,rheight);
				Object pixels=ip.getPixels();
				return pixels;
			}else{
				return null;
			}
		} catch(Exception e){
			return null;
		}
	}
	
	public int getType(){
		Object pix=getSubImage(0,0,1,1,0,0,0);
		if(pix instanceof float[]) return 2;
		else if(pix instanceof short[]) return 1;
		else if(pix instanceof byte[]) return 0;
		else if(pix instanceof double[]) return 3;
		else return 4;
	}
	
	public float[] getLine(int y,int channel,int slice,int frame){
		Object pix=getSubImage(0,y,width,1,channel,slice,frame);
		return algutils.convert_arr_float(pix);
	}
	
	public float getPixel(int x,int y,int channel,int slice,int frame){
		Object pix=getSubImage(x,y,1,1,channel,slice,frame);
		if(pix!=null){
			if(pix instanceof float[]){
				return ((float[])pix)[0];
			} else if(pix instanceof short[]){
				return (float)(((short[])pix)[0]&0xffff);
			} else {
				return (float)(((byte[])pix)[0]&0xff);
			}
		}
		return Float.NaN;
	}
	
	public float[] getCArray(int x,int y,int slice,int frame){
		float[] carray=new float[channels];
		for(int i=0;i<channels;i++){
			carray[i]=getPixel(x,y,i,slice,frame);
		}
		return carray;
	}
	
	public float[][] getCCarpet(int y,int slice,int frame){
		float[][] carray=new float[channels][];
		for(int i=0;i<channels;i++){
			carray[i]=getLine(y,i,slice,frame);
		}
		return carray;
	}

	public int getNSeries(String directory,String fname){
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			r.close();
			return nseries;
		}catch(FormatException e){
			return 0;
		}catch(IOException e){
			return 0;
		}
	}

	public int get_stack_index(int channel,int slice,int frame,int channels,int slices,int frames,String order){
		// note that channel, slice, and frame must be 0 based
		if(order.equals("XYZCT")){
			return slice+slices*channel+channels*slices*frame;
		}else if(order.equals("XYCZT")){
			return channel+slice*channels+frame*channels*slices;
		}else if(order.equals("XYCTZ")){
			return channel+frame*channels+slice*frames*channels;
		}else if(order.equals("XYZTC")){
			return slice+slices*frame+slices*frames*channels;
		}else if(order.equals("XYTCZ")){
			return frame+frames*channel+frames*channels*slice;
		}else if(order.equals("XYTZC")){
			return frame+frames*slice+frames*slices*channel;
		}else{
			return 0;
		}
	}

	public void dumpMetaData(Hashtable<String,Object> metadata){
		ArrayList<String> keys=new ArrayList<String>(metadata.keySet());
		Collections.sort(keys);
		StringBuffer sb=new StringBuffer();
		for(String key:keys){
			sb.append(key);
			sb.append(" , ");
			sb.append(metadata.get(key));
			sb.append("\n");
		}
		IJ.log(sb.toString());
	}

}
