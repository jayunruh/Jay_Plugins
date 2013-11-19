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

public class LOCI_file_reader{
	// this plugin simply uses the loci library to open files
	public int nseries;
	public boolean nometa=false;

	public ImagePlus get_loci_imp(String directory,String fname){
		return get_loci_imp(directory,fname,false);
	}

	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta){
		return get_loci_imp(directory,fname,outmeta,0);
	}

	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta,int series){
		return get_loci_imp(directory,fname,outmeta,series,false,"");
	}

	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta,int series,boolean proj,String projstat){
		IMetadata omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			if(series>=nseries)
				series=0;
			r.setSeries(series);
			int num=r.getImageCount();
			int width=r.getSizeX();
			int height=r.getSizeY();
			int channels=r.getSizeC();
			int slices=r.getSizeZ();
			int frames=r.getSizeT();
			String order=r.getDimensionOrder();

			if(outmeta&&!nometa){
				Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
				if(globalMeta!=null)
					dumpMetaData(globalMeta);
				Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
				if(seriesMeta!=null)
					dumpMetaData(seriesMeta);
			}
			String name=""+fname;
			if(nseries>1&&!nometa)
				name=omexmlMetadata.getImageName(series);
			else if(nseries>1)
				name=name+series;
			float psize=1.0f;
			float zsize=1.0f;
			float tsize=1.0f;
			if(!nometa){
				if(omexmlMetadata.getPixelsPhysicalSizeX(series)!=null)
					psize=omexmlMetadata.getPixelsPhysicalSizeX(series).getValue().floatValue();
				if(omexmlMetadata.getPixelsPhysicalSizeZ(series)!=null)
					zsize=omexmlMetadata.getPixelsPhysicalSizeZ(series).getValue().floatValue();
				if(omexmlMetadata.getPixelsTimeIncrement(series)!=null)
					tsize=omexmlMetadata.getPixelsTimeIncrement(series).floatValue();
			}
			ImageStack stack=new ImageStack(width,height);
			if(proj){
				int counter=0;
				for(int i=0;i<frames;i++){
					for(int j=0;j<channels;j++){
						Object[] tempzstack=new Object[slices];
						for(int k=0;k<slices;k++){
							int index=get_stack_index(j,k,i,channels,slices,frames,order);
							ImageProcessor ip=r.openProcessors(index)[0];
							tempzstack[k]=ip.getPixels();
						}
						Object projslice=algutils.get_stack_proj_stat(projstat,tempzstack,width,height,slices,null);
						stack.addSlice("",projslice);
						IJ.showProgress(counter,frames*channels);
						counter++;
					}
				}
				slices=1;
			}else{
				int counter=0;
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							int index=get_stack_index(k,j,i,channels,slices,frames,order);
							ImageProcessor ip=r.openProcessors(index)[0];
							stack.addSlice(ip);
							IJ.showProgress(counter,frames*slices*channels);
							counter++;
						}
					}
				}
			}
			r.close();
			// shuffle(stack,channels,slices,frames,order);
			ImagePlus imp=new ImagePlus(name,stack);
			if(!nometa){
				jutils.set_psize(imp,psize);
				jutils.set_pdepth(imp,zsize);
				jutils.set_pinterval(imp,tsize);
			}
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(channels,slices,frames);
			if(channels>1){
				return new CompositeImage(imp,CompositeImage.COLOR);
			}else{
				return imp;
			}
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}

	public ImagePlus get_loci_imp_tile(String directory,String fname,boolean outmeta,int series,int xtiles,int ytiles){
		IMetadata omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			if(series>=nseries)
				series=0;
			r.setSeries(series);
			int num=r.getImageCount();
			int width=r.getSizeX();
			int height=r.getSizeY();
			int channels=r.getSizeC();
			int slices=r.getSizeZ();
			int frames=r.getSizeT();
			String order=r.getDimensionOrder();
			int newwidth=(int)((float)width/(float)xtiles);
			int newheight=(int)((float)height/(float)ytiles);

			if(outmeta&&!nometa){
				Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
				if(globalMeta!=null)
					IJ.log("Global Metadata");
					dumpMetaData(globalMeta);
				Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
				if(seriesMeta!=null)
					IJ.log("Series Metadata");
					dumpMetaData(seriesMeta);
			}
			String name=""+fname;
			if(nseries>1&&!nometa)
				name=omexmlMetadata.getImageName(series);
			else if(nseries>1)
				name=name+series;
			float psize=1.0f;
			float zsize=1.0f;
			float tsize=1.0f;
			if(!nometa){
				if(omexmlMetadata.getPixelsPhysicalSizeX(series)!=null)
					psize=omexmlMetadata.getPixelsPhysicalSizeX(series).getValue().floatValue();
				if(omexmlMetadata.getPixelsPhysicalSizeZ(series)!=null)
					zsize=omexmlMetadata.getPixelsPhysicalSizeZ(series).getValue().floatValue();
				if(omexmlMetadata.getPixelsTimeIncrement(series)!=null)
					tsize=omexmlMetadata.getPixelsTimeIncrement(series).floatValue();
			}
			ImageStack stack=new ImageStack(newwidth,newheight);
			int counter=0;
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int y=0;y<ytiles;y++){
						for(int x=0;x<xtiles;x++){
							for(int k=0;k<channels;k++){
								int index=get_stack_index(k,j,i,channels,slices,frames,order);
								ImageProcessor ip=r.openProcessors(index,x*newwidth,y*newheight,newwidth,newheight)[0];
								stack.addSlice(ip);
								IJ.showProgress(counter,frames*slices*xtiles*ytiles*channels);
								counter++;
							}
						}
					}
				}
			}
			r.close();
			// shuffle(stack,channels,slices,frames,order);
			ImagePlus imp=new ImagePlus(name,stack);
			if(!nometa){
				jutils.set_psize(imp,psize);
				jutils.set_pdepth(imp,zsize);
				jutils.set_pinterval(imp,tsize);
			}
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(channels,slices*xtiles*ytiles,frames);
			if(channels>1){
				return new CompositeImage(imp,CompositeImage.COLOR);
			}else{
				return imp;
			}
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}
	
	public Hashtable<String,Object> getGlobalMetaData(String directory,String fname){
		IMetadata omexmlMetadata=omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			//if(series>=nseries) series=0;
			//r.setSeries(series);
			Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
			r.close();
			return globalMeta;
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}
	
	public String get_metadata_value(String directory,String fname,int series,String key,boolean dump){
		IMetadata omexmlMetadata=omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			if(series>=nseries) series=0;
			r.setSeries(series);
			Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
			Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
			if(dump){
				IJ.log("Global MetaData");
				dumpMetaData(globalMeta);
				IJ.log("Series MetaData");
				dumpMetaData(seriesMeta);
			}
			if(globalMeta!=null){
				String temp=findMetaDataKeyVal(globalMeta,key);
				if(temp!=null){
					r.close();
					return temp;
				}
			}
			if(seriesMeta!=null){
				String temp=findMetaDataKeyVal(seriesMeta,key);
				if(temp!=null){
					r.close();
					return temp;
				}
			}
			r.close();
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
		return null;
	}

	public String[] getSeriesNames(String directory,String fname){
		IMetadata omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageReader r=new ImageReader();
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			String[] names=new String[nseries];
			if(!nometa){
				for(int i=0;i<nseries;i++){
					r.setSeries(i);
					names[i]=omexmlMetadata.getImageName(i);
					if(names[i]==null||names[i]=="")
						names[i]="Series"+i;
				}
			}else{
				for(int i=0;i<nseries;i++){
					names[i]=fname+(i+1);
				}
			}
			r.close();
			return names;
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
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

	public void shuffle(ImageStack stack,int channels,int slices,int frames,String order){
		// here we convert any stack to xyczt order
		if(order.equals("XYZCT")){
			Object[] pix1=stack.getImageArray();
			String[] labels1=stack.getSliceLabels();
			Object[] pix2=new Object[pix1.length];
			System.arraycopy(pix1,0,pix2,0,pix1.length);
			String[] labels2=new String[labels1.length];
			System.arraycopy(labels1,0,labels2,0,labels1.length);
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<channels;k++){
						pix1[k+j*channels+i*slices*channels]=pix2[j+k*slices+i*channels*slices];
						labels1[k+j*channels+i*slices*channels]=labels2[j+k*slices+i*channels*slices];
					}
				}
			}
		}else{
			return;
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
	
	public static String findMetaDataKeyVal(Hashtable<String,Object> metadata,String fkey){
		ArrayList<String> keys=new ArrayList<String>(metadata.keySet());
		Collections.sort(keys);
		for(String key:keys){
			if(key.indexOf(fkey)>=0){
				return metadata.get(key).toString();
			}
		}
		return null;
	}

}
