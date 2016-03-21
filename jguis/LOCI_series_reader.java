/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.process.ImageProcessor;
import jalgs.FrameInterface;
import jalgs.algutils;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import ome.units.UNITS;
import ome.units.quantity.Length;
//import org.slf4j.LoggerFactory;
//import ch.qos.logback.classic.Level;
//import ch.qos.logback.classic.Logger;
import loci.formats.ChannelSeparator;
import loci.formats.FormatException;
import loci.formats.ImageReader;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;

public class LOCI_series_reader implements FrameInterface{
	// this plugin interactively loads series' from a multiseries loci supported file
	public int nseries;
	public boolean nometa=false;
	public ImageProcessorReader r;
	public IMetadata omexmlMetadata;
	public int s;
	public String[] names;

	public LOCI_series_reader(String directory,String fname,boolean nometa){
		this.nometa=nometa;
		this.s=-1;
		omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
			r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			int maxdigits=Integer.toString(nseries).length();
			names=new String[nseries];
			if(!nometa){
				for(int i=0;i<nseries;i++){
					r.setSeries(i);
					names[i]=omexmlMetadata.getImageName(i);
					if(names[i]==null||names[i]=="")
						names[i]="Series"+pad_number(i+1,maxdigits);
				}
			}else{
				for(int i=0;i<nseries;i++){
					names[i]=fname+pad_number(i+1,maxdigits);
				}
			}
		}catch(FormatException e){
			IJ.log("Format Exception");
			dispose();
		}catch(IOException e){
			IJ.log("IO Exception");
			dispose();
		}
	}
	
	public int getNSeries(){return nseries;}
	
	public String[] getNames(){return names;}

	public Object getNextFrame(){
		s+=1;
		if(s>=nseries){
			dispose();
			return null;
		}
		if(r==null) return null;
		try{
			r.setSeries(s);
			//int num=r.getImageCount();
			int width=r.getSizeX();
			int height=r.getSizeY();
			int channels=r.getSizeC();
			int slices=r.getSizeZ();
			int frames=r.getSizeT();
			String order=r.getDimensionOrder();
			float psize=1.0f;
			float zsize=1.0f;
			float tsize=1.0f;
			if(!nometa){
				//Length temp=new Length(1.0,UNITS.MICROM);
				if(omexmlMetadata.getPixelsPhysicalSizeX(s)!=null)
					psize=omexmlMetadata.getPixelsPhysicalSizeX(s).value().floatValue();
				if(omexmlMetadata.getPixelsPhysicalSizeZ(s)!=null)
					zsize=omexmlMetadata.getPixelsPhysicalSizeZ(s).value().floatValue();
				if(omexmlMetadata.getPixelsTimeIncrement(s)!=null)
					tsize=omexmlMetadata.getPixelsTimeIncrement(s).value().floatValue();
			}
			ImageStack stack=new ImageStack(width,height);
			//int counter=0;
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<channels;k++){
						int index=get_stack_index(k,j,i,channels,slices,frames,order,null);
						//IJ.log(""+k+"\t "+j+"\t "+i+"\t "+index);
						ImageProcessor ip=r.openProcessors(index)[0];
						stack.addSlice(ip);
						//IJ.showProgress(counter,frames*slices*channels);
						//counter++;
					}
				}
			}
			ImagePlus imp=new ImagePlus(names[s],stack);
			if(!nometa){
				jutils.set_psize(imp,psize);
				jutils.set_pdepth(imp,zsize);
				jutils.set_pinterval(imp,tsize);
			}
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(channels,slices,frames);
			if(channels>1){
				imp=new CompositeImage(imp,CompositeImage.COLOR);
			}
			return imp;
		} catch(FormatException e){
			IJ.log("Format Exception");
			return null;
		}catch(IOException e){
			IJ.log("IO Exception");
			return null;
		}
	}

	public void dispose(){
		try{
			r.close();
			r=null;
		}catch(IOException e){
			r=null;
			return;
		}
	}

	public int get_stack_index(int channel,int slice,int frame,int channels,int slices,int frames,String order,int[] chlengths){
		// note that channel, slice, and frame must be 0 based
		if(order.equals("XYZCT")){
			//here we have the possibility of unequal slice numbers for channels
			if(chlengths==null) return slice+channel*slices+frame*channels*slices;
			else{
    			int temp=0;
    			for(int i=0;i<channel;i++) temp+=chlengths[i];
    			if(slice<chlengths[channel]) temp+=slice;
    			else temp+=chlengths[channel]-1;
    			if(frames==1 || frame==0) return temp;
    			else{
    				int fsize=0;
    				for(int i=0;i<channels;i++) fsize+=chlengths[i];
    				return fsize*frame+temp;
    			}
			}
		}else if(order.equals("XYCZT")){
			return channel+slice*channels+frame*channels*slices;
		}else{
			return 0;
		}
	}
	
	public String pad_number(int num,int len){
		String temp=Integer.toString(num);
		if(temp.length()==len) return temp;
		if(temp.length()>len) return temp.substring(0,len);
		for(int i=temp.length();i<len;i++){
			temp="0"+temp;
		}
		return temp;
	}

}
