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
//import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.process.ImageProcessor;
import ij.process.LUT;
import jalgs.algutils;
import jalgs.jdataio;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;




//import ome.units.UNITS;
//import ome.units.quantity.Length;
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

public class LOCI_file_reader{
	// this plugin simply uses the loci library to open files
	public int nseries;
	public boolean nometa=false;
	
	public static void main(String[] args){
		//the first argument should be the incoming file name
		File f=new File(args[0]);
		//System.out.println(args[0]);
		LOCI_file_reader lfr=new LOCI_file_reader();
		System.out.println(f.getParent());
		System.out.println(f.getName());
		ImagePlus imp=lfr.get_loci_imp_simple(f.getParent()+File.separator,f.getName(),0);
		if(imp==null){
			System.out.println("image did not load");
			return;
		}
		FileSaver fs=new FileSaver(imp);
		fs.saveAsTiffStack(args[1]);
		imp.close();
		System.out.println(args[0]+"=>"+args[1]+" complete");
	}

	public ImagePlus get_loci_imp(String directory,String fname){
		return get_loci_imp(directory,fname,false);
	}
	
	public ImagePlus get_loci_imp(String path){
		return get_loci_imp(path,false,0,false,"");
	}

	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta){
		return get_loci_imp(directory,fname,outmeta,0);
	}

	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta,int series){
		return get_loci_imp(directory,fname,outmeta,series,false,"");
	}
	
	public ImagePlus get_loci_imp_simple(String directory,String fname,int series){
		return get_loci_imp_simple(directory+fname,series);
	}
	
	public ImagePlus get_loci_imp_simple(String path,int series){
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			String fname=(new File(path)).getName();
			r.setId(path);
			nseries=r.getSeriesCount();
			if(series>=nseries)
				series=0;
			r.setSeries(series);
			int num=r.getImageCount(); //IJ.log("n images = "+num);
			int width=r.getSizeX();
			int height=r.getSizeY();
			String name=""+fname+series;
			ImageStack stack=new ImageStack(width,height);
			for(int i=0;i<num;i++){
				stack.addSlice(r.openProcessors(i)[0]);
				IJ.showProgress(i,num);
			}
			r.close();
			// shuffle(stack,channels,slices,frames,order);
			ImagePlus imp=new ImagePlus(name,stack);
			return imp;
		}catch(FormatException e){
			jdataio jdio=new jdataio();
			System.out.println(jdio.getExceptionTrace(e));
			return null;
		}catch(IOException e){
			jdataio jdio=new jdataio();
			System.out.println(jdio.getExceptionTrace(e));
			return null;
		}
	}
	
	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta,int series,boolean proj,String projstat,int refchan){
		return get_loci_imp(directory+fname,outmeta,series,proj,projstat,refchan);
	}
	
	public ImagePlus get_loci_imp(String path,boolean outmeta,int series,boolean proj,String projstat,int refchan){
		return get_loci_subimp(path,outmeta,series,proj,projstat,refchan,null);
	}
	
	/************
	 * this version collects a subset of channels, slices, and frames
	 * @param directory
	 * @param fname
	 * @param outmeta
	 * @param series
	 * @param proj
	 * @param projstat
	 * @param refchan
	 * @param limits: an int array with minch,maxch,minz,maxz,mint,maxt (0 based)
	 * @return
	 */
	public ImagePlus get_loci_subimp(String directory,String fname,boolean outmeta,int series,boolean proj,String projstat,int refchan,int[] limits){
		return get_loci_subimp(directory+fname,outmeta,series,proj,projstat,refchan,limits);
	}
	
	/************
	 * this version collects a subset of channels, slices, and frames
	 * @param directory
	 * @param fname
	 * @param outmeta
	 * @param series
	 * @param proj
	 * @param projstat
	 * @param refchan
	 * @param limits: an int array with minch,maxch,minz,maxz,mint,maxt (0 based)
	 * @return
	 */
	public ImagePlus get_loci_subimp(String path,boolean outmeta,int series,boolean proj,String projstat,int refchan,int[] limits){
		IMetadata omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			String fname=(new File(path)).getName();
			r.setId(path);
			nseries=r.getSeriesCount();
			if(series>=nseries)
				series=0;
			r.setSeries(series);
			int num=r.getImageCount(); //IJ.log("n images = "+num);
			int width=r.getSizeX();
			int height=r.getSizeY();
			int channels=r.getSizeC();
			int slices=r.getSizeZ();
			int frames=r.getSizeT();
			String order=r.getDimensionOrder();
			int[] chsizes=null;
			if(refchan>=0){
				//here the slices is actually max slices
				//int totslices=channels*slices*frames;
				//int framesize=slices*(channels-1)+1;
				chsizes=new int[channels];
				for(int i=0;i<channels;i++) chsizes[i]=slices;
				chsizes[refchan]=1;
				//float newframes=(float)num/(float)framesize;
				//frames=(int)(newframes+0.0001f);
			}

			if(outmeta&&!nometa){
				Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
				if(globalMeta!=null){
					IJ.log("Global Metadata");
					dumpMetaData(globalMeta);
				}
				Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
				if(seriesMeta!=null){
					IJ.log("Series Metadata");
					dumpMetaData(seriesMeta);
				}
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
				//Length temp=new Length(1.0,UNITS.MICROM);
				if(omexmlMetadata.getPixelsPhysicalSizeX(series)!=null)
					psize=omexmlMetadata.getPixelsPhysicalSizeX(series).value().floatValue();
				if(omexmlMetadata.getPixelsPhysicalSizeZ(series)!=null)
					zsize=omexmlMetadata.getPixelsPhysicalSizeZ(series).value().floatValue();
				if(omexmlMetadata.getPixelsTimeIncrement(series)!=null)
					tsize=omexmlMetadata.getPixelsTimeIncrement(series).value().floatValue();
			}
			if(limits==null){
				limits=new int[]{-1,-1,-1,-1,-1,-1};
			}
			if(limits[0]<0) limits[0]=0;
			if(limits[1]<0 || limits[1]>=channels) limits[1]=channels-1;
			if(limits[2]<0) limits[2]=0;
			if(limits[3]<0 || limits[3]>=slices) limits[3]=slices-1;
			if(limits[4]<0) limits[4]=0;
			if(limits[5]<0 || limits[5]>=frames) limits[5]=frames-1;
			int tempslices=limits[3]-limits[2]+1;
			int tempframes=limits[5]-limits[4]+1;
			int tempchannels=limits[1]-limits[0]+1;
			Object[] luts=new Object[channels];
			ImageStack stack=new ImageStack(width,height);
			if(proj){
				int counter=0;
				for(int i=limits[4];i<=limits[5];i++){
					for(int j=limits[0];j<=limits[1];j++){
						Object[] tempzstack=new Object[tempslices];
						for(int k=limits[2];k<=limits[3];k++){
							int index=get_stack_index(j,k,i,channels,slices,frames,order,chsizes);
							ImageProcessor ip=r.openProcessors(index)[0];
							tempzstack[k-limits[2]]=ip.getPixels();
						}
						Object projslice=algutils.get_stack_proj_stat(projstat,tempzstack,width,height,tempslices,null);
						stack.addSlice("",projslice);
						if(luts[j]==null){
							luts[j]=r.get8BitLookupTable();
							if(luts[j]==null) luts[j]=r.get16BitLookupTable();
						}
						IJ.showProgress(counter,tempframes*tempchannels);
						counter++;
					}
				}
				tempslices=1;
			}else{
				int counter=0;
				for(int i=limits[4];i<=limits[5];i++){
					for(int j=limits[2];j<=limits[3];j++){
						for(int k=limits[0];k<=limits[1];k++){
							int index=get_stack_index(k,j,i,channels,slices,frames,order,chsizes);
							//IJ.log(""+k+"\t "+j+"\t "+i+"\t "+index);
							ImageProcessor ip=r.openProcessors(index)[0];
							stack.addSlice(ip);
							if(luts[k]==null){
								luts[k]=r.get8BitLookupTable();
								if(luts[k]==null) luts[k]=r.get16BitLookupTable();
							}
							IJ.showProgress(counter,tempframes*tempslices*tempchannels);
							counter++;
						}
					}
				}
			}
			r.close();
			LUT[] luts2=new LUT[channels];
			for(int i=0;i<channels;i++){
				if(luts[i]==null){
					luts2=null;
					break;
				}
				//here assuming luts is a byte or a short
				if(luts[i] instanceof byte[][]){
					luts2[i]=new LUT(((byte[][])luts[i])[0],((byte[][])luts[i])[1],((byte[][])luts[i])[2]);
				} else {
					byte[][] temp=new byte[3][];
					temp[0]=algutils.convert_arr_byte(((Object[])luts[i])[0]);
					temp[1]=algutils.convert_arr_byte(((Object[])luts[i])[1]);
					temp[2]=algutils.convert_arr_byte(((Object[])luts[i])[2]);
					luts2[i]=new LUT(temp[0],temp[1],temp[2]);
				}
			}
			// shuffle(stack,channels,slices,frames,order);
			ImagePlus imp=new ImagePlus(name,stack);
			if(!nometa){
				jutils.set_psize(imp,psize);
				jutils.set_pdepth(imp,zsize);
				jutils.set_pinterval(imp,tsize);
			}
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(tempchannels,tempslices,tempframes);
			if(tempchannels>1){
				CompositeImage ci=new CompositeImage(imp,CompositeImage.COLOR);
				if(luts2!=null) ci.setLuts(luts2);
				return ci;
			}else{
				return imp;
			}
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}

	public ImagePlus get_loci_imp(String directory,String fname,boolean outmeta,int series,boolean proj,String projstat){
		return get_loci_imp(directory,fname,outmeta,series,proj,projstat,-1);
	}
	
	public ImagePlus get_loci_imp(String path,boolean outmeta,int series,boolean proj,String projstat){
		return get_loci_imp(path,outmeta,series,proj,projstat,-1);
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
			//omexmlMetadata.getObjectiveImmersion(0,0);
			if(!nometa){
				if(omexmlMetadata.getPixelsPhysicalSizeX(series)!=null)
					psize=omexmlMetadata.getPixelsPhysicalSizeX(series).value().floatValue();
				if(omexmlMetadata.getPixelsPhysicalSizeZ(series)!=null)
					zsize=omexmlMetadata.getPixelsPhysicalSizeZ(series).value().floatValue();
				if(omexmlMetadata.getPixelsTimeIncrement(series)!=null)
					tsize=omexmlMetadata.getPixelsTimeIncrement(series).value().floatValue();
			}
			ImageStack stack=new ImageStack(newwidth,newheight);
			int counter=0;
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int y=0;y<ytiles;y++){
						for(int x=0;x<xtiles;x++){
							for(int k=0;k<channels;k++){
								int index=get_stack_index(k,j,i,channels,slices,frames,order,null);
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
	
	public int[] get_imp_sizes(String directory,String fname,int series){
		IMetadata omexmlMetadata=null;
		//if(!nometa) omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			//if(!nometa) r.setMetadataStore(omexmlMetadata);
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
			r.close();
			return new int[]{width,height,channels,slices,frames};
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}
	
	public Hashtable<String,Object> getGlobalMetaData(String directory,String fname){
		IMetadata omexmlMetadata=MetadataTools.createOMEXMLMetadata();
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
	
	public String[][] getOMEXMLObjectiveInfo(String directory,String fname,int series){
		//here we get selected info about the objective from the omexmlmetadata
		IMetadata omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			if(series>=nseries) series=0;
			r.setSeries(series);
			String[][] objinfo=new String[7][2];
			objinfo[0][0]="Correction";
			objinfo[0][1]=omexmlMetadata.getObjectiveCorrection(0,0).toString();
			objinfo[1][0]="Immersion";
			objinfo[1][1]=omexmlMetadata.getObjectiveImmersion(0,0).toString();
			objinfo[2][0]="LensNA";
			objinfo[2][1]=omexmlMetadata.getObjectiveLensNA(0,0).toString();
			objinfo[3][0]="NominalMagnification";
			objinfo[3][1]=omexmlMetadata.getObjectiveNominalMagnification(0,0).toString();
			objinfo[4][0]="Manufacturer";
			objinfo[4][1]=omexmlMetadata.getObjectiveManufacturer(0,0).toString();
			objinfo[5][0]="Model";
			objinfo[5][1]=omexmlMetadata.getObjectiveModel(0,0).toString();
			objinfo[6][0]="WorkingDistance";
			objinfo[6][1]=omexmlMetadata.getObjectiveWorkingDistance(0,0).toString();
			r.close();
			return objinfo;
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}
	
	public Hashtable<String,Object> getSeriesMetaData(String directory,String fname,int series){
		IMetadata omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			if(series>=nseries) series=0;
			r.setSeries(series);
			//Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
			Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
			r.close();
			return seriesMeta;
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
	}
	
	public String get_metadata_value(String directory,String fname,int series,String key,boolean dump){
		IMetadata omexmlMetadata=MetadataTools.createOMEXMLMetadata();
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
	
	public String[] batch_get_series_metadata_value(String directory,String fname,String key){
		String[] keys={key};
		return batch_get_series_metadata_value(directory,fname,keys)[0];
		//return null;
	}
	
	public String[][] batch_get_series_metadata_value(String directory,String fname,String[] key){
		IMetadata omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			String[][] vals=new String[key.length][nseries];
			Hashtable<String,Object> globalMeta=r.getGlobalMetadata();
			for(int i=0;i<nseries;i++){
				r.setSeries(i);
				Hashtable<String,Object> seriesMeta=r.getSeriesMetadata();
				for(int j=0;j<key.length;j++){
					Object temp=seriesMeta.get(key[j]);
					if(temp==null) globalMeta.get(key[j]);
					if(temp!=null) vals[j][i]=temp.toString();
				}
			}
			r.close();
			return vals;
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
		//return null;
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
	
	public String pad_number(int num,int len){
		String temp=Integer.toString(num);
		if(temp.length()==len) return temp;
		if(temp.length()>len) return temp.substring(0,len);
		for(int i=temp.length();i<len;i++){
			temp="0"+temp;
		}
		return temp;
	}
	
	public boolean batchExportSeries(String directory,String fname,String outdir,boolean nometa){
		this.nometa=nometa;
		IMetadata omexmlMetadata=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			int maxdigits=Integer.toString(nseries).length();
			String[] names=new String[nseries];
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
			for(int s=0;s<nseries;s++){
				r.setSeries(s);
				int num=r.getImageCount();
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
				int counter=0;
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							int index=get_stack_index(k,j,i,channels,slices,frames,order,null);
							//IJ.log(""+k+"\t "+j+"\t "+i+"\t "+index);
							ImageProcessor ip=r.openProcessors(index)[0];
							stack.addSlice(ip);
							IJ.showProgress(counter,frames*slices*channels);
							counter++;
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
				FileSaver fs=new FileSaver(imp);
				fs.saveAsTiffStack(outdir+names[s]+".tif");
				imp.close();
				IJ.showStatus(names[s]+" exported");
				IJ.showProgress(s,nseries);
				if(IJ.escapePressed()) break;
			}
			r.close();
		}catch(FormatException e){
			return false;
		}catch(IOException e){
			return false;
		}
		return true;
	}
	
	public ImagePlus[] batchOpenSeries(String directory,String fname,boolean nometa){
		this.nometa=nometa;
		IMetadata omexmlMetadata=null;
		ImagePlus[] imps=null;
		if(!nometa)
			omexmlMetadata=MetadataTools.createOMEXMLMetadata();
		ImageProcessorReader r=new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		try{
			if(!nometa)
				r.setMetadataStore(omexmlMetadata);
			r.setId(directory+fname);
			nseries=r.getSeriesCount();
			imps=new ImagePlus[nseries];
			int maxdigits=Integer.toString(nseries).length();
			String[] names=new String[nseries];
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
			for(int s=0;s<nseries;s++){
				r.setSeries(s);
				int num=r.getImageCount();
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
				int counter=0;
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							int index=get_stack_index(k,j,i,channels,slices,frames,order,null);
							//IJ.log(""+k+"\t "+j+"\t "+i+"\t "+index);
							ImageProcessor ip=r.openProcessors(index)[0];
							stack.addSlice(ip);
							IJ.showProgress(counter,frames*slices*channels);
							counter++;
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
				imps[s]=imp;
				/*FileSaver fs=new FileSaver(imp);
				fs.saveAsTiffStack(outdir+names[s]+".tif");
				imp.close();*/
				IJ.showStatus(names[s]+" exported");
				IJ.showProgress(s,nseries);
				if(IJ.escapePressed()) break;
			}
			r.close();
		}catch(FormatException e){
			return null;
		}catch(IOException e){
			return null;
		}
		return imps;
	}

}
