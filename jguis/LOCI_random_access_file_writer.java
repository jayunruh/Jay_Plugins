/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.ImageJ;
import ij.ImagePlus;
import jalgs.jdataio;

import java.io.IOException;

import loci.common.DataTools;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.FormatTools;
import loci.formats.IFormatWriter;
import loci.formats.ImageWriter;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import ome.units.quantity.Length;
import ome.units.quantity.Time;
import ome.xml.model.enums.DimensionOrder;
import ome.xml.model.enums.EnumerationException;
import ome.xml.model.enums.PixelType;
import ome.xml.model.primitives.PositiveFloat;
import ome.xml.model.primitives.PositiveInteger;

public class LOCI_random_access_file_writer{
	// this plugin simply uses the loci library to write tif files by random access
	public int nseries,width,height,channels,num,slices,frames;
	public float psize,tsize,zsize;
	public String order,name;
	public boolean littleEndian;
	public IFormatWriter w;
	
	public LOCI_random_access_file_writer(String directory,String fname,int width,int height,int channels,int slices,int frames){
		this.width=width;
		this.height=height;
		this.channels=channels;
		this.slices=slices;
		this.frames=frames;
		String fname2=fname.substring(0);
		if(!fname2.endsWith(".tif")) fname2+=".tif";
		//assume the data is 32 bit
		int ptype=FormatTools.FLOAT;
		try{
			w=new ImageWriter().getWriter(directory+fname);
		}catch(FormatException e){
			return;
		}
		w.setWriteSequentially(false);
		//FileInfo fi=imp.getOriginalFileInfo();
		//String xml = fi == null ? null : fi.description == null ? null :
	    //    fi.description.indexOf("xml") == -1 ? null : fi.description;
		OMEXMLService service=null;
		IMetadata store=null;
        try {
            ServiceFactory factory = new ServiceFactory();
            service = factory.getInstance(OMEXMLService.class);
            store = service.createOMEXMLMetadata(null);
        }
        catch (DependencyException de) { }
        catch (ServiceException se) { }
        store.createRoot();
        store.setPixelsSizeX(new PositiveInteger(width),0);
        store.setPixelsSizeY(new PositiveInteger(height),0);
        store.setPixelsSizeZ(new PositiveInteger(slices),0);
        store.setPixelsSizeC(new PositiveInteger(channels),0);
        store.setPixelsSizeT(new PositiveInteger(frames),0);
        if(store.getImageID(0)==null){
        	store.setImageID(MetadataTools.createLSID("Image",0),0);
        }
        if(store.getPixelsID(0)==null){
        	store.setPixelsID(MetadataTools.createLSID("Pixels",0),0);
        }
     // always reset the pixel type
        // this prevents problems if the user changed the bit depth of the image
        try {
          store.setPixelsType(PixelType.fromString(
            FormatTools.getPixelTypeString(ptype)), 0);
        }
        catch (EnumerationException e) { }

        if (store.getPixelsBinDataCount(0) == 0 ||
          store.getPixelsBinDataBigEndian(0, 0) == null)
        {
          store.setPixelsBinDataBigEndian(Boolean.FALSE, 0, 0);
        }
        String ORDER="XYCZT";
        if (store.getPixelsDimensionOrder(0) == null) {
          try {
            store.setPixelsDimensionOrder(DimensionOrder.fromString(ORDER), 0);
          }
          catch (EnumerationException e) { }
        }

        for (int c=0; c<channels; c++) {
          if (c >= store.getChannelCount(0) || store.getChannelID(0, c) == null) {
            String lsid = MetadataTools.createLSID("Channel", 0, c);
            store.setChannelID(lsid, 0, c);
          }
          store.setChannelSamplesPerPixel(new PositiveInteger(channels), 0, 0);
        }

       // Calibration cal = imp.getCalibration();

        store.setPixelsPhysicalSizeX(new Length(1.0,null), 0);
        store.setPixelsPhysicalSizeY(new Length(1.0,null), 0);
        store.setPixelsPhysicalSizeZ(new Length(1.0,null), 0);
        store.setPixelsTimeIncrement(new Time(1.0,null), 0);
        
        //Object info = imp.getProperty("Info");
        //if (info != null) {
          //String imageInfo = info.toString();
        	String imageInfo=getDescriptionString();
          if (imageInfo != null) {
            String[] lines = imageInfo.split("\n");
            for (String line : lines) {
              int eq = line.lastIndexOf("=");
              if (eq > 0) {
                String key = line.substring(0, eq).trim();
                String value = line.substring(eq + 1).trim();

                if (key.endsWith("BitsPerPixel")) {
                  w.setValidBitsPerPixel(Integer.parseInt(value));
                  break;
               }
              }
            }
          }
        //}
        
        w.setMetadataRetrieve(store);
        
        littleEndian =
            !w.getMetadataRetrieve().getPixelsBinDataBigEndian(0, 0).booleanValue();
        w.setInterleaved(false);
	}
	
	/** Returns a string containing information about the specified image. */
	public String getDescriptionString(){
		//Calibration cal=imp.getCalibration();
		StringBuffer sb=new StringBuffer(100);
		sb.append("ImageJ="+ImageJ.VERSION+"\n");
		sb.append("images="+channels*slices*frames+"\n");
		//int channels=imp.getNChannels();
		if(channels>1)
			sb.append("channels="+channels+"\n");
		//int slices=imp.getNSlices();
		if(slices>1)
			sb.append("slices="+slices+"\n");
		//int frames=stackframes/(channels*slices);
		sb.append("frames="+frames+"\n");
		//if(imp.isHyperStack()){
		//	sb.append("hyperstack=true\n");
		//}else{
			if(slices>1||channels>1)
				sb.append("hyperstack=true\n");
		//}
		//if(imp.isComposite()){
		if(channels>1){
			//String mode=((CompositeImage)imp).getModeAsString();
			String mode="COLOR";
			sb.append("mode="+mode+"\n");
		}
		//if(cal.pixelWidth>0)
		//	sb.append("unit="+(cal.getUnit().equals("\u00B5m")?"um":cal.getUnit())+"\n");

		// get min and max display values
		//double min=ip.getMin();
		//double max=ip.getMax();
		double min=0.0;
		double max=100.0;
		//int type=imp.getType();
		int type=ImagePlus.GRAY32;
		boolean enhancedLut=(type==ImagePlus.GRAY8||type==ImagePlus.COLOR_256)&&(min!=0.0||max!=255.0);
		if(enhancedLut||type==ImagePlus.GRAY16||type==ImagePlus.GRAY32){
			sb.append("min="+min+"\n");
			sb.append("max="+max+"\n");
		}

		// get non-zero origins
		/*if(cal.xOrigin!=0.0)
			sb.append("xorigin="+cal.xOrigin+"\n");
		if(cal.yOrigin!=0.0)
			sb.append("yorigin="+cal.yOrigin+"\n");
		if(cal.zOrigin!=0.0)
			sb.append("zorigin="+cal.zOrigin+"\n");
		if(cal.info!=null&&cal.info.length()<=64&&cal.info.indexOf('=')==-1&&cal.info.indexOf('\n')==-1)
			sb.append("info="+cal.info+"\n");*/
		sb.append((char)0);
		return new String(sb);
	}
	
	public void dispose(){
		try{
			w.close();
		}catch(IOException e){
			w=null;
		}
	}
	
	public boolean writeSubImage(float[] subimage,int x,int y,int rwidth,int rheight,int channel,int slice,int frame){
		byte[] buf=DataTools.floatsToBytes(subimage,littleEndian);
		try{
			w.saveBytes(frame*slices*channels+slice*channels+channel,buf,x,y,rwidth,rheight);
		}catch(FormatException e){
			return false;
		}catch(IOException e){
			return false;
		}
		return true;
	}
	
	public boolean writeLine(float[] line,int y,int channel,int slice,int frame){
		return writeSubImage(line,0,y,width,1,channel,slice,frame);
	}
	
	public boolean writePixel(float val,int x,int y,int channel,int slice,int frame){
		return writeSubImage(new float[]{val},x,y,1,1,channel,slice,frame);
	}
	
	public boolean writeCArray(float[] carray,int x,int y,int slice,int frame){
		for(int i=0;i<channels;i++){
			if(!writePixel(carray[i],x,y,i,slice,frame)){
				return false;
			}
		}
		return true;
	}
	
	public boolean writeCCarpet(float[][] carray,int y,int slice,int frame){
		for(int i=0;i<channels;i++){
			if(!writeLine(carray[i],y,i,slice,frame)) return false;
		}
		return true;
	}

}
