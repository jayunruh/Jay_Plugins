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
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.TiffDecoder;
import ij.measure.Calibration;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.ShortProcessor;
import jalgs.FrameInterface;
import jalgs.jdataio;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class Tiff_Writer{
	// many parts of this code have been adapted from the ImageJ TiffEncoder
	double[] displayRanges;
	byte[][] channelLuts;
	int bitsPerSample,nEntries,samplesPerPixel,photoInterp,imageSize;
	int nSliceLabels,nMetaDataEntries,extraMetaDataEntries,ifdSize,imageOffset;
	int nMetaDataTypes,scaleSize,metaDataSize,width,height,stackframes;
	private boolean littleEndian=true;
	private byte buffer[]=new byte[8];
	byte[] description;
	public BufferedOutputStream outg;
	public long nextIFDg;

	private ImagePlus imp;
	private ImageProcessor ip;
	private Calibration cal;
	private FrameInterface finterface;

	long stackSize;
	
	public Tiff_Writer(){} //this is only if saving as single plane tif

	public Tiff_Writer(ImagePlus imp,int nframes,FrameInterface finterface){
		this.imp=imp;
		this.finterface=finterface;
		this.ip=imp.getProcessor();
		this.width=imp.getWidth();
		this.height=imp.getHeight();
		this.stackframes=nframes;
		this.cal=imp.getCalibration();
	}
	
	/**********************
	 * Use this version if you want to change the width and height
	 * @param imp
	 * @param width
	 * @param height
	 * @param nframes
	 * @param finterface
	 */
	public Tiff_Writer(ImagePlus imp,int width,int height,int nframes,FrameInterface finterface){
		this.imp=imp;
		this.finterface=finterface;
		this.ip=imp.getProcessor();
		this.width=width;
		this.height=height;
		this.stackframes=nframes;
		this.cal=imp.getCalibration();
	}
	
	public boolean saveAsTiff(String path,ImagePlus imp){
		this.imp=imp;
		ip=imp.getProcessor();
		width=imp.getWidth();
		height=imp.getHeight();
		stackframes=1;
		cal=imp.getCalibration();
		saveSizes();
		try{
			DataOutputStream out=new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)));
			byte[] hdr={73,73,42,0,8,0,0,0};
			// write the header
			out.write(hdr);
			long nextIFD=0L;
			writeIFD(out,imageOffset,(int)nextIFD);
			if(ip instanceof ColorProcessor)
				writeBitsPerPixel(out);
			writeDescription(out);
			if(scaleSize>0)
				writeScale(out);
			if(metaDataSize>0)
				writeMetaData(out);
			// write the image data
			jdataio jdio=new jdataio();
			if(ip instanceof ByteProcessor){
				jdio.writebytearray(out,(byte[])imp.getProcessor().getPixels());
			}else{
				if(ip instanceof ShortProcessor){
					jdio.writeintelshortarray(out,(short[])imp.getProcessor().getPixels());
				}else{
					if(ip instanceof FloatProcessor){
						jdio.writeintelfloatarray(out,(float[])imp.getProcessor().getPixels());
					}else{
						jdio.writeintelintarray(out,(int[])imp.getProcessor().getPixels());
					}
				}
			}
			out.close();
		}catch(IOException e){
			IJ.showMessage(e.getMessage());
			return false;
		}
		return true;
	}

	public boolean saveAsTiffStack(String path){
		if(imp.isComposite()){
			saveDisplayRangesAndLuts();
		}
		saveSizes();
		try{
			BufferedOutputStream out=new BufferedOutputStream(new FileOutputStream(path));
			byte[] hdr={73,73,42,0,8,0,0,0};
			// write the header
			out.write(hdr);
			long nextIFD=imageOffset+stackframes;
			if(nextIFD+stackframes*ifdSize>0xffffffffL)
				nextIFD=0L;
			writeIFD(out,imageOffset,(int)nextIFD);
			if(ip instanceof ColorProcessor)
				writeBitsPerPixel(out);
			writeDescription(out);
			if(scaleSize>0)
				writeScale(out);
			if(metaDataSize>0)
				writeMetaData(out);
			// write the image data
			jdataio jdio=new jdataio();
			for(int i=0;i<stackframes;i++){
				if(ip instanceof ByteProcessor){
					jdio.writebytearray(out,(byte[])finterface.getNextFrame());
				}else{
					if(ip instanceof ShortProcessor){
						jdio.writeintelshortarray(out,(short[])finterface.getNextFrame());
					}else{
						if(ip instanceof FloatProcessor){
							jdio.writeintelfloatarray(out,(float[])finterface.getNextFrame());
						}else{
							jdio.writeintelintarray(out,(int[])finterface.getNextFrame());
						}
					}
				}
			}
			if(nextIFD>0L){
				int ifdSize2=ifdSize;
				if(metaDataSize>0){
					metaDataSize=0;
					nEntries-=2;
					ifdSize2-=2*12;
				}
				for(int i=2;i<=stackframes;i++){
					if(i==stackframes)
						nextIFD=0;
					else
						nextIFD+=ifdSize2;
					imageOffset+=imageSize;
					writeIFD(out,imageOffset,(int)nextIFD);
				}
			}
			out.close();
		}catch(IOException e){
			IJ.showMessage(e.getMessage());
			return false;
		}
		return true;
	}
	
	/************
	 * this version uses the saveFrame method to write each frame instead of the frame interface
	 * @param path
	 * @return
	 */
	public boolean saveAsTiffStack2(String path){
		if(imp.isComposite()){
			saveDisplayRangesAndLuts();
		}
		saveSizes();
		try{
			outg=new BufferedOutputStream(new FileOutputStream(path));
			byte[] hdr={73,73,42,0,8,0,0,0};
			// write the header
			outg.write(hdr);
			nextIFDg=imageOffset+stackframes;
			if(nextIFDg+stackframes*ifdSize>0xffffffffL)
				nextIFDg=0L;
			writeIFD(outg,imageOffset,(int)nextIFDg);
			if(ip instanceof ColorProcessor)
				writeBitsPerPixel(outg);
			writeDescription(outg);
			if(scaleSize>0)
				writeScale(outg);
			if(metaDataSize>0)
				writeMetaData(outg);
			// write the image data frame by frame using the saveFrame method
		}catch(IOException e){
			IJ.showMessage(e.getMessage());
			return false;
		}
		return true;
	}
	
	public void saveFrame(Object pix){
		jdataio jdio=new jdataio();
		if(pix instanceof byte[]){
			jdio.writebytearray(outg,(byte[])pix);
		}else{
			if(pix instanceof short[]){
				jdio.writeintelshortarray(outg,(short[])pix);
			}else{
				if(pix instanceof float[]){
					jdio.writeintelfloatarray(outg,(float[])pix);
				}else{
					jdio.writeintelintarray(outg,(int[])pix);
				}
			}
		}
	}
	
	public boolean endWrite(){
		try{
    		if(nextIFDg>0L){
    			int ifdSize2=ifdSize;
    			if(metaDataSize>0){
    				metaDataSize=0;
    				nEntries-=2;
    				ifdSize2-=2*12;
    			}
    			for(int i=2;i<=stackframes;i++){
    				if(i==stackframes)
    					nextIFDg=0;
    				else
    					nextIFDg+=ifdSize2;
    				imageOffset+=imageSize;
    				writeIFD(outg,imageOffset,(int)nextIFDg);
    			}
    		}
    		outg.close();
		}catch(IOException e){
			IJ.showMessage(e.getMessage());
			return false;
		}
		return true;
	}

	void saveSizes(){
		bitsPerSample=8;
		samplesPerPixel=1;
		nEntries=9;
		photoInterp=1;
		int bytesPerPixel=1;
		int bpsSize=0;
		int colorMapSize=0;
		if(ip instanceof ShortProcessor){
			bitsPerSample=16;
			bytesPerPixel=2;
		}else{
			if(ip instanceof FloatProcessor){
				bitsPerSample=32;
				bytesPerPixel=4;
				nEntries++;
			}else{
				if(ip instanceof ColorProcessor){
					photoInterp=2;
					samplesPerPixel=3;
					bytesPerPixel=3;
					bpsSize=6;
				}
			}
		}
		scaleSize=0;
		if(cal.pixelWidth!=0){
			nEntries+=3;
			scaleSize=16;
		}
		makeDescriptionString();
		nEntries++; // for description entry
		imageSize=width*height*bytesPerPixel;
		stackSize=(long)imageSize*stackframes;
		metaDataSize=getMetaDataSize();
		if(metaDataSize>0)
			nEntries+=2;
		ifdSize=2+nEntries*12+4;
		int descriptionSize=description.length;
		imageOffset=8+ifdSize+bpsSize+descriptionSize+scaleSize+colorMapSize+nMetaDataEntries*4+metaDataSize;
	}

	int getMetaDataSize(){
		// if (stackSize+IMAGE_START>0xffffffffL) return 0;
		nSliceLabels=0;
		int size=0;
		int nTypes=0;
		if(displayRanges!=null){
			nMetaDataEntries++;
			size+=displayRanges.length*8;
			nTypes++;
		}
		if(channelLuts!=null){
			for(int i=0;i<channelLuts.length;i++){
				if(channelLuts[i]!=null)
					size+=channelLuts[i].length;
			}
			nTypes++;
			nMetaDataEntries+=channelLuts.length;
		}
		if(nMetaDataEntries>0)
			nMetaDataEntries++; // add entry for header
		int hdrSize=4+nTypes*8;
		if(size>0)
			size+=hdrSize;
		nMetaDataTypes=nTypes;
		return size;
	}

	void saveDisplayRangesAndLuts(){
		CompositeImage ci=(CompositeImage)imp;
		int channels=imp.getNChannels();
		displayRanges=new double[channels*2];
		for(int i=1;i<=channels;i++){
			LUT lut=ci.getChannelLut(i);
			displayRanges[(i-1)*2]=lut.min;
			displayRanges[(i-1)*2+1]=lut.max;
		}
		if(ci.hasCustomLuts()){
			channelLuts=new byte[channels][];
			for(int i=0;i<channels;i++){
				LUT lut=ci.getChannelLut(i+1);
				byte[] bytes=lut.getBytes();
				if(bytes==null){
					channelLuts=null;
					break;
				}
				channelLuts[i]=bytes;
			}
		}
	}

	/** Returns a string containing information about the specified image. */
	public String getDescriptionString(){
		Calibration cal=imp.getCalibration();
		StringBuffer sb=new StringBuffer(100);
		sb.append("ImageJ="+ImageJ.VERSION+"\n");
		sb.append("images="+stackframes+"\n");
		int channels=imp.getNChannels();
		if(channels>1)
			sb.append("channels="+channels+"\n");
		int slices=imp.getNSlices();
		if(slices>1)
			sb.append("slices="+slices+"\n");
		int frames=stackframes/(channels*slices);
		sb.append("frames="+frames+"\n");
		if(imp.isHyperStack()){
			sb.append("hyperstack=true\n");
		}else{
			if(slices>1||channels>1)
				sb.append("hyperstack=true\n");
		}
		if(imp.isComposite()){
			String mode=((CompositeImage)imp).getModeAsString();
			sb.append("mode="+mode+"\n");
		}
		if(cal.pixelWidth>0)
			sb.append("unit="+(cal.getUnit().equals("\u00B5m")?"um":cal.getUnit())+"\n");

		// get min and max display values
		double min=ip.getMin();
		double max=ip.getMax();
		int type=imp.getType();
		boolean enhancedLut=(type==ImagePlus.GRAY8||type==ImagePlus.COLOR_256)&&(min!=0.0||max!=255.0);
		if(enhancedLut||type==ImagePlus.GRAY16||type==ImagePlus.GRAY32){
			sb.append("min="+min+"\n");
			sb.append("max="+max+"\n");
		}

		// get non-zero origins
		if(cal.xOrigin!=0.0)
			sb.append("xorigin="+cal.xOrigin+"\n");
		if(cal.yOrigin!=0.0)
			sb.append("yorigin="+cal.yOrigin+"\n");
		if(cal.zOrigin!=0.0)
			sb.append("zorigin="+cal.zOrigin+"\n");
		if(cal.info!=null&&cal.info.length()<=64&&cal.info.indexOf('=')==-1&&cal.info.indexOf('\n')==-1)
			sb.append("info="+cal.info+"\n");
		sb.append((char)0);
		return new String(sb);
	}

	/**
	 * Creates an optional image description string for saving calibration data.
	 * For stacks, also saves the stack size so ImageJ can open the stack
	 * without decoding an IFD for each slice.
	 */
	void makeDescriptionString(){
		String sdescription=getDescriptionString();
		if(sdescription.charAt(sdescription.length()-1)!=(char)0)
			sdescription+=" ";
		description=sdescription.getBytes();
		description[description.length-1]=(byte)0;
	}

	/** Writes one 12-byte IFD entry. */
	void writeEntry(OutputStream out,int tag,int fieldType,int count,int value) throws IOException{
		writeShort(out,tag);
		writeShort(out,fieldType);
		writeInt(out,count);
		if(count==1&&fieldType==3){
			writeShort(out,value);
			writeShort(out,0);
		}else
			writeInt(out,value); // may be an offset
	}

	/** Writes one IFD (Image File Directory). */
	void writeIFD(OutputStream out,int imageOffset,int nextIFD) throws IOException{
		int tagDataOffset=8+ifdSize;
		writeShort(out,nEntries);
		writeEntry(out,TiffDecoder.NEW_SUBFILE_TYPE,4,1,0);
		writeEntry(out,TiffDecoder.IMAGE_WIDTH,4,1,width);
		writeEntry(out,TiffDecoder.IMAGE_LENGTH,4,1,height);
		if(ip instanceof ColorProcessor){
			writeEntry(out,TiffDecoder.BITS_PER_SAMPLE,3,3,tagDataOffset);
			tagDataOffset+=6;
		}else
			writeEntry(out,TiffDecoder.BITS_PER_SAMPLE,3,1,bitsPerSample);
		writeEntry(out,TiffDecoder.PHOTO_INTERP,3,1,photoInterp);
		writeEntry(out,TiffDecoder.IMAGE_DESCRIPTION,2,description.length,tagDataOffset);
		tagDataOffset+=description.length;
		writeEntry(out,TiffDecoder.STRIP_OFFSETS,4,1,imageOffset);
		writeEntry(out,TiffDecoder.SAMPLES_PER_PIXEL,3,1,samplesPerPixel);
		writeEntry(out,TiffDecoder.ROWS_PER_STRIP,3,1,height);
		writeEntry(out,TiffDecoder.STRIP_BYTE_COUNT,4,1,imageSize);
		if(cal.pixelWidth!=0){
			writeEntry(out,TiffDecoder.X_RESOLUTION,5,1,tagDataOffset);
			writeEntry(out,TiffDecoder.Y_RESOLUTION,5,1,tagDataOffset+8);
			tagDataOffset+=16;
			int unit=1;
			if(cal.getUnit().equals("inch"))
				unit=2;
			else if(cal.getUnit().equals("cm"))
				unit=3;
			writeEntry(out,TiffDecoder.RESOLUTION_UNIT,3,1,unit);
		}
		if(ip instanceof FloatProcessor){
			int format=3;
			writeEntry(out,TiffDecoder.SAMPLE_FORMAT,3,1,format);
		}
		if(metaDataSize>0){
			writeEntry(out,TiffDecoder.META_DATA_BYTE_COUNTS,4,nMetaDataEntries,tagDataOffset);
			writeEntry(out,TiffDecoder.META_DATA,1,metaDataSize,tagDataOffset+4*nMetaDataEntries);
			tagDataOffset+=nMetaDataEntries*4+metaDataSize;
		}
		writeInt(out,nextIFD);
	}

	/** Writes the 6 bytes of data required by RGB BitsPerSample tag. */
	void writeBitsPerPixel(OutputStream out) throws IOException{
		int bitsPerPixel=8;
		writeShort(out,bitsPerPixel);
		writeShort(out,bitsPerPixel);
		writeShort(out,bitsPerPixel);
	}

	/** Writes the variable length ImageDescription string. */
	void writeDescription(OutputStream out) throws IOException{
		out.write(description,0,description.length);
	}

	/**
	 * Writes the 16 bytes of data required by the XResolution and YResolution
	 * tags.
	 */
	void writeScale(OutputStream out) throws IOException{
		double xscale=1.0/cal.pixelWidth;
		double yscale=xscale;
		if(cal.pixelHeight>0&&cal.pixelHeight!=cal.pixelWidth){
			yscale=1.0/cal.pixelHeight;
		}
		double scale=1000000.0;
		if(xscale>1000.0)
			scale=1000.0;
		writeInt(out,(int)(xscale*scale));
		writeInt(out,(int)scale);
		writeInt(out,(int)(yscale*scale));
		writeInt(out,(int)scale);
	}

	/**
	 * Writes image metadata ("info" image propery, stack slice labels, channel
	 * display ranges, luts, ROIs, overlays and extra metadata).
	 */
	void writeMetaData(OutputStream out) throws IOException{

		// write byte counts (META_DATA_BYTE_COUNTS tag)
		writeInt(out,4+nMetaDataTypes*8); // header size
		String[] slicelabels=imp.getStack().getSliceLabels();
		for(int i=0;i<nSliceLabels;i++){
			if(slicelabels[i]==null)
				writeInt(out,2);
			else
				writeInt(out,slicelabels[i].length()*2);
		}
		if(displayRanges!=null)
			writeInt(out,displayRanges.length*8);
		if(channelLuts!=null){
			for(int i=0;i<channelLuts.length;i++)
				writeInt(out,channelLuts[i].length);
		}

		// write header (META_DATA tag header)
		writeInt(out,0x494a494a); // "IJIJ"
		if(nSliceLabels>0){
			writeInt(out,0x6c61626c); // type="labl"
			writeInt(out,nSliceLabels); // count
		}
		if(displayRanges!=null){
			writeInt(out,0x72616e67); // type="rang"
			writeInt(out,1); // count
		}
		if(channelLuts!=null){
			writeInt(out,0x6c757473); // type="luts"
			writeInt(out,channelLuts.length); // count
		}

		// write data (META_DATA tag body)
		for(int i=0;i<nSliceLabels;i++){
			if(slicelabels[i]!=null){
				writeChars(out,slicelabels[i]);
			}else{
				writeChars(out," ");
			}
		}
		if(displayRanges!=null){
			for(int i=0;i<displayRanges.length;i++)
				writeDouble(out,displayRanges[i]);
		}
		if(channelLuts!=null){
			for(int i=0;i<channelLuts.length;i++)
				out.write(channelLuts[i]);
		}
	}

	final void writeShort(OutputStream out,int v) throws IOException{
		if(littleEndian){
			out.write(v&255);
			out.write((v>>>8)&255);
		}else{
			out.write((v>>>8)&255);
			out.write(v&255);
		}
	}

	final void writeInt(OutputStream out,int v) throws IOException{
		if(littleEndian){
			out.write(v&255);
			out.write((v>>>8)&255);
			out.write((v>>>16)&255);
			out.write((v>>>24)&255);
		}else{
			out.write((v>>>24)&255);
			out.write((v>>>16)&255);
			out.write((v>>>8)&255);
			out.write(v&255);
		}
	}

	final void writeLong(OutputStream out,long v) throws IOException{
		if(littleEndian){
			buffer[7]=(byte)(v>>>56);
			buffer[6]=(byte)(v>>>48);
			buffer[5]=(byte)(v>>>40);
			buffer[4]=(byte)(v>>>32);
			buffer[3]=(byte)(v>>>24);
			buffer[2]=(byte)(v>>>16);
			buffer[1]=(byte)(v>>>8);
			buffer[0]=(byte)v;
			out.write(buffer,0,8);
		}else{
			buffer[0]=(byte)(v>>>56);
			buffer[1]=(byte)(v>>>48);
			buffer[2]=(byte)(v>>>40);
			buffer[3]=(byte)(v>>>32);
			buffer[4]=(byte)(v>>>24);
			buffer[5]=(byte)(v>>>16);
			buffer[6]=(byte)(v>>>8);
			buffer[7]=(byte)v;
			out.write(buffer,0,8);
		}
	}

	final void writeDouble(OutputStream out,double v) throws IOException{
		writeLong(out,Double.doubleToLongBits(v));
	}

	final void writeChars(OutputStream out,String s) throws IOException{
		int len=s.length();
		if(littleEndian){
			for(int i=0;i<len;i++){
				int v=s.charAt(i);
				out.write(v&255);
				out.write((v>>>8)&255);
			}
		}else{
			for(int i=0;i<len;i++){
				int v=s.charAt(i);
				out.write((v>>>8)&255);
				out.write(v&255);
			}
		}
	}

}
