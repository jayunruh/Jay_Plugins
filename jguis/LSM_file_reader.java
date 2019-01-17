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
import ij.LookUpTable;
import ij.io.FileInfo;
import ij.io.FileSaver;
import ij.io.ImageReader;
import ij.io.RandomAccessStream;
import ij.measure.Calibration;

import java.awt.Color;
import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;

public class LSM_file_reader{
	// this code was released under the GNU:
	/*
	 * Copyright (C) 2002-2009 Patrick Pirrotte, Jerome Mutterer, Yannick Krempp
	 * 
	 * This program is free software; you can redistribute it and/or modify it
	 * under the terms of the GNU General Public License as published by the
	 * Free Software Foundation; either version 2 of the License, or (at your
	 * option) any later version.
	 * 
	 * This program is distributed in the hope that it will be useful, but
	 * WITHOUT ANY WARRANTY; without even the implied warranty of
	 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
	 * Public License for more details.
	 * 
	 * You should have received a copy of the GNU General Public License along
	 * with this program; if not, write to the Free Software Foundation, Inc.,
	 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
	 * 
	 * Modified by Jay Unruh
	 */
	public static char micro='\u00b5';

	public static String micrometer=micro+"m";
	
	public static void main(String[] args){
		//the first argument should be the incoming file name, the second the outgoing
		File f=new File(args[0]);
		LSM_file_reader lfr=new LSM_file_reader();
		ImagePlus imp=lfr.open(f);
		if(imp.getNFrames()>1){
			
		}
		FileSaver fs=new FileSaver(imp);
		fs.saveAsTiffStack(args[1]);
		imp.close();
		System.out.println(args[0]+"=>"+args[1]+" complete");
	}

	public ImagePlus open(String directory,String filename){
		return open(new File(directory,filename));
	}

	public ImagePlus open(File fimg){
		ImagePlus imp=null;
		RandomAccessFile file;
		LsmFileInfo lsm;
		try{
			file=new RandomAccessFile(fimg,"r");
			RandomAccessStream stream=new RandomAccessStream(file);
			lsm=new LsmFileInfo();
			lsm.fileName=fimg.getName();
			lsm.directory=fimg.getParent();
			if(isLSMfile(stream)){
				ImageDirectory imDir=readImageDirectoy(stream,8);
				// IJ.log("Time Interval = "+imDir.TIF_CZ_LSMINFO.TimeInterval);
				lsm.imageDirectories.add(imDir);
				while(imDir.OFFSET_NEXT_DIRECTORY!=0){
					imDir=readImageDirectoy(stream,imDir.OFFSET_NEXT_DIRECTORY);
					lsm.imageDirectories.add(imDir);
				}
				imp=open(stream,lsm);
				stream.close();
			}else
				return null;
		}catch(FileNotFoundException e){
			return null;
		}catch(IOException e){
			return null;
		}
		return imp;
	}

	public boolean isLSMfile(RandomAccessStream stream){
		boolean identifier=false;
		long ID=0;
		try{
			stream.seek(2);
			ID=swap(stream.readShort());
			if(ID==42){
				identifier=true;
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		return identifier;
	}

	private long getTagCount(RandomAccessStream stream,long position){
		long tags=0;
		try{
			stream.seek((int)position);
			tags=swap(stream.readShort());
		}catch(IOException e){
			e.printStackTrace();
		}
		return tags;
	}

	private ImageDirectory readImageDirectoy(RandomAccessStream stream,long startPosition){
		ImageDirectory imDir=new ImageDirectory();
		long tags=getTagCount(stream,startPosition);
		byte[] tag;
		int tagtype=0;
		int MASK=0x00ff;
		long MASK2=0x000000ff;
		long currentTagPosition=0;
		// needed because sometimes offset do not fit in
		// the imDir structure and are placed elsewhere
		long stripOffset=0,stripByteOffset=0;
		for(int i=0;i<tags;i++){
			currentTagPosition=startPosition+2+i*12;
			tag=readTag(stream,(int)currentTagPosition);
			tagtype=((tag[1]&MASK)<<8)|((tag[0]&MASK)<<0);

			switch(tagtype){
			case 254:
				imDir.TIF_NEWSUBFILETYPE=((tag[11]&MASK2)<<24)|((tag[10]&MASK2)<<16)|((tag[9]&MASK2)<<8)|(tag[8]&MASK2);
				break;
			case 256:
				imDir.TIF_IMAGEWIDTH=((tag[11]&MASK2)<<24)|((tag[10]&MASK2)<<16)|((tag[9]&MASK2)<<8)|(tag[8]&MASK2);
				break;
			case 257:
				imDir.TIF_IMAGELENGTH=((tag[11]&MASK2)<<24)|((tag[10]&MASK2)<<16)|((tag[9]&MASK2)<<8)|(tag[8]&MASK2);
				break;
			case 258:
				imDir.TIF_BITSPERSAMPLE_LENGTH=((tag[7]&MASK2)<<24)|((tag[6]&MASK2)<<16)|((tag[5]&MASK2)<<8)|(tag[4]&MASK2);
				imDir.TIF_BITSPERSAMPLE_CHANNEL[0]=((tag[8]&MASK2)<<0);
				imDir.TIF_BITSPERSAMPLE_CHANNEL[1]=((tag[9]&MASK2)<<0);
				imDir.TIF_BITSPERSAMPLE_CHANNEL[2]=((tag[10]&MASK2)<<0);
				break;
			case 259:
				imDir.TIF_COMPRESSION=((tag[8]&MASK2)<<0);
				break;
			case 262:
				imDir.TIF_PHOTOMETRICINTERPRETATION=((tag[8]&MASK2)<<0);
				break;
			case 273:
				imDir.TIF_STRIPOFFSETS_LENGTH=((tag[7]&MASK2)<<24)|((tag[6]&MASK2)<<16)|((tag[5]&MASK2)<<8)|(tag[4]&MASK2);
				stripOffset=((tag[11]&MASK2)<<24)|((tag[10]&MASK2)<<16)|((tag[9]&MASK2)<<8)|(tag[8]&MASK2);
				break;
			case 277:
				imDir.TIF_SAMPLESPERPIXEL=((tag[8]&MASK2)<<0);
				break;
			case 279:
				imDir.TIF_STRIPBYTECOUNTS_LENGTH=((tag[7]&MASK2)<<24)|((tag[6]&MASK2)<<16)|((tag[5]&MASK2)<<8)|(tag[4]&MASK2);
				stripByteOffset=((tag[11]&MASK2)<<24)|((tag[10]&MASK2)<<16)|((tag[9]&MASK2)<<8)|(tag[8]&MASK2);
				break;
			case 317:
				imDir.TIF_PREDICTOR=((tag[8]&MASK2)<<0);
				break;
			case 34412:
				imDir.TIF_CZ_LSMINFO_OFFSET=((tag[11]&MASK2)<<24)|((tag[10]&MASK2)<<16)|((tag[9]&MASK2)<<8)|(tag[8]&MASK2);
				break;
			default:
				break;
			}
		}
		imDir.TIF_STRIPOFFSETS=new long[(int)imDir.TIF_STRIPOFFSETS_LENGTH];
		if(imDir.TIF_STRIPOFFSETS_LENGTH==1){
			imDir.TIF_STRIPOFFSETS[0]=stripOffset;
		}else{
			imDir.TIF_STRIPOFFSETS=getIntTable(stream,stripOffset,(int)imDir.TIF_STRIPOFFSETS_LENGTH);
		}
		imDir.TIF_STRIPBYTECOUNTS=new long[(int)imDir.TIF_STRIPBYTECOUNTS_LENGTH];
		if(imDir.TIF_STRIPBYTECOUNTS_LENGTH==1){
			imDir.TIF_STRIPBYTECOUNTS[0]=stripByteOffset;
		}else{
			imDir.TIF_STRIPBYTECOUNTS=getIntTable(stream,stripByteOffset,(int)imDir.TIF_STRIPBYTECOUNTS_LENGTH);
		}
		try{
			stream.seek((int)(currentTagPosition+12));
			int offset_next_directory=swap(stream.readInt());
			imDir.OFFSET_NEXT_DIRECTORY=offset_next_directory;
		}catch(IOException e){
			e.printStackTrace();
		}
		if(imDir.TIF_CZ_LSMINFO_OFFSET!=0){
			imDir.TIF_CZ_LSMINFO=getCZ_LSMINFO(stream,imDir.TIF_CZ_LSMINFO_OFFSET);
		}
		return imDir;
	}

	private byte[] readTag(RandomAccessStream stream,int position){
		byte[] tag=new byte[12];
		try{
			stream.seek(position);
			stream.readFully(tag);
		}catch(IOException e){
			e.printStackTrace();
		}
		return tag;
	}

	private long[] getIntTable(RandomAccessStream stream,long position,int count){
		long[] offsets=new long[count];
		try{
			stream.seek((int)position);
			for(int i=0;i<count;i++)
				offsets[i]=swap(stream.readInt());
		}catch(IOException e){
			e.printStackTrace();
		}
		return offsets;
	}

	private CZ_LSMInfo getCZ_LSMINFO(RandomAccessStream stream,long position){
		CZ_LSMInfo cz=new CZ_LSMInfo();
		try{
			if(position==0)
				return cz;
			stream.seek((int)position+8);
			cz.DimensionX=swap(stream.readInt());
			cz.DimensionY=swap(stream.readInt());
			cz.DimensionZ=swap(stream.readInt());

			// number of channels
			cz.DimensionChannels=swap(stream.readInt());
			// Timestack size
			cz.DimensionTime=swap(stream.readInt());

			cz.IntensityDataType=swap(stream.readInt());

			cz.ThumbnailX=swap(stream.readInt());
			cz.ThumbnailY=swap(stream.readInt());
			cz.VoxelSizeX=swap(stream.readDouble());
			cz.VoxelSizeY=swap(stream.readDouble());
			cz.VoxelSizeZ=swap(stream.readDouble());

			stream.seek((int)position+88);
			cz.ScanType=swap(stream.readShort());
			stream.seek((int)position+108);
			cz.OffsetChannelColors=swap(stream.readInt());
			cz.TimeInterval=swap(stream.readDouble());
			// stream.seek((int) position + 120);
			cz.OffsetChannelDataTypes=swap(stream.readInt());
			stream.seek((int)position+132);
			cz.OffsetTimeStamps=swap(stream.readInt());

			if(cz.OffsetChannelDataTypes!=0){
				cz.OffsetChannelDataTypesValues=getOffsetChannelDataTypesValues(stream,cz.OffsetChannelDataTypes,cz.DimensionChannels);
			}
			if(cz.OffsetChannelColors!=0){
				ChannelNamesAndColors channelNamesAndColors=getChannelNamesAndColors(stream,cz.OffsetChannelColors,cz.DimensionChannels);
				cz.channelNamesAndColors=channelNamesAndColors;
			}
		}catch(IOException getCZ_LSMINFO_exception){
			getCZ_LSMINFO_exception.printStackTrace();
		}
		return cz;
	}

	private int[] getOffsetChannelDataTypesValues(RandomAccessStream stream,long position,long channelCount){
		int[] OffsetChannelDataTypesValues=new int[(int)channelCount];
		try{
			stream.seek((int)position);
			for(int i=0;i<channelCount;i++){
				OffsetChannelDataTypesValues[i]=swap(stream.readInt());
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		return OffsetChannelDataTypesValues;
	}

	private ChannelNamesAndColors getChannelNamesAndColors(RandomAccessStream stream,long position,long channelCount){
		ChannelNamesAndColors channelNamesAndColors=new ChannelNamesAndColors();
		try{
			stream.seek((int)position);
			channelNamesAndColors.BlockSize=swap(stream.readInt());
			channelNamesAndColors.NumberColors=swap(stream.readInt());
			channelNamesAndColors.NumberNames=swap(stream.readInt());
			channelNamesAndColors.ColorsOffset=swap(stream.readInt());
			channelNamesAndColors.NamesOffset=swap(stream.readInt());
			channelNamesAndColors.Mono=swap(stream.readInt());
			// reserved 4 words
			stream.seek((int)channelNamesAndColors.NamesOffset+(int)position);
			channelNamesAndColors.ChannelNames=new String[(int)channelCount];
			// long Namesize = channelNamesAndColors.BlockSize-
			// channelNamesAndColors.NamesOffset;
			for(int j=0;j<channelCount;j++){
				long size=swap(stream.readInt());
				channelNamesAndColors.ChannelNames[j]=readSizedNULLASCII(stream,size);
			}
			stream.seek((int)channelNamesAndColors.ColorsOffset+(int)position);
			channelNamesAndColors.Colors=new long[(int)(channelNamesAndColors.NumberColors)];
			for(int j=0;j<(int)(channelNamesAndColors.NumberColors);j++){
				channelNamesAndColors.Colors[j]=swap(stream.readInt());
			}
		}catch(IOException e){
			e.printStackTrace();
		}
		return channelNamesAndColors;
	}

	public static String readSizedNULLASCII(RandomAccessStream stream,long s){
		int offset=0;
		String tempstr=new String("");
		int in=0;
		char ch;
		boolean addchar=true;
		try{
			while(offset<s){
				in=stream.read();
				if(in==-1){
					break;
				}
				ch=(char)in;
				if(addchar==true){
					String achar=new Character(ch).toString();
					if(ch!=0x00){
						tempstr+=achar;
					}else{
						addchar=false;
					}
				}
				offset++;
			}
		}catch(IOException Read_ASCII_exception){
			Read_ASCII_exception.printStackTrace();
		}
		return tempstr;
	}

	/*
	 * apply_colors, applies color gradient; function taken out from Lut_Panel
	 * plugin
	 */
	private void applyColors(ImagePlus imp,int channel,Color[] gc,int i){
		FileInfo fi=new FileInfo();
		int size=256;
		fi.reds=new byte[size];
		fi.greens=new byte[size];
		fi.blues=new byte[size];
		fi.lutSize=size;
		float nColorsfl=size;
		float interval=size;
		float iR=gc[0].getRed();
		float iG=gc[0].getGreen();
		float iB=gc[0].getBlue();
		float idR=gc[1].getRed()-gc[0].getRed();
		float idG=gc[1].getGreen()-gc[0].getGreen();
		float idB=gc[1].getBlue()-gc[0].getBlue();
		idR=(idR/interval);
		idG=(idG/interval);
		idB=(idB/interval);
		int a=0;
		for(a=(int)(interval*0);a<(int)(interval*(0)+interval);a++,iR+=idR,iG+=idG,iB+=idB){
			fi.reds[a]=(byte)(iR);
			fi.greens[a]=(byte)(iG);
			fi.blues[a]=(byte)(iB);
		}
		int b=(int)(interval*0+interval)-1;
		fi.reds[b]=(byte)(gc[1].getRed());
		fi.greens[b]=(byte)(gc[1].getGreen());
		fi.blues[b]=(byte)(gc[1].getBlue());
		nColorsfl=size;
		if(nColorsfl>0){
			if(nColorsfl<size){
				interpolate(size,fi.reds,fi.greens,fi.blues,(int)nColorsfl);
			}
			showLut(imp,channel,fi,true);
			return;
		}
	}

	/*
	 * interpolate, modified from the ImageJ method by Wayne Rasband.
	 */
	private void interpolate(int size,byte[] reds,byte[] greens,byte[] blues,int nColors){
		byte[] r=new byte[nColors];
		byte[] g=new byte[nColors];
		byte[] b=new byte[nColors];
		System.arraycopy(reds,0,r,0,nColors);
		System.arraycopy(greens,0,g,0,nColors);
		System.arraycopy(blues,0,b,0,nColors);
		double scale=nColors/(float)size;
		int i1,i2;
		double fraction;
		for(int i=0;i<size;i++){
			i1=(int)(i*scale);
			i2=i1+1;
			if(i2==nColors)
				i2=nColors-1;
			fraction=i*scale-i1;
			reds[i]=(byte)((1.0-fraction)*(r[i1]&255)+fraction*(r[i2]&255));
			greens[i]=(byte)((1.0-fraction)*(g[i1]&255)+fraction*(g[i2]&255));
			blues[i]=(byte)((1.0-fraction)*(b[i1]&255)+fraction*(b[i2]&255));
		}
	}

	/*
	 * showLut, applies the new Lut on the actual image
	 */

	private void showLut(ImagePlus imp,int channel,FileInfo fi,boolean showImage){
		if(imp!=null){
			if(imp.getType()==ImagePlus.COLOR_RGB)
				IJ.error("Color tables cannot be assiged to RGB Images.");
			else{
				IndexColorModel cm=null;
				cm=new IndexColorModel(8,256,fi.reds,fi.greens,fi.blues);
				imp.setPosition(channel+1,imp.getSlice(),imp.getFrame());
				if(imp.isComposite()){
					((CompositeImage)imp).setChannelColorModel(cm);
					((CompositeImage)imp).updateChannelAndDraw();
				}else{
					imp.getProcessor().setColorModel(cm);
					imp.updateAndDraw();
				}
			}
		}
	}

	public ImagePlus open(RandomAccessStream stream,LsmFileInfo lsmFi){
		ImageDirectory firstImDir=(ImageDirectory)lsmFi.imageDirectories.get(0);
		if(firstImDir==null){
			return null;
		} // should not be if it is a true LSM file
		CZ_LSMInfo cz=firstImDir.TIF_CZ_LSMINFO;
		if(cz==null){
			return null;
		} // should not be, first Directory should have a CZ...
		ImagePlus imps=null;
		if(cz.ScanType<7){
			imps=readStack(stream,lsmFi,cz);
		}
		return imps;
	}

	private ImagePlus readStack(RandomAccessStream stream,LsmFileInfo lsmFi,CZ_LSMInfo cz){
		ImageDirectory firstImDir=(ImageDirectory)lsmFi.imageDirectories.get(0);
		lsmFi.url="";
		lsmFi.fileFormat=FileInfo.TIFF;
		lsmFi.pixelDepth=cz.VoxelSizeZ*1000000;
		lsmFi.pixelHeight=cz.VoxelSizeY*1000000;
		lsmFi.pixelWidth=cz.VoxelSizeX*1000000;
		lsmFi.unit=micrometer;
		lsmFi.valueUnit=micrometer;
		//lsmFi.nImages=1;
		lsmFi.intelByteOrder=true;
		ImageStack st=null;
		int datatype=(int)cz.IntensityDataType;
		if(datatype==0)
			datatype=cz.OffsetChannelDataTypesValues[0];
		switch(datatype){
		case 1:
			lsmFi.fileType=FileInfo.GRAY8;
			break;
		case 2:
			lsmFi.fileType=FileInfo.GRAY16_UNSIGNED;
			break;
		case 3:
			lsmFi.fileType=FileInfo.GRAY16_UNSIGNED;
			break;
		case 5:
			lsmFi.fileType=FileInfo.GRAY32_FLOAT;
			break;
		default:
			lsmFi.fileType=FileInfo.GRAY8;
			break;
		}
		ColorModel cm=null;
		if(lsmFi.fileType==FileInfo.COLOR8&&lsmFi.lutSize>0){
			cm=new IndexColorModel(8,lsmFi.lutSize,lsmFi.reds,lsmFi.greens,lsmFi.blues);
		}else{
			cm=LookUpTable.createGrayscaleColorModel(lsmFi.whiteIsZero);
		}
		st=new ImageStack((int)firstImDir.TIF_IMAGEWIDTH,(int)firstImDir.TIF_IMAGELENGTH,cm);
		//st=new ImageStack((int)firstImDir.TIF_IMAGEWIDTH,(int)firstImDir.TIF_IMAGELENGTH,null);
		firstImDir=null;
		ImageReader reader=null;
		int flength=0;
		lsmFi.stripOffsets=new int[1];
		lsmFi.stripLengths=new int[1];
		for(int imageCounter=0;imageCounter<lsmFi.imageDirectories.size();imageCounter++){
			ImageDirectory imDir=(ImageDirectory)lsmFi.imageDirectories.get(imageCounter);
			for(int i=0;i<imDir.TIF_STRIPBYTECOUNTS.length;i++){
				if(imDir.TIF_COMPRESSION==5){
					lsmFi.compression=FileInfo.LZW;
					flength=(int)new File(lsmFi.directory+System.getProperty("file.separator")+lsmFi.fileName).length();
					if(imDir.TIF_PREDICTOR==2)
						lsmFi.compression=FileInfo.LZW_WITH_DIFFERENCING;
				}else{
					lsmFi.compression=0;
				}
			}
			if(imDir.TIF_NEWSUBFILETYPE==0){
				lsmFi.width=(int)imDir.TIF_IMAGEWIDTH;
				lsmFi.height=(int)imDir.TIF_IMAGELENGTH;
				Object pixels;
				for(int channelCount=0;channelCount<(int)(cz.DimensionChannels);channelCount++){
					datatype=(int)cz.IntensityDataType;
					if(datatype==0){
						datatype=cz.OffsetChannelDataTypesValues[channelCount];
					}
					// IJ.log(""+datatype);
					switch(datatype){
					case 1:
						lsmFi.fileType=FileInfo.GRAY8;
						break;
					case 2:
						lsmFi.fileType=FileInfo.GRAY16_UNSIGNED;
						break;
					case 3:
						lsmFi.fileType=FileInfo.GRAY16_UNSIGNED;
						break;
					case 5:
						lsmFi.fileType=FileInfo.GRAY32_FLOAT;
						break;
					default:
						lsmFi.fileType=FileInfo.GRAY8;
						break;
					}
					lsmFi.stripLengths[0]=(int)imDir.TIF_STRIPBYTECOUNTS[channelCount];
					lsmFi.stripOffsets[0]=(int)imDir.TIF_STRIPOFFSETS[channelCount];
					reader=new ImageReader(lsmFi);
					if(channelCount<imDir.TIF_STRIPOFFSETS_LENGTH){
						if(lsmFi.stripLengths[0]+lsmFi.stripOffsets[0]>flength){
							lsmFi.stripLengths[0]=flength-lsmFi.stripOffsets[0];
						}
						try{
							stream.seek(lsmFi.stripOffsets[0]);
						}catch(IOException e){
							e.printStackTrace();
						}
						pixels=reader.readPixels(stream);
						st.addSlice("",pixels);
					}
				}
			}
		}
		ImagePlus imp=new ImagePlus(lsmFi.fileName,st);
		//IJ.log(""+(int)cz.DimensionChannels);
		int tempframes=(int)cz.DimensionTime;
		if(cz.ScanType==2) tempframes=1;
		imp.setDimensions((int)cz.DimensionChannels,(int)cz.DimensionZ,tempframes);
		//IJ.log(""+(int)cz.DimensionChannels+" , "+(int)cz.DimensionZ+" , "+(int)cz.DimensionTime);
		if((int)cz.DimensionChannels>1){
			imp=new CompositeImage(imp,CompositeImage.COLOR);
		}
		lsmFi.nImages=st.getSize();
		imp.setFileInfo(lsmFi);
		Calibration cal=new Calibration();
		cal.setUnit(lsmFi.unit);
		cal.pixelDepth=lsmFi.pixelDepth;
		cal.pixelHeight=lsmFi.pixelHeight;
		cal.pixelWidth=lsmFi.pixelWidth;
		imp.setCalibration(cal);
		Color[] color=new Color[2];
		color[0]=new Color(0,0,0);
		for(int channel=0;channel<(int)cz.DimensionChannels;channel++){
			int r=(int)(cz.channelNamesAndColors.Colors[channel]&255);
			int g=(int)((cz.channelNamesAndColors.Colors[channel]>>8)&255);
			int b=(int)((cz.channelNamesAndColors.Colors[channel]>>16)&255);
			color[1]=new Color(r,g,b);
			if(r==0&&g==0&&b==0){
				color[1]=Color.white;
			}
			applyColors(imp,channel,color,2);
		}
		if(imp.getOriginalFileInfo().fileType==FileInfo.GRAY16_UNSIGNED){
			double min=imp.getProcessor().getMin();
			double max=imp.getProcessor().getMax();
			imp.getProcessor().setMinAndMax(min,max);
		}
		int stackPosition=1;
		for(int i=1;i<=cz.DimensionTime;i++){
			for(int j=1;j<=cz.DimensionZ;j++){
				for(int k=1;k<=cz.DimensionChannels;k++){
					if(stackPosition<=imp.getStackSize()){
						String label=cz.channelNamesAndColors.ChannelNames[k-1];
						st.setSliceLabel(label,stackPosition++);
					}
				}
			}
		}
		setInfo(imp,lsmFi);
		return imp;
	}

	public ImagePlus setInfo(ImagePlus imp,LsmFileInfo lsm){
		ImageDirectory imDir=(ImageDirectory)lsm.imageDirectories.get(0);
		if(imDir==null)
			return null;
		CZ_LSMInfo cz=imDir.TIF_CZ_LSMINFO;
		String infos=new String();
		String stacksize=IJ.d2s(cz.DimensionZ,0);
		String width=IJ.d2s(lsm.width,0);
		String height=IJ.d2s(lsm.height,0);
		String channels=IJ.d2s(cz.DimensionChannels,0);
		String scantype="";
		switch(cz.ScanType){
		case 0:
			scantype="Normal X-Y-Z scan";
			break;
		case 1:
			scantype="Z scan";
			break;
		case 2:
			scantype="Line scan";
			break;
		case 3:
			scantype="Time series X-Y";
			break;
		case 4:
			scantype="Time series X-Y";
			break;
		case 5:
			scantype="Time series - Means of ROIs";
			break;
		case 6:
			scantype="Time series - X-Y-Z";
			break;
		default:
			scantype="UNKNOWN !";
			break;
		}
		String voxelsize_x=IJ.d2s(cz.VoxelSizeX*1000000,2)+" "+micrometer;
		String voxelsize_y=IJ.d2s(cz.VoxelSizeY*1000000,2)+" "+micrometer;
		String voxelsize_z=IJ.d2s(cz.VoxelSizeZ*1000000,2)+" "+micrometer;
		String timestacksize=IJ.d2s(cz.DimensionTime,0);
		if(cz.ScanType==2) timestacksize=IJ.d2s(1,0);
		String plane_width=IJ.d2s(cz.DimensionX*cz.VoxelSizeX,2)+" "+micrometer;
		String plane_height=IJ.d2s(cz.DimensionY*cz.VoxelSizeY,2)+" "+micrometer;
		String volume_depth=IJ.d2s(cz.DimensionZ*cz.VoxelSizeZ,2)+" "+micrometer;
		String time_interval=IJ.d2s(cz.TimeInterval,3)+" seconds";
		infos="Filename: "+lsm.fileName+"\n";
		infos+="Width: "+width+"\n";
		infos+="Height: "+height+"\n";
		infos+="Channels: "+channels+"\n";
		infos+="Z-stack size:"+stacksize+"\n";
		infos+="T-stack size: "+timestacksize+"\n";
		infos+="Scan type: "+scantype+"\n";
		infos+="Voxel size X: "+voxelsize_x+"\n";
		infos+="Voxel size Y: "+voxelsize_y+"\n";
		infos+="Voxel size Z: "+voxelsize_z+"\n";
		infos+="Plane width: "+plane_width+"\n";
		infos+="Plane height: "+plane_height+"\n";
		infos+="Plane depth: "+volume_depth+"\n";
		infos+="Time interval: "+time_interval+"\n";
		imp.setProperty("Info",infos);
		return imp;
	}

	private short swap(short x){
		return (short)((x<<8)|((x>>8)&0xff));
	}

	private int swap(int x){
		return (swap((short)x)<<16)|(swap((short)(x>>16))&0xffff);
	}

	private long swap(long x){
		return ((long)swap((int)(x))<<32)|(swap((int)(x>>32))&0xffffffffL);
	}

	private double swap(double x){
		return Double.longBitsToDouble(swap(Double.doubleToLongBits(x)));
	}

	protected class ImageDirectory{
		public long TIF_NEWSUBFILETYPE=0;
		public long TIF_IMAGEWIDTH=0;
		public long TIF_IMAGELENGTH=0;
		public long TIF_BITSPERSAMPLE_LENGTH=0;
		public long[] TIF_BITSPERSAMPLE_CHANNEL=new long[3];
		public long TIF_COMPRESSION=0;
		public long TIF_PHOTOMETRICINTERPRETATION=0;
		public long TIF_STRIPOFFSETS_LENGTH=0;
		public long[] TIF_STRIPOFFSETS;
		public long TIF_SAMPLESPERPIXEL=0;
		public long TIF_STRIPBYTECOUNTS_LENGTH=0;
		public long[] TIF_STRIPBYTECOUNTS;
		public long TIF_PLANARCONFIGURATION=0;
		public long TIF_PREDICTOR=0;
		public long TIF_COLORMAP=0;
		public long TIF_CZ_LSMINFO_OFFSET=0; // OFFSET
		public CZ_LSMInfo TIF_CZ_LSMINFO; // STRUCT
		public long OFFSET_NEXT_DIRECTORY=0;
	}

	protected class LsmFileInfo extends FileInfo{
		public ArrayList imageDirectories=new ArrayList();
	}

	private class CZ_LSMInfo{
		public long DimensionX=0;
		public long DimensionY=0;
		public long DimensionZ=0;
		public long DimensionChannels=0;
		public long DimensionTime=0;
		public long IntensityDataType=0;
		public long ThumbnailX=0;
		public long ThumbnailY=0;
		public double VoxelSizeX=0;
		public double VoxelSizeY=0;
		public double VoxelSizeZ=0;
		public int ScanType=0;
		public long OffsetChannelColors=0;
		public ChannelNamesAndColors channelNamesAndColors;
		public long OffsetChannelDataTypes=0;
		public int[] OffsetChannelDataTypesValues;
		public long OffsetTimeStamps;
		public double TimeInterval;
	}

	private class ChannelNamesAndColors{
		public long BlockSize=0;
		public long NumberColors=0;
		public long NumberNames=0;
		public long ColorsOffset=0;
		public long NamesOffset=0;
		public long Mono=0;
		public long[] Reserved;
		public long[] Colors;
		public String[] ChannelNames;
	}

}
