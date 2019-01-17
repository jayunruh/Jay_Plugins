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
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.IOException;

public class DV_file_reader{
	public DVFile file;
	public FileInfo fi;
	public boolean complex;

	public ImagePlus open(String dir,String fname){
		try{
			file=new DVFile(dir,fname);
		}catch(IOException e){
			IJ.error(e.getMessage());
			return null;
		}
		fi=new FileInfo();
		fi.fileFormat=FileInfo.RAW;
		fi.fileName=fname;
		fi.directory=dir;
		fi.width=file.getImageWidth();
		fi.height=file.getImageHeight();
		fi.nImages=file.getNumOfImages();
		fi.intelByteOrder=file.intelByteOrder();
		fi.pixelWidth=file.getPixelWidth();
		fi.pixelHeight=file.getPixelHeight();
		fi.pixelDepth=file.getPixelDepth();
		fi.offset=file.getImageDataOffset();
		int pixType=file.getPixelType();
		switch(pixType){
		case 0:
			fi.fileType=FileInfo.GRAY8;
			break;
		case 1:
			fi.fileType=FileInfo.GRAY16_SIGNED;
			break;
		case 2:
			fi.fileType=FileInfo.GRAY32_FLOAT;
			break;
		case 5:
			fi.fileType=FileInfo.GRAY16_SIGNED;
			break;
		case 6:
			fi.fileType=FileInfo.GRAY16_UNSIGNED;
			break;
		case 7:
			fi.fileType=FileInfo.GRAY32_INT;
			break;
		case 4:
			//this is the complex data type
			complex=true;
			fi.nImages=2*file.getNumOfImages();
			fi.fileType=FileInfo.GRAY32_FLOAT;
			break;
		default: {
			IJ.error("unsupported pixel type:"+file.getPixelTypeString());
			return null;
		}
		}
		fi.unit="um";
		fi.info=file.getMetaDataString();
		ImageStack stack=new ImageStack(fi.width,fi.height);
		if(complex){
			for(int i=0;i<fi.nImages;i+=2){
				ImageProcessor[] temp=getComplexProcessor(i+1);
    			stack.addSlice("",temp[0]);
    			stack.addSlice("",temp[1]);
    		}
		} else {
    		for(int i=0;i<fi.nImages;i++){
    			stack.addSlice("",getProcessor(i+1));
    		}
		}
		ImagePlus imp=new ImagePlus(fname,stack);
		imp.setProperty("Info",fi.info);
		imp.setFileInfo(fi);
		Calibration cal=imp.getCalibration();
		cal.pixelWidth=fi.pixelWidth;
		cal.pixelHeight=fi.pixelHeight;
		cal.pixelDepth=fi.pixelDepth;
		cal.setUnit(fi.unit);
		imp.setOpenAsHyperStack(true);
		int nslices=file.getNSlices(); if(complex) nslices*=2;
		imp.setDimensions(file.getNChannels(),file.getNSlices(),file.getNFrames());
		if(file.getNChannels()>1){
			imp=new CompositeImage(imp,CompositeImage.COLOR);
		}
		return imp;
	}

	public int getWidth(){
		return fi.width;
	}

	public int getHeight(){
		return fi.height;
	}

	public int getSize(){
		return fi.nImages;
	}
	
	public ImageProcessor[] getComplexProcessor(int stackslice){
		ImageProcessor ip1=getProcessor(stackslice);
		ImageProcessor ip2=getProcessor(stackslice+1);
		float[] pix1=(float[])ip1.getPixels();
		float[] pix2=(float[])ip2.getPixels();
		float[] combined=new float[pix1.length+pix2.length];
		System.arraycopy(pix1,0,combined,0,pix1.length);
		System.arraycopy(pix2,0,combined,pix1.length,pix2.length);
		float[] rpix=new float[pix1.length];
		float[] ipix=new float[pix1.length];
		for(int i=0;i<combined.length;i+=2){
			rpix[i/2]=combined[i];
			ipix[i/2]=combined[i+1];
		}
		ImageProcessor ip3=new FloatProcessor(getWidth(),getHeight(),rpix,null);
		ip3.flipVertical();
		ImageProcessor ip4=new FloatProcessor(getWidth(),getHeight(),ipix,null);
		ip4.flipVertical();
		return new ImageProcessor[]{ip3,ip4};
	}

	public ImageProcessor getProcessor(int stackslice){
		long size=fi.width*fi.height*fi.getBytesPerPixel();
		long oldlongOffset=fi.longOffset;
		int n=map_slice_number(stackslice)-1;
		fi.longOffset=fi.getOffset()+n*(size+fi.gapBetweenImages);
		int oldoffset=fi.offset;
		fi.offset=0;
		int oldnImages=fi.nImages;
		fi.nImages=1;
		FileOpener fo=new FileOpener(fi);
		ImagePlus timp=fo.open(false);
		fi.offset=oldoffset;
		fi.nImages=oldnImages;
		fi.longOffset=oldlongOffset;
		ImageProcessor tip=timp.getProcessor();
		tip.flipVertical();
		return tip;
	}

	public int map_slice_number(int stackslice){
		// this maps our ZTC data to the input CZT slice
		int channels=file.getNChannels();
		int frames=file.getNFrames();
		int slices=file.getNSlices();
		// first find out where we are in the CZT stack
		int mslice=stackslice-1;
		int frame=(int)((float)mslice/(float)(channels*slices));
		int temp=mslice-frame*slices*channels;
		int slice=(int)((float)temp/(float)channels);
		int channel=temp-slice*channels;
		// now map to the new space
		return channel*frames*slices+frame*slices+slice+1;
	}

}
