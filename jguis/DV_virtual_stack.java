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
import ij.VirtualStack;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.measure.Calibration;
import ij.process.ImageProcessor;

import java.io.IOException;

public class DV_virtual_stack extends VirtualStack{
	public DVFile file;
	public FileInfo fi;

	public DV_virtual_stack(String dir,String fname){
		try{
			file=new DVFile(dir,fname);
		}catch(IOException e){
			IJ.error(e.getMessage());
			return;
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
		default: {
			IJ.error("unsupported pixel type:"+file.getPixelTypeString());
			return;
		}
		}
		fi.unit="um";
		fi.info=file.getMetaDataString();
		ImagePlus imp=new ImagePlus(fname,this);
		imp.setProperty("Info",fi.info);
		imp.setFileInfo(fi);
		Calibration cal=imp.getCalibration();
		cal.pixelWidth=fi.pixelWidth;
		cal.pixelHeight=fi.pixelHeight;
		cal.pixelDepth=fi.pixelDepth;
		cal.setUnit(fi.unit);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(file.getNChannels(),file.getNSlices(),file.getNFrames());
		if(file.getNChannels()>1){
			imp=new CompositeImage(imp,CompositeImage.COLOR);
		}
		imp.show();
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

	public String getSliceLabel(int slice){
		return "";
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
