/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;
import jalgs.jfit.*;

public class KISS_unmixing_jru_v1 implements PlugIn {

	public void run(String arg) {
		//basic procedure for each image: subback,bin,unmix
		//at the end merge the unmixed images
		ImagePlus[] imps=jutils.selectImages(false,6,new String[]{"488_Image","561_Image","633_Image","488_spectra","561_spectra","633_spectra"});
		if(imps==null) return;
		//now subtract the background
		Roi roi=null;
		for(int i=0;i<3;i++){
			roi=imps[i].getRoi();
			if(roi!=null) break;
		}
		if(roi==null){
			IJ.error("Select a background roi on at least one image");
			return;
		}
		ImagePlus[] subimps=new ImagePlus[3];
		for(int i=0;i<3;i++){
			subimps[i]=roiavgsub(imps[i],roi);
		}
		//next spatially bin the images 2 x 2
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Spatial_Bin?",false);
		gd2.addNumericField("Bin_By",2,0);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		boolean bin=gd2.getNextBoolean();
		int binby=(int)gd2.getNextNumber();
		if(bin){
			for(int i=0;i<3;i++){
				subimps[i]=bin2D(subimps[i],binby);
			}
		}
		int[] nchan={imps[0].getNChannels(),imps[1].getNChannels(),imps[2].getNChannels()};
		GenericDialog gd5=new GenericDialog("Select Spectral Regions");
		int start1=0; int end1=12;  if(end1>=nchan[0]) end1=nchan[0]-1;
		int start2=0; int end2=22; if(end2>=nchan[1]) end2=nchan[1]-1;
		int start3=9; int end3=22; if(end3>=nchan[2]) end3=nchan[2]-1;
		gd5.addNumericField("488_Start",start1,0);
		gd5.addNumericField("488_End",end1,0);
		gd5.addNumericField("561_Start",start2,0);
		gd5.addNumericField("561_End",end2,0);
		gd5.addNumericField("633_Start",start3,0);
		gd5.addNumericField("633_End",end3,0);
		gd5.showDialog(); if(gd5.wasCanceled()) return;
		start1=(int)gd5.getNextNumber(); if(start1<0) start1=0;
		end1=(int)gd5.getNextNumber(); if(end1>=nchan[0]) end1=nchan[0]-1;
		start2=(int)gd5.getNextNumber(); if(start2<0) start2=0;
		end2=(int)gd5.getNextNumber(); if(end2>=nchan[1]) end2=nchan[1]-1;
		start3=(int)gd5.getNextNumber(); if(start3<0) start3=0;
		end3=(int)gd5.getNextNumber(); if(end3>=nchan[2]) end3=nchan[2]-1;
		//now unmix
		subimps[0]=unmix(subimps[0],imps[3].getWindow(),start1,end1);
		subimps[1]=unmix(subimps[1],imps[4].getWindow(),start2,end2);
		subimps[2]=unmix(subimps[2],imps[5].getWindow(),start3,end3);
		ImageStack blankstack=new ImageStack(subimps[0].getWidth(),subimps[1].getHeight());
		blankstack.addSlice("",new float[subimps[0].getWidth()*subimps[1].getHeight()]);
		//and merge the windows
		ImageStack[] stacks={subimps[0].getStack(),subimps[1].getStack(),subimps[2].getStack(),blankstack};
		int totchannels=stacks[0].getSize()+stacks[1].getSize()+stacks[2].getSize();
		String[] channeloptions=new String[totchannels+1];
		int[][] optionindices=new int[totchannels+1][2];
		int counter=0;
		for(int i=0;i<stacks[0].getSize();i++){
			channeloptions[counter]="488_image"+(i+1); optionindices[counter][0]=0; optionindices[counter][1]=i+1; counter++;
		}
		for(int i=0;i<stacks[1].getSize();i++){
			channeloptions[counter]="561_image"+(i+1); optionindices[counter][0]=1; optionindices[counter][1]=i+1; counter++;
		}
		for(int i=0;i<stacks[2].getSize();i++){
			channeloptions[counter]="633_image"+(i+1); optionindices[counter][0]=2; optionindices[counter][1]=i+1; counter++;
		}
		channeloptions[totchannels]="null";
		optionindices[totchannels][0]=3; optionindices[totchannels][1]=0;
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Ch1",channeloptions,channeloptions[0]);
		//gd.addChoice("Ch2",channeloptions,channeloptions[1]);
		gd.addChoice("Ch2",channeloptions,channeloptions[stacks[0].getSize()]);
		gd.addChoice("Ch3",channeloptions,channeloptions[stacks[0].getSize()+1]);
		gd.addChoice("Ch4",channeloptions,channeloptions[stacks[0].getSize()+stacks[1].getSize()]);
		gd.addChoice("Ch5",channeloptions,channeloptions[stacks[0].getSize()+stacks[1].getSize()+1]);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int[][] indices=new int[5][];
		for(int i=0;i<5;i++){
			int temp=gd.getNextChoiceIndex();
			indices[i]=optionindices[temp];
		}
		
		//Object[] newpix={stacks[0].getPixels(1),stacks[0].getPixels(2),stacks[1].getPixels(1),stacks[1].getPixels(2),stacks[2].getPixels(1),stacks[2].getPixels(2)};
		//Object[] newpix={stacks[indices[0][0]].getPixels(indices[0][1]),stacks[indices[1][0]].getPixels(indices[1][1]),stacks[indices[2][0]].getPixels(indices[2][1]),stacks[indices[3][0]].getPixels(indices[3][1]),stacks[indices[4][0]].getPixels(indices[4][1]),stacks[indices[5][0]].getPixels(indices[5][1])};
		Object[] newpix={stacks[indices[0][0]].getPixels(indices[0][1]),stacks[indices[1][0]].getPixels(indices[1][1]),stacks[indices[2][0]].getPixels(indices[2][1]),stacks[indices[3][0]].getPixels(indices[3][1]),stacks[indices[4][0]].getPixels(indices[4][1])};
		ImageStack finstack=jutils.array2stack(newpix,subimps[0].getWidth(),subimps[1].getHeight());
		ImagePlus imp5=new ImagePlus("Merged Image",finstack);
		imp5.copyScale(subimps[0]);
		imp5.setOpenAsHyperStack(true);
		imp5.setDimensions(newpix.length,1,1);
		(new CompositeImage(imp5,CompositeImage.COMPOSITE)).show();
		//IJ.run("set lut colors jru v1", "color1=blue color2=blue color3=green color4=green color5=red color6=red");
		IJ.run("set lut colors jru v1", "color1=blue color2=green color3=green color4=red color5=red");
	}

	public ImagePlus unmix(ImagePlus imp,ImageWindow plot,int start,int end){
		int npix=imp.getWidth()*imp.getHeight();
		float[][] spectra=(float[][])jutils.runPW4VoidMethod(plot,"getYValues");
		float[][] contr=(new linear_unmix(spectra,start,end)).unmix(jutils.stack2array(imp.getStack()),npix,true);
		ImageStack stack=jutils.array2stack(contr,imp.getWidth(),imp.getHeight());
		ImagePlus retimp=new ImagePlus("Unmixed",stack);
		retimp.setOpenAsHyperStack(true);
		retimp.setDimensions(imp.getNChannels(),imp.getNSlices(),imp.getNFrames());
		retimp.copyScale(imp);
		return retimp;
	}

	public ImagePlus roiavgsub(ImagePlus imp,Roi roi){
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		int width=imp.getWidth(); int height=imp.getHeight();
		boolean[] mask=jstatistics.poly2mask(roi.getPolygon(),width,height);
		Rectangle r=roi.getPolygon().getBounds();
		int[] lims={r.x,r.x+r.width,r.y,r.y+r.height};
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<size;i++){
			Object temp2=stack.getPixels(i+1);
			float avg=jstatistics.getstatistic("Avg",temp2,width,height,mask,lims,null);
			float[] pixels=algutils.convert_arr_float(temp2);
			float[] temp=new float[width*height];
			for(int j=0;j<width*height;j++) temp[j]=pixels[j]-avg;
			retstack.addSlice("",(Object)temp);
		}
		ImagePlus retimp=new ImagePlus("Binned",retstack);
		retimp.setOpenAsHyperStack(true);
		retimp.setDimensions(imp.getNChannels(),imp.getNSlices(),imp.getNFrames());
		retimp.copyScale(imp);
		return retimp;
	}

	public ImagePlus bin2D(ImagePlus imp,int binby){
		ImageStack stack=imp.getStack();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack retstack=new ImageStack(width/binby,height/binby);
		for(int i=0;i<stack.getSize();i++){
			float[] image=algutils.convert_arr_float2(stack.getPixels(i+1));
			float[] binned=jsmooth.bin2D(image,width,height,binby,true);
			retstack.addSlice("",binned);
		}
		ImagePlus retimp=new ImagePlus("Binned",retstack);
		retimp.setOpenAsHyperStack(true);
		retimp.setDimensions(imp.getNChannels(),imp.getNSlices(),imp.getNFrames());
		retimp.copyScale(imp);
		jutils.set_psize(retimp,jutils.get_psize(imp)*2.0);
		return retimp;
	}

}
