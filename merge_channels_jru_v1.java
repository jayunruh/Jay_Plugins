/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.*;
import jguis.*;
import java.awt.image.*;

public class merge_channels_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Number of Images to Merge",5,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int nimages=(int)gd.getNextNumber();
		ImagePlus[] images=jutils.selectImages(false,nimages);
		if(images==null){return;}
		int width=images[0].getWidth(); int height=images[0].getHeight();
		int frames=images[0].getNFrames(); int slices=images[0].getNSlices();
		int totchannels=0;
		for(int i=0;i<nimages;i++){
			totchannels+=images[i].getNChannels();
			if(images[i].getNFrames()>frames) frames=images[i].getNFrames();
			if(images[i].getNSlices()>slices) slices=images[i].getNSlices();
		}
		LUT[] luts=new LUT[totchannels];
		double[] mins=new double[totchannels];
		double[] maxs=new double[totchannels];
		ImageStack outstack=new ImageStack(width,height);
		int counter=0;
		for(int i=0;i<frames*slices;i++){
			for(int j=0;j<images.length;j++){
				int channels=images[j].getNChannels();
				int tslices=images[j].getNFrames()*images[j].getNSlices();
				ImageStack stack=images[j].getStack();
				for(int k=0;k<channels;k++){
					outstack.addSlice("",get_proc(stack,k,i,channels,tslices));
					//outstack.addSlice("",stack.getProcessor(1+k+j*channels));
				}
				if(i==0){
					for(int k=0;k<channels;k++){
						mins[counter]=stack.getProcessor(1+k).getMin();
						maxs[counter]=stack.getProcessor(1+k).getMax();
						luts[counter]=new LUT((IndexColorModel)stack.getProcessor(1+k).getColorModel(),mins[counter],maxs[counter]);
						counter++;
					}
				}
			}

			//images[i].changes=false;
			//images[i].close();
		}
		for(int i=0;i<totchannels;i++){
			outstack.getProcessor(i+1).setMinAndMax(mins[i],maxs[i]);
			outstack.getProcessor(i+1).setColorModel(luts[i]);
		}
		ImagePlus outimp=new ImagePlus("Merged Image",outstack);
		outimp.copyScale(images[0]);
		outimp.setOpenAsHyperStack(true);
		outimp.setDimensions(totchannels,slices,frames);
		CompositeImage ci=new CompositeImage(outimp,CompositeImage.COMPOSITE);
		ci.setLuts(luts);
		ci.show();
	}

	public ImageProcessor get_proc(ImageStack stack,int channel,int slice,int channels,int slices){
		int tchannel=channel;
		if(tchannel>=channels) tchannel=channels-1;
		int tslice=slice;
		if(tslice>=slices) tslice=slices-1;
		return stack.getProcessor(tslice*channels+tchannel+1);
	}

}
