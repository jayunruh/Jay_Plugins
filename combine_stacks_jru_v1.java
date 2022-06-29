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

public class combine_stacks_jru_v1 implements PlugIn {
	//here we can combine stacks in multiple ways (merge channels, merge slices, and merge frames)

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Number of Images to Merge",5,0);
		String[] dims={"c","z","t"};
		gd.addChoice("Merge_Dimension",dims,dims[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int nimages=(int)gd.getNextNumber();
		int mergeindex=gd.getNextChoiceIndex();
		ImagePlus[] images=jutils.selectImages(false,nimages);
		if(images==null){return;}
		int width=images[0].getWidth(); int height=images[0].getHeight();
		int frames=images[0].getNFrames(); int slices=images[0].getNSlices(); int channels=images[0].getNChannels();
		if(mergeindex==0){
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
			for(int i=0;i<frames;i++){
				for(int s=0;s<slices;s++){
					for(int j=0;j<images.length;j++){
						int tchannels=images[j].getNChannels();
						int tslices=images[j].getNSlices();
						int tframes=images[j].getNFrames();
						ImageStack stack=images[j].getStack();
						for(int k=0;k<tchannels;k++){
							outstack.addSlice("",get_proc(stack,k,s,i,tchannels,tslices,tframes));
						}
						if(i==0){
							for(int k=0;k<tchannels;k++){
								mins[counter]=stack.getProcessor(1+k).getMin();
								maxs[counter]=stack.getProcessor(1+k).getMax();
								luts[counter]=new LUT((IndexColorModel)stack.getProcessor(1+k).getColorModel(),mins[counter],maxs[counter]);
								counter++;
							}
						}
					}
				}
			}
			for(int i=0;i<totchannels;i++){
				outstack.getProcessor(i+1).setMinAndMax(mins[i],maxs[i]);
				outstack.getProcessor(i+1).setColorModel(luts[i]);
			}
			ImagePlus outimp=new ImagePlus("Merged Image",outstack);
			outimp.setOpenAsHyperStack(true);
			outimp.setDimensions(totchannels,slices,frames);
			CompositeImage ci=new CompositeImage(outimp,CompositeImage.COMPOSITE);
			ci.setLuts(luts);
			ci.show();
		}
		if(mergeindex==1){
			//here we merge z slices
			int totslices=0;
			for(int i=0;i<nimages;i++){
				totslices+=images[i].getNSlices();
				if(images[i].getNFrames()>frames) frames=images[i].getNFrames();
				if(images[i].getNChannels()>channels) channels=images[i].getNChannels();
			}
			ImageStack outstack=new ImageStack(width,height);
			for(int i=0;i<frames;i++){
				for(int j=0;j<nimages;j++){
					int tchannels=images[j].getNChannels();
					int tslices=images[j].getNSlices();
					int tframes=images[j].getNFrames();
					ImageStack stack=images[j].getStack();
					for(int s=0;s<tslices;s++){
						for(int k=0;k<channels;k++){
							outstack.addSlice("",get_proc(stack,k,s,i,tchannels,tslices,tframes));
						}
					}
				}
			}
			ImagePlus outimp=new ImagePlus("Merged Image",outstack);
			outimp.setOpenAsHyperStack(true);
			outimp.setDimensions(channels,totslices,frames);
			if(images[0].isComposite()){
				CompositeImage ci=new CompositeImage(outimp,CompositeImage.COMPOSITE);
				ci.copyLuts(images[0]);
				ci.show();
			} else {
				outimp.show();
			}
		}
		if(mergeindex==2){
			//here we merge time points
			int totframes=0;
			for(int i=0;i<nimages;i++){
				totframes+=images[i].getNFrames();
				if(images[i].getNSlices()>slices) slices=images[i].getNSlices();
				if(images[i].getNChannels()>channels) channels=images[i].getNChannels();
			}
			ImageStack outstack=new ImageStack(width,height);
			for(int j=0;j<nimages;j++){
				int tchannels=images[j].getNChannels();
				int tslices=images[j].getNSlices();
				int tframes=images[j].getNFrames();
				ImageStack stack=images[j].getStack();
				for(int i=0;i<tframes;i++){
					for(int s=0;s<slices;s++){
						for(int k=0;k<channels;k++){
							outstack.addSlice("",get_proc(stack,k,s,i,tchannels,tslices,tframes));
						}
					}
				}
			}
			ImagePlus outimp=new ImagePlus("Merged Image",outstack);
			outimp.setOpenAsHyperStack(true);
			outimp.setDimensions(channels,slices,totframes);
			if(images[0].isComposite()){
				CompositeImage ci=new CompositeImage(outimp,CompositeImage.COMPOSITE);
				ci.copyLuts(images[0]);
				ci.show();
			} else {
				outimp.show();
			}
		}
	}

	public ImageProcessor get_proc(ImageStack stack,int channel,int slice,int frame,int channels,int slices,int frames){
		//here channel, slice, and frame are zero based
		int tchannel=channel;
		if(tchannel>=channels) tchannel=channels-1;
		int tslice=slice;
		if(tslice>=slices) tslice=slices-1;
		int tframe=frame;
		if(tframe>=frames) tframe=frames-1;
		return stack.getProcessor(tframe*slices*channels+tslice*channels+tchannel+1);
	}

}
