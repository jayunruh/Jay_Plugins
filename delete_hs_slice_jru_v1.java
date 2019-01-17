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

public class delete_hs_slice_jru_v1 implements PlugIn {
	int slices,frames,channels;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] choices={"Channel","Z Slice","Time Point"};
		gd.addChoice("Delete Current",choices,choices[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int index=gd.getNextChoiceIndex();
		ImagePlus imp=WindowManager.getCurrentImage();
		channels=imp.getNChannels();
		slices=imp.getNSlices();
		frames=imp.getNFrames();
		int currchan=imp.getChannel()-1;
		int currslice=imp.getSlice()-1;
		int currframe=imp.getFrame()-1;
		if(!imp.lock()){return;}
		if(index==0){
			delete_channel(imp,currchan);
		}
		if(index==1){
			delete_slice(imp,currslice);
		}
		if(index==2){
			delete_frame(imp,currframe);
		}
		imp.unlock();
		imp.updateAndDraw();
	}

	private void delete_channel(ImagePlus imp,int channel){
		ImageStack stack=imp.getStack();
		LUT[] luts=null;
		if(imp.isComposite()){
			luts=new LUT[channels-1];
			LUT[] oldluts=((CompositeImage)imp).getLuts();
			int counter=0;
			for(int i=0;i<channels;i++){
				if(i!=channel){
					luts[counter]=oldluts[i];
					counter++;
				}
			}
		}
		int counter=0;
		for(int i=0;i<frames;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					if(k==channel){
						stack.deleteSlice(k+j*channels+i*slices*channels+1-counter);
						counter++;
					}
				}
			}
		}
		boolean composite=imp.isComposite();
		imp.setStack(stack);
		imp.setDimensions(channels-1,slices,frames);
		if(composite){
			((CompositeImage)imp).setLuts(luts);
		} else {
			if(luts!=null && (channels-1)==1){
				imp.getProcessor().setColorModel(luts[0]);
			}
		}
		imp.setSlice(1);
	}

	private void delete_slice(ImagePlus imp,int slice){
		ImageStack stack=imp.getStack();
		int counter=0;
		for(int i=0;i<frames;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					if(j==slice){
						stack.deleteSlice(k+j*channels+i*slices*channels+1-counter);
						counter++;
					}
				}
			}
		}
		imp.setDimensions(channels,slices-1,frames);
	}

	private void delete_frame(ImagePlus imp,int frame){
		ImageStack stack=imp.getStack();
		int counter=0;
		for(int i=0;i<frames;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					if(i==frame){
						stack.deleteSlice(k+j*channels+i*slices*channels+1-counter);
						counter++;
					}
				}
			}
		}
		imp.setDimensions(channels,slices,frames-1);
	}

}
