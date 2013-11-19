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

public class duplicate_hs_slice_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] choices={"Channel","Z Slice","Time Point"};
		gd.addChoice("Duplicate_Current",choices,choices[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int index2=gd.getNextChoiceIndex();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int channels=imp.getNChannels();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		ImageStack stack=imp.getStack();
		int currchan=imp.getChannel()-1;
		int currslice=imp.getSlice()-1;
		int currframe=imp.getFrame()-1;
		ImageStack retstack=new ImageStack(width,height);
		if(index2==0){
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<channels;k++){
						if(k==currchan){
							int index=k+j*channels+i*slices*channels+1;
							retstack.addSlice(stack.getSliceLabel(index),stack.getPixels(index));
						}
					}
				}
			}
			ImagePlus imp2=new ImagePlus("Duplicate",retstack);
			imp2.copyScale(imp);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(1,slices,frames);
			imp2.show();
		}
		if(index2==1){
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<channels;k++){
						if(j==currslice){
							int index=k+j*channels+i*slices*channels+1;
							retstack.addSlice(stack.getSliceLabel(index),stack.getPixels(index));
						}
					}
				}
			}
			ImagePlus imp2=new ImagePlus("Duplicate",retstack);
			imp2.copyScale(imp);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(channels,1,frames);
			if(channels>1 && imp.isComposite()){
				CompositeImage ci=new CompositeImage(imp2,((CompositeImage)imp).getMode());
				ci.copyLuts(imp);
				ci.show();
			} else {
				imp2.show();
			}
		}
		if(index2==2){
			int counter=0;
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<channels;k++){
						if(i==currframe){
							int index=k+j*channels+i*slices*channels+1;
							retstack.addSlice(stack.getSliceLabel(index),stack.getPixels(index));
						}
					}
				}
			}
			ImagePlus imp2=new ImagePlus("Duplicate",retstack);
			imp2.copyScale(imp);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(channels,slices,1);
			if(channels>1 && imp.isComposite()){
				CompositeImage ci=new CompositeImage(imp2,((CompositeImage)imp).getMode());
				ci.copyLuts(imp);
				ci.show();
			} else {
				imp2.show();
			}
		}
	}

}
