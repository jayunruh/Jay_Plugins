/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;
import jalgs.jsim.*;

public class add_random_rois_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int currframe=imp.getT();
		int currchan=imp.getC();
		if(slices==1){
			slices=frames;
			frames=1;
			currframe=1;
		}
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Numb_of_rois",1000,0);
		gd.addCheckbox("Limit to mask",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int nrois=(int)gd.getNextNumber();
		boolean limit=gd.getNextBoolean();
		rngs random=new rngs();
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		for(int i=0;i<nrois;i++){
			int x=(int)random.unidev(0.0,(double)width-0.00001);
			int y=(int)random.unidev(0.0,(double)height-0.00001);
			int z=(int)random.unidev(0.0,(double)slices-0.00001);
			if(limit){
				byte[] mask=(byte[])stack.getPixels(currchan-1+z*channels+(currframe-1)*slices*channels+1);
				while(mask[x+y*width]==(byte)0){
					x=(int)random.unidev(0.0,(double)width-0.00001);
					y=(int)random.unidev(0.0,(double)height-0.00001);
					z=(int)random.unidev(0.0,(double)slices-0.00001);
					mask=(byte[])stack.getPixels(currchan-1+z*channels+(currframe-1)*slices*channels+1);
				}
			}
			PointRoi roi=new PointRoi(x,y);
			roi.setPosition(z+1);
			rman.addRoi(roi);
			rman.getRoi(rman.getCount()-1).setPosition(currchan,z+1,currframe);
			IJ.showProgress(i,nrois);
		}
	}

}
