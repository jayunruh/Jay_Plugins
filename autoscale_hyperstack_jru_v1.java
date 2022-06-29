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
import jalgs.*;

public class autoscale_hyperstack_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int nchans=imp.getNChannels();
		int nslices=imp.getNSlices();
		int nframes=imp.getNFrames();
		float[][] stats=new float[nchans][2];
		int maxslice=0;
		for(int i=0;i<nchans;i++){
			stats[i][0]=jstatistics.getstatistic("Min",stack.getPixels(i+1),null);
			stats[i][1]=jstatistics.getstatistic("Max",stack.getPixels(i+1),null);
		}
		int counter=0;
		float overmax=stats[0][1];
		for(int i=0;i<nframes;i++){
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
					Object pix=stack.getPixels(counter+1);
					counter++;
					float min=jstatistics.getstatistic("Min",pix,null);
					float max=jstatistics.getstatistic("Max",pix,null);
					if(min<stats[k][0]) stats[k][0]=min;
					if(max>stats[k][1]) stats[k][1]=max;
					if(max>overmax){
						overmax=max;
						maxslice=j;
					}
				}
			}
		}
		int currframe=imp.getFrame();
		int currslice=imp.getSlice();	
		int currch=imp.getChannel();
		for(int i=0;i<nchans;i++){
			imp.setPositionWithoutUpdate(i+1,currslice,currframe);
			imp.setDisplayRange(stats[i][0],stats[i][1]);
		}
		if(imp instanceof CompositeImage) ((CompositeImage)imp).reset();
		imp.setPosition(currch,maxslice,currframe);
		imp.updateAndDraw();
	}

}
