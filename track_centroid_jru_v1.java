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
import jalgs.jseg.*;

public class track_centroid_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] centeroptions={"use_upper_left","use_center","use_initial"};
		gd.addChoice("Traj_Origin",centeroptions,centeroptions[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int origindex=gd.getNextChoiceIndex();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		Rectangle rect=new Rectangle(0,0,width,height);
		if(imp.getRoi()!=null) rect=imp.getRoi().getBounds();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		float[][] traj=new float[2][slices];
		if(stack.getPixels(1) instanceof byte[]){
			for(int i=0;i<slices;i++){
				byte[] image=(byte[])stack.getPixels(i+1);
				float[] centroid=measure_object.centroid(image,width,height,rect);
				traj[0][i]=centroid[0];
				traj[1][i]=centroid[1];
			}
		} else {
			for(int i=0;i<slices;i++){
				float[] image=(float[])stack.getPixels(i+1);
				float[] centroid=measure_object.centroid(image,width,height,rect);
				traj[0][i]=centroid[0];
				traj[1][i]=centroid[1];
			}
		}
		float subx=0.0f; float suby=0.0f;
		if(origindex==1){
			subx=(float)(width/2); suby=(float)(height/2);
		}
		if(origindex==2){
			subx=traj[0][0]; suby=traj[1][0];
		}
		if(slices>1){
			for(int i=0;i<slices;i++){
				traj[0][i]-=subx;
				traj[1][i]-=suby;
			}
			new PlotWindow4("Max Track","x","y",traj[0],traj[1]).draw();
		} else {
			IJ.log("maxx = "+(traj[0][0]-subx)+" , maxy = "+(traj[1][0]-suby));
		}
	}

}
