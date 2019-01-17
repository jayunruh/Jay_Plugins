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
import jalgs.jseg.*;

public class filter_objects_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		boolean clearedges=true;
		gd.addCheckbox("Clear_Edges?",clearedges);
		int minarea=10;
		gd.addNumericField("Min Size (pixels)",minarea,0);
		int maxarea=100;
		gd.addNumericField("Max Size (pixels)",maxarea,0);
		boolean newimage=true;
		gd.addCheckbox("Create_New_Image?",newimage);
		gd.addCheckbox("Output_Indexed?",false);
		gd.addCheckbox("Create_Border (circ)",false);
		gd.addNumericField("Border_Thickness (pixels)",3,0);
		gd.addNumericField("Border_Gap (pixels)",0,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		clearedges=gd.getNextBoolean();
		minarea=(int)gd.getNextNumber();
		maxarea=(int)gd.getNextNumber();
		newimage=gd.getNextBoolean();
		boolean outindexed=gd.getNextBoolean();
		boolean border=gd.getNextBoolean();
		int bordrad=(int)gd.getNextNumber();
		int bordgap=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		findblobs3 fb=new findblobs3(width,height);
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<size;i++){
			byte[] data=(byte[])stack.getProcessor(i+1).getPixels();
			float[] objects=fb.dofindblobs(data);
			if(clearedges) fb.clear_edges(objects);
			int[] filter={minarea,maxarea};
			fb.filter_area(objects,filter);
			if(border) objects=fb.get_circ(objects,bordrad,bordgap);
			if(newimage || outindexed){
				if(!outindexed) retstack.addSlice("",fb.tobinary(objects,true));
				else retstack.addSlice("",objects);
			} else {
				stack.setPixels(fb.tobinary(objects,true),i+1);
				//data=fb.tobinary(objects,true);
			}
		}
		if(newimage || outindexed){
			new ImagePlus(imp.getTitle()+" filtered",retstack).show();
		} else {
			imp.setStack(null,stack);
			imp.updateAndDraw();
		}
	}

}
