/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
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

public class hist_2_overlay_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Histogram","Image"});
		if(imps==null) return;
		ImageWindow iw=imps[0].getWindow();
		int[] indices=(int[])jutils.runReflectionMethod(iw,"getroiindices",null);
		if(indices==null){
			IJ.error("Select Roi First");
			return;
		}
		int width=imps[1].getWidth(); int height=imps[1].getHeight();
		int[] mask=new int[width*height];
		//for(int i=0;i<mask.length;i++) mask[i]=0xff000000;
		for(int i=0;i<indices.length;i++) mask[indices[i]]=0xffff0000;
		ColorProcessor cp=new ColorProcessor(width,height,mask);
		ImageRoi roi=new ImageRoi(0,0,cp);
		roi.setZeroTransparent(true);
		imps[1].setOverlay(new Overlay(roi));
	}

}
