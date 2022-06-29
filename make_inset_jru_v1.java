/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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

public class make_inset_jru_v1 implements PlugIn, DialogListener {
	int xpos,ypos,inwidth,inheight;
	ImagePlus imp;

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"destination","inset"});
		if(imps==null) return;
		imp=imps[0];
		inwidth=imps[1].getWidth();
		inheight=imps[1].getHeight();
		xpos=10;
		ypos=10;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X position",xpos,0);
		gd.addNumericField("Y position",ypos,0);
		gd.addDialogListener(this);
		gd.showDialog(); if(gd.wasCanceled()){imp.setHideOverlay(true); return;}
		xpos=(int)gd.getNextNumber();
		ypos=(int)gd.getNextNumber();
		imp.setHideOverlay(true);
		//IJ.run(imp, "Insert...", "source="+imps[1].getTitle()+" destination="+imps[0].getTitle()+" x=10 y=10");
		insert(imps[1],imps[0],xpos,ypos,true);
	}

	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e){
		//here we update the current frame's overlay with the inset square
		if(IJ.escapePressed()){return false;}
		xpos=(int)gd.getNextNumber();
		ypos=(int)gd.getNextNumber();
		Roi roi=new Roi(xpos,ypos,inwidth,inheight);
		Overlay overlay=new Overlay();
		roi.setStrokeColor(Color.red);
		overlay.add(roi);
		imp.setOverlay(overlay);
		return true;
	}

	public void insert(ImagePlus imp1, ImagePlus imp2, int x, int y,boolean border) {
		ImageStack stack1 = imp1.getStack();
		ImageStack stack2 = imp2.getStack();
		int size1 = stack1.getSize();
		int size2 = stack2.getSize();
		ImageProcessor ip1, ip2;
		for (int i=1; i<=size2; i++) {
			ip1 = stack1.getProcessor(i<=size1?i:size1);
			ip2 = stack2.getProcessor(i);
			ip2.setColor(Toolbar.getForegroundColor());
			ip2.insert(ip1, x, y);
			if(border){
				ip2.drawRect(x,y,inwidth,inheight);
			}
			stack2.setPixels(ip2.getPixels(), i);
		}
		imp2.setStack(null, stack2);
	}

}
