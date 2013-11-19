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

public class edit_objects_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp2=WindowManager.getCurrentImage();
		int width=imp2.getWidth(); int height=imp2.getHeight();
		int slices=imp2.getStack().getSize();
		if(slices!=2){
			float[] pix=(float[])imp2.getProcessor().convertToFloat().getPixels();
			if(!(imp2.getProcessor() instanceof ByteProcessor)){
				IJ.showMessage("Need thresholded object image");
				return;
			}
			ImageStack dispstack=new ImageStack(width,height);
			dispstack.addSlice("",pix);
			dispstack.addSlice("",pix.clone());
			imp2=new ImagePlus("Outlined Objects",dispstack);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(2,1,1);
			imp2.copyScale(WindowManager.getCurrentImage());
			imp2=new CompositeImage(imp2,CompositeImage.COMPOSITE);
			imp2.show();
		}
		CompositeImage imp=(CompositeImage)imp2;
		imp.setPosition(1,1,1);
		LUT graylut=jutils.get_lut_for_color(Color.white);
		imp.setChannelColorModel(graylut);
		imp.setPosition(2,1,1);
		LUT redlut=jutils.get_lut_for_color(Color.red);
		imp.setChannelColorModel(redlut);
		imp.setPosition(1,1,1);
		imp.updateAndRepaintWindow();
		threshold_panel tp=new threshold_panel();
		findblobs3 fb=new findblobs3(width,height);
		float[] temp=(float[])imp.getStack().getPixels(2);
		float[] objects=fb.dofindblobs(temp,254.0f);
		tp.init(imp,objects,false);
		threshold_panel.launch_frame(tp);
	}

}
