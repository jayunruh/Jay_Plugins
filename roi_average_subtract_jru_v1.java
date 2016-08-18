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
import jalgs.*;

public class roi_average_subtract_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		Object[] temp=exec(imp);
		ImagePlus imp2 =(ImagePlus)exec(imp)[0];
		imp2.copyScale(imp);
		imp2.setOpenAsHyperStack(true);
		imp2.setDimensions(imp.getNChannels(),imp.getNSlices(),imp.getNFrames());
		if(imp.getNChannels()>1){
			CompositeImage ci=new CompositeImage(imp2);
			if(imp.isComposite()){
				LUT[] lut=((CompositeImage)imp).getLuts();
				ci.setLuts(lut);
				ci.resetDisplayRanges();
				ci.show();
			} else {
				ci.show();
			}
		} else {
			imp2.show();
		}
	}

	public Object[] exec(ImagePlus imp){
		return exec(imp,imp.getRoi());
	}

	public Object[] exec(ImagePlus imp,Roi roi){
		ImageStack stack = imp.getStack();
		int slices=stack.getSize();
		int width=stack.getWidth();
		int height=stack.getHeight();
		ImageStack result_stack=new ImageStack(width,height);
		boolean[] mask=jstatistics.poly2mask(roi.getPolygon(),width,height);
		Rectangle r=roi.getPolygon().getBounds();
		int[] lims={r.x,r.x+r.width,r.y,r.y+r.height};
		for(int i=0;i<slices;i++){
			Object temp2=stack.getPixels(i+1);
			float avg=jstatistics.getstatistic("Avg",temp2,width,height,mask,lims,null);
			float[] pixels=algutils.convert_arr_float(temp2);
			float[] temp=new float[width*height];
			for(int j=0;j<width*height;j++) temp[j]=pixels[j]-avg;
			result_stack.addSlice("",(Object)temp);
			IJ.showProgress(i,slices);
		}
		return new Object[]{new ImagePlus("Subtracted",result_stack)};
	}

}
