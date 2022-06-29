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

public class fill_roi_float_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Fill Value",0.0,5,15,null);
		gd.addCheckbox("Single Slice",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float fval=(float)gd.getNextNumber();
		boolean single=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int type=0;
		if(stack.getPixels(1) instanceof short[]) type=1;
		if(stack.getPixels(1) instanceof byte[]) type=2;
		short sval=(short)0;
		byte bval=(byte)0;
		if(type==1){
			if(fval>0 && fval<65535.0f) sval=(short)((int)fval);
			else if(fval>=65535.0f) sval=(short)65535;
		}
		if(type==2){
			if(fval>0 && fval<255.0f) bval=(byte)((int)fval);
			else if(fval>=255.0f) bval=(byte)255;
		}
		Roi roi=imp.getRoi();
		if(single){
			ImageProcessor ip=imp.getProcessor();
			Object temp=ip.getPixels();
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					if(roi.contains(j,i)){
						if(type==0) ((float[])temp)[j+i*width]=fval;
						if(type==1) ((short[])temp)[j+i*width]=sval;
						if(type==2) ((byte[])temp)[j+i*width]=bval;
					}
				}
			}
		} else {
		for(int k=0;k<stack.getSize();k++){
			Object temp=stack.getPixels(k+1);
			for(int i=0;i<height;i++){
				for(int j=0;j<width;j++){
					if(roi.contains(j,i)){
						if(type==0) ((float[])temp)[j+i*width]=fval;
						if(type==1) ((short[])temp)[j+i*width]=sval;
						if(type==2) ((byte[])temp)[j+i*width]=bval;
					}
				}
			}
		}
		}
		imp.updateAndDraw();
		/*ImagePlus imp2=new ImagePlus("Filled Image",new FloatProcessor(width,height,temp,null));
		imp2.copyScale(imp);
		imp2.show();*/
	}

}
