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
import jguis.*;

public class mask_hist_jru_v1 implements PlugIn {

	public void run(String arg) {
		String[] labels={"Image","Mask"};
		ImagePlus imp[]=jutils.selectImages(true,2,labels);
		if(imp==null){return;}
		ImagePlus imp1 = imp[0];
		int height = imp1.getHeight();
		int width = imp1.getWidth();
		ImagePlus imp2=imp[1];
		float[] mask=null;
		if(imp2!=null){
			mask=(float[])((imp2.getProcessor().convertToFloat()).getPixels());
		} else {
			mask=new float[width*height];
			for(int i=0;i<width*height;i++){mask[i]=1.0f;}
		}

		ImageStack datastack=imp1.getStack();
		int slices=datastack.getSize();
		boolean isfloat=(datastack.getPixels(1) instanceof float[]);
		int nis=0;
		int nnot=0;
		for(int i=0;i<width*height;i++){
			if(mask[i]>0.0f){nis++;}
			else{nnot++;}
		}
		float[] valsis=new float[slices*nis];
		float[] valsnot=new float[slices*nnot];
		for(int i=0;i<slices;i++){
			Object temppixels=datastack.getPixels(i+1);
			int counteris=0;
			int counternot=0;
			for(int j=0;j<width*height;j++){
				if(mask[j]>0.0f){
					if(isfloat){valsis[counteris+i*nis]=((float[])temppixels)[j];}
					else{valsis[counteris+i*nis]=(float)(((short[])temppixels)[j]&0xffff);}
					if(Float.isNaN(valsis[counteris+i*nis]) || Float.isInfinite(valsis[counteris+i*nis])){
						valsis[counteris+i*nis]=0.0f;
					}
					counteris++;
				} else {
					if(isfloat){valsnot[counternot+i*nnot]=((float[])temppixels)[j];}
					else{valsnot[counternot+i*nnot]=(float)(((short[])temppixels)[j]&0xffff);}
					if(Float.isNaN(valsnot[counternot+i*nnot]) || Float.isInfinite(valsnot[counternot+i*nnot])){
						valsnot[counternot+i*nnot]=0.0f;
					}
					counternot++;
				}
			}
		}
		if(nis>0)
			new PlotWindowHist("Mask Histogram","Value","Frequency",valsis,3).draw();
		if(nnot>0)
			new PlotWindowHist("Mask Not Histogram","Value","Frequency",valsnot,3).draw();
	}

}
