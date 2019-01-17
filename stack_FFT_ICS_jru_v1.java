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
import jalgs.jfft.*;

public class stack_FFT_ICS_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean avg_sections,interpolate;
		int i,height,width,start,end;
		float avg,variance,overall_avg,g0;
		ImagePlus imp = WindowManager.getCurrentImage();
		height = imp.getHeight();
		width = imp.getWidth();
		ImageStack stack = imp.getStack();
		int size = stack.getSize();
		GenericDialog gd = new GenericDialog("Options");
		gd.addCheckbox("Avg Quadrants?",true);
		gd.addCheckbox("Interpolate G(0)?",true);
		gd.addNumericField("Beginning Slice",1.0,0);
		gd.addNumericField("End Slice",(double)size,0);
		boolean brightcorr=false;
		gd.addCheckbox("Brightcorr?",brightcorr);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		avg_sections = gd.getNextBoolean();
		interpolate = gd.getNextBoolean();
		start = (int)gd.getNextNumber();
		if(start<1){start=1;}
		end = (int)gd.getNextNumber();
		if(end>size){end=size;}
		size = end-start+1;
		brightcorr=gd.getNextBoolean();
		int[] windex=fftutils.get_best_index(width,false,19);
		int[] hindex=fftutils.get_best_index(height,false,19);
		//IJ.log(""+windex[1]+" , "+hindex[1]);
		autocorr2D ac2D=new autocorr2D(windex[1],hindex[1],windex[0],hindex[0]);
		float[] ac=new float[windex[1]*hindex[1]];
		for(i=start;i<=end;i++){
			//get the appropriate stack image
			float[] pixels=(float[])get_center_crop(stack.getPixels(i),width,height,windex[1],hindex[1]);
			float[] tempac=ac2D.doautocorr2D(pixels,true,true,avg_sections,brightcorr);
			for(int j=0;j<windex[1]*hindex[1];j++){
				ac[j]+=tempac[j]/(float)size;
			}
			IJ.showProgress(i-start,size);
		}
		if(interpolate){
			ac[(hindex[1]/2)*windex[1]+windex[1]/2]=(ac[(hindex[1]/2)*windex[1]+windex[1]/2-1]+ac[(hindex[1]/2)*windex[1]+windex[1]/2+1])/2.0f;
		}
		//output the autocorrelated files
		ImagePlus imp2=new ImagePlus("Autocorrelation",new FloatProcessor(windex[1],hindex[1],ac,null));
		imp2.copyScale(imp);
		imp2.show();
	}

	public Object get_center_crop(Object pixels,int width,int height,int newwidth,int newheight){
		if(newwidth!=width || newheight!=height){
			return algutils.get_region(pixels,width/2,height/2,newwidth,newheight,width,height);
		} else {
			return algutils.convert_arr_float(pixels);
		}
	}

}
