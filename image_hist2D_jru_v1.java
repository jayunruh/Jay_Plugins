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

public class image_hist2D_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X_Channel",1,0);
		gd.addNumericField("Y_Channel",2,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int xchan=(int)gd.getNextNumber();
		int ychan=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int currchan=imp.getChannel();
		int nchannels=imp.getNChannels();
		if(xchan>nchannels || ychan>nchannels || xchan<1 || ychan<1){
			IJ.showMessage("Channel out of bounds");
			return;
		}
		int currslice=imp.getSlice();
		int nslices=imp.getNSlices();
		int currframe=imp.getFrame();
		int nframes=imp.getNFrames();
		ImageStack stack=imp.getStack();
		Object pixx=jutils.get3DSlice(stack,currframe-1,currslice-1,xchan-1,nframes,nslices,nchannels);
		Object pixy=jutils.get3DSlice(stack,currframe-1,currslice-1,ychan-1,nframes,nslices,nchannels);
		float[] pixxf=jutils.convert_arr_float(pixx);
		float[] pixyf=jutils.convert_arr_float(pixy);
		Roi roi=imp.getRoi();
		if(roi!=null){
			boolean[] mask=jutils.roi2mask(roi,width,height);
			int npix=0;
			for(int i=0;i<width*height;i++) if(mask[i]) npix++;
			float[] temppixxf=new float[npix];
			float[] temppixyf=new float[npix];
			int counter=0;
			for(int i=0;i<width*height;i++){
				if(mask[i]){
					temppixxf[counter]=pixxf[i];
					temppixyf[counter]=pixyf[i];
					counter++;
				}
			}
			pixxf=temppixxf;
			pixyf=temppixyf;
		}
		new PlotWindow2DHist(imp.getTitle()+" 2D histogram","Channel "+xchan,"Channel "+ychan,pixxf,pixyf,null).draw();
	}

}
