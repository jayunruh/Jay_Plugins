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
import jguis.*;

public class find_carpet_max_column_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Interpolate?",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean interp=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		float[] data=(float[])((FloatProcessor)imp.getProcessor()).getPixels();
		float[] trajectory=new float[height];
		for(int i=0;i<height;i++){
			float[] temp=new float[width];
			System.arraycopy(data,i*width,temp,0,width);
			float maxval=temp[0];
			int maxx=0;
			for(int j=1;j<width;j++){
				if(temp[j]>maxval){
					maxx=j;
					maxval=temp[j];
				}
			}
			if(!interp) trajectory[i]=(float)maxx;
			else trajectory[i]=interpolation.get_local_max1D(temp,maxx);
		}
		new PlotWindow4("Max Column","row","Max Column",trajectory).draw();
	}

}
