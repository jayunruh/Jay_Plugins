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
import jalgs.*;

public class difference_movie_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Frame_Delay (>0)",1,0);
		gd.addCheckbox("Correct_for_avg_change",true);
		gd.addCheckbox("Abs_value",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int delay=(int)gd.getNextNumber();
		boolean correct=gd.getNextBoolean();
		boolean abs=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		if(frames==1){
			frames=slices;
			slices=1;
		}
		int width=imp.getWidth();
		int height=imp.getHeight();
		int newframes=frames-delay;
		Object[][][] dpix=new Object[slices][channels][];
		for(int i=0;i<slices;i++){
			for(int j=0;j<channels;j++){
				Object[] tseries=jutils.get3DTSeries(stack,i,j,frames,slices,channels);
				float[] avgprofile=null;
				if(correct){
					avgprofile=new float[frames];
					for(int k=0;k<frames;k++) avgprofile[k]=jstatistics.getstatistic("Avg",tseries[k],null);
				}
				Object[] newseries=new Object[newframes];
				float[] source=algutils.convert_arr_float(tseries[0]);
				for(int k=0;k<newframes;k++){
					float[] target=algutils.convert_arr_float(tseries[k+delay]);
					if(correct){
						float mult=avgprofile[k+delay]/avgprofile[k];
						for(int l=0;l<width*height;l++){
							source[l]=(target[l]-source[l]*mult);
							if(abs) source[l]=(float)Math.abs(source[l]);
						}
					} else {
						for(int l=0;l<width*height;l++){
							source[l]=(target[l]-source[l]);
							if(abs) source[l]=(float)Math.abs(source[l]);
						}
					}
					newseries[k]=source;
					if(delay==1) source=target;
					else source=algutils.convert_arr_float(tseries[k+1]);
				}
				dpix[i][j]=newseries;
			}
		}
		ImageStack dstack=new ImageStack(width,height);
		for(int i=0;i<newframes;i++){
			for(int j=0;j<slices;j++){
				for(int k=0;k<channels;k++){
					dstack.addSlice("",dpix[j][k][i]);
				}
			}
		}
		ImagePlus imp2=jutils.create_hyperstack("Difference Image",dstack,imp,newframes,slices,channels);
		imp2.show();
	}

}
