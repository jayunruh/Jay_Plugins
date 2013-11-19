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
import ij.text.*;
import ij.util.*;
import jguis.*;
import ij.measure.*;

public class plot_line_profiles_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		int width = imp.getWidth();
		int height=imp.getHeight();
		FloatProcessor fp = (FloatProcessor)imp.getProcessor().convertToFloat();
		float[] pixels = (float[])fp.getPixels();
		float[] minmax=getminmax(pixels);
		GenericDialog gd = new GenericDialog("Options");
		boolean vertline=false;
		gd.addCheckbox("Vertical?",vertline);
		boolean movie=false;
		gd.addCheckbox("Plot_Stack?",movie);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		vertline=gd.getNextBoolean();
		movie=gd.getNextBoolean();
		float pwidth=(float)imp.getCalibration().pixelWidth;
		float pheight=(float)imp.getCalibration().pixelHeight;
		if(!vertline){
			if(movie){
				float[][][] linevals=new float[height][1][width];
				float[][][] xvals=new float[height][1][width];
				for(int i=0;i<height;i++){
					for(int j=0;j<width;j++){
						xvals[i][0][j]=pwidth*j;
						linevals[i][0][j]=pixels[i*width+j];
					}
				}
				PlotStack4 plot = new PlotStack4("Profiles","x","Intensity",xvals,linevals);
				plot.draw();
				plot.setLimitsAll(0,pwidth*(width-1),minmax[0],minmax[1]);	
			} else {
				float[][] linevals=new float[height][width];
				float[][] xvals=new float[height][width];
				for(int i=0;i<height;i++){
					for(int j=0;j<width;j++){
						xvals[i][j]=pwidth*j;
						linevals[i][j]=pixels[i*width+j];
					}
				}
				PlotWindow4 plot = new PlotWindow4("Profiles","x","Intensity",xvals,linevals,null);
				plot.draw();
			}
		}
		else{
			if(movie){
				float[][][] linevals=new float[width][1][height];
				float[][][] xvals=new float[width][1][height];
				for(int i=0;i<width;i++){
					for(int j=0;j<height;j++){
						xvals[i][0][j]=pheight*j;
						linevals[i][0][j]=pixels[j*width+i];
					}
				}
				PlotStack4 plot = new PlotStack4("Profiles","intensity","pixel",xvals,linevals);
				plot.draw();
				plot.setLimitsAll(0,pheight*(height-1),minmax[0],minmax[1]);
			} else {
				float[][] linevals=new float[width][height];
				float[][] xvals=new float[width][height];
				for(int i=0;i<width;i++){
					for(int j=0;j<height;j++){
						xvals[i][j]=pheight*j;
						linevals[i][j]=pixels[j*width+i];
					}
				}
				PlotWindow4 plot = new PlotWindow4("Profiles","intensity","pixel",xvals,linevals,null);
				plot.draw();
			}
		}
	}

	float[] getminmax(float[] data){
		float[] temp={data[0],data[0]};
		for(int i=0;i<data.length;i++){
			if(data[i]<temp[0]){temp[0]=data[i];}
			if(data[i]>temp[1]){temp[1]=data[i];}
		}
		return temp;
	}

	void floatarray2text(float[] data){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<data.length;i++){
			sb.append(""+data[i]+"\n");
		}
		TextWindow tw = new TextWindow("Profile Values",sb.toString(),200,400);
	}

}
