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

public class plot_line_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		int width = imp.getWidth();
		int height=imp.getHeight();
		FloatProcessor fp = (FloatProcessor)imp.getProcessor();
		float[] pixels = (float[])fp.getPixels();
		GenericDialog gd = new GenericDialog("Options");
		boolean vertline=false;
		gd.addCheckbox("Vertical?",vertline);
		int linenum=(int)((float)height/2.0f);
		gd.addNumericField("Line number",linenum,0);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		vertline=gd.getNextBoolean();
		linenum=(int)gd.getNextNumber();
		Calibration cal=imp.getCalibration();
		float psize=(float)cal.pixelWidth;
		if(vertline){
			psize=(float)cal.pixelHeight;
		}
		if(!vertline){
			float[] linevals=new float[width];
			float[] xvals=new float[width];
			for(int i=0;i<width;i++){
				linevals[i]=pixels[linenum*width+i];
				xvals[i]=psize*(float)i;
			}
			floatarray2text(linevals);
			PlotWindow4 plot = new PlotWindow4("Profile "+linenum,"x","Intensity",xvals,linevals);
			plot.draw();
		}
		else{
			float[] linevals=new float[height];
			float[] xvals=new float[height];
			for(int i=0;i<height;i++){
				linevals[i]=pixels[linenum+i*width];
				xvals[i]=psize*(float)i;
			}
			floatarray2text(linevals);
			PlotWindow4 plot = new PlotWindow4("Profile "+linenum,"x","Intensity",xvals,linevals);
			plot.draw();
		}
	}

	void floatarray2text(float[] data){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<data.length;i++){
			sb.append(""+data[i]+"\n");
		}
		TextWindow tw = new TextWindow("Profile Values",sb.toString(),200,400);
	}

}
