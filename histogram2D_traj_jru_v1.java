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

public class histogram2D_traj_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		String ylabel=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("x_axis_series",1,0);
		gd.addNumericField("y_axis_series",2,0);
		gd.addCheckbox("all_series_interleaved",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int ser1=(int)gd.getNextNumber()-1;
		int ser2=(int)gd.getNextNumber()-1;
		boolean allser=gd.getNextBoolean();
		if(allser){
			int totpts=0;
			for(int i=0;i<npts.length;i+=2){
				totpts+=npts[i];
			}
			float[] data1=new float[totpts];
			float[] data2=new float[totpts];
			int offset=0;
			for(int i=0;i<npts.length;i+=2){
				System.arraycopy(yvals[i],0,data1,offset,npts[i]);
				System.arraycopy(yvals[i+1],0,data2,offset,npts[i]);
				offset+=npts[i];
			}
			new PlotWindow2DHist("Histogram","Channel1","Channel2",data1,data2,null).draw();
		} else {
			new PlotWindow2DHist("Histogram","Channel"+(ser1+1),"Channel"+(ser2+1),yvals[ser1],yvals[ser2],null).draw();
		}
	}

}
