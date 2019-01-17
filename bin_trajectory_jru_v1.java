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

public class bin_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		int binby=2;
		gd.addNumericField("Bin by",binby,0);
		boolean cum=true;
		gd.addCheckbox("Cumulative",cum);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		binby=(int)gd.getNextNumber();
		cum=gd.getNextBoolean();
		ImageWindow iw=WindowManager.getCurrentWindow();
		String title=(String)jutils.runPW4VoidMethod(iw,"getPlotTitle");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		int length=yvals[0].length;
		int newlength=(int)((float)length/(float)binby);
		float[][] newyvals=new float[yvals.length][newlength];
		float[][] newxvals=new float[yvals.length][newlength];
		for(int k=0;k<yvals.length;k++){
			//double xinc=((double)xvals[k][1]-(double)xvals[k][0])*(double)binby;
			//double xstart=(double)xvals[k][0];
			for(int i=0;i<newlength;i++){
				for(int j=0;j<binby;j++){
					newyvals[k][i]+=yvals[k][i*binby+j];
					newxvals[k][i]+=xvals[k][i*binby+j];
				}
				 if(!cum){newyvals[k][i]/=(float)binby;}
				newxvals[k][i]/=(float)binby;
				//newxvals[k][i]=(float)(xinc*(double)i+xstart);
			}
		}
		new PlotWindow4(title+" Binned by "+binby,"x","y",newxvals,newyvals,null).draw();
	}

}
