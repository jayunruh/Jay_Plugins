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

public class smooth_trajectory_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Smooth_x",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean smoothx=gd.getNextBoolean();
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		for(int i=0;i<yvals.length;i++){
			int templength=yvals[i].length;
			float[] temp=new float[templength];
			temp[0]=0.5f*(yvals[i][0]+yvals[i][1]);
			for(int j=1;j<(templength-1);j++){
				temp[j]=(yvals[i][j-1]+yvals[i][j]+yvals[i][j+1])/3.0f;
			}
			temp[templength-1]=0.5f*(yvals[i][templength-2]+yvals[i][templength-1]);
			System.arraycopy(temp,0,yvals[i],0,templength);
			if(smoothx){
				float[] temp2=new float[templength];
				temp2[0]=0.5f*(xvals[i][0]+xvals[i][1]);
				for(int j=1;j<(templength-1);j++){
					temp2[j]=(xvals[i][j-1]+xvals[i][j]+xvals[i][j+1])/3.0f;
				}
				temp2[templength-1]=0.5f*(xvals[i][templength-2]+xvals[i][templength-1]);
				System.arraycopy(temp2,0,xvals[i],0,templength);
			}
		}
		jutils.runPW4VoidMethod(iw,"yautoscale");
		if(smoothx) jutils.runPW4VoidMethod(iw,"xautoscale");
	}

}
