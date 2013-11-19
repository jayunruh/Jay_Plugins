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

public class normalize_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		String[] normoptions={"Max","Integral","Zero_Point","Min_Max"};
		gd.addChoice("Normalization",normoptions,normoptions[0]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int normindex=gd.getNextChoiceIndex();

		ImageWindow iw=WindowManager.getCurrentWindow();
		PlotWindow4 pw=jutils.getPW4Copy(iw);
		float[][] yvals=pw.getYValues();
		float[][] xvals=pw.getXValues();
		int[] npts=pw.getNpts();
		String title=pw.getPlotTitle();
		String xlabel=pw.getxLabel();
		String ylabel=pw.getyLabel();
		int nseries=pw.getNSeries();
		int maxpts=yvals[0].length;
		float[][] newyvals=new float[nseries][maxpts];

		if(normindex==0){
			for(int i=0;i<nseries;i++){
				float max=yvals[i][0];
				for(int j=1;j<maxpts;j++){
					if(yvals[i][j]>max){max=yvals[i][j];}
				}
				for(int j=0;j<maxpts;j++){
					yvals[i][j]/=max;
				}
			}
		} else {
			if(normindex==1){
				for(int i=0;i<nseries;i++){
					float integral=0.0f;
					for(int j=0;j<maxpts;j++){
						integral+=yvals[i][j];
					}
					for(int j=0;j<maxpts;j++){
						yvals[i][j]/=integral;
					}
				}
			} else {
				if(normindex==2){
					for(int i=0;i<nseries;i++){
						float minx=(float)Math.abs(xvals[i][0]);
						float minxy=yvals[i][0];
						for(int j=0;j<npts[i];j++){
							if(Math.abs(xvals[i][j])<minx){
								minx=(float)Math.abs(xvals[i][j]);
								minxy=yvals[i][j];
							}
						}
						for(int j=0;j<maxpts;j++){
							yvals[i][j]/=minxy;
						}
					}
				} else {
					for(int i=0;i<nseries;i++){
						float max=yvals[i][0];
						float min=yvals[i][0];
						for(int j=1;j<npts[i];j++){
							if(yvals[i][j]>max) max=yvals[i][j];
							if(yvals[i][j]<min) min=yvals[i][j];
						}
						for(int j=0;j<npts[i];j++){
							yvals[i][j]=(yvals[i][j]-min)/(max-min);
						}
					}
				}
			}
		}
		pw.yautoscale();
	}
}
