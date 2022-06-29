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
import java.io.*;
import jalgs.*;
import jguis.*;

public class trajectory_statistics_jru_v2 implements PlugIn {
	//this plugin calculates some statistic for each image in a stack (or region of interest) and plots it
	int index,histbins;
	float histstart,histend;

	public void run(String arg) {
		String[] options=jstatistics.stats;
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Spectrum Statistic?",options,options[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		index=gd.getNextChoiceIndex();
		String statistic=options[index];
		//get the image and its info
		ImageWindow iw=WindowManager.getCurrentWindow();
		int selected=-1;
		float[][] yvals=null;
		int[] npts=null;
		if(jutils.isPlotHist(iw)){
			selected=((Integer)jutils.runPW4VoidMethod(iw,"getSelected")).intValue();
			yvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			npts=new int[]{yvals[0].length};
		} else {
			selected=((Integer)jutils.runPW4VoidMethod(iw,"getSelected")).intValue();
			yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		}
		float[] extras=jutils.getStatsOptions(statistic);
		/*if(statistic=="Mode"){
			extras=gethistoptions();
			if(extras==null){return;}
		}
		if(statistic=="ConditionalAvg"){
			extras=getcondavgoptions();
			if(extras==null){return;}
		}*/
		if(selected>=0 && selected<yvals.length){
			float[] temp=new float[npts[selected]];
			System.arraycopy(yvals[selected],0,temp,0,npts[selected]);
			float stat=jstatistics.getstatistic(statistic,temp,extras);
			IJ.log(statistic+" = \t"+stat);
			if(statistic.equals("Percentile")){
				for(int i=1;i<extras.length;i++) IJ.log(statistic+" = \t"+extras[i]);
			}
		} else {
			String collabels="Trajectory\t"+statistic;
			if(statistic.equals("Percentile")){
				for(int j=1;j<extras.length;j++) collabels=collabels+"\t"+statistic+(j+1);
			}
			TextWindow tw=new TextWindow("Trajectory_Statistics",collabels,"",400,200);
			for(int i=0;i<yvals.length;i++){
				float[] temp=new float[npts[i]];
				System.arraycopy(yvals[i],0,temp,0,npts[i]);
				float stat=jstatistics.getstatistic(statistic,temp,extras);
				StringBuffer sb=new StringBuffer();
				sb.append(""+(i+1)+"\t"+stat);
				//IJ.log(statistic+""+(i+1)+" = \t"+stat);
				if(statistic.equals("Percentile")){
					//for(int j=1;j<extras.length;j++) IJ.log(statistic+""+(i+1)+" = \t"+extras[j]);
					for(int j=1;j<extras.length;j++) sb.append("\t"+extras[j]);
				}
				tw.append(sb.toString());
			}
		}
	}

	private float[] gethistoptions(){
		GenericDialog gd=new GenericDialog("Histogram Options");
		gd.addNumericField("Histogram Bins",histbins,0);
		gd.addNumericField("Histogram Start",histstart,5,10,null);
		gd.addNumericField("Histogram End",histend,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		float[] temp=new float[3];
		histbins=(int)gd.getNextNumber();
		histstart=(float)gd.getNextNumber();
		histend=(float)gd.getNextNumber();
		temp[0]=(float)histbins; temp[1]=histstart; temp[2]=histend;
		return temp;
	}

	private float[] getcondavgoptions(){
		GenericDialog gd=new GenericDialog("Conditional Avg Options");
		gd.addNumericField("Max",100.0,5,10,null);
		gd.addNumericField("Min",0.0,5,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return null;}
		float[] temp=new float[2];
		temp[0]=(float)gd.getNextNumber();
		temp[1]=(float)gd.getNextNumber();
		return temp;
	}

}
