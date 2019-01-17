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
import jalgs.jfit.*;
import ij.text.*;

public class traj_MAD_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Start_tau",0,0);
		gd.addNumericField("End_tau",5,0);
		gd.addNumericField("Min_length",10,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int start=(int)gd.getNextNumber();
		int end=(int)gd.getNextNumber();
		int minlength=(int)gd.getNextNumber();
		boolean[] valid=new boolean[yvals.length];
		int nvalid=0;
		for(int i=0;i<yvals.length;i++){
			if(npts[i]>=minlength){nvalid++; valid[i]=true;}
		}
		float[][] newyvals=new float[nvalid][end-start+1];
		float[][] newxvals=new float[nvalid][end-start+1];
		int counter=0;
		for(int j=0;j<yvals.length;j++){
			if(valid[j]){
				for(int i=start;i<=end;i++){
					for(int k=0;k<(npts[j]-i-1);k++){
						//calculate the angle between the vectors
						float[] vec1={xvals[j][k+1]-xvals[j][k],yvals[j][k+1]-yvals[j][k]}; normvec(vec1);
						float[] vec2={xvals[j][k+i+1]-xvals[j][k+i],yvals[j][k+i+1]-yvals[j][k+i]}; normvec(vec2);
						float angle=vecangle(vec1,vec2);
						newyvals[counter][i-start]+=angle*angle;
						//newyvals[counter][i-start]+=angle;
					}
					newxvals[counter][i-start]=(float)i;
					newyvals[counter][i-start]/=(float)(npts[j]-i-1);
				}
				counter++;
			}
		}
		new PlotWindow4("MAD Plot","tau","MAD",newxvals,newyvals,null).draw();
	}

	public void normvec(float[] vec){
		float length=(float)Math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]);
		vec[0]/=length;
		vec[1]/=length;
	}

	public float vecangle(float[] vec1,float[] vec2){
		float angle=(float)Math.acos(vec1[0]*vec2[0]+vec1[1]*vec2[1]);
		if(Float.isNaN(angle) || Float.isInfinite(angle)){return 0.0f;}
		return angle;
		//return (float)(vec1[0]*vec2[0]+vec1[1]*vec2[1]);
	}

}
