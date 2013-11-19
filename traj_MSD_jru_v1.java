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

public class traj_MSD_jru_v1 implements PlugIn {
	float startalpha,endalpha;

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][][] zvals=null;
		boolean threeD=iw.getClass().getName().equals("jguis.PlotWindow3D");
		if(threeD) zvals=(float[][][])jutils.runPW4VoidMethod(iw,"getZValues");
		int[] npts;
		if(threeD) npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
		else npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Start_tau",0,0);
		gd.addNumericField("End_tau",5,0);
		gd.addNumericField("Min_length",10,0);
		gd.addCheckbox("Fit_MSD",true);
		gd.addCheckbox("Cross_Corr",false);
		gd.addCheckbox("One_D",false);
		startalpha=0.01f;
		gd.addNumericField("Min_Alpha",startalpha,5,15,null);
		endalpha=2.0f;
		gd.addNumericField("Max_Alpha",endalpha,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int start=(int)gd.getNextNumber();
		int end=(int)gd.getNextNumber();
		int minlength=(int)gd.getNextNumber();
		boolean fit=gd.getNextBoolean();
		boolean cc=gd.getNextBoolean();
		boolean oneD=gd.getNextBoolean();
		startalpha=(float)gd.getNextNumber();
		endalpha=(float)gd.getNextNumber();
		boolean[] valid=new boolean[yvals.length];
		int nvalid=0;
		if(!cc){
			for(int i=0;i<yvals.length;i++){
				if(npts[i]>=minlength){nvalid++; valid[i]=true;}
			}
		} else {
			for(int i=0;i<yvals.length;i+=2){
				if(npts[i]>=minlength && npts[i+1]>=minlength){nvalid++; valid[i/2]=true;}
			}
		}
		float[][] newyvals=new float[nvalid][end-start+1];
		float[][] newxvals=new float[nvalid][end-start+1];
		int counter=0;
		if(!cc){
			if(!oneD){
				for(int j=0;j<yvals.length;j++){
					if(valid[j]){
						for(int i=start;i<=end;i++){
							for(int k=0;k<(npts[j]-i);k++){
								newyvals[counter][i-start]+=(xvals[j][k+i]-xvals[j][k])*(xvals[j][k+i]-xvals[j][k])+(yvals[j][k+i]-yvals[j][k])*(yvals[j][k+i]-yvals[j][k]);
								if(zvals!=null) newyvals[counter][i-start]+=(zvals[0][j][k+i]-zvals[0][j][k])*(zvals[0][j][k+i]-zvals[0][j][k]);
							}
							newxvals[counter][i-start]=(float)i;
							newyvals[counter][i-start]/=(float)(npts[j]-i);
						}
						counter++;
					}
				}
			} else {
				for(int j=0;j<yvals.length;j++){
					if(valid[j]){
						for(int i=start;i<=end;i++){
							for(int k=0;k<(npts[j]-i);k++){
								newyvals[counter][i-start]+=(yvals[j][k+i]-yvals[j][k])*(yvals[j][k+i]-yvals[j][k]);
							}
							newxvals[counter][i-start]=(float)i;
							newyvals[counter][i-start]/=(float)(npts[j]-i);
						}
						counter++;
					}
				}
			}
		} else {
			for(int j=0;j<yvals.length/2;j++){
				if(valid[j]){
					for(int i=start;i<=end;i++){
						for(int k=0;k<(npts[2*j]-i);k++){
							newyvals[counter][i-start]+=(xvals[2*j][k+i]-xvals[2*j][k])*(xvals[2*j+1][k+i]-xvals[2*j+1][k])+(yvals[2*j][k+i]-yvals[2*j][k])*(yvals[2*j][k+i]-yvals[2*j+1][k]);
							if(zvals!=null) newyvals[counter][i-start]+=(zvals[0][2*j][k+i]-zvals[0][2*j][k])*(zvals[0][2*j+1][k+i]-zvals[0][2*j+1][k]);
						}
						newxvals[counter][i-start]=(float)i;
						newyvals[counter][i-start]/=(float)(npts[j]-i);
					}
					counter++;
				}
			}
		}
		new PlotWindow4("MSD Plot","tau","MSD",newxvals,newyvals,null).draw();
		if(fit){
			StringBuffer sb=new StringBuffer();
			for(int i=0;i<newyvals.length;i++){
				float[] params=fit_MSD2(newxvals[i],newyvals[i]);
				sb.append(""+params[0]+"\t"+params[1]+"\n");
				IJ.showProgress(i,newyvals.length);
			}
			new TextWindow("MSD Fits","alpha\tD",sb.toString(),400,200);
		}
	}

	public float[] fit_MSD(float[] xvals,float[] yvals){
		//start by removing zero point if necessary
		float[] newxvals=null;
		float[] newyvals=null;
		if(xvals[0]==0.0f){
			newxvals=new float[xvals.length-1];
			System.arraycopy(xvals,1,newxvals,0,xvals.length-1);
			newyvals=new float[xvals.length-1];
			System.arraycopy(yvals,1,newyvals,0,xvals.length-1);
		} else {
			newxvals=xvals.clone();
			newyvals=yvals.clone();
		}
		for(int i=0;i<newxvals.length;i++){
			newxvals[i]=(float)Math.log(newxvals[i]);
			newyvals[i]=(float)Math.log(newyvals[i]);
		}
		float[] params=(new linleastsquares()).get_amp_offset(newxvals,newyvals,true);
		params[1]=(float)Math.exp(params[1]);
		params[1]/=4.0f;
		return params;
	}

	public float[] fit_MSD2(float[] xvals,float[] yvals){
		boolean starting=true;
		float minalpha=0.0f; float minamp=0.0f;
		double minc2=0.0;
		for(float alpha=startalpha;alpha<=endalpha;alpha+=0.01f){
			float[] func=new float[xvals.length];
				for(int i=0;i<xvals.length;i++){
					func[i]=(float)Math.pow(xvals[i],alpha);
				}
			float[] params=(new linleastsquares()).get_amp_offset(func,yvals,false);
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,params);
			if(starting){
				minc2=c2; minalpha=alpha; minamp=params[0];
				starting=false;
			} else {
				if(c2<minc2){
					minc2=c2; minalpha=alpha; minamp=params[0];
				}
			}
		}
		float[] params={minalpha,minamp/4.0f,(float)minc2};
		return params;
	}

}
