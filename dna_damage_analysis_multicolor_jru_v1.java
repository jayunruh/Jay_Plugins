/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.jseg.*;
import ij.text.*;

public class dna_damage_analysis_multicolor_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		int maxpts=(int)jstatistics.getstatistic("Max",npts,null);
		int bleachchan=0;
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addNumericField("Bleach Channel",bleachchan+1,0);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		bleachchan=(int)gd2.getNextNumber()-1;
		int nsets=npts.length/6;
		float[][] rtxvals1=new float[nsets][maxpts];
		float[][] rtvals1=new float[nsets][maxpts];
		int[] rtnpts1=new int[nsets];
		float[][] rtxvals2=new float[nsets][maxpts];
		float[][] rtvals2=new float[nsets][maxpts];
		int[] rtnpts2=new int[nsets];
		TextWindow tw=jutils.selectTable("DNA Damage Recruitment");
		if(tw==null) tw=new TextWindow("DNA Damage Recruitment","title\tch1initial\tch1maxrt\tch2initial\tch2maxrt","",400,200);
		//note that st stands for stripe and off is for areas off the damage stripe
		for(int i=0;i<npts.length/6;i++){
			int damageindex=findbleach(yvals[i*6+bleachchan],npts[i*6+bleachchan]);
			int predamagestart=damageindex-4;
			IJ.log("set "+i+" damage pos = "+damageindex);
			float stpredam1=getavg(yvals[i*6],npts[i*6],predamagestart,damageindex-1);
			float nucpredam1=getavg(yvals[i*6+4],npts[i*6+4],predamagestart,damageindex-1);
			rtnpts1[i]=npts[i*6];
			float stpredam2=getavg(yvals[i*6+1],npts[i*6+1],predamagestart,damageindex-1);
			float nucpredam2=getavg(yvals[i*6+5],npts[i*6+5],predamagestart,damageindex-1);
			rtnpts2[i]=npts[i*6+1];
			for(int j=0;j<npts[i*6];j++){
				rtxvals1[i][j]=j-damageindex-1;
				rtvals1[i][j]=(yvals[i*6][j]/stpredam1)/(yvals[i*6+4][j]/nucpredam1);
				rtxvals2[i][j]=j-damageindex-1;
				rtvals2[i][j]=(yvals[i*6+1][j]/stpredam2)/(yvals[i*6+5][j]/nucpredam2);
			}
			float maxrt1=get_3x_smoothed_max(rtvals1[i],rtnpts1[i],damageindex+1);
			float maxrt2=get_3x_smoothed_max(rtvals2[i],rtnpts2[i],damageindex+1);
			tw.append(iw.getTitle()+"-"+(i+1)+"\t"+stpredam1+"\t"+maxrt1+"\t"+stpredam2+"\t"+maxrt2);
		}
		new PlotWindow4("Rt_profiles_ch1","time","R(t)",rtxvals1,rtvals1,rtnpts1).draw();
		new PlotWindow4("Rt_profiles_ch2","time","R(t)",rtxvals2,rtvals2,rtnpts2).draw();
	}

	public float get_3x_smoothed_max(float[] traj,int npts,int startpt){
		float[] temp=(float[])algutils.get_subarray(traj,0,npts);
		temp=jsmooth.smooth1D(temp);
		temp=jsmooth.smooth1D(temp);
		temp=jsmooth.smooth1D(temp);
		float max=temp[startpt];
		for(int i=startpt;i<temp.length;i++) if(temp[i]>max) max=temp[i];
		return max;
	}

	public float[] getregion(float[] yvals,int npts,int start,int length){
		float[] temp=new float[length];
		for(int i=start;i<(start+length);i++){
			float val=yvals[0];
			if(i>=0 && i<npts) val=yvals[i];
			if(i>=npts) val=yvals[npts-1];
			temp[i-start]=val;
		}
		return temp;
	}

	public float getavg(float[] yvals,int npts,int start,int end){
		int start2=start;
		if(start2<0) start2=0;
		int end2=end;
		if(end2>=npts) end2=(npts-1);
		float sum=0.0f;
		int count=0;
		for(int i=start2;i<=end2;i++){
			sum+=yvals[i];
			count++;
		}
		return sum/(float)count;
	}

	public int findbleach(float[] yvals,int npts){
		float[] temp=(float[])algutils.get_subarray(yvals,0,npts);
		//temp=jsmooth.smooth1D(temp);
		float minder=0.0f;
		int minindex=0;
		float prevval=temp[0];
		float[] ders=new float[npts-1];
		for(int i=1;i<npts;i++){
			float derivative=temp[i]-prevval;
			prevval=temp[i];
			ders[i-1]=derivative;
			if(derivative<minder){
				minder=derivative;
				minindex=i;
			}
		}
		//new PlotWindow4("derivative","x","y",ders).draw();
		return minindex;
	}

}
