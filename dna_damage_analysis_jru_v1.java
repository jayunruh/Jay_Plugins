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

public class dna_damage_analysis_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		int maxpts=(int)jstatistics.getstatistic("Max",npts,null);
		int nsets=npts.length/3;
		float[][] rtxvals=new float[nsets][maxpts];
		float[][] rtvals=new float[nsets][maxpts];
		int[] rtnpts=new int[nsets];
		int fretlength=10;
		float[][] stfretvals=new float[nsets][fretlength];
		float[][] offfretvals=new float[nsets][fretlength];
		TextWindow tw=jutils.selectTable("DNA Damage Recruitment");
		if(tw==null) tw=new TextWindow("DNA Damage Recruitment","title\tinitial\tmaxrt","",400,200);
		//note that st stands for stripe and off is for areas off the damage stripe
		for(int i=0;i<npts.length/3;i++){
			int damageindex=findbleach(yvals[i*3],npts[i*3]);
			int predamagestart=damageindex-4;
			IJ.log("set "+i+" damage pos = "+damageindex);
			float stpredam=getavg(yvals[i*3],npts[i*3],predamagestart,damageindex-1);
			float nucpredam=getavg(yvals[i*3+2],npts[i*3+2],predamagestart,damageindex-1);
			rtnpts[i]=npts[i*3];
			for(int j=0;j<npts[i*3];j++){
				rtxvals[i][j]=j-damageindex-1;
				rtvals[i][j]=(yvals[i*3][j]/stpredam)/(yvals[i*3+2][j]/nucpredam);
			}
			float maxrt=get_3x_smoothed_max(rtvals[i],rtnpts[i],damageindex+1);
			tw.append(iw.getTitle()+"-"+(i+1)+"\t"+stpredam+"\t"+maxrt);
		}
		new PlotWindow4("Rt_profiles","time","intensity",rtxvals,rtvals,rtnpts).draw();
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
		float minder=0.0f;
		int minindex=0;
		float prevval=yvals[0];
		for(int i=1;i<npts;i++){
			float derivative=yvals[i]-prevval;
			if(derivative<minder){
				minder=derivative;
				minindex=i;
			}
		}
		return minindex;
	}

}
