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
import jalgs.*;
import jalgs.jfft.*;

public class traj_pearson_spatial_cc_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float psize=xvals2[0][1]-xvals2[0][0];
		//int length=yvals[0].length;
		//crosscorr cc=new crosscorr(length,false);
		int npairs=yvals.length/2;
		float[][] tics=new float[npairs][];
		float[][] txvals=new float[npairs][];
		int[] tnpts=new int[npairs];
		for(int j=0;j<npairs;j++){
			float[] temp1=(float[])algutils.get_subarray(yvals[j*2],0,npts[j*2]);
			float[] temp2=(float[])algutils.get_subarray(yvals[j*2+1],0,npts[j*2+1]);
			float[][] temp=(new crosscorr(npts[j*2],false)).docrosscorrnofft(temp1,temp2,false);
			//float[][] temp=cc.docrosscorrnofft(yvals[j*2],yvals[j*2+1],false);
			float stdev1=jstatistics.getstatistic("StDev",temp1,null);
			float stdev2=jstatistics.getstatistic("StDev",temp2,null);
			float[] ics=new float[npts[j*2]];
			float[] xvals=new float[npts[j*2]];
			for(int i=0;i<npts[j*2];i++){
				int position=i+npts[j*2]/2;
				if(position>=npts[j*2]) position-=npts[j*2];
				ics[i]=temp[0][position]*temp[1][0]*temp[1][1]/(stdev1*stdev2);
				xvals[i]=psize*(float)(i-npts[j*2]/2);
			}
			tics[j]=ics;
			txvals[j]=xvals;
			tnpts[j]=npts[j*2];
		}
		if(npairs==1) new PlotWindow4("Spatial Correlation","shift","Pearson",txvals[0],tics[0]).draw();
		else{
			//clip these plots so they all have the same x axis
			float maxmin=txvals[0][0];
			int minpts=txvals[0].length;
			for(int i=1;i<npairs;i++){
				if(txvals[i][0]>maxmin){
					maxmin=txvals[i][0];
					minpts=txvals[i].length;
				}
			}
			//IJ.log(""+maxmin+" , "+minpts);
			float[][] tics2=new float[npairs][minpts];
			float[][] txvals2=new float[npairs][minpts];
			for(int i=0;i<npairs;i++){
				int off=Math.round((maxmin-txvals[i][0])/psize);
				System.arraycopy(tics[i],off,tics2[i],0,minpts);
				System.arraycopy(txvals[i],off,txvals2[i],0,minpts);
			}
			new PlotWindow4("Spatial Correlation","shift","Pearson",txvals2,tics2,null).draw();
		}
	}

}
