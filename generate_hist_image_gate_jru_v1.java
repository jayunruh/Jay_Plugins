/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;

public class generate_hist_image_gate_jru_v1 implements PlugIn {
	//this plugin works on a histogram image
	//the y axis is linear
	//the x axis is logarithmic
	//bin size is assumed to be 4 x 4
	//for every bin, we calculated the 99th? percentile in y
	//then we make a polyline selection that follows that profile
	//below a lower limit and above an upper limit, profile is straight

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_X",1.0,5,15,null);
		gd.addNumericField("Max_X",1000.0,5,15,null);
		gd.addNumericField("Min_Y",-0.1,5,15,null);
		gd.addNumericField("Max_Y",0.5,5,15,null);
		gd.addNumericField("Gate_Percentile",99.0,5,15,null);
		gd.addNumericField("Lower_Limit",2.0,5,15,null);
		gd.addNumericField("Upper_Limit",104.0,5,15,null);
		gd.addNumericField("Shift (bins)",2,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float xmin=(float)gd.getNextNumber();
		float xmax=(float)gd.getNextNumber();
		float ymin=(float)gd.getNextNumber();
		float ymax=(float)gd.getNextNumber();
		float gateper=(float)gd.getNextNumber();
		float lowlim=(float)gd.getNextNumber();
		float uplim=(float)gd.getNextNumber();
		int shift=(int)gd.getNextNumber();
		float[] pix=(float[])imp.getProcessor().getPixels();
		//generate the x axis (logarithmic)
		float logxmin=(float)Math.log(xmin);
		float logxmax=(float)Math.log(xmax);
		float logxinc=(logxmax-logxmin)/(float)(width-1);
		float[] xvals=new float[width];
		for(int i=0;i<width;i++){
			float logval=logxmin+logxinc*(float)i;
			xvals[i]=(float)Math.exp(logval);
		}
		//now bin by 4
		float[] xvalsbinned=new float[width/4];
		for(int i=0;i<width;i+=4){
			for(int j=0;j<4;j++){
				xvalsbinned[i/4]+=xvals[i+j];
			}
			xvalsbinned[i/4]*=0.25f;
		}
		//now get the y vals
		float[] yvals=new float[height];
		float yinc=(ymax-ymin)/(float)(height-1);
		for(int i=0;i<height;i++){
			//yvals[i]=ymin+yinc*(float)i;
			yvals[i]=(float)i;
		}
		//and bin them by 4
		float[] yvalsbinned=new float[height/4];
		for(int i=0;i<height;i+=4){
			for(int j=0;j<4;j++){
				yvalsbinned[i/4]+=yvals[i+j];
			}
			yvalsbinned[i/4]*=0.25f;
		}
		//now start on the first full bin (first above lowlim) and end on the last full bin (last-1 above uplim)
		int xpos=0;
		while(xvalsbinned[xpos]<lowlim){
			xpos++;
		}
		int firstbin=xpos;
		while(xpos<xvalsbinned.length && xvalsbinned[xpos]<uplim){
			xpos++;
		}
		int lastbin=xpos-1;
		//now go through and get the gate value for each bin
		float[][] hists=new float[width/4][height/4];
		float[][] dupyvals=new float[width/4][height/4];
		for(int i=0;i<width;i+=4){
			dupyvals[i/4]=yvalsbinned;
			for(int j=(height-1);j>=0;j-=4){
				hists[i/4][(height-j+1)/4]=pix[i+j*width];
			}
		}
		//optionally plot the histograms
		//new PlotWindow4("Histograms","bin","hist",dupyvals,hists,null).draw();
		//and finally find the appropriate percentile for each histogram between firstbin and lastbin
		float[] gatevals=new float[width/4];
		for(int i=firstbin;i<=lastbin;i++){
			gatevals[i]=getHistPercentile(hists[i],yvalsbinned,gateper);
		}
		//now fill in the lower and upper bins
		for(int i=0;i<firstbin;i++) gatevals[i]=gatevals[firstbin];
		for(int i=(lastbin+1);i<gatevals.length;i++) gatevals[i]=gatevals[lastbin];
		//should we smooth this?--yes
		smooth(gatevals,5);
		//now plot the gate profile for reference
		new PlotWindow4("Gate Profile","x","y",xvalsbinned,gatevals).draw();
		//now create the roi
		//need to map the gate values back to original bins
		int nroipoints=width/4+4;
		int[][] coords=new int[2][nroipoints];
		for(int i=0;i<width/4;i++){
			coords[0][i+4]=i*4+2;
			//coords[1][i+4]=height-2-(int)((gatevals[i]-ymin)/yinc);
			coords[1][i+4]=height-7-(int)gatevals[i]-shift;
		}
		coords[0][0]=width-1; coords[1][0]=coords[1][nroipoints-1];
		coords[0][1]=width-1; coords[1][1]=height-1;
		coords[0][2]=0; coords[1][2]=height-1;
		coords[0][3]=0; coords[1][3]=coords[1][4];
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) rman=new RoiManager();
		rman.addRoi(new PolygonRoi(coords[0],coords[1],nroipoints,Roi.POLYGON));
	}

	public void smooth(float[] traj,int width){
		//does a y smooth on the trajectory
		//width should be an odd number
		float[] traj2=traj.clone();
		int shift=width/2;
		for(int i=shift;i<(traj.length-shift-1);i++){
			traj2[i]=0.0f;
			for(int j=(i-shift);j<(i-shift+width);j++){
				traj2[i]+=traj[j];
			}
			traj2[i]/=(float)width;
		}
		for(int i=0;i<traj.length;i++) traj[i]=traj2[i];
	}

	public float getHistPercentile(float[] hist,float[] axis,float percentile){
		//this returns a hist percentile axis units
		//note that percentile is a percentile, not a fraction
		//get the normalized cumulative histogram
		float sum=jstatistics.getstatistic("Sum",hist,null);
		float[] cum=new float[hist.length];
		cum[0]=hist[0]/sum;
		for(int i=1;i<hist.length;i++){
			cum[i]=cum[i-1]+hist[i]/sum;
		}
		float frac=percentile/100.0f;
		//finally get the percentile
		for(int i=0;i<hist.length;i++){
			if(cum[i]==frac){
				return axis[i];
			} else if(cum[i]>frac){
				float rem=cum[i]-frac;
				float interp=axis[i-1]+rem*(axis[i]-axis[i-1]);
				return interp;
			}
		}
		return axis[hist.length-1];
	}

}
