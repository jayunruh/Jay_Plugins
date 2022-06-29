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

public class trajectories_2_flowmap_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[] limits=(float[])jutils.runPW4VoidMethod(iw,"getLimits");
		float[][] xdata=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] ydata=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		String[] annotations=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[] framestarts=new float[npts.length];
		float[] frameends=new float[npts.length];
		float maxend=0f;
		for(int i=0;i<annotations.length;i++){
			if(annotations!=null) framestarts[i]=Float.parseFloat(annotations[i]);
			frameends[i]=framestarts[i]+(float)(npts[i]-1);
			if(frameends[i]>maxend) maxend=frameends[i];
		}
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X_Start",limits[0],5,10,null);
		gd.addNumericField("X_End",limits[1],5,10,null);
		gd.addNumericField("Y_Start",limits[2],5,10,null);
		gd.addNumericField("Y_End",limits[3],5,10,null);
		gd.addNumericField("Frame_Start",0,5,15,null);
		gd.addNumericField("Frame_End",maxend-1,5,15,null);
		gd.addNumericField("Box_Size",32,5,10,null);
		gd.addNumericField("Box_Spacing",16,5,10,null);
		gd.addNumericField("Min_Trajectories (per box)",5,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		float xmin=(float)gd.getNextNumber();
		float xmax=(float)gd.getNextNumber();
		float ymin=(float)gd.getNextNumber();
		float ymax=(float)gd.getNextNumber();
		float tmin=(float)gd.getNextNumber();
		float tmax=(float)gd.getNextNumber();
		float boxsize=(float)gd.getNextNumber();
		float stepsize=(float)gd.getNextNumber();
		int minnum=(int)gd.getNextNumber();
		int xsteps=(int)(Math.abs(xmax-xmin)/stepsize);
		int ysteps=(int)(Math.abs(ymax-ymin)/stepsize);
		float[][] dx=new float[xsteps][ysteps];
		float[][] dy=new float[xsteps][ysteps];
		FloatProcessor cp=new FloatProcessor((int)stepsize*xsteps,(int)stepsize*ysteps);
		cp.setLineWidth(2);
		float magthresh=0.0f;
		for(int i=0;i<ysteps;i++){
			float ystart=(float)i*stepsize;
			float yend=ystart+boxsize;
			for(int j=0;j<xsteps;j++){
				float xstart=(float)j*stepsize;
				float xend=xstart+boxsize;
				int count=0;
				for(int k=0;k<npts.length;k++){
					for(int l=1;l<npts[k];l++){
						float t=(float)(framestarts[k]+l);
						if(xdata[k][l]>=xstart && xdata[k][l]<xend && ydata[k][l]>=ystart && ydata[k][l]<yend && t>=tmin && t<=tmax){
							count++;
							dx[j][i]+=(xdata[k][l]-xdata[k][l-1]);
							dy[j][i]+=(ydata[k][l]-ydata[k][l-1]);
						}
					}
				}
				if(count>minnum){
					dx[j][i]/=(float)count;
					dy[j][i]/=(float)count;
				} else {
					dx[j][i]=0.0f;
					dy[j][i]=0.0f;
				}
				float mag=(float)Math.sqrt(dx[j][i]*dx[j][i]+dy[j][i]*dy[j][i]);
				if(mag>magthresh){
					cp.setValue(mag);
					float xnorm=boxsize*0.5f*dx[j][i]/mag;
					float ynorm=boxsize*0.5f*dy[j][i]/mag;
					jutils.draw_arrow(cp,(int)(xstart+0.5f*stepsize-0.5f*xnorm),(int)(ystart+0.5f*stepsize-0.5f*ynorm),(int)(xstart+0.5f*stepsize+0.5f*xnorm),(int)(ystart+0.5f*stepsize+0.5f*ynorm));
				}
			}
		}
		new ImagePlus("Flow Map",cp).show();
		ImageStack stack=new ImageStack(xsteps,ysteps);
		float[] tempdx=new float[xsteps*ysteps];
		float[] tempdy=new float[xsteps*ysteps];
		for(int i=0;i<ysteps;i++){
			for(int j=0;j<xsteps;j++){
				tempdx[j+i*xsteps]=dx[j][i];
				tempdy[j+i*xsteps]=dy[j][i];
			}
		}
		stack.addSlice("",tempdx);
		stack.addSlice("",tempdy);
		new ImagePlus("Velocities",stack).show();
		IJ.log("dx");
		IJ.log(table_tools.print_float_array(dx));
		IJ.log("dy");
		IJ.log(table_tools.print_float_array(dy));
	}

}
