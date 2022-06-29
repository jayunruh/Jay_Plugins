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
import jguis.*;

public class pair_correlation_function_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("raw histogram",false);
		gd.addCheckbox("2D_normalization (otherwise 3D)",true);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean rawhist=gd.getNextBoolean();
		boolean norm2d=gd.getNextBoolean();
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals1=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals1=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[] xvals=null; float[] yvals=null; float[] zvals=null;
		if(iw.getClass().getName().equals("jguis.PlotWindow3D")){
			int[] npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
			xvals=make1D(xvals1,npts);
			yvals=make1D(yvals1,npts);
			float[][][] zvals1=(float[][][])jutils.runPW4VoidMethod(iw,"getZValues");
			zvals=make1D(zvals1[0],npts);
		} else {
			int[] npts=((int[])jutils.runPW4VoidMethod(iw,"getNpts"));
			xvals=make1D(xvals1,npts);
			yvals=make1D(yvals1,npts);
			zvals=new float[yvals.length];
		}
		int length=yvals.length;
		float maxx=maxval(xvals);
		float minx=minval(xvals);
		float maxy=maxval(yvals);
		float miny=minval(yvals);
		float maxz=maxval(zvals);
		float minz=minval(zvals);
		
		int maxr=(int)getdist(minx,miny,minz,maxx,maxy,maxz);
		int rbinsize=1;
		int histlength=maxr/rbinsize+1;
		float[] rhist=new float[histlength];
		for(int i=0;i<length;i++){
			float x=xvals[i];
			float y=yvals[i];
			float z=zvals[i];
			for(int j=(i+1);j<length;j++){
				float r=getdist(x,y,z,xvals[j],yvals[j],zvals[j]);
				//float r=(float)Math.sqrt((x-xvals[j])*(x-xvals[j])+(y-yvals[j])*(y-yvals[j]));
				int histval=(int)(r/(float)rbinsize);
				if(histval<histlength){
					rhist[histval]+=2.0f; //count 2 for the donor and acceptor pair
				}
			}
		}
		float[] rvals=new float[histlength];
		for(int i=0;i<histlength;i++){
			if(!rawhist){
				if(norm2d){
					//normalize by the volume of each bin band (pcf will be proportional to particle density)
					rhist[i]/=(float)(length-1)*2.0f*(float)Math.PI*((float)i+0.5f)*(float)rbinsize; 
				} else {
					//normalize by the volume of each bin shell (pcf will be proportional to particle density)
					rhist[i]/=(float)(length-1)*4.0f*(float)Math.PI*(((float)i+0.5f)*((float)i+0.5f)*(float)rbinsize+0.08333f*(float)rbinsize*(float)rbinsize*(float)rbinsize); 
				}
			}
			rvals[i]=(float)(rbinsize*i);
		}
		new PlotWindow4("Pair Correlation Function","r","G(r)",rvals,rhist).draw();
	}

	public float getdist(float x1,float y1,float z1,float x2,float y2,float z2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	}

	public float[] make1D(float[][] data,int[] npts){
		int totpts=0;
		for(int i=0;i<npts.length;i++) totpts+=npts[i];
		float[] newdata=new float[totpts];
		int counter=0;
		for(int i=0;i<npts.length;i++){
			System.arraycopy(data[i],0,newdata,counter,npts[i]);
			counter+=npts[i];
		}
		return newdata;
	}

	public float minval(float[] data){
		float min=data[0];
		for(int i=1;i<data.length;i++){
			if(data[i]<min){min=data[i];}
		}
		return min;
	}

	public float maxval(float[] data){
		float max=data[0];
		for(int i=1;i<data.length;i++){
			if(data[i]>max){max=data[i];}
		}
		return max;
	}

}
