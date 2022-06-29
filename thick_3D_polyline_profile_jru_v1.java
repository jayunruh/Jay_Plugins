/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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

public class thick_3D_polyline_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"3D Profile","Z Stack"});
		if(imps==null) return;
		ImageWindow iw=imps[0].getWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		float[][][] zvals=(float[][][])jutils.runPW4VoidMethod(iw,"getZValues");
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
		if(sel<0) sel=0;
		int[][] npts=(int[][])jutils.runPW4VoidMethod(iw,"getNpts");
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Thickness",5,0);
		gd.addNumericField("Z_Ratio",jutils.get_zratio(imps[1]),5,15,null);
		gd.addNumericField("End Extension (frac of len)",0.3f,5,15,null);
		gd.addCheckbox("Straighten",false);
		gd.addCheckbox("Single Slice",false);
		gd.addCheckbox("Ignore Calibration",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int thickness=(int)gd.getNextNumber();
		float zratio=(float)gd.getNextNumber();
		float ext=(float)gd.getNextNumber();
		boolean straighten=gd.getNextBoolean();
		boolean single=gd.getNextBoolean();
		boolean ignore=gd.getNextBoolean();
		float psize=(float)jutils.get_psize(imps[1]);
		if(ignore) psize=1.0f;
		float[] pxvals=new float[npts[0][sel]];
		float[] pyvals=new float[npts[0][sel]];
		float[] pzvals=new float[npts[0][sel]];
		for(int i=0;i<npts[0][sel];i++){
			pxvals[i]=xvals[sel][i]/psize;
			pyvals[i]=yvals[sel][i]/psize;
			pzvals[i]=zvals[0][sel][i]/psize;
		}
		if(pxvals.length<2){
			//if we have a single point, just draw a horizontal line at 1/3 of the thickness
			float temp=pxvals[0]-(float)thickness/6.0f;
			pxvals=new float[]{temp,temp+(float)thickness/3.0f};
			pyvals=new float[]{pyvals[0],pyvals[0]};
			pzvals=new float[]{pzvals[0],pzvals[0]};
		}
		int width=imps[1].getWidth(); int height=imps[1].getHeight();
		int length=profiler.get3DPolygonLength(pxvals,pyvals,pzvals,false);
		float extlen=ext*(float)length; //extension length in pixel units
		if(extlen>0){
			float len1=profiler.get3DLength(new float[]{pxvals[0],pyvals[0],pzvals[0],pxvals[1],pyvals[1],pzvals[1]});
			float xinc=(pxvals[1]-pxvals[0])/len1;
			float yinc=(pyvals[1]-pyvals[0])/len1;
			float zinc=(pzvals[1]-pzvals[0])/len1;
			pxvals[0]-=xinc*extlen;
			pyvals[0]-=yinc*extlen;
			pzvals[0]-=zinc*extlen;
			int e=pxvals.length-1;
			float len2=profiler.get3DLength(new float[]{pxvals[e],pyvals[e],pzvals[e],pxvals[e-1],pyvals[e-1],pzvals[e-1]});
			xinc=(pxvals[e]-pxvals[e-1])/len2;
			yinc=(pyvals[e]-pyvals[e-1])/len2;
			zinc=(pzvals[e]-pzvals[e-1])/len2;
			pxvals[e]+=xinc*extlen;
			pyvals[e]+=yinc*extlen;
			pzvals[e]+=zinc*extlen;
		}
		Object[] stack=jutils.stack2array(imps[1].getStack());
		int chans=imps[1].getNChannels();
		if(chans==1){
			if(straighten){
				int length10=profiler.get3DPolygonLength(pxvals,pyvals,pzvals,false);
				if(single){
					float[] straightened=profiler.get3DStraightened(stack,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
					new ImagePlus("Thick 3D Straightened",new FloatProcessor(thickness,length10,straightened,null)).show();
				} else {
					float[][] straightened=profiler.get3DThickStraightened(stack,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
					new ImagePlus("Thick 3D Straightened",jutils.array2stack(straightened,thickness,length10)).show();
				}
			} else {
				float[] profile=profiler.get3DThickProfile(stack,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
				float[] xvals2=new float[profile.length];
				float startx=-extlen;
				for(int i=0;i<xvals2.length;i++){
					xvals2[i]=psize*((float)i+startx);
				}
				new PlotWindow4("Thick 3D Profile","position (um)","Intensity",xvals2,profile).draw();
			}
		} else {
			int slices=stack.length/chans;
			if(straighten){
				int length10=profiler.get3DPolygonLength(pxvals,pyvals,pzvals,false);
				if(single){
					float[][] straightened=new float[chans][];
					for(int i=0;i<chans;i++){
						Object[] stack2=algutils.get3DZSeries(stack,0,i,1,slices,chans);
						straightened[i]=profiler.get3DStraightened(stack2,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
					}
					ImagePlus straightimp=new ImagePlus("Thick 3D Straightened",jutils.array2stack(straightened,thickness,length10));
					straightimp.setOpenAsHyperStack(true);
					straightimp.setDimensions(chans,1,1);
					(new CompositeImage(straightimp,CompositeImage.COMPOSITE)).show();
				} else {
					float[][] straightened=new float[chans*thickness][];
					for(int i=0;i<chans;i++){
						Object[] stack2=algutils.get3DZSeries(stack,0,i,1,slices,chans);
						float[][] temp=profiler.get3DThickStraightened(stack2,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
						for(int j=0;j<thickness;j++) straightened[j*chans+i]=temp[j];
					}
					ImagePlus straightimp=new ImagePlus("Thick 3D Straightened",jutils.array2stack(straightened,thickness,length10));
					straightimp.setOpenAsHyperStack(true);
					straightimp.setDimensions(chans,thickness,1);
					(new CompositeImage(straightimp,CompositeImage.COMPOSITE)).show();
				}
			} else {
				Object[] stack2=algutils.get3DZSeries(stack,0,0,1,slices,chans);
				float[] profile=profiler.get3DThickProfile(stack,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
				float[] xvals2=new float[profile.length];
				float startx=-extlen;
				for(int i=0;i<xvals2.length;i++){
					xvals2[i]=psize*((float)i+startx);
				}
				PlotWindow4 pw=new PlotWindow4("Thick 3D Profile","position (um)","Intensity",xvals2,profile);
				pw.draw();
				for(int i=1;i<chans;i++){
					stack2=algutils.get3DZSeries(stack,0,i,1,slices,chans);
					float[] temp=profiler.get3DThickProfile(stack2,width,height,pxvals,pyvals,pzvals,false,thickness,0,zratio);
					pw.addPoints(temp,true);
				}
			}
		}
	}

}
