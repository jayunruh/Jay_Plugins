/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
package jguis;

import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.LookUpTable;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.filter.BackgroundSubtracter;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.ShortProcessor;
import ij.text.TextWindow;
import jalgs.algutils;
import jalgs.interpolation;
import jalgs.jdataio;
import jalgs.jstatistics;
import jalgs.jseg.findblobs3;
import jalgs.jseg.measure_object;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.image.IndexColorModel;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.Vector;

public class jutils{

	public static String formatted_string(double number){
		double absnumber=Math.abs(number);
		if(absnumber>=1000.0||absnumber<0.01){
			DecimalFormat expformat=new DecimalFormat("0.00E0");
			return expformat.format(number);
		}else{
			if(absnumber>=100.0){
				DecimalFormat tempformat=new DecimalFormat("000.0");
				return tempformat.format(number);
			}else{
				if(absnumber>=10.0){
					DecimalFormat tempformat=new DecimalFormat("00.00");
					return tempformat.format(number);
				}else{
					DecimalFormat tempformat=new DecimalFormat("0.000");
					return tempformat.format(number);
				}
			}
		}
	}

	public static boolean isPlotFamily(ImageWindow iw){
		return(isPlot(iw)||isPlotHist(iw)||iw.getClass().getName().equals("jguis.PlotWindow3D") || iw.getClass().getName().equals("jguis.PlotWindowColumn"));
	}
	
	public static boolean isPlotHist(ImageWindow iw){
		if(iw==null) return false;
		return (iw.getClass().getName().equals("jguis.PlotWindowHist") || iw.getClass().getName().equals("jguis.PlotWindow2DHist"));
	}
	
	public static boolean is3DPlot(ImageWindow iw){
		if(iw==null) return false;
		return (iw.getClass().getName().equals("jguis.PlotWindow3D"));
	}

	public static boolean isPlot(ImageWindow iw){
		if(iw==null) return false;
		return(isPW4(iw)||isPW(iw));
	}

	public static boolean isPW4(ImageWindow iw){
		if(iw==null) return false;
		return iw.getClass().getName().equals("jguis.PlotWindow4");
	}

	public static boolean isPW(ImageWindow iw){
		if(iw==null) return false;
		return iw.getClass().getName().equals("ij.gui.PlotWindow");
	}
	
	public static Plot4 plot2Plot4(Plot plot){
		try{
			//need to update this for new data model use allPlotObjects?
			Class<?> temp=plot.getClass();
			//Field datafield=temp.getDeclaredField("storedData");
			//datafield.setAccessible(true);
			//ArrayList data=(ArrayList)datafield.get(plot);
			/*Field plotobjfield=temp.getDeclaredField("allPlotObjects");
			plotobjfield.setAccessible(true);
			Vector<?> allPlotObjects=(Vector<?>)plotobjfield.get(plot);*/
			/*Field maincurvefield=temp.getDeclaredField("getMainCurveObject");
			maincurvefield.setAccessible(true);
			Object po=maincurvefield.get(plot);*/
			/*Field xlabelfield=temp.getDeclaredField("xLabel");
			xlabelfield.setAccessible(true);
			String xlabel=(String)xlabelfield.get(plot);
			Field ylabelfield=temp.getDeclaredField("yLabel");
			ylabelfield.setAccessible(true);
			String ylabel=(String)ylabelfield.get(plot);*/
			String xlabel=plot.getLabel('x');
			String ylabel=plot.getLabel('y');
			if(xlabel.equals("")) xlabel="x";
			if(ylabel.equals("")) xlabel="y";
			double[] limits=plot.getLimits();
			double xmin=limits[0]; double xmax=limits[1]; double ymin=limits[2]; double ymax=limits[3];
			/*Field xminfield=temp.getDeclaredField("xMin");
			xminfield.setAccessible(true);
			double xmin=(Double)xminfield.get(plot);
			Field yminfield=temp.getDeclaredField("yMin");
			yminfield.setAccessible(true);
			double ymin=(Double)yminfield.get(plot);
			Field xmaxfield=temp.getDeclaredField("xMax");
			xmaxfield.setAccessible(true);
			double xmax=(Double)xmaxfield.get(plot);
			Field ymaxfield=temp.getDeclaredField("yMax");
			ymaxfield.setAccessible(true);
			double ymax=(Double)ymaxfield.get(plot);*/
			/*Field ebarsfield=temp.getDeclaredField("errorBars");
			ebarsfield.setAccessible(true);
			float[] ebars=(float[])ebarsfield.get(plot);*/
			//int nseries=data.size()/2;
			float[] xvals=plot.getXValues();
			float[] yvals=plot.getYValues();
			//Plot4 p4=new Plot4(xlabel,ylabel,(float[])data.get(0),(float[])data.get(1));
			Plot4 p4=new Plot4(xlabel,ylabel,xvals,yvals);
			/*for(int i=2;i<data.size();i+=2){
				p4.addPoints((float[])data.get(i),(float[])data.get(i+1),true);
			}
			int maxpts=p4.getmaxpts();
			if(ebars!=null){
				float[][] ebars2=new float[nseries][maxpts];
				System.arraycopy(ebars,0,ebars2[0],0,ebars.length);
				p4.addErrors(ebars2);
			}*/
			p4.setLimits(new float[]{(float)xmin,(float)xmax,(float)ymin,(float)ymax});
			return p4;
		//}catch(NoSuchFieldException e){
		//	IJ.log("no such field exception");
		}catch(IllegalArgumentException e){
			IJ.log("illegal argument exception");
		//}catch(IllegalAccessException e){
		//	IJ.log("illegal access exception");
		}
		return null;
	}

	public static PlotWindow pw42pw(ImageWindow iw){
		if(iw.getClass().getName().equals("jguis.PlotWindow4")){
			float[][] yvals=(float[][])runPW4VoidMethod(iw,"getYValues");
			float[][] xvals=(float[][])runPW4VoidMethod(iw,"getXValues");
			int[] npts=(int[])runPW4VoidMethod(iw,"getNpts");
			String[] labels=(String[])runPW4VoidMethod(iw,"getAllLabels");
			float[] limits=(float[])runPW4VoidMethod(iw,"getLimits");
			int[] shapes=(int[])runPW4VoidMethod(iw,"getShapes");
			int[] colors=(int[])runPW4VoidMethod(iw,"getColors");
			Plot plot=new Plot(labels[0],labels[1],labels[2],(float[])algutils.get_subarray(xvals[0],0,npts[0]),(float[])algutils.get_subarray(yvals[0],0,npts[0]));
			plot.setLimits(limits[0],limits[1],limits[2],limits[3]);
			for(int i=1;i<yvals.length;i++){
				plot.setColor(Plot4.java_colors[i]);
				plot.addPoints((float[])algutils.get_subarray(xvals[i],0,npts[i]),(float[])algutils.get_subarray(yvals[i],0,npts[i]),PlotWindow.LINE);
			}
			plot.setColor(Plot4.java_colors[0]);
			return plot.show();
		}else{
			if(iw.getClass().getName().equals("ij.gui.PlotWindow")){
				return (PlotWindow)iw;
			}else{
				return null;
			}
		}
	}

	public static Plot getPWPlot(PlotWindow pw){
		//return (Plot)getReflectionField(pw,"plot");
		return pw.getPlot();
	}

	public static PlotWindow4 pw2pw4(PlotWindow pw){
		Plot4 p4=plot2Plot4(getPWPlot(pw));
		return new PlotWindow4(pw.getTitle(),p4);
	}

	public static Object runPW4VoidMethod(ImageWindow iw,String method){
		Class<?> temp=iw.getClass();
		Object data=null;
		if(temp.getName().equals("jguis.PlotWindow4")||temp.getName().equals("jguis.PlotWindow3D")||temp.getName().equals("jguis.PlotWindowHist") || temp.getName().equals("jguis.PlotWindow2DHist")){
			data=runReflectionMethod(iw,method,null,null);
		}else if(temp.getName().equals("ij.gui.PlotWindow")){
			Plot4 plot=plot2Plot4(getPWPlot((PlotWindow)iw)); // first convert into a Plot
			if(method.equals("getAllLabels")){
				data=new String[]{iw.getTitle(),plot.getxLabel(),plot.getyLabel()};
			}else{
				data=runReflectionMethod(plot,method,null,null);
			}
		}else{
			;
		}
		return data;
	}

	public static PlotWindow4 getPW4Copy(ImageWindow iw){
		if(iw.getClass().getName().equals("jguis.PlotWindow4")){
			float[][] yvals=(float[][])runPW4VoidMethod(iw,"getYValues");
			float[][] newyvals=new float[yvals.length][yvals[0].length];
			copyfloat2dimvector(yvals,newyvals);
			float[][] xvals=(float[][])runPW4VoidMethod(iw,"getXValues");
			float[][] newxvals=new float[yvals.length][yvals[0].length];
			copyfloat2dimvector(xvals,newxvals);
			int[] npts=(int[])runPW4VoidMethod(iw,"getNpts");
			int[] newnpts=new int[npts.length];
			System.arraycopy(npts,0,newnpts,0,npts.length);
			String[] labels=(String[])runPW4VoidMethod(iw,"getAllLabels");
			String[] newlabels=new String[labels.length];
			copystringarray(labels,newlabels);
			float[] limits=(float[])runPW4VoidMethod(iw,"getLimits");
			int[] shapes=(int[])runPW4VoidMethod(iw,"getShapes");
			int[] colors=(int[])runPW4VoidMethod(iw,"getColors");
			boolean[] logaxes=(boolean[])runPW4VoidMethod(iw,"getLogAxes");
			PlotWindow4 pw=new PlotWindow4(labels[0],labels[1],labels[2],newxvals,newyvals,newnpts);
			pw.draw();
			pw.setLogAxes(logaxes[0],logaxes[1]);
			int[] newshapes=pw.getShapes();
			System.arraycopy(shapes,0,newshapes,0,shapes.length);
			int[] newcolors=pw.getColors();
			System.arraycopy(colors,0,newcolors,0,colors.length);
			pw.setLimits(limits);
			boolean showerrors=(Boolean)runPW4VoidMethod(iw,"getShowErrors");
			if(showerrors){
				float[][][] errors=(float[][][])runPW4VoidMethod(iw,"getErrors");
				if(errors!=null) pw.addErrors(errors);
			}
			String[] annot=(String[])runPW4VoidMethod(iw,"getAnnotations");
			if(annot!=null) pw.getPlot().setAnnotations(annot);
			return pw;
		}else{
			if(iw.getClass().getName().equals("ij.gui.PlotWindow")){
				PlotWindow4 pw2=pw2pw4((PlotWindow)iw);
				pw2.draw();
				return pw2;
			}else{
				if(iw.getClass().getName().equals("jguis.PlotWindowHist")){
					return getPW4SelCopy(iw,0);
				} else {
					return null;
				}
			}
		}
	}

	public static Plot4 getPW4PlotCopy(ImageWindow iw){
		if(iw.getClass().getName().equals("jguis.PlotWindow4")){
			float[][] yvals=(float[][])runPW4VoidMethod(iw,"getYValues");
			float[][] newyvals=new float[yvals.length][yvals[0].length];
			copyfloat2dimvector(yvals,newyvals);
			float[][] xvals=(float[][])runPW4VoidMethod(iw,"getXValues");
			float[][] newxvals=new float[yvals.length][yvals[0].length];
			copyfloat2dimvector(xvals,newxvals);
			int[] npts=(int[])runPW4VoidMethod(iw,"getNpts");
			int[] newnpts=new int[npts.length];
			System.arraycopy(npts,0,newnpts,0,npts.length);
			String[] labels=(String[])runPW4VoidMethod(iw,"getAllLabels");
			String[] newlabels=new String[labels.length];
			copystringarray(labels,newlabels);
			float[] limits=(float[])runPW4VoidMethod(iw,"getLimits");
			int[] shapes=(int[])runPW4VoidMethod(iw,"getShapes");
			int[] colors=(int[])runPW4VoidMethod(iw,"getColors");
			boolean[] logaxes=(boolean[])runPW4VoidMethod(iw,"getLogAxes");
			Plot4 p4=new Plot4(newlabels[1],newlabels[2],newxvals,newyvals,newnpts);
			p4.setLogAxes(logaxes[0],logaxes[1]);
			int[] newshapes=p4.getShapes();
			System.arraycopy(shapes,0,newshapes,0,shapes.length);
			int[] newcolors=p4.getColors();
			System.arraycopy(colors,0,newcolors,0,colors.length);
			p4.setLimits(limits);
			boolean showerrors=(Boolean)runPW4VoidMethod(iw,"getShowErrors");
			if(showerrors){
				float[][][] errors=(float[][][])runPW4VoidMethod(iw,"getErrors");
				if(errors!=null) p4.addErrors(errors);
			}
			String[] annot=(String[])runPW4VoidMethod(iw,"getAnnotations");
			if(annot!=null) p4.setAnnotations(annot);
			return p4;
		}else{
			if(iw.getClass().getName().equals("ij.gui.PlotWindow")){
				PlotWindow pw=(PlotWindow)iw;
				return plot2Plot4(getPWPlot(pw));
			}else{
				return null;
			}
		}
	}

	public static PlotWindow4 getPW4SelCopy(ImageWindow iw){
		int selected=0;
		if(!iw.getClass().getName().equals("ij.gui.PlotWindow")) selected=((Integer)runPW4VoidMethod(iw,"getSelected")).intValue();
		return getPW4SelCopy(iw,selected);
	}

	public static PlotWindow4 getPW4SelCopy(ImageWindow iw,int selected){
		if(iw.getClass().getName().equals("jguis.PlotWindow4")){
			float[][] yvals=(float[][])runPW4VoidMethod(iw,"getYValues");
			if(selected<=0) selected=0;
			if(selected>=yvals.length) selected=0;
			int[] npts=(int[])runPW4VoidMethod(iw,"getNpts");
			float[] newyvals=new float[npts[selected]];
			System.arraycopy(yvals[selected],0,newyvals,0,npts[selected]);
			float[][] xvals=(float[][])runPW4VoidMethod(iw,"getXValues");
			float[] newxvals=new float[npts[selected]];
			System.arraycopy(xvals[selected],0,newxvals,0,npts[selected]);
			String[] labels=(String[])runPW4VoidMethod(iw,"getAllLabels");
			String[] newlabels=new String[labels.length];
			copystringarray(labels,newlabels);
			float[] limits=(float[])runPW4VoidMethod(iw,"getLimits");
			int[] shapes=(int[])runPW4VoidMethod(iw,"getShapes");
			int[] colors=(int[])runPW4VoidMethod(iw,"getColors");
			boolean[] logaxes=(boolean[])runPW4VoidMethod(iw,"getLogAxes");
			PlotWindow4 pw=new PlotWindow4(labels[0]+"-series"+(selected+1),labels[1],labels[2],newxvals,newyvals);
			pw.draw();
			pw.setLogAxes(logaxes[0],logaxes[1]);
			int[] newshapes=pw.getShapes();
			newshapes[0]=shapes[selected];
			int[] newcolors=pw.getColors();
			newcolors[0]=colors[selected];
			pw.setLimits(limits);
			boolean showerrors=(Boolean)runPW4VoidMethod(iw,"getShowErrors");
			if(showerrors){
				float[][][] errors=(float[][][])runPW4VoidMethod(iw,"getErrors");
				if(errors!=null){
					float[][][] errs2=new float[2][1][npts[selected]];
					System.arraycopy(errors[0][selected],0,errs2[0][0],0,npts[selected]);
					System.arraycopy(errors[1][selected],0,errs2[1][0],0,npts[selected]);
					pw.addErrors(errs2);
				}
			}
			return pw;
		}else if(iw.getClass().getName().equals("ij.gui.PlotWindow")){
			Plot4 p4=getPW4PlotCopy(iw);
			PlotWindow pw=(PlotWindow)iw;
			float[] yvals=pw.getYValues();
			float[] xvals=pw.getXValues();
			PlotWindow4 pw2=new PlotWindow4(iw.getTitle()+"-series"+(selected+1),p4.getxLabel(),p4.getyLabel(),xvals,yvals);
			pw2.draw();
			float[] errors=p4.getErrors(0,false);
			if(errors!=null){
				pw2.addErrors(new float[][]{errors});
			}
			return pw2;
		}else if(iw.getClass().getName().equals("jguis.PlotWindowHist")){
			Object plot=runReflectionMethod(iw,"getPlot",null,null);
			float[][] hist=(float[][])runReflectionMethod(plot,"getHistogram",null,null);
			boolean[] logaxes=(boolean[])runReflectionMethod(plot,"getLogAxes",null,null);
			float[] limits=(float[])runReflectionMethod(plot,"getLimits",null,null);
			String[] labels=(String[])runReflectionMethod(iw,"getAllLabels",null,null);
			PlotWindow4 pw=new PlotWindow4(labels[0]+"-series"+0,labels[1],labels[2],hist[0],hist[1]);
			pw.draw();
			pw.setLogAxes(logaxes[0],logaxes[1]);
			pw.setLimits(limits);
			return pw;
		}else{
			return null;
		}
	}
	
	public static PlotWindow3D getPW3DSelCopy(ImageWindow iw,int selected) {
		float[][] yvals=(float[][])runPW4VoidMethod(iw,"getYValues");
		if(selected<=0) selected=0;
		if(selected>=yvals.length) selected=0;
		int[][] npts=(int[][])runPW4VoidMethod(iw,"getNpts");
		float[][][] zvals=(float[][][])runPW4VoidMethod(iw,"getZValues");
		float[][] xvals=(float[][])runPW4VoidMethod(iw,"getXValues");
		Object plot=runPW4VoidMethod(iw,"getPlot");
		String[] labels=(String[])runPW4VoidMethod(iw,"getAllLabels");
		String[] newlabels=new String[labels.length];
		copystringarray(labels,newlabels);
		float[] limits=(float[])runPW4VoidMethod(iw,"getLimits");
		int[] shapes=(int[])runPW4VoidMethod(iw,"getShapes");
		int[] colors=(int[])runPW4VoidMethod(iw,"getColors");
		boolean[] logaxes=(boolean[])runPW4VoidMethod(iw,"getLogAxes");
		PlotWindow3D pw=null;
		if(plot.getClass().getName().equals("jguis.Traj3D")) {
			float[] newyvals=new float[npts[0][selected]];
			System.arraycopy(yvals[selected],0,newyvals,0,npts[0][selected]);
			float[] newxvals=new float[npts[0][selected]];
			System.arraycopy(xvals[selected],0,newxvals,0,npts[0][selected]);
			float[] newzvals=new float[npts[0][selected]];
			System.arraycopy(zvals[0][selected],0,newzvals,0,npts[0][selected]);
			Traj3D plot2=new Traj3D(labels[1],labels[2],labels[3],newxvals,newyvals,newzvals);
			pw=new PlotWindow3D(labels[0]+"-series"+(selected+1),plot2);
		} else {
			float[] newyvals=new float[npts[1][selected]];
			System.arraycopy(yvals[selected],0,newyvals,0,npts[1][selected]);
			float[] newxvals=new float[npts[selected][0]];
			System.arraycopy(xvals[selected],0,newxvals,0,npts[0][selected]);
			float[][] newzvals=new float[npts[0][selected]][npts[1][selected]];
			for(int j=0;j<npts[0][selected];j++) {
				System.arraycopy(zvals[selected][j],0,newzvals[j],0,npts[1][selected]);
			}
			Plot3D plot2=new Plot3D(labels[1],labels[2],labels[3],newxvals,newyvals,newzvals);
			pw=new PlotWindow3D(labels[0]+"-series"+(selected+1),plot2);
		}
		pw.draw();
		pw.setLogAxes(logaxes[0],logaxes[1],logaxes[2]);
		int[] newshapes=pw.getShapes();
		newshapes[0]=shapes[selected];
		int[] newcolors=pw.getColors();
		newcolors[0]=colors[selected];
		pw.setLimits(limits);
		return pw;
	}
	
	public static Plot3D getPW3DPlotCopy(ImageWindow iw) {
		float[][] yvals=(float[][])runPW4VoidMethod(iw,"getYValues");
		int[][] npts=(int[][])runPW4VoidMethod(iw,"getNpts");
		float[][][] zvals=(float[][][])runPW4VoidMethod(iw,"getZValues");
		float[][] xvals=(float[][])runPW4VoidMethod(iw,"getXValues");
		Object plot=runPW4VoidMethod(iw,"getPlot");
		String[] labels=(String[])runPW4VoidMethod(iw,"getAllLabels");
		String[] newlabels=new String[labels.length];
		copystringarray(labels,newlabels);
		float[] limits=(float[])runPW4VoidMethod(iw,"getLimits");
		int[] shapes=(int[])runPW4VoidMethod(iw,"getShapes");
		int[] colors=(int[])runPW4VoidMethod(iw,"getColors");
		boolean[] logaxes=(boolean[])runPW4VoidMethod(iw,"getLogAxes");
		if(plot.getClass().getName().equals("jguis.Traj3D")) {
			float[][] newyvals=algutils.clone_multidim_array(yvals);
			float[][] newxvals=algutils.clone_multidim_array(xvals);
			float[][] newzvals=algutils.clone_multidim_array(zvals[0]);
			int[] newnpts=npts[0].clone();
			Traj3D plot2=new Traj3D(labels[1],labels[2],labels[3],newxvals,newyvals,newzvals,newnpts);
			plot2.setLogAxes(logaxes[0],logaxes[1],logaxes[2]);
			int[] newshapes=plot2.getShapes();
			for(int i=0;i<newshapes.length;i++) newshapes[i]=shapes[i];
			int[] newcolors=plot2.getColors();
			for(int i=0;i<newcolors.length;i++) newcolors[i]=colors[i];
			plot2.setLimits(limits);
			return plot2;
		} else {
			float[][] newyvals=algutils.clone_multidim_array(yvals);
			float[][] newxvals=algutils.clone_multidim_array(xvals);
			float[][][] newzvals=new float[zvals.length][][];
			int[][] newnpts=algutils.clone_multidim_array(npts);
			for(int j=0;j<zvals.length;j++) {
				newzvals[j]=algutils.clone_multidim_array(zvals[j]);
			}
			Plot3D plot2=new Plot3D(labels[1],labels[2],labels[3],newxvals,newyvals,newzvals,newnpts);
			plot2.setLogAxes(logaxes[0],logaxes[1],logaxes[2]);
			int[] newshapes=plot2.getShapes();
			for(int i=0;i<newshapes.length;i++) newshapes[i]=shapes[i];
			int[] newcolors=plot2.getColors();
			for(int i=0;i<newcolors.length;i++) newcolors[i]=colors[i];
			return plot2;
		}
	}

	public static void savePW4(ImageWindow iw,String filename){
		Class<?> temp=iw.getClass();
		if(temp.getName().equals("ij.gui.PlotWindow")){
			savePW(filename,(PlotWindow)iw);
		} else if(isPlotFamily(iw)){
			runReflectionMethod(iw,"saveAsObject",new Object[]{filename});
		}
	}

	public static void savePW4(ImageWindow iw,String filename,int outtype){
		Class<?> temp=iw.getClass();
		if(temp.getName().equals("jguis.PlotWindow4")){
			runReflectionMethod(iw,"saveAs",new Object[]{filename,new Integer(outtype)});
		}
	}

	public static void savePW(String filename,PlotWindow pw){
		Plot4 p4=plot2Plot4(getPWPlot(pw));
		p4.saveplot2file(filename);
	}
	
	public static Object constructReflectionObject(Object tclass,Object[] args){
		// here we automatically assume that number types are primitive
		if(args==null)
			return constructReflectionObject(tclass,null,null);
		Class[] argcs=new Class[args.length];
		for(int i=0;i<args.length;i++) argcs[i]=args[i].getClass();
		transformClasses(argcs);
		return constructReflectionObject(tclass,args,argcs);
	}
	
	public static Object constructReflectionObject(Object tclass,Object[] args,Class[] argcs){
		try{
			Class<?> temp=tclass.getClass();
			Constructor cons=temp.getDeclaredConstructor(argcs);
			cons.setAccessible(true);
			try{
				return cons.newInstance(args);
			}catch(IllegalArgumentException e){
				IJ.log("illegal argument exception");
			}catch(InstantiationException e){
				IJ.log("instantiation exception");
			}catch(IllegalAccessException e){
				IJ.log("illegal access exception");
			}catch(InvocationTargetException e){
				IJ.log("invocation target exception");
			}
		}catch(SecurityException e){
			IJ.log("security exception");
		}catch(NoSuchMethodException e){
			IJ.log("no such method exception");
		}
		return null;
	}
	
	public static void transformClasses(Class[] argcs){
		for(int i=0;i<argcs.length;i++){
			try{
				if(argcs[i]==Class.forName("ij.CompositeImage")) argcs[i]=Class.forName("ij.ImagePlus");
			}catch(ClassNotFoundException e){}
			if(argcs[i]==Integer.class) argcs[i]=Integer.TYPE;
			if(argcs[i]==Float.class) argcs[i]=Float.TYPE;
			if(argcs[i]==Double.class) argcs[i]=Double.TYPE;
			if(argcs[i]==Short.class) argcs[i]=Short.TYPE;
			if(argcs[i]==Byte.class) argcs[i]=Byte.TYPE;
			if(argcs[i]==Boolean.class) argcs[i]=Boolean.TYPE;
		}
	}

	public static Object runReflectionMethod(Object obj,String method,Object[] args){
		// here we automatically assume that number types are primitive
		if(args==null)
			return runReflectionMethod(obj,method,null,null);
		Class[] argcs=new Class[args.length];
		for(int i=0;i<args.length;i++) argcs[i]=args[i].getClass();
		transformClasses(argcs);
		return runReflectionMethod(obj,method,args,argcs);
	}

	public static Object runReflectionMethod(Object obj,String method,Object[] args,Class[] argcs){
		try{
			Class<?> temp=obj.getClass();
			Method meth=temp.getDeclaredMethod(method,argcs);
			meth.setAccessible(true);
			try{
				Object data=meth.invoke(obj,args);
				return data;
			}catch(IllegalAccessException e){
				IJ.log("illegal access exception");
			}catch(InvocationTargetException e){
				IJ.log("invocation target exception");
				e.printStackTrace();
			}catch(ClassCastException e){
				IJ.log(e.getMessage());
			}
		}catch(NoSuchMethodException e){
			IJ.log("no such method exception");
		}
		return null;
	}
	
	public static Object[] getReflectionMethods(Object obj){
		Class<?> temp=obj.getClass();
		Method[] meth=temp.getDeclaredMethods();
		Class[][] argcs=new Class[meth.length][];
		for(int i=0;i<meth.length;i++) argcs[i]=meth[i].getParameterTypes();
		return new Object[]{meth,argcs};
	}

	public static Object getReflectionField(Object obj,String fieldname){
		try{
			Class<?> temp=obj.getClass();
			Field field=temp.getDeclaredField(fieldname);
			field.setAccessible(true);
			return field.get(obj);
		}catch(NoSuchFieldException e){
			IJ.log("no such field exception");
		}catch(IllegalArgumentException e){
			IJ.log("illegal argument exception");
		}catch(IllegalAccessException e){
			IJ.log("illegal access exception");
		}
		return null;
	}
	
	public static void setReflectionField(Object obj,String fieldname,Object value) {
		try{
			Class<?> temp=obj.getClass();
			Field field=temp.getDeclaredField(fieldname);
			field.setAccessible(true);
			if(value instanceof Boolean) field.setBoolean(obj,((Boolean)value).booleanValue());
			if(value instanceof Byte) field.setByte(obj,((Byte)value).byteValue());
			if(value instanceof Short) field.setShort(obj,((Short)value).shortValue());
			if(value instanceof Float) field.setFloat(obj,((Float)value).floatValue());
			if(value instanceof Double) field.setDouble(obj,((Double)value).doubleValue());
			if(value instanceof Integer) field.setInt(obj,((Integer)value).intValue());
		}catch(NoSuchFieldException e){
			IJ.log("no such field exception");
		}catch(IllegalArgumentException e){
			IJ.log("illegal argument exception");
		}catch(IllegalAccessException e){
			IJ.log("illegal access exception");
		}
	}

	public static void copyfloat2dimvector(float[][] data,float[][] dest){
		for(int i=0;i<data.length;i++){
			System.arraycopy(data[i],0,dest[i],0,data[i].length);
		}
	}

	public static void copystringarray(String[] data,String[] dest){
		for(int i=0;i<data.length;i++){
			dest[i]=data[i].substring(0);
		}
	}

	public static void draw_dashed_line(ImageProcessor ip,int x1,int y1,int x2,int y2,int dash_length){
		int length=(int)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		int nunits=(int)((float)length/(float)(2*dash_length));
		int halfunits=(int)((float)length/(float)dash_length);
		float xinc=(float)(x2-x1)/(float)length;
		float yinc=(float)(x2-x1)/(float)length;
		float x=x1;
		float y=y1;
		ip.moveTo((int)x,(int)y);
		for(int i=0;i<nunits;i++){
			x+=xinc*dash_length;
			y+=yinc*dash_length;
			ip.lineTo((int)x,(int)y);
			x+=xinc*dash_length;
			y+=yinc*dash_length;
			ip.moveTo((int)x,(int)y);
		}
		if(halfunits>2*nunits){
			ip.lineTo(x2,y2);
		}
	}

	public static void draw_arrow(ImageProcessor ip,int x1,int y1,int x2,int y2){
		float length=(float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		if(length==0.0f){
			return;
		}
		float xinc=(x2-x1)/length;
		float yinc=(y2-y1)/length;
		float crossx=x2-xinc*5.0f; //go back 5 pixels from the second point
		float crossy=y2-yinc*5.0f;
		float x3=crossx-yinc*3.5f;
		float x4=crossx+yinc*3.5f;
		float y3=crossy+xinc*3.5f;
		float y4=crossy-xinc*3.5f;
		ip.drawLine(x1,y1,x2,y2);
		ip.drawLine(x2,y2,Math.round(x3),Math.round(y3));
		ip.drawLine(x2,y2,Math.round(x4),Math.round(y4));
	}

	public static void draw_arrow2(ImageProcessor ip,int x1,int y1,int x2,int y2){
		//this isn't really an arrow but rather a line with a dot at the end
		float length=(float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		if(length==0.0f){
			return;
		}
		ip.drawLine(x1,y1,x2,y2);
		ip.drawPixel(x2-1,y2-1);
		ip.drawPixel(x2,y2-1);
		ip.drawPixel(x2+1,y2-1);
		ip.drawPixel(x2-1,y2);
		ip.drawPixel(x2+1,y2);
		ip.drawPixel(x2-1,y2+1);
		ip.drawPixel(x2,y2+1);
		ip.drawPixel(x2+1,y2+1);
	}

	public static Polygon circle2poly(float x,float y,float radius){
		float dphi=2.0f/radius;
		int nangles=(int)((2.0f*(float)Math.PI)/dphi);
		int[] xpts=new int[nangles];
		int[] ypts=new int[nangles];
		for(int i=0;i<nangles;i++){
			float phi=dphi*i;
			xpts[i]=(int)(radius*(float)Math.cos(phi)+x);
			ypts[i]=(int)(radius*(float)Math.sin(phi)+y);
		}
		return new Polygon(xpts,ypts,nangles);
	}

	public static Polygon ellipse2poly(float x1,float y1,float x2,float y2,float aspectRatio){
		// code adapted from the EllipseRoi code in ImageJ
		int[][] coords=new int[2][72];
		double centerX=(x1+x2)/2.0;
		double centerY=(y1+y2)/2.0;
		double dx=x2-x1;
		double dy=y2-y1;
		double major=Math.sqrt(dx*dx+dy*dy);
		double minor=major*aspectRatio;
		double phiB=Math.atan2(dy,dx);
		double alpha=phiB*180.0/Math.PI;
		int nPoints=0;
		for(int i=0;i<72;i++){
			double degrees=i*5.0;
			double beta1=degrees/180.0*Math.PI;
			dx=Math.cos(beta1)*major/2.0;
			dy=Math.sin(beta1)*minor/2.0;
			double beta2=Math.atan2(dy,dx);
			double rad=Math.sqrt(dx*dx+dy*dy);
			double beta3=beta2+alpha/180.0*Math.PI;
			double dx2=Math.cos(beta3)*rad;
			double dy2=Math.sin(beta3)*rad;
			coords[0][nPoints]=(int)(centerX+dx2);
			coords[1][nPoints]=(int)(centerY+dy2);
			nPoints++;
		}
		return new Polygon(coords[0],coords[1],nPoints);
	}

	public static void draw_ellipse(ImageProcessor ip,float x1,float y1,float x2,float y2,float aspectRatio){
		// start by creating the converting the ellipse into a polygon
		Polygon poly=ellipse2poly(x1,y1,x2,y2,aspectRatio);
		draw_polygon(ip,poly,true);
	}

	public static void fill_ellipse(ImageProcessor ip,float x1,float y1,float x2,float y2,float aspectRatio){
		// start by creating the converting the ellipse into a polygon
		Polygon poly=ellipse2poly(x1,y1,x2,y2,aspectRatio);
		fill_polygon(ip,poly);
	}

	public static void draw_polygon(ImageProcessor ip,Polygon poly,boolean closed){
		if(poly==null){
			return;
		}
		int length=poly.npoints;
		int[] xpoints=poly.xpoints;
		int[] ypoints=poly.ypoints;
		for(int i=1;i<length;i++){
			ip.drawLine(xpoints[i-1],ypoints[i-1],xpoints[i],ypoints[i]);
		}
		if(closed){
			ip.drawLine(xpoints[length-1],ypoints[length-1],xpoints[0],ypoints[0]);
		}
	}

	public static void fill_polygon(ImageProcessor ip,Polygon poly){
		if(poly==null){
			return;
		}
		Rectangle r=poly.getBounds();
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					ip.drawPixel(j,i);
				}
			}
		}
	}

	public static int[] intval2rgb(int value){
		int[] temp=new int[3];
		temp[0]=(value&0xff0000)>>16;
		temp[1]=(value&0xff00)>>8;
		temp[2]=value&0xff;
		return temp;
	}

	public static int rgb2intval(int r,int g,int b){
		int temp=0xff000000|(r<<16)|(g<<8)|b;
		return temp;
	}
	
	public static int argb2intval(int a,int r,int g,int b){
		int temp= (a<<24)|(r<<16)|(g<<8)|b;
		return temp;
	}

	public static int rgb2intval(float r,float g,float b){
		int r2=(int)r;
		if(r2>255)
			r2=255;
		if(r2<0)
			r2=0;
		int g2=(int)g;
		if(g2>255)
			g2=255;
		if(g2<0)
			g2=0;
		int b2=(int)b;
		if(r2>255)
			b2=255;
		if(b2<0)
			b2=0;
		int temp=0xff000000|(r2<<16)|(g2<<8)|b2;
		return temp;
	}

	public static byte[][] intval2rgb(int[] values){
		byte[][] temp=new byte[3][values.length];
		for(int i=0;i<values.length;i++){
			int[] temp2=intval2rgb(values[i]);
			temp[0][i]=(byte)temp2[0];
			temp[1][i]=(byte)temp2[1];
			temp[2][i]=(byte)temp2[2];
		}
		return temp;
	}

	public static int rgb2intval(byte r,byte g,byte b){
		int temp=0xff000000|((r&0xff)<<16)|((g&0xff)<<8)|(b&0xff);
		return temp;
	}

	public static boolean isgray(int value){
		return isgray(value,2);
	}

	public static boolean isgray(int value,int tolerance){
		int[] temp=jutils.intval2rgb(value);
		if(Math.abs(temp[1]-temp[0])>(double)tolerance){
			return false;
		}
		if(Math.abs(temp[2]-temp[0])>(double)tolerance){
			return false;
		}
		if(Math.abs(temp[0]-255)<(double)tolerance){
			// this is white
			return false;
		}
		if(temp[0]<tolerance){
			// this is black
			return false;
		}
		return true;
	}

	public static boolean isblack(int value,int tolerance){
		int[] temp=jutils.intval2rgb(value);
		if(temp[0]<tolerance&&temp[1]<tolerance&&temp[2]<tolerance){
			return true;
		}else{
			return false;
		}
	}

	public static boolean iswhite(int value,int tolerance){
		int[] temp=jutils.intval2rgb(value);
		int limit=255-tolerance;
		if(temp[0]>limit&&temp[1]>limit&&temp[2]>limit){
			return true;
		}else{
			return false;
		}
	}

	public static Color get_closest_color(int value){
		return colors[get_closest_color_index(value)];
	}
	
	public static String get_closest_color_name(int value){
		return colornames[get_closest_color_index(value)];
	}
	
	public static String get_closest_color_name(Color color){
		int value=color.getRGB();
		return colornames[get_closest_color_index(value)];
	}

	public static Color[] colors={Color.black,Color.blue,Color.green,Color.red,Color.magenta,Color.cyan,Color.yellow,Color.orange};

	public static String[] colornames={"black","blue","green","red","magenta","cyan","yellow","orange"};

	public static int[] get_color_RGB(Color color){
		int[] temp={color.getRed(),color.getGreen(),color.getBlue()};
		return temp;
	}

	public static int get_closest_color_index(int value){
		int[] temp=jutils.intval2rgb(value);
		int[][] cd=new int[8][];
		for(int i=0;i<8;i++){
			cd[i]=get_color_RGB(colors[i]);
		}
		cd[0][0]=255;
		cd[0][1]=255;
		cd[0][2]=255;
		float shortest=(float)Math.sqrt((temp[0]-cd[0][0])*(temp[0]-cd[0][0])+(temp[1]-cd[0][1])*(temp[1]-cd[0][1])+(temp[2]-cd[0][2])*(temp[2]-cd[0][2]));
		int sindex=0;
		for(int i=1;i<8;i++){
			float distance=(float)Math.sqrt((temp[0]-cd[i][0])*(temp[0]-cd[i][0])+(temp[1]-cd[i][1])*(temp[1]-cd[i][1])+(temp[2]-cd[i][2])*(temp[2]-cd[i][2]));
			if(distance<shortest){
				sindex=i;
				shortest=distance;
			}
		}
		return sindex;
	}

	public static int get_LUT_color_index(LookUpTable lut){
		byte[] greens=lut.getGreens();
		byte[] reds=lut.getReds();
		byte[] blues=lut.getBlues();
		int mapsize=lut.getMapSize();
		int[] max={reds[mapsize-1]&0xff,greens[mapsize-1]&0xff,blues[mapsize-1]&0xff};
		int max2=max[0];
		if(max[1]>max[0]){
			max2=max[1];
		}
		if(max[2]>max[1]){
			max2=max[2];
		}
		float multiplier=255.0f/max2;
		int[] maxnorm={(int)(max[0]*multiplier),(int)(max[1]*multiplier),(int)(max[2]*multiplier)};
		return get_closest_color_index(rgb2intval(maxnorm[0],maxnorm[1],maxnorm[2]));
	}

	public static int get_LUT_color_index(LUT lut){
		byte[] greens=new byte[256];
		byte[] reds=new byte[256];
		byte[] blues=new byte[256];
		lut.getGreens(greens);
		lut.getReds(reds);
		lut.getBlues(blues);
		int mapsize=lut.getMapSize();
		int[] max={reds[mapsize-1]&0xff,greens[mapsize-1]&0xff,blues[mapsize-1]&0xff};
		int max2=max[0];
		if(max[1]>max[0]){
			max2=max[1];
		}
		if(max[2]>max[1]){
			max2=max[2];
		}
		float multiplier=255.0f/max2;
		int[] maxnorm={(int)(max[0]*multiplier),(int)(max[1]*multiplier),(int)(max[2]*multiplier)};
		return get_closest_color_index(rgb2intval(maxnorm[0],maxnorm[1],maxnorm[2]));
	}

	public static LUT lookuptable2LUT(LookUpTable lut){
		return new LUT(lut.getReds(),lut.getGreens(),lut.getBlues());
	}

	public static float get3DSliceStat(ImageStack is,int frame,int slice,int channel,int frames,int slices,int channels,String stat){
		Object pixels=get3DSlice(is,frame,slice,channel,frames,slices,channels);
		return jstatistics.getstatistic(stat,pixels,null);
	}

	public static Object get3DSlice(ImageStack is,int frame,int slice,int channel,int frames,int slices,int channels){
		return is.getPixels(1+channel+slice*channels+frame*channels*slices);
	}
	
	public static Object get3DSliceInterp(ImageStack is,int frame,float slice,int channel,int frames,int slices,int channels){
		int prev=(int)slice;
		Object sl1=is.getPixels(1+channel+prev*channels+frame*channels*slices);
		if(prev==slice) return sl1;
		int next=prev+1;
		Object sl2=is.getPixels(1+channel+next*channels+frame*channels*slices);
		return interpolation.interpz(sl1,sl2,is.getWidth(),is.getHeight(),slice-prev);
	}

	public static void set3DSlice(ImageStack is,Object pixels,int frame,int slice,int channel,int frames,int slices,int channels){
		is.setPixels(pixels,1+channel+slice*channels+frame*channels*slices);
	}

	public static Object[] get3DTSeries(ImageStack is,int slice,int channel,int frames,int slices,int channels){
		Object[] temp=new Object[frames];
		for(int i=0;i<frames;i++){
			temp[i]=is.getPixels(i*slices*channels+slice*channels+channel+1);
		}
		return temp;
	}

	public static void set3DTSeries(Object[] source,ImageStack is,int channel,int slice,int frames,int slices,int channels){
		for(int i=0;i<frames;i++){
			is.setPixels(source[i],i*slices*channels+slice*channels+channel+1);
		}
	}

	public static Object[] get3DZSeries(ImageStack is,int channel,int frame,int frames,int slices,int channels){
		Object[] temp=new Object[slices];
		for(int i=0;i<slices;i++){
			temp[i]=is.getPixels(frame*slices*channels+i*channels+channel+1);
		}
		return temp;
	}

	public static void set3DZSeries(Object[] source,ImageStack is,int channel,int frame,int frames,int slices,int channels){
		for(int i=0;i<slices;i++){
			is.setPixels(source[i],frame*slices*channels+i*channels+channel+1);
		}
	}

	public static Object[] get3DCSeries(ImageStack is,int slice,int frame,int frames,int slices,int channels){
		Object[] temp=new Object[channels];
		for(int i=0;i<channels;i++){
			temp[i]=is.getPixels(frame*slices*channels+slice*channels+i+1);
		}
		return temp;
	}

	public static void set3DCSeries(Object[] source,ImageStack is,int slice,int frame,int frames,int slices,int channels){
		for(int i=0;i<channels;i++){
			is.setPixels(source[i],frame*slices*channels+slice*channels+i+1);
		}
	}

	public static float[] get3DProjZStat(ImageStack is,int frame,int channel,int frames,int slices,int channels,String stat){
		int npixels=is.getWidth()*is.getHeight();
		boolean isfloat=(is.getPixels(1) instanceof float[]);
		boolean isshort=(is.getPixels(1) instanceof short[]);
		Object[] zstack=new Object[slices];
		for(int i=0;i<slices;i++){
			zstack[i]=get3DSlice(is,frame,i,channel,frames,slices,channels);
		}
		float[] projection=new float[npixels];
		for(int i=0;i<npixels;i++){
			if(isfloat){
				float[] temp=new float[slices];
				for(int j=0;j<slices;j++){
					temp[j]=((float[])zstack[j])[i];
				}
				projection[i]=jstatistics.getstatistic(stat,temp,null);
			}else{
				if(isshort){
					short[] temp=new short[slices];
					for(int j=0;j<slices;j++){
						temp[j]=((short[])zstack[j])[i];
					}
					projection[i]=jstatistics.getstatistic(stat,temp,null);
				}else{
					byte[] temp=new byte[slices];
					for(int j=0;j<slices;j++){
						temp[j]=((byte[])zstack[j])[i];
					}
					projection[i]=jstatistics.getstatistic(stat,temp,null);
				}
			}
		}
		return projection;
	}

	public static float[][] get3DProjZStat(ImageStack is,String stat,int startslice,int slices,float[] options){
		int npixels=is.getWidth()*is.getHeight();
		boolean isfloat=(is.getPixels(1) instanceof float[]);
		boolean isshort=(is.getPixels(1) instanceof short[]);
		float[][] projection=null;
		boolean ispercentile=false;
		if(stat=="Percentile"){
			projection=new float[options.length][npixels];
			ispercentile=true;
		}else
			projection=new float[1][npixels];
		Object[] zstack=is.getImageArray();
		if(!ispercentile){
			for(int i=0;i<npixels;i++){
				if(isfloat){
					float[] temp=new float[slices];
					for(int j=startslice;j<(startslice+slices);j++){
						temp[j]=((float[])zstack[j])[i];
					}
					projection[0][i]=jstatistics.getstatistic(stat,temp,options);
				}else{
					if(isshort){
						short[] temp=new short[slices];
						for(int j=startslice;j<(startslice+slices);j++){
							temp[j]=((short[])zstack[j])[i];
						}
						projection[0][i]=jstatistics.getstatistic(stat,temp,options);
					}else{
						byte[] temp=new byte[slices];
						for(int j=startslice;j<(startslice+slices);j++){
							temp[j]=((byte[])zstack[j])[i];
						}
						projection[0][i]=jstatistics.getstatistic(stat,temp,options);
					}
				}
			}
		}else{
			for(int i=0;i<npixels;i++){
				float[] tempoptions=options.clone();
				if(isfloat){
					float[] temp=new float[slices];
					for(int j=startslice;j<(startslice+slices);j++){
						temp[j]=((float[])zstack[j])[i];
					}
					jstatistics.getstatistic(stat,temp,tempoptions);
					for(int j=0;j<tempoptions.length;j++)
						projection[j][i]=tempoptions[j];
				}else{
					if(isshort){
						short[] temp=new short[slices];
						for(int j=startslice;j<(startslice+slices);j++){
							temp[j]=((short[])zstack[j])[i];
						}
						jstatistics.getstatistic(stat,temp,tempoptions);
						for(int j=0;j<tempoptions.length;j++)
							projection[j][i]=tempoptions[j];
					}else{
						byte[] temp=new byte[slices];
						for(int j=startslice;j<(startslice+slices);j++){
							temp[j]=((byte[])zstack[j])[i];
						}
						jstatistics.getstatistic(stat,temp,tempoptions);
						for(int j=0;j<tempoptions.length;j++)
							projection[j][i]=tempoptions[j];
					}
				}
			}
		}
		return projection;
	}

	public static Object get3DProjZStat2(ImageStack is,int frame,int channel,int frames,int slices,int channels,String stat){
		// here we convert back to original data type
		int npixels=is.getWidth()*is.getHeight();
		boolean isfloat=(is.getPixels(1) instanceof float[]);
		boolean isshort=(is.getPixels(1) instanceof short[]);
		Object[] zstack=new Object[slices];
		for(int i=0;i<slices;i++){
			zstack[i]=get3DSlice(is,frame,i,channel,frames,slices,channels);
		}
		float[] projection=new float[npixels];
		for(int i=0;i<npixels;i++){
			if(isfloat){
				float[] temp=new float[slices];
				for(int j=0;j<slices;j++){
					temp[j]=((float[])zstack[j])[i];
				}
				projection[i]=jstatistics.getstatistic(stat,temp,null);
			}else{
				if(isshort){
					short[] temp=new short[slices];
					for(int j=0;j<slices;j++){
						temp[j]=((short[])zstack[j])[i];
					}
					projection[i]=jstatistics.getstatistic(stat,temp,null);
				}else{
					byte[] temp=new byte[slices];
					for(int j=0;j<slices;j++){
						temp[j]=((byte[])zstack[j])[i];
					}
					projection[i]=jstatistics.getstatistic(stat,temp,null);
				}
			}
		}
		if(isfloat)
			return projection;
		if(isshort)
			return convert_arr_short(projection);
		return convert_arr_byte(projection);
	}

	public static float[] get3DProjTStat(ImageStack is,int slice,int channel,int frames,int slices,int channels,String stat){
		int npixels=is.getWidth()*is.getHeight();
		boolean isfloat=(is.getPixels(1) instanceof float[]);
		boolean isbyte=(is.getPixels(1) instanceof byte[]);
		Object[] zstack=new Object[frames];
		for(int i=0;i<frames;i++){
			zstack[i]=get3DSlice(is,i,slice,channel,frames,slices,channels);
		}
		float[] projection=new float[npixels];
		for(int i=0;i<npixels;i++){
			if(isfloat){
				float[] temp=new float[frames];
				for(int j=0;j<frames;j++){
					temp[j]=((float[])zstack[j])[i];
				}
				projection[i]=jstatistics.getstatistic(stat,temp,null);
			}else{
				if(isbyte){
					byte[] temp=new byte[slices];
					for(int j=0;j<frames;j++){
						temp[j]=((byte[])zstack[j])[i];
					}
					projection[i]=jstatistics.getstatistic(stat,temp,null);
				}else{
					short[] temp=new short[slices];
					for(int j=0;j<frames;j++){
						temp[j]=((short[])zstack[j])[i];
					}
					projection[i]=jstatistics.getstatistic(stat,temp,null);
				}
			}
		}
		return projection;
	}

	public static float[] get_slice_copy(ImageStack stack,int slice){
		Object pix=stack.getPixels(slice);
		return algutils.convert_arr_float(pix);
	}

	public static float[] getStatsOptions(String stat){
		if(stat=="Mode"){
			GenericDialog gd=new GenericDialog("Stat Options");
			gd.addNumericField("#_of_Hist_Bins",20,0);
			gd.addNumericField("Hist_Start",0.0,5,15,null);
			gd.addNumericField("Hist_End",1.0,5,15,null);
			gd.showDialog();
			if(gd.wasCanceled()){
				return null;
			}
			float[] temp=new float[3];
			temp[0]=(int)gd.getNextNumber();
			temp[1]=(float)gd.getNextNumber();
			temp[2]=(float)gd.getNextNumber();
			return temp;
		}
		if(stat=="ConditionalAvg"){
			GenericDialog gd=new GenericDialog("Stat Options");
			gd.addNumericField("Lower_Limit",0.0,5,15,null);
			gd.addNumericField("Upper_Limit",1.0,5,15,null);
			gd.showDialog();
			if(gd.wasCanceled()){
				return null;
			}
			float[] temp=new float[2];
			temp[0]=(float)gd.getNextNumber();
			temp[1]=(float)gd.getNextNumber();
			return temp;
		}
		if(stat=="Count"){
			GenericDialog gd=new GenericDialog("Stat Options");
			gd.addNumericField("Lower_Limit",0.0,5,15,null);
			gd.addNumericField("Upper_Limit",1.0,5,15,null);
			gd.showDialog();
			if(gd.wasCanceled()){
				return null;
			}
			float[] temp=new float[2];
			temp[0]=(float)gd.getNextNumber();
			temp[1]=(float)gd.getNextNumber();
			return temp;
		}
		if(stat=="Percentile"){
			GenericDialog gd=new GenericDialog("Stat Options");
			gd.addNumericField("Percentile",95.0,5,15,null);
			gd.addNumericField("#_of_Percentiles?",1,0);
			gd.showDialog();
			if(gd.wasCanceled()){
				return null;
			}
			float temp=(float)gd.getNextNumber();
			int number=(int)gd.getNextNumber();
			if(number==1){
				float[] temp2={temp};
				return temp2;
			}
			float[] temp2=new float[number];
			temp2[0]=temp;
			GenericDialog gd2=new GenericDialog("More Stat Options");
			for(int i=0;i<(number-1);i++)
				gd2.addNumericField("Percentile"+(i+2),95.0,5,15,null);
			gd2.showDialog();
			if(gd2.wasCanceled()){
				return null;
			}
			for(int i=0;i<(number-1);i++)
				temp2[i+1]=(float)gd2.getNextNumber();
			return temp2;
		}
		return null;
	}

	public static Object[] stack2array(ImageStack is){
		Object[] temp=new Object[is.getSize()];
		for(int i=0;i<is.getSize();i++){
			temp[i]=is.getPixels(i+1);
		}
		return temp;
	}
	
	public static ImageStack array2stack(Object[] arr,int width,int height){
		ImageStack stack=new ImageStack(width,height);
		for(int i=0;i<arr.length;i++) stack.addSlice("",arr[i]);
		return stack;
	}
	
	public static ImageStack array2stack(float[][] arr,int width,int height){
		ImageStack stack=new ImageStack(width,height);
		for(int i=0;i<arr.length;i++) stack.addSlice("",arr[i]);
		return stack;
	}

	public static boolean[] roi2mask(Roi roi,int width,int height){
		if(roi==null)
			return null;
		boolean[] retmask=new boolean[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				retmask[j+i*width]=roi.contains(j,i);
			}
		}
		return retmask;
	}

	public static ImageProcessor custom_8_bit(ColorProcessor cp){
		// here we convert an RGB image to 8 bit using 17 gray levels and 15
		// green and red levels as well as combinations
		int[] pixels=(int[])cp.getPixels();
		byte[] pixels8=new byte[pixels.length];
		byte[][] LUT=new byte[3][256];
		for(int i=0;i<16;i++){
			for(int j=0;j<3;j++){
				LUT[j][i]=(byte)(i*16);
			}
		}
		LUT[0][16]=(byte)255;
		LUT[1][16]=(byte)255;
		LUT[2][16]=(byte)255;
		// here j will index over red and i will index over green
		for(int i=0;i<15;i++){
			for(int j=0;j<15;j++){
				LUT[0][j+i*15+17]=(byte)((int)((j)*18.215));
				LUT[1][j+i*15+17]=(byte)((int)((i)*18.215));
			}
		}
		IndexColorModel cm=new IndexColorModel(8,256,LUT[0],LUT[1],LUT[2]);
		// now parse the original colors into the predefined channels
		for(int i=0;i<pixels.length;i++){
			if(isblack(pixels[i],9)){
				pixels8[i]=(byte)0;
			}else{
				if(iswhite(pixels[i],15)){
					pixels8[i]=(byte)16;
				}else{
					int[] rgb=intval2rgb(pixels[i]);
					if(isgray(pixels[i],15)){
						int graylevel=(int)((rgb[0]+rgb[1]+rgb[2])/(3.0*16.0));
						pixels8[i]=(byte)graylevel;
					}else{
						int redlevel=(int)((rgb[0])/18.214);
						int greenlevel=(int)((rgb[1])/18.214);
						pixels8[i]=(byte)(redlevel+greenlevel*15+17);
					}
				}
			}
			// pixels8[i]=(byte)31;
		}
		return new ByteProcessor(cp.getWidth(),cp.getHeight(),pixels8,cm);
	}

	public static void set_plugin_options(String plugname,String[] vals){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){
				a.mkdir();
			}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+plugname+".jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			for(int i=0;i<vals.length;i++){
				d.write(vals[i]+"\n");
			}
			d.close();
		}catch(IOException e){
			IJ.showMessage("error writing file");
		}
		return;
	}

	public static String[] get_plugin_options(String plugname,int noptions){
		String dir=System.getProperty("user.home");
		String[] vals=new String[noptions];
		try{
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+plugname+".jrn");
			BufferedReader d=new BufferedReader(new FileReader(b));
			for(int i=0;i<noptions;i++){
				vals[i]=d.readLine();
				if(vals[i]==null){
					d.close();
					return null;
				}
			}
			d.close();
		}catch(IOException e){
			return null;
		}
		return vals;
	}

	public static LUT get_lut_for_color(Color color){
		float maxblue=color.getBlue();
		float maxred=color.getRed();
		float maxgreen=color.getGreen();
		if(maxblue<=0.0f&&maxred<=0.0f&&maxgreen<=0.0f){
			maxblue=255.0f;
			maxred=255.0f;
			maxgreen=255.0f;
		}
		byte[] r=new byte[256];
		byte[] g=new byte[256];
		byte[] b=new byte[256];
		for(int i=0;i<256;i++){
			float fraction=(i)/255.0f;
			r[i]=(byte)((int)(maxred*fraction));
			g[i]=(byte)((int)(maxgreen*fraction));
			b[i]=(byte)((int)(maxblue*fraction));
		}
		return new LUT(r,g,b);
	}

	public static Object convert_array(Object oldarr,int type){
		// can convert from and to byte, short, or float
		switch(type){
		case 0:
			return convert_arr_byte(oldarr);
		case 1:
			return convert_arr_short(oldarr);
		case 2:
			return convert_arr_float(oldarr);
		}
		return null;
	}

	public static byte[] convert_arr_byte(Object oldarr){
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			byte[] newarr=new byte[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=(int)temparr[i];
				if(temp>=0){
					if(temp<256){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)255;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			byte[] newarr=new byte[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xffff;
				if(temp>=0){
					if(temp<256){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)255;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			byte[] newarr=new byte[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		return null;
	}

	public static float[] convert_arr_float(Object oldarr){
		if(oldarr instanceof double[]){
			double[] temparr=(double[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				newarr[i]=(float)temparr[i];
			}
			return newarr;
		}
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			float[] newarr=new float[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xffff;
				newarr[i]=temp;
			}
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			float[] newarr=new float[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xff;
				newarr[i]=temp;
			}
			return newarr;
		}
		return null;
	}

	public static short[] convert_arr_short(Object oldarr){
		if(oldarr instanceof float[]){
			float[] temparr=(float[])oldarr;
			short[] newarr=new short[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=(int)temparr[i];
				if(temp>=0){
					if(temp<65536){
						newarr[i]=(byte)temp;
					}else{
						newarr[i]=(byte)65536;
					}
				}else{
					newarr[i]=(byte)0;
				}
			}
			return newarr;
		}
		if(oldarr instanceof short[]){
			short[] temparr=(short[])oldarr;
			short[] newarr=new short[temparr.length];
			System.arraycopy(temparr,0,newarr,0,temparr.length);
			return newarr;
		}
		if(oldarr instanceof byte[]){
			byte[] temparr=(byte[])oldarr;
			short[] newarr=new short[temparr.length];
			for(int i=0;i<temparr.length;i++){
				int temp=temparr[i]&0xff;
				newarr[i]=(short)temp;
			}
			return newarr;
		}
		return null;
	}

	public static Object[] getImageWindowList(boolean addnull){
		int[] ids=WindowManager.getIDList();
		String[] titles=null;
		if(addnull){
			titles=new String[ids.length+1];
		}else{
			titles=new String[ids.length];
		}
		for(int i=0;i<ids.length;i++){
			ImagePlus imp=WindowManager.getImage(ids[i]);
			if(imp!=null){
				titles[i]=imp.getTitle();
			}else{
				titles[i]="";
			}
		}
		if(addnull){
			titles[ids.length]="null";
		}
		Object[] retvals={ids,titles};
		return retvals;
	}

	public static Object[] getPlotWindowList(boolean addnull){
		int[] ids=WindowManager.getIDList();
		String[] titles=null;
		if(addnull){
			titles=new String[ids.length+1];
		}else{
			titles=new String[ids.length];
		}
		for(int i=0;i<ids.length;i++){
			ImagePlus imp=WindowManager.getImage(ids[i]);
			ImageWindow iw=imp.getWindow();
			if(isPlot(iw)){
				titles[i]=imp.getTitle();
			}else{
				titles[i]="NA";
			}
		}
		if(addnull){
			titles[ids.length]="null";
		}
		Object[] retvals={ids,titles};
		return retvals;
	}
	
	public static Object[] getPlotFamilyWindowList(boolean addnull){
		int[] ids=WindowManager.getIDList();
		String[] titles=null;
		if(addnull){
			titles=new String[ids.length+1];
		}else{
			titles=new String[ids.length];
		}
		for(int i=0;i<ids.length;i++){
			ImagePlus imp=WindowManager.getImage(ids[i]);
			ImageWindow iw=imp.getWindow();
			if(isPlotFamily(iw)){
				titles[i]=imp.getTitle();
			}else{
				titles[i]="NA";
			}
		}
		if(addnull){
			titles[ids.length]="null";
		}
		Object[] retvals={ids,titles};
		return retvals;
	}

	public static Object[] getTableWindowList(boolean addnull){
		Frame[] niframes=WindowManager.getNonImageWindows();
		String[] titles=new String[niframes.length];
		Frame[] niframes2=new Frame[niframes.length];
		int counter=0;
		for(int i=0;i<niframes.length;i++){
			if(niframes[i] instanceof TextWindow && !niframes[i].getTitle().equals("Log")){
				titles[counter]=niframes[i].getTitle();
				niframes2[counter]=niframes[i];
				counter++;
			}
		}
		int ntitles=counter;
		if(addnull) ntitles++;
		String[] titles2=new String[ntitles];
		Frame[] niframes3=new Frame[ntitles];
		for(int i=0;i<counter;i++){titles2[i]=titles[i]; niframes3[i]=niframes2[i];}
		if(addnull){
			titles[counter]="null";
			niframes3[counter]=null;
		}
		Object[] retvals={niframes3,titles2};
		return retvals;
	}

	public static ImagePlus[] selectImages(boolean addnull,int nimages,String[] labels){
		Object[] windowlist=jutils.getImageWindowList(addnull);
		String[] titles=(String[])windowlist[1];
		int[] ids=(int[])windowlist[0];
		GenericDialog gd=new GenericDialog("Select Images");
		for(int i=0;i<nimages;i++){
			gd.addChoice(labels[i],titles,titles[0]);
		}
		gd.showDialog();
		if(gd.wasCanceled()){
			return null;
		}
		ImagePlus[] windows=new ImagePlus[nimages];
		for(int i=0;i<nimages;i++){
			int index=gd.getNextChoiceIndex();
			if(index==ids.length){
				windows[i]=null;
			}else{
				windows[i]=WindowManager.getImage(ids[index]);
			}
		}
		return windows;
	}

	public static ImagePlus[] selectImages(boolean addnull,int nimages){
		String[] labels=new String[nimages];
		for(int i=0;i<nimages;i++){
			labels[i]="Image"+(i+1);
		}
		return selectImages(addnull,nimages,labels);
	}

	public static TextWindow[] selectTables(boolean addnull,int ntables,String[] labels){
		Object[] windowlist=jutils.getTableWindowList(addnull);
		String[] titles=(String[])windowlist[1];
		Frame[] frames=(Frame[])windowlist[0];
		GenericDialog gd=new GenericDialog("Select Tables");
		for(int i=0;i<ntables;i++){
			gd.addChoice(labels[i],titles,titles[0]);
		}
		gd.showDialog();
		if(gd.wasCanceled()){
			return null;
		}
		TextWindow[] windows=new TextWindow[ntables];
		for(int i=0;i<ntables;i++){
			int index=gd.getNextChoiceIndex();
			if(index==frames.length){
				windows[i]=null;
			}else{
				windows[i]=(TextWindow)frames[index];
			}
		}
		return windows;
	}

	public static TextWindow[] selectTables(boolean addnull,int ntables){
		String[] labels=new String[ntables];
		for(int i=0;i<ntables;i++){
			labels[i]="Table"+(i+1);
		}
		return selectTables(addnull,ntables,labels);
	}

	/************************************
	 * selects a currently open table in ImageJ
	 * @param title: the title of the table
	 * @return
	 */
	public static TextWindow selectTable(String title){
		Object[] windowlist=jutils.getTableWindowList(false);
		String[] titles=(String[])windowlist[1];
		Frame[] frames=(Frame[])windowlist[0];
		for(int i=0;i<titles.length;i++){
			if(titles[i].equals(title)){
				return (TextWindow)frames[i];
			}
		}
		return null;
	}
	
	public static String getUniqueTableName(String title){
		String temp=title.substring(0);
		TextWindow tw=selectTable(temp);
		int increment=1;
		while(tw!=null){
			temp=temp+"-"+increment;
			tw=selectTable(temp);
			increment++;
		}
		return temp;
	}
	
	public static ImageWindow[] selectPlotFamily(boolean addnull,int nimages){
		String[] labels=new String[nimages];
		for(int i=0;i<nimages;i++){
			labels[i]="Plot"+(i+1);
		}
		return selectPlotFamily(addnull,nimages,labels);
	}
	
	public static ImageWindow[] selectPlotFamily(boolean addnull,int nimages,String[] labels){
		Object[] windowlist=jutils.getPlotFamilyWindowList(addnull);
		String[] titles=(String[])windowlist[1];
		int[] ids=(int[])windowlist[0];
		GenericDialog gd=new GenericDialog("Select Plots");
		for(int i=0;i<nimages;i++){
			gd.addChoice(labels[i],titles,titles[0]);
		}
		gd.showDialog();
		if(gd.wasCanceled()){
			return null;
		}
		ImageWindow[] windows=new ImageWindow[nimages];
		for(int i=0;i<nimages;i++){
			int index=gd.getNextChoiceIndex();
			if(index==ids.length){
				windows[i]=null;
			}else{
				windows[i]=WindowManager.getImage(ids[index]).getWindow();
			}
		}
		return windows;
	}

	public static ImageWindow[] selectPlots(boolean addnull,int nimages,String[] labels){
		Object[] windowlist=jutils.getPlotWindowList(addnull);
		String[] titles=(String[])windowlist[1];
		int[] ids=(int[])windowlist[0];
		GenericDialog gd=new GenericDialog("Select Plots");
		for(int i=0;i<nimages;i++){
			gd.addChoice(labels[i],titles,titles[0]);
		}
		gd.showDialog();
		if(gd.wasCanceled()){
			return null;
		}
		ImageWindow[] windows=new ImageWindow[nimages];
		for(int i=0;i<nimages;i++){
			int index=gd.getNextChoiceIndex();
			if(index==ids.length){
				windows[i]=null;
			}else{
				windows[i]=WindowManager.getImage(ids[index]).getWindow();
			}
		}
		return windows;
	}

	public static ImageWindow[] selectPlots(boolean addnull,int nimages){
		String[] labels=new String[nimages];
		for(int i=0;i<nimages;i++){
			labels[i]="Plot"+(i+1);
		}
		return selectPlots(addnull,nimages,labels);
	}
	
	public static Object[] selectMixedWindows(boolean addnull,int nwindows,String[] labels,int[] types){
		//here we select multiple windows with different types: 0=image, 1=plot family, 2=table
		Object[] plotlist=jutils.getPlotFamilyWindowList(addnull);
		Object[] imagelist=jutils.getImageWindowList(addnull);
		Object[] tablelist=getTableWindowList(addnull);
		String[] plottitles=(String[])plotlist[1];
		int[] plotids=(int[])plotlist[0];
		String[] imagetitles=(String[])imagelist[1];
		int[] imageids=(int[])imagelist[0];
		String[] tabletitles=(String[])tablelist[1];
		Frame[] tableids=(Frame[])tablelist[0];
		GenericDialog gd=new GenericDialog("Select Windows");
		for(int i=0;i<nwindows;i++){
			if(types[i]==0) gd.addChoice(labels[i],imagetitles,imagetitles[0]);
			if(types[i]==1) gd.addChoice(labels[i],plottitles,plottitles[0]);
			if(types[i]==2) gd.addChoice(labels[i],tabletitles,tabletitles[0]);
		}
		gd.showDialog();
		if(gd.wasCanceled()){
			return null;
		}
		Object[] windows=new Object[nwindows];
		for(int i=0;i<nwindows;i++){
			int index=gd.getNextChoiceIndex();
			if(types[i]==0){
    			if(index==imageids.length) windows[i]=null;
    			else windows[i]=WindowManager.getImage(imageids[index]);
			}
			if(types[i]==1){
    			if(index==plotids.length) windows[i]=null;
    			else windows[i]=WindowManager.getImage(imageids[index]).getWindow();
			}
			if(types[i]==2){
    			if(index==tableids.length) windows[i]=null;
    			else windows[i]=tableids[index];
			}
		}
		return windows;
	}

	public static void copyColorRoi(Polygon roi,ImagePlus source,ImagePlus destination,int x,int y,Rectangle clip){
		int[] srcpix=(int[])source.getProcessor().getPixels();
		copyColorRoi(roi,srcpix,source.getWidth(),source.getHeight(),destination,x,y,clip);
	}

	public static void copyColorRoi(Polygon roi,int[] source,int srcwidth,int srcheight,ImagePlus destination,int x,int y,Rectangle clip){
		Rectangle r=roi.getBounds();
		int[] srcpix=source;
		int[] dstpix=(int[])destination.getProcessor().getPixels();
		Rectangle clip2=null;
		if(clip==null){
			clip2=new Rectangle(destination.getWidth()-1,destination.getHeight()-1);
		}else{
			clip2=new Rectangle(clip.x,clip.y,clip.width,clip.height);
			if(clip2.x<0) clip2.x=0;
			if(clip2.y<0) clip2.y=0;
			if((clip2.x+clip2.width)>=destination.getWidth()) clip2.width=(destination.getWidth()-clip2.x-1);
			if((clip2.y+clip2.height)>=destination.getHeight()) clip2.width=(destination.getHeight()-clip2.y-1);
		}
		int dstwidth=destination.getWidth();
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(clip2.contains(x+j-r.x,y+i-r.y)){
					if(roi.contains(j,i)){
						if(j>=0 && i>=0 && j<srcwidth && i<srcheight) dstpix[x+j-r.x+(y+i-r.y)*dstwidth]=srcpix[j+i*srcwidth];
					}
				}
			}
		}
	}

	public static void clearColorImp(ImagePlus imp){
		int[] pixels=(int[])imp.getProcessor().getPixels();
		for(int i=0;i<pixels.length;i++){
			pixels[i]=0xff000000;
		}
	}

	public static ImageProcessor binProcessor(ImageProcessor ip,int binx,int biny){
		int width=ip.getWidth();
		int height=ip.getHeight();
		int newwidth=(int)((float)width/(float)binx);
		int newheight=(int)((float)height/(float)biny);
		float totbin=binx*biny;
		Object pixels=ip.getPixels();
		if(pixels instanceof float[]){
			float[] pix2=new float[newwidth*newheight];
			for(int i=0;i<newheight;i++){
				for(int j=0;j<newwidth;j++){
					for(int k=0;k<biny;k++){
						for(int l=0;l<binx;l++){
							pix2[j+i*newwidth]+=((float[])pixels)[l+j*binx+(k+i*biny)*width]/totbin;
						}
					}
				}
			}
			return new FloatProcessor(newwidth,newheight,pix2,ip.getColorModel());
		}else{
			if(pixels instanceof short[]){
				short[] pix2=new short[newwidth*newheight];
				for(int i=0;i<newheight;i++){
					for(int j=0;j<newwidth;j++){
						float avg=0.0f;
						for(int k=0;k<biny;k++){
							for(int l=0;l<binx;l++){
								avg+=(((short[])pixels)[l+j*binx+(k+i*biny)*width]&0xffff)/totbin;
							}
						}
						pix2[j+i*newwidth]=(short)avg;
					}
				}
				return new ShortProcessor(newwidth,newheight,pix2,ip.getColorModel());
			}else{
				if(pixels instanceof byte[]){
					byte[] pix2=new byte[newwidth*newheight];
					for(int i=0;i<newheight;i++){
						for(int j=0;j<newwidth;j++){
							float avg=0.0f;
							for(int k=0;k<biny;k++){
								for(int l=0;l<binx;l++){
									avg+=(((byte[])pixels)[l+j*binx+(k+i*biny)*width]&0xff)/totbin;
								}
							}
							pix2[j+i*newwidth]=(byte)avg;
						}
					}
					return new ByteProcessor(newwidth,newheight,pix2,ip.getColorModel());
				}else{
					int[] pix2=new int[newwidth*newheight];
					for(int i=0;i<newheight;i++){
						for(int j=0;j<newwidth;j++){
							float avgr=0.0f;
							float avgg=0.0f;
							float avgb=0.0f;
							for(int k=0;k<biny;k++){
								for(int l=0;l<binx;l++){
									int[] temp=intval2rgb(((int[])pixels)[l+j*binx+(k+i*biny)*width]);
									avgr+=temp[0]/totbin;
									avgg+=temp[1]/totbin;
									avgb+=temp[2]/totbin;
								}
							}
							pix2[j+i*newwidth]=rgb2intval(avgr,avgg,avgb);
						}
					}
					return new ColorProcessor(newwidth,newheight,pix2);
				}
			}
		}
	}

	public static ImageProcessor cropProcessor(ImageProcessor ip,Rectangle r){
		int width=ip.getWidth();
		int height=ip.getHeight();
		if(width==r.width&&height==r.height){
			return ip;
		}
		if(r.x<0||r.y<0||(r.x+r.width)>=width||(r.y+r.height)>=height){
			return null;
		}
		Object pixels=ip.getPixels();
		if(pixels instanceof float[]){
			float[] pix2=new float[r.width*r.height];
			for(int i=0;i<r.height;i++){
				System.arraycopy(pixels,r.x+(r.y+i)*width,pix2,i*r.width,r.width);
			}
			return new FloatProcessor(r.width,r.height,pix2,ip.getColorModel());
		}else{
			if(pixels instanceof short[]){
				short[] pix2=new short[r.width*r.height];
				for(int i=0;i<r.height;i++){
					System.arraycopy(pixels,r.x+(r.y+i)*width,pix2,i*r.width,r.width);
				}
				return new ShortProcessor(r.width,r.height,pix2,ip.getColorModel());
			}else{
				if(pixels instanceof byte[]){
					byte[] pix2=new byte[r.width*r.height];
					for(int i=0;i<r.height;i++){
						System.arraycopy(pixels,r.x+(r.y+i)*width,pix2,i*r.width,r.width);
					}
					return new ByteProcessor(r.width,r.height,pix2,ip.getColorModel());
				}else{
					int[] pix2=new int[r.width*r.height];
					for(int i=0;i<r.height;i++){
						System.arraycopy(pixels,r.x+(r.y+i)*width,pix2,i*r.width,r.width);
					}
					return new ColorProcessor(r.width,r.height,pix2);
				}
			}
		}
	}

	public static void set_psize(ImagePlus imp,double psize,String unit){
		Calibration cal=imp.getCalibration();
		cal.pixelWidth=psize;
		cal.pixelHeight=psize;
		cal.setUnit(unit);
	}

	public static void set_psize(ImagePlus imp,double psize){
		set_psize(imp,psize,"um");
	}

	public static double get_psize(ImagePlus imp){
		return imp.getCalibration().pixelWidth;
	}

	public static void set_pdepth(ImagePlus imp,double depth){
		imp.getCalibration().pixelDepth=depth;
	}

	public static double get_pdepth(ImagePlus imp){
		return imp.getCalibration().pixelDepth;
	}
	
	public static double get_zratio(ImagePlus imp){
		return get_pdepth(imp)/get_psize(imp);
	}

	public static void set_pinterval(ImagePlus imp,double interval){
		imp.getCalibration().frameInterval=interval;
	}

	public static double get_pinterval(ImagePlus imp){
		return imp.getCalibration().frameInterval;
	}

	/*****************************
	 * here we provide the luts as well
	 * @param title
	 * @param stack
	 * @param frames
	 * @param slices
	 * @param channels
	 * @param composite
	 * @param luts
	 * @return
	 */
	public static ImagePlus create_hyperstack(String title,ImageStack stack,int frames,int slices,int channels,boolean composite,LUT[] luts){
		ImagePlus imp=new ImagePlus(title,stack);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(channels,slices,frames);
		if(channels>1&&composite){
			int model=CompositeImage.COLOR;
			if(composite)
				model=CompositeImage.COMPOSITE;
			CompositeImage ci=new CompositeImage(imp,model);
			if(luts!=null){
				ci.setLuts(luts);
			}
			return ci;
		}else{
			return imp;
		}
	}

	/**********************************
	 * here the hyperstack has the exact same dimensions as the template
	 * @param title
	 * @param stack
	 * @param template
	 * @return
	 */
	public static ImagePlus create_hyperstack(String title,ImageStack stack,ImagePlus template){
		ImagePlus imp=new ImagePlus(title,stack);
		imp.copyScale(template);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(template.getNChannels(),template.getNSlices(),template.getNFrames());
		if(template.isComposite()){
			int model=((CompositeImage)template).getMode();
			CompositeImage ci=new CompositeImage(imp,model);
			LUT[] lut=((CompositeImage)template).getLuts();
			ci.setLuts(lut);
			return ci;
		}else{
			return imp;
		}
	}

	/************************
	 * this creates a hyperstack with different number of frames, slices, and channels than the template
	 * @param title
	 * @param stack
	 * @param template
	 * @param frames
	 * @param slices
	 * @param channels
	 * @return
	 */
	public static ImagePlus create_hyperstack(String title,ImageStack stack,ImagePlus template,int frames,int slices,int channels){
		ImagePlus imp=new ImagePlus(title,stack);
		imp.copyScale(template);
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(channels,slices,frames);
		if(template.isComposite()){
			int model=((CompositeImage)template).getMode();
			CompositeImage ci=new CompositeImage(imp,model);
			if(template.getNChannels()==channels){
				LUT[] lut=((CompositeImage)template).getLuts();
				ci.setLuts(lut);
			}
			return ci;
		}else{
			return imp;
		}
	}

	public static ImagePlus create_threshold_image(findblobs3 fb,float[] objects,float[] image,boolean showmask,int[] ids,boolean showids){
		ImageStack dispstack=new ImageStack(fb.width,fb.height);
		dispstack.addSlice("",image);
		dispstack.addSlice("",new float[fb.width*fb.height]);
		ImagePlus display=new ImagePlus("Outlined Objects",dispstack);
		display.setOpenAsHyperStack(true);
		display.setDimensions(2,1,1);
		display=new CompositeImage(display,CompositeImage.COMPOSITE);
		LUT graylut=jutils.get_lut_for_color(Color.white);
		double dispmin=jstatistics.getstatistic("Min",image,null);
		double dispmax=jstatistics.getstatistic("Max",image,null);
		graylut.min=dispmin;
		graylut.max=dispmax;
		((CompositeImage)display).setChannelLut(graylut,1);
		((CompositeImage)display).setDisplayRange(dispmin,dispmax);
		LUT redlut=jutils.get_lut_for_color(Color.red);
		redlut.min=0.0;
		redlut.max=255.0;
		((CompositeImage)display).setChannelLut(redlut,2);
		display.show();
		update_threshold_image(display,showmask,fb,objects,ids,showids);
		return display;
	}

	public static ImagePlus create_threshold_image(findblobs3 fb,float[] objects,float[] image,boolean showmask,int[] ids,boolean showids,Polygon[] outlines,float[][] coords){
		ImageStack dispstack=new ImageStack(fb.width,fb.height);
		dispstack.addSlice("",image);
		dispstack.addSlice("",new float[fb.width*fb.height]);
		ImagePlus display=new ImagePlus("Outlined Objects",dispstack);
		display.setOpenAsHyperStack(true);
		display.setDimensions(2,1,1);
		display=new CompositeImage(display,CompositeImage.COMPOSITE);
		LUT graylut=jutils.get_lut_for_color(Color.white);
		double dispmin=jstatistics.getstatistic("Min",image,null);
		double dispmax=jstatistics.getstatistic("Max",image,null);
		graylut.min=dispmin;
		graylut.max=dispmax;
		((CompositeImage)display).setChannelLut(graylut,1);
		((CompositeImage)display).setDisplayRange(dispmin,dispmax);
		LUT redlut=jutils.get_lut_for_color(Color.red);
		redlut.min=0.0;
		redlut.max=255.0;
		((CompositeImage)display).setChannelLut(redlut,2);
		display.show();
		update_threshold_image(display,showmask,fb,objects,ids,showids,outlines,coords);
		return display;
	}

	public static void update_threshold_image(ImagePlus imp,boolean showmask,findblobs3 fb,float[] objects,int[] ids,boolean showids){
		ImageProcessor ip=imp.getStack().getProcessor(2);
		float[] disppix=(float[])ip.getPixels();
		if(showmask){
			float[] temp=(float[])jutils.convert_array(fb.tobinary(objects,true),2);
			System.arraycopy(temp,0,disppix,0,disppix.length);
		}else{
			for(int i=0;i<disppix.length;i++){
				disppix[i]=0.0f;
			}
			Polygon[] objects2=fb.get_object_outlines(objects);
			ip.setValue(255);
			for(int i=0;i<objects2.length;i++){
				jutils.draw_polygon(ip,objects2[i],true);
			}
		}
		if(showids){
			if(ids==null){
				ids=new int[fb.nobjects];
				for(int i=0;i<fb.nobjects;i++)
					ids[i]=i+1;
			}
			ip.setFont(new Font("SansSerif",Font.PLAIN,10));
			ip.setJustification(ImageProcessor.CENTER_JUSTIFY);
			FontMetrics fm=ip.getFontMetrics();
			int stringheight=(int)(fm.getHeight()/2.0f);
			if(showmask)
				ip.setValue(0);
			float[][] coords=measure_object.centroids(objects,fb.width,fb.height);
			for(int i=0;i<fb.nobjects;i++){
				ip.drawString(""+ids[i],(int)coords[i][0],(int)coords[i][1]+stringheight);
			}
		}
		imp.updateAndDraw();
	}

	public static void update_threshold_image(ImagePlus imp,boolean showmask,findblobs3 fb,float[] objects,int[] ids,boolean showids,Polygon[] outlines,float[][] coords){
		ImageProcessor ip=imp.getStack().getProcessor(2);
		float[] disppix=(float[])ip.getPixels();
		if(showmask){
			float[] temp=(float[])jutils.convert_array(fb.tobinary(objects,true),2);
			System.arraycopy(temp,0,disppix,0,disppix.length);
		}else{
			for(int i=0;i<disppix.length;i++){
				disppix[i]=0.0f;
			}
			Polygon[] objects2=outlines;
			ip.setValue(255);
			for(int i=0;i<objects2.length;i++){
				jutils.draw_polygon(ip,objects2[i],true);
			}
		}
		if(showids){
			if(ids==null){
				ids=new int[fb.nobjects];
				for(int i=0;i<fb.nobjects;i++)
					ids[i]=i+1;
			}
			ip.setFont(new Font("SansSerif",Font.PLAIN,10));
			ip.setJustification(ImageProcessor.CENTER_JUSTIFY);
			FontMetrics fm=ip.getFontMetrics();
			int stringheight=(int)(fm.getHeight()/2.0f);
			if(showmask)
				ip.setValue(0);
			for(int i=0;i<ids.length;i++){
				ip.drawString(""+ids[i],(int)coords[i][0],(int)coords[i][1]+stringheight);
			}
		}
		imp.updateAndDraw();
	}
	
	/*************
	 * this is a simple implementation of the ImageJ rolling ball background subtraction
	 * @param image
	 * @param ballrad
	 * @param width
	 * @param height
	 * @return
	 */
	public static float[] sub_roll_ball_back(float[] image,float ballrad,int width,int height){
		FloatProcessor fp2=new FloatProcessor(width,height,image.clone(),null);
		fp2.snapshot();
		BackgroundSubtracter bs=new BackgroundSubtracter();
		bs.rollingBallBackground(fp2,ballrad,false,false,false,true,true);
		return (float[])fp2.getPixels();
	}
	
	/************
	 * does rolling ball for stacks
	 * @param stack
	 * @param ballrad
	 * @param width
	 * @param height
	 * @return
	 */
	public static float[][] sub_roll_ball_back(Object[] stack,float ballrad,int width,int height){
		float[][] retstack=new float[stack.length][];
		for(int i=0;i<stack.length;i++){
			retstack[i]=sub_roll_ball_back(algutils.convert_arr_float(stack[i]),ballrad,width,height);
		}
		return retstack;
	}
	
	public static float[] gaussian_blur(float[] image,float blurstdev,int width,int height) {
		FloatProcessor fp=new FloatProcessor(width,height,image.clone(),null);
		(new GaussianBlur()).blurFloat(fp,blurstdev,blurstdev,0.0002);
		return (float[])fp.getPixels();
	}
	
	public static void load_plugins_list(){
		ClassLoader cl=IJ.getClassLoader();
		/*IJ.log(cl.toString());
		try{
			Enumeration<URL> resource=cl.getResources("Plugins_List.xls");
			IJ.log(resource.nextElement().toString());
		}catch(IOException e){
			return;
		}*/
		InputStream is=cl.getResourceAsStream("Plugins_List.xls");
		jdataio jdio=new jdataio();
		ByteArrayOutputStream buffer=new ByteArrayOutputStream();
		int len=0;
		byte[] data=new byte[16384];
		try{
			while((len=is.read(data,0,data.length))!=-1){
				buffer.write(data,0,len);
			}
			buffer.flush();
		}catch(IOException e){
			IJ.log(jdio.getExceptionTrace(e));
			return;
		}
		List<List<String>> listtable=table_tools.table2listtable(buffer.toString(),"\t");
		String[] col_labels=table_tools.list2stringarray(listtable.get(0));
		listtable.remove(0);
		table_tools.create_table("Jay_Plugins_List",listtable,col_labels);
	}

	public static void run_command_in_IJ_thread(String command,String args){
		pluginrunner pr=new pluginrunner(command,args);
		Thread thread=new Thread(pr);
		thread.start();
	}

}

class pluginrunner implements Runnable{
	private final String command;
	private final String args;

	public pluginrunner(String command,String args){
		this.command=command;
		this.args=args;
	}

	public void run(){
		IJ.run(command,args);
	}
}
