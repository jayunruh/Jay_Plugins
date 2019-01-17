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

public class combine_all_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Combine_all_series",true);
		gd.addNumericField("Series_to_combine",1,0);
		gd.addCheckbox("Delete_Originals",false);
		gd.addCheckbox("Select_Plots",false);
		gd.addNumericField("Number_of_Plots (if selecting)",2,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean combineall=gd.getNextBoolean();
		int combseries=(int)gd.getNextNumber()-1;
		boolean delor=gd.getNextBoolean();
		boolean selplots=gd.getNextBoolean();
		int nplots=(int)gd.getNextNumber();
		ImageWindow[] windows=null;
		if(selplots){
			windows=jutils.selectPlotFamily(false,nplots);
			if(windows==null) return;
		} else {
			int[] wList = WindowManager.getIDList();
			windows=new ImageWindow[wList.length];
			nplots=0;
			for(int i=0;i<wList.length;i++){
				ImagePlus imp = WindowManager.getImage(wList[i]);
				ImageWindow iw=imp.getWindow();
				if(jutils.isPlotFamily(iw)){windows[nplots]=iw; nplots++;}
			}
		}
		if(combineall){
			if(windows[0].getClass().getName().equals("jguis.PlotWindow3D")){
				Plot3D plot=jutils.getPW3DPlotCopy(windows[0]);
				for(int j=1;j<nplots;j++){
					Object[] pdata=getPlotValues(windows[j]); 
					float[][] xvals2=(float[][])pdata[0]; 
					float[][] yvals2=(float[][])pdata[1];
					float[][][] zvals2=(float[][][])pdata[2];
					int[][] npts=(int[][])pdata[3];
					if(plot.getClass().getName().equals("jguis.Traj3D")){
						for(int i=0;i<yvals2.length;i++){
							float[] newxvals=new float[npts[0][i]];
							System.arraycopy(xvals2[i],0,newxvals,0,npts[0][i]);
							float[] newyvals=new float[npts[0][i]];
							System.arraycopy(yvals2[i],0,newyvals,0,npts[0][i]);
							float[] newzvals=new float[npts[0][i]];
							System.arraycopy(zvals2[0][i],0,newzvals,0,npts[0][i]);
							((Traj3D)plot).addPoints(newxvals,newyvals,newzvals,true);
						}
					} else {
						for(int i=0;i<yvals2.length;i++){
							float[] newxvals=new float[npts[i][0]];
							System.arraycopy(xvals2[i],0,newxvals,0,npts[i][0]);
							float[] newyvals=new float[npts[i][1]];
							System.arraycopy(yvals2[i],0,newyvals,0,npts[i][1]);
							float[][] newzvals=new float[npts[i][0]][npts[i][1]];
							for(int k=0;k<npts[i][0];k++) System.arraycopy(zvals2[i][k],0,newzvals[k],0,npts[i][1]);
							plot.addPoints(newxvals,newyvals,newzvals,true);
						}
					}
				}
				(new PlotWindow3D("Plot Copy",plot)).draw();
			} else {	
				PlotWindow4 pw=jutils.getPW4Copy(windows[0]);
				for(int j=1;j<nplots;j++){
					Object[] pdata=getPlotValues(windows[j]); 
					float[][] xvals2=(float[][])pdata[0]; 
					float[][] yvals2=(float[][])pdata[1];
					int[] npts=(int[])pdata[2];
					boolean showerrs=(Boolean)jutils.runPW4VoidMethod(windows[j],"getShowErrors");
					float[][][] errs=null;
					if(showerrs) errs=(float[][][])jutils.runPW4VoidMethod(windows[j],"getErrors");
					for(int i=0;i<yvals2.length;i++){
						float[] newxvals=new float[npts[i]];
						System.arraycopy(xvals2[i],0,newxvals,0,npts[i]);
						float[] newyvals=new float[npts[i]];
						System.arraycopy(yvals2[i],0,newyvals,0,npts[i]);
						pw.addPoints(newxvals,newyvals,true);
						if(showerrs && errs!=null){
							pw.addSeriesErrors(pw.getNSeries()-1,new float[][]{errs[0][i],errs[1][i]});
						}
					}
				}
			}
		} else {
			if(windows[0].getClass().getName().equals("jguis.PlotWindow3D")){
				Plot3D plot=jutils.getPW3DPlotCopy(windows[0]);
				for(int j=1;j<nplots;j++){
					Object[] pdata=getPlotValues(windows[j]); 
					float[][] xvals2=(float[][])pdata[0]; 
					float[][] yvals2=(float[][])pdata[1];
					float[][][] zvals2=(float[][][])pdata[2];
					int[][] npts=(int[][])pdata[3];
					if(plot.getClass().getName().equals("jguis.Traj3D")){
						for(int i=0;i<yvals2.length;i++){
							float[] newxvals=new float[npts[0][i]];
							System.arraycopy(xvals2[i],0,newxvals,0,npts[0][i]);
							float[] newyvals=new float[npts[0][i]];
							System.arraycopy(yvals2[i],0,newyvals,0,npts[0][i]);
							float[] newzvals=new float[npts[0][i]];
							System.arraycopy(zvals2[0][i],0,newzvals,0,npts[0][i]);
							((Traj3D)plot).addPoints(newxvals,newyvals,newzvals,true);
						}
					} else {
						for(int i=0;i<yvals2.length;i++){
							float[] newxvals=new float[npts[i][0]];
							System.arraycopy(xvals2[i],0,newxvals,0,npts[i][0]);
							float[] newyvals=new float[npts[i][1]];
							System.arraycopy(yvals2[i],0,newyvals,0,npts[i][1]);
							float[][] newzvals=new float[npts[i][0]][npts[i][1]];
							for(int k=0;k<npts[i][0];k++) System.arraycopy(zvals2[i][k],0,newzvals[k],0,npts[i][1]);
							plot.addPoints(newxvals,newyvals,newzvals,true);
						}
					}
				}
				(new PlotWindow3D("Plot Copy",plot)).draw();
			} else {
				PlotWindow4 pw=jutils.getPW4SelCopy(windows[0],combseries);
				for(int j=1;j<nplots;j++){
					Object[] pdata=getPlotValues(windows[j]); 
					float[][] xvals2=(float[][])pdata[0]; 
					float[][] yvals2=(float[][])pdata[1];
					int[] npts=(int[])pdata[2];
					float[] newxvals=new float[npts[combseries]];
					System.arraycopy(xvals2[combseries],0,newxvals,0,npts[combseries]);
					float[] newyvals=new float[npts[combseries]];
					System.arraycopy(yvals2[combseries],0,newyvals,0,npts[combseries]);
					pw.addPoints(newxvals,newyvals,true);
					boolean showerrs=(Boolean)jutils.runPW4VoidMethod(windows[j],"getShowErrors");
					if(showerrs){
						float[][][] errs=(float[][][])jutils.runPW4VoidMethod(windows[j],"getErrors");
						if(errs!=null) pw.addSeriesErrors(pw.getNSeries()-1,new float[][]{errs[0][combseries],errs[1][combseries]});
					}
				}
			}
		}
		if(delor){
			for(int i=0;i<nplots;i++){
				windows[i].close();
			}
		}
	}

	public Object[] getPlotValues(ImageWindow iw){
		if(iw.getClass().getName().equals("jguis.PlotWindowHist")){
			Object plot=jutils.runReflectionMethod(iw,"getPlot",null,null);
			float[][] hist=(float[][])jutils.runReflectionMethod(plot,"getHistogram",null,null);
			int[] npts={hist[0].length};
			return new Object[]{new float[][]{hist[0]},new float[][]{hist[1]},npts};
		} else if(iw.getClass().getName().equals("jguis.PlotWindow3D")){
			float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			float[][][] zvals2=(float[][][])jutils.runPW4VoidMethod(iw,"getZValues");
			int[][] npts=(int[][])jutils.runPW4VoidMethod(iw,"getNpts");
			return new Object[]{xvals2,yvals2,zvals2,npts};
		} else {
			float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
			return new Object[]{xvals2,yvals2,npts};
		}
	}
}
