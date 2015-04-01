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
		gd.addNumericField("Series_to_combine",0,0);
		gd.addCheckbox("Delete_Originals",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean combineall=gd.getNextBoolean();
		int combseries=(int)gd.getNextNumber();
		boolean delor=gd.getNextBoolean();
		int[] wList = WindowManager.getIDList();
		ImageWindow[] windows=new ImageWindow[wList.length];
		int nplots=0;
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			ImageWindow iw=imp.getWindow();
			if(iw.getClass().getName().equals("jguis.PlotWindow4") || iw.getClass().getName().equals("ij.gui.PlotWindow")){windows[nplots]=iw; nplots++;}
		}
		if(combineall){
			PlotWindow4 pw=jutils.getPW4Copy(windows[0]);
			for(int j=1;j<nplots;j++){
				float[][] xvals2=(float[][])jutils.runPW4VoidMethod(windows[j],"getXValues");
				float[][] yvals2=(float[][])jutils.runPW4VoidMethod(windows[j],"getYValues");
				int[] npts=(int[])jutils.runPW4VoidMethod(windows[j],"getNpts");
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
		} else {
			PlotWindow4 pw=jutils.getPW4SelCopy(windows[0],combseries);
			for(int j=1;j<nplots;j++){
				float[][] xvals2=(float[][])jutils.runPW4VoidMethod(windows[j],"getXValues");
				float[][] yvals2=(float[][])jutils.runPW4VoidMethod(windows[j],"getYValues");
				int[] npts=(int[])jutils.runPW4VoidMethod(windows[j],"getNpts");
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
		if(delor){
			for(int i=0;i<nplots;i++){
				windows[i].close();
			}
		}
	}
}
