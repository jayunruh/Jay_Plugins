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

public class combine_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			ImageWindow iw=imp.getWindow();
			if(iw.getClass().getName().equals("jguis.PlotWindow4")){titles[i]=imp.getTitle();}
			else{titles[i]="NA";}
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Destination",titles,titles[0]);
		gd.addChoice("Source",titles,titles[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		int index2 = gd.getNextChoiceIndex();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);
		ImageWindow iw1=imp.getWindow();
		ImageWindow iw2=imp2.getWindow();
		PlotWindow4 pw=jutils.getPW4Copy(iw1);
		float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iw2,"getXValues");
		float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iw2,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw2,"getNpts");
		for(int i=0;i<yvals2.length;i++){
			float[] newxvals=new float[npts[i]];
			System.arraycopy(xvals2[i],0,newxvals,0,npts[i]);
			float[] newyvals=new float[npts[i]];
			System.arraycopy(yvals2[i],0,newyvals,0,npts[i]);
			pw.addPoints(newxvals,newyvals,true);
		}
	}
}
