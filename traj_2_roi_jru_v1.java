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

public class traj_2_roi_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we create a segmented line roi from a plot
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[] xvals=((float[][])jutils.runPW4VoidMethod(iw,"getXValues"))[0];
		float[] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0];
		int[] intxvals=new int[xvals.length];
		int[] intyvals=new int[yvals.length];
		for(int i=0;i<xvals.length;i++){
			intxvals[i]=(int)xvals[i];
			intyvals[i]=(int)yvals[i];
		}
		PolygonRoi roi=new PolygonRoi(intxvals,intyvals,xvals.length,Roi.FREELINE);
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			titles[i]=imp.getTitle();
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Image",titles,titles[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		imp.setRoi(roi);
		imp.updateAndDraw();
	}

}
