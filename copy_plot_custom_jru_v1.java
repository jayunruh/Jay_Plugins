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
import ij.plugin.PlugIn;
import jguis.*;
import java.io.*;
import java.awt.datatransfer.*;

public class copy_plot_custom_jru_v1 implements PlugIn, ClipboardOwner {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("copy_column_titles",true);
		gd.addCheckbox("copy_first_xvals",true);
		gd.addCheckbox("copy_other_xvals",false);
		gd.addCheckbox("copy_first_yvals",true);
		gd.addCheckbox("copy_other_yvals",true);
		gd.addCheckbox("copy_errors",true);
		gd.addCheckbox("copy_upper_and_lower_errors",false);
		gd.addCheckbox("padding_as_NaN",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean copytitles=gd.getNextBoolean();
		boolean copyfirstx=gd.getNextBoolean();
		boolean copyotherx=gd.getNextBoolean();
		boolean copyfirsty=gd.getNextBoolean();
		boolean copyothery=gd.getNextBoolean();
		boolean copyerrs=gd.getNextBoolean();
		boolean copybotherrs=gd.getNextBoolean();
		boolean padnan=gd.getNextBoolean();
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[][][] errs=(float[][][])jutils.runPW4VoidMethod(iw,"getErrors");
		Clipboard systemClipboard = null;
		try {systemClipboard = IJ.getInstance().getToolkit().getSystemClipboard();}
		catch (Exception e) {systemClipboard = null; }
		if (systemClipboard==null)
			{IJ.error("Unable to copy to Clipboard."); return;}
		IJ.showStatus("Copying plot values...");
		StringBuffer sb = new StringBuffer();
		int length=yvals[0].length;
		int nseries=yvals.length;
		if(copytitles){
			for (int j=0; j<nseries; j++) {
				if(j==0 && copyfirstx) sb.append("x"+(j+1)+"\t");
				if(j==0 && copyfirsty) sb.append("y"+(j+1)+"\t");
				if(j>0 && copyotherx) sb.append("x"+(j+1)+"\t");
				if(j>0 && copyothery) sb.append("y"+(j+1)+"\t");
				if(copyerrs && errs!=null){
					if(copybotherrs) sb.append("err1\terr2\t");
					else sb.append("err1\t");
				}
			}
		}
		for (int i=0; i<length; i++) {
			for (int j=0; j<nseries; j++) {
				String xval=""+xvals[j][i];
				String yval=""+yvals[j][i];
				if(padnan){
					if(i>=npts[j]){xval="NaN"; yval="NaN";}
				}
				if(j==0 && copyfirstx) sb.append(xval+"\t");
				if(j==0 && copyfirsty) sb.append(yval+"\t");
				if(j>0 && copyotherx) sb.append(xval+"\t");
				if(j>0 && copyothery) sb.append(yval+"\t");
				if(copyerrs && errs!=null){
					String err1=""+errs[0][j][i];
					String err2=""+errs[1][j][i];
					if(padnan){
						if(i>=npts[j]){err1="NaN"; err2="NaN";}
					}
					if(copybotherrs) sb.append(err1+"\t"+err2+"\t");
					else sb.append(err1+"\t");
				}
				if(j==(nseries-1)){sb.deleteCharAt(sb.length()-1);}
			}
			sb.append("\n");
		}
		String text = sb.toString();
		StringSelection contents = new StringSelection(text);
		systemClipboard.setContents(contents, this);
		IJ.showStatus(text.length() + " characters copied to Clipboard");
	}

	public void lostOwnership(Clipboard clipboard, Transferable contents) {}

}
