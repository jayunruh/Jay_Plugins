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

public class edit_plot_colors_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		int[] colors=(int[])jutils.runPW4VoidMethod(iw,"getColors");
		int selected=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
		String[] clabels={"black","blue","green","red","magenta","cyan","yellow","orange"};
		GenericDialog gd=new GenericDialog("Options");
		int length=colors.length;
		if(length>clabels.length){
			length=clabels.length;
		}
		if(selected>=0){
			gd.addChoice("Series "+(selected+1)+" color",clabels,clabels[colors[selected]%8]);
			gd.showDialog(); if(gd.wasCanceled()){return;}
			colors[selected]=gd.getNextChoiceIndex();
			jutils.runPW4VoidMethod(iw,"updatePlot");
		} else {
			for(int i=0;i<length;i++){
				gd.addChoice("Series"+(i+1)+" color",clabels,clabels[colors[i]]);
			}
			gd.showDialog(); if(gd.wasCanceled()){return;}
			for(int i=0;i<colors.length;i++){
				colors[i]=gd.getNextChoiceIndex();
			}
			jutils.runPW4VoidMethod(iw,"updatePlot");
		}
	}

}
