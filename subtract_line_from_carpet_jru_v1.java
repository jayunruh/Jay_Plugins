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

public class subtract_line_from_carpet_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){
				titles[i]=imp.getTitle();
			} else {
				titles[i]="NA";
			}
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Plot Window",titles,titles[0]);
		gd.addChoice("Carpet",titles,titles[0]);
		gd.addCheckbox("Horizontal_Line",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		int index2 = gd.getNextChoiceIndex();
		boolean horizontal=gd.getNextBoolean();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);

		ImageWindow iw=imp.getWindow();
		float[][] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"));
		int[] npts=((int[])jutils.runPW4VoidMethod(iw,"getNpts"));
		int length=imp2.getHeight();
		int width=imp2.getWidth();
		ImageStack stack=imp2.getStack();
		int slices=stack.getSize();
		ImageStack outstack=new ImageStack(width,length);
		if(!horizontal){
			for(int k=0;k<slices;k++){
				float[] subimage=new float[width*length];
				float[] orig=jutils.convert_arr_float(stack.getPixels(k+1));
				float[] line=yvals[0];
				if(yvals.length>1) line=yvals[k];
				for(int i=0;i<length;i++){
					for(int j=0;j<width;j++){
						subimage[j+i*width]=orig[j+i*width]-line[i];
					}
				}
				outstack.addSlice("",subimage);
			}
		} else {
			for(int k=0;k<slices;k++){
				float[] subimage=new float[width*length];
				float[] orig=jutils.convert_arr_float(stack.getPixels(k+1));
				float[] line=yvals[0];
				if(yvals.length>1) line=yvals[k];
				for(int i=0;i<length;i++){
					for(int j=0;j<width;j++){
						subimage[j+i*width]=orig[j+i*width]-line[j];
					}
				}
				outstack.addSlice("",subimage);
			}
		}
		(new ImagePlus("Subtracted",outstack)).show();
	}

}
