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

public class stack_editor_jru_v1 implements PlugIn {
//this plugin deletes slices from a stack
	public void run(String arg) {
		ImagePlus imp = WindowManager.getCurrentImage();
		ImageStack stack = imp.getStack();
		int slices = stack.getSize();
		GenericDialog gd = new GenericDialog("Edit Stack");
		String[] temp = {"Delete_<=_Slice_a","Delete_>=_Slice_b","a_<=_Delete_=>_b"};
		gd.addChoice("Action",temp,temp[0]);
		gd.addNumericField("Slice_a",1.0,0);
		gd.addNumericField("Slice_b",(double)slices,0);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int action = gd.getNextChoiceIndex();
		int a = (int)gd.getNextNumber();
		int b = (int)gd.getNextNumber();
		if(b>slices){b=slices;}
		if(a<1){a=1;}
		if(a>b){a=b;}
		if(!imp.lock()){return;}
		if(action==0) {for(int i=1; i<=a; i++){stack.deleteSlice(1); IJ.showProgress(i,a);}}
		if(action==1) {for(int i=b;i<=slices;i++){stack.deleteLastSlice(); IJ.showProgress(i,b);}}
		if(action==2) {for(int i=a;i<=b;i++){stack.deleteSlice(a); IJ.showProgress(i,a);}}
		imp.setStack(null,stack);
		imp.setSlice(1);
		imp.unlock();
	}

}
