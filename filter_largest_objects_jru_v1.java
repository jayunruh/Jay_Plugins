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
import jalgs.*;
import jalgs.jseg.*;

public class filter_largest_objects_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		int nlargest=1;
		gd.addNumericField("#_of_largest_objects",nlargest,0);
		boolean newimage=true;
		gd.addCheckbox("Create_New_Image?",newimage);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		nlargest=(int)gd.getNextNumber();
		newimage=gd.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		findblobs3 fb=new findblobs3(width,height);
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		ImageStack retstack=new ImageStack(width,height);
		for(int i=0;i<size;i++){
			byte[] data=(byte[])stack.getProcessor(i+1).getPixels();
			float[] objects=fb.dofindblobs(data);
			int[] areas=fb.get_areas(objects);
			int nobjects=areas.length;
			int[] arearank=jsort.get_javasort_order(areas);
			float[] newobjects=new float[width*height];
			for(int j=0;j<nlargest;j++){
				if(j<nobjects){
					fb.copy_object(objects,(float)(arearank[nobjects-j-1]+1),newobjects,(float)(j+1));
				}
			}
			if(newimage){
				retstack.addSlice("",fb.tobinary(newobjects,true));
			} else {
				stack.setPixels(fb.tobinary(newobjects,true),i+1);
				//data=fb.tobinary(objects,true);
			}
		}
		if(newimage){
			new ImagePlus(imp.getTitle()+" filtered",retstack).show();
		} else {
			imp.setStack(null,stack);
			imp.updateAndDraw();
		}
	}

}
