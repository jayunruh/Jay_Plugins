/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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

public class get_imaris_dataset_jru_v1 implements PlugIn {

	public void run(String arg) {
			float[] dimensions=ImarisXT_utils.getDimensions();
			if(dimensions==null){
				IJ.log("Couldn't load stack dimensions");
				return;
			}
			int width=(int)dimensions[0];
			int height=(int)dimensions[1];
			int chans=(int)dimensions[2];
			int slices=(int)dimensions[3];
			int frames=(int)dimensions[4];
			float psize=dimensions[5];
			float pdepth=dimensions[6];
			Object[] dataset=ImarisXT_utils.getDataSet();
			if(dataset==null){
				IJ.log("Couldn't load dataset");
				return;
			}
			ImageStack stack=jutils.array2stack(dataset,width,height);
			ImagePlus imp=new ImagePlus("Imaris Dataset",stack);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(chans,slices,frames);
			jutils.set_psize(imp,psize);
			jutils.set_pdepth(imp,pdepth);
			if(chans>1) (new CompositeImage(imp,CompositeImage.COMPOSITE)).show();
			else imp.show();

	}

}
