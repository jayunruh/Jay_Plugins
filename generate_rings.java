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

public class generate_rings implements PlugIn {

	public void run(String arg) {
		//here we output a sphere
		GenericDialog gd=new GenericDialog("Options");
		double radius=17.0;
		gd.addNumericField("Radius",radius,5,15,null);
		double spacing=5.0;
		gd.addNumericField("Spacing",spacing,5,15,null);
		int size=256;
		gd.addNumericField("Image Size",size,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		radius=gd.getNextNumber();
		spacing=gd.getNextNumber();
		size=(int)gd.getNextNumber();
		Object[] images=new Object[size];
		for(int i=0;i<size;i++){images[i]=new float[size*size];}
		double dphi=1.0/radius;
		int y1=(int)((double)size/2.0-spacing/2.0);
		int y2=(int)((double)y1+spacing);
		for(double phi=0.0;phi<(2.0*Math.PI);phi+=dphi){
			int x=(int)(radius*Math.cos(phi)+(double)size/2.0);
			int z=(int)(radius*Math.sin(phi)+(double)size/2.0);
			((float[])images[z])[x+y1*size]=1.0f;
			((float[])images[z])[x+y2*size]=1.0f;
		}
		ImageStack retstack=new ImageStack(size,size);
		for(int i=0;i<size;i++){
			retstack.addSlice("",images[i]);
		}
		new ImagePlus("Rings",retstack).show();
	}

}
