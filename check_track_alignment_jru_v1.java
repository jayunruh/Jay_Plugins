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

public class check_track_alignment_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(true,2,new String[]{"Movie","Track"});
		ImagePlus imp = imps[0];
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		Roi roi=imp.getRoi();
		float[] xvals=null;
		float[] yvals=null;
		boolean hastrack;
		if(imps[1]!=null){
			ImageWindow iw=imps[1].getWindow();
			xvals=((float[][])jutils.runPW4VoidMethod(iw,"getXValues"))[0].clone();
			yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0].clone();
			hastrack=true;
		} else {
			Rectangle r=roi.getBounds();
			xvals=new float[slices]; yvals=new float[slices];
			for(int i=0;i<slices;i++){xvals[i]=r.x; yvals[i]=r.y;}
			hastrack=false;
		}
		for(int x=0;x<slices;x++){
			imp.setRoi(new PointRoi((int)xvals[x],(int)yvals[x]));
			imp.setSlice(x+1);
			(new WaitForUserDialog("Update point if necessary")).show();
			Rectangle tempr=imp.getRoi().getBounds();
			xvals[x]=tempr.x;
			yvals[x]=tempr.y;
			if(!hastrack && x<(slices-1)){
				xvals[x+1]=xvals[x];
				yvals[x+1]=yvals[x];
			}
		}
		new PlotWindow4("Updated_Track","x","y",xvals,yvals).draw();
	}

}
