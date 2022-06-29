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
import java.awt.Frame;
import java.awt.Rectangle;
import ij.plugin.*;
import jguis.*;
import java.util.*;

public class man_track_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp =WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int frames=imp.getNFrames();
		boolean zstack=(frames==1);
		if(zstack) frames=imp.getNSlices();
		int startframe=imp.getFrame()-1;
		int currslice=imp.getSlice();
		int currchan=imp.getChannel();
		if(zstack){
			int temp=startframe+1;
			startframe=currslice-1;
			currslice=temp;
		}
		Roi roi=imp.getRoi();
		if(roi==null){
			IJ.error("put a point roi at the starting point");
			return;
		}
		Rectangle r=roi.getBounds();
		List<float[]> traj=new ArrayList<float[]>();
		int currx=r.x;
		int curry=r.y;
		for(int x=startframe;x<frames;x++){
			imp.setRoi(new PointRoi(currx,curry));
			if(zstack) imp.setPosition(currchan,x+1,currslice);
			else imp.setPosition(currchan,currslice,x+1);
			WaitForUserDialog dg=new WaitForUserDialog("Update point if necessary (esc to abort)");
			dg.show();
			if(dg.escPressed()) break;
			Rectangle tempr=imp.getRoi().getBounds();
			traj.add(new float[]{tempr.x,tempr.y});
			currx=tempr.x; curry=tempr.y;
		}
		float[] xvals=new float[traj.size()];
		float[] yvals=new float[traj.size()];
		for(int i=0;i<traj.size();i++){
			xvals[i]=traj.get(i)[0];
			yvals[i]=traj.get(i)[1];
		}
		PlotWindow4 pw=new PlotWindow4("Track","x","y",xvals,yvals);
		pw.p3.setAnnotations(new String[]{""+startframe});
		pw.draw();
	}

}
