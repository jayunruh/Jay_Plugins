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

public class combine_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow[] iw=jutils.selectPlotFamily(false,2,new String[]{"Destination","Source"});
		if(iw==null) return;
		if(iw[0].getClass().getName().equals("jguis.PlotWindow4")){
			PlotWindow4 pw=jutils.getPW4Copy(iw[0]);
			int origser=(Integer)jutils.runPW4VoidMethod(iw[0],"getNSeries");
			float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iw[1],"getXValues");
			float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iw[1],"getYValues");
			int[] colors=(int[])jutils.runPW4VoidMethod(iw[1],"getColors");
			int[] shapes=(int[])jutils.runPW4VoidMethod(iw[1],"getShapes");
			int[] npts=(int[])jutils.runPW4VoidMethod(iw[1],"getNpts");
			for(int i=0;i<yvals2.length;i++){
				float[] newxvals=new float[npts[i]];
				System.arraycopy(xvals2[i],0,newxvals,0,npts[i]);
				float[] newyvals=new float[npts[i]];
				System.arraycopy(yvals2[i],0,newyvals,0,npts[i]);
				pw.addPoints(newxvals,newyvals,true);
			}
			int[] newcolors=pw.getColors();
			int[] newshapes=pw.getShapes();
			for(int i=origser;i<newcolors.length;i++){
				newcolors[i]=colors[i-origser];
				newshapes[i]=shapes[i-origser];
			}
			pw.updatePlot();
			String newtitle=pw.getTitle()+"-1";
			pw.getImagePlus().setTitle(newtitle);
		} else {
			//assume for now that this is a 3D trajectory (not a 3D mesh)
			float[][] xvals1=(float[][])jutils.runPW4VoidMethod(iw[0],"getXValues");
			float[][] yvals1=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
			float[][][] zvals1=(float[][][])jutils.runPW4VoidMethod(iw[0],"getZValues");
			int[][] npts1=(int[][])jutils.runPW4VoidMethod(iw[0],"getNpts");
			float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iw[1],"getXValues");
			float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iw[1],"getYValues");
			float[][][] zvals2=(float[][][])jutils.runPW4VoidMethod(iw[1],"getZValues");
			int[][] npts2=(int[][])jutils.runPW4VoidMethod(iw[1],"getNpts");
			int maxlen=xvals1[0].length;
			if(xvals2[0].length>maxlen) maxlen=xvals2[0].length;
			float[][] newxvals=new float[xvals1.length+xvals2.length][maxlen];
			float[][] newyvals=new float[xvals1.length+xvals2.length][maxlen];
			float[][] newzvals=new float[xvals1.length+xvals2.length][maxlen];
			int[] newnpts=new int[xvals1.length+xvals2.length];
			for(int i=0;i<xvals1.length;i++){
				newnpts[i]=npts1[0][i];
				System.arraycopy(xvals1[i],0,newxvals[i],0,newnpts[i]);
				System.arraycopy(yvals1[i],0,newyvals[i],0,newnpts[i]);
				System.arraycopy(zvals1[0][i],0,newzvals[i],0,newnpts[i]);
			}
			int temp=xvals1.length;
			for(int i=0;i<xvals2.length;i++){
				newnpts[temp+i]=npts2[0][i];
				System.arraycopy(xvals2[i],0,newxvals[temp+i],0,newnpts[temp+i]);
				System.arraycopy(yvals2[i],0,newyvals[temp+i],0,newnpts[temp+i]);
				System.arraycopy(zvals2[0][i],0,newzvals[temp+i],0,newnpts[temp+i]);
			}
			String xlab=(String)jutils.runPW4VoidMethod(iw[0],"getxLabel");
			String ylab=(String)jutils.runPW4VoidMethod(iw[0],"getyLabel");
			String zlab=(String)jutils.runPW4VoidMethod(iw[0],"getzLabel");
			Traj3D t3d=new Traj3D(xlab,ylab,zlab,newxvals,newyvals,newzvals,newnpts);
			new PlotWindow3D("Combined Plot",t3d).draw();
		}
	}
}
