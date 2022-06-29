/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.text.*;

public class filter_point_pairs_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow[] iws=jutils.selectPlotFamily(true,2,new String[]{"source","filter"});
		if(iws==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Distance Threshold",10.0f,5,15,null);
		gd.addCheckbox("Output Destination Points",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float distthresh=(float)gd.getNextNumber();
		boolean outdest=gd.getNextBoolean();
		float[][] xvals1=(float[][])jutils.runPW4VoidMethod(iws[0],"getXValues");
		float[][] yvals1=(float[][])jutils.runPW4VoidMethod(iws[0],"getYValues");
		if(iws[1]==null){
		if(iws[0].getClass().getName().equals("jguis.PlotWindow3D")){
			int[] npts1=((int[][])jutils.runPW4VoidMethod(iws[0],"getNpts"))[0];
			float[][][] zvals1=(float[][][])jutils.runPW4VoidMethod(iws[0],"getZValues");
			if(npts1.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals1,yvals1,zvals1,npts1,0);
				xvals1=(float[][])temp[0]; yvals1=(float[][])temp[1]; zvals1=(float[][][])temp[2]; npts1=(int[])temp[3];
			}
			Object[] temp=getClosestPoints(xvals1,yvals1,zvals1);
			float[] mindist=(float[])temp[0];
			int[] dest=(int[])temp[1];
			StringBuffer sb=new StringBuffer();
			int npairs=0;
			for(int i=0;i<xvals1.length;i++){
				sb.append(mindist[i]+"\t"+dest[i]+"\n");
				if(mindist[i]<=distthresh) npairs++;
			}
			sb.deleteCharAt(sb.length()-1);
			new TextWindow("Min_Distances","Dist\tDest_Index",sb.toString(),200,200);
			showPairs(xvals1,yvals1,zvals1,mindist,distthresh,npairs);
			if(outdest){
				showDest(xvals1,yvals1,zvals1,mindist,distthresh,dest);
			}
		} else {
			int[] npts1=(int[])jutils.runPW4VoidMethod(iws[0],"getNpts");
			if(npts1.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals1,yvals1,npts1,0);
				xvals1=(float[][])temp[0]; yvals1=(float[][])temp[1]; npts1=(int[])temp[2];
			}
			Object[] temp=getClosestPoints(xvals1,yvals1);
			float[] mindist=(float[])temp[0];
			int[] dest=(int[])temp[1];
			StringBuffer sb=new StringBuffer();
			int npairs=0;
			for(int i=0;i<xvals1.length;i++){
				sb.append(mindist[i]+"\t"+dest[i]+"\n");
				if(mindist[i]<=distthresh) npairs++;
			}
			sb.deleteCharAt(sb.length()-1);
			new TextWindow("Min_Distances","Dist\tDest_Index",sb.toString(),200,200);
			showPairs(xvals1,yvals1,mindist,distthresh,npairs);
			if(outdest){
				showDest(xvals1,yvals1,mindist,distthresh,dest);
			}
		}
		} else {
		float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iws[1],"getXValues");
		float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iws[1],"getYValues");
		if(iws[0].getClass().getName().equals("jguis.PlotWindow3D")){
			int[] npts1=((int[][])jutils.runPW4VoidMethod(iws[0],"getNpts"))[0];
			int[] npts2=((int[][])jutils.runPW4VoidMethod(iws[1],"getNpts"))[0];
			float[][][] zvals1=(float[][][])jutils.runPW4VoidMethod(iws[0],"getZValues");
			float[][][] zvals2=(float[][][])jutils.runPW4VoidMethod(iws[1],"getZValues");
			if(npts1.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals1,yvals1,zvals1,npts1,0);
				xvals1=(float[][])temp[0]; yvals1=(float[][])temp[1]; zvals1=(float[][][])temp[2]; npts1=(int[])temp[3];
			}
			if(npts2.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals2,yvals2,zvals2,npts2,0);
				xvals2=(float[][])temp[0]; yvals2=(float[][])temp[1]; zvals2=(float[][][])temp[2]; npts2=(int[])temp[3];
			}
			Object[] temp=getClosestPoints(xvals1,yvals1,zvals1,xvals2,yvals2,zvals2);
			float[] mindist=(float[])temp[0];
			int[] dest=(int[])temp[1];
			StringBuffer sb=new StringBuffer();
			int npairs=0;
			for(int i=0;i<xvals1.length;i++){
				sb.append(mindist[i]+"\t"+dest[i]+"\n");
				if(mindist[i]<=distthresh) npairs++;
			}
			sb.deleteCharAt(sb.length()-1);
			new TextWindow("Min_Distances","Dist\tDest_Index",sb.toString(),200,200);
			showPairs(xvals1,yvals1,zvals1,mindist,distthresh,npairs);
			if(outdest){
				showDest(xvals2,yvals2,zvals2,mindist,distthresh,dest);
			}
		} else {
			int[] npts1=(int[])jutils.runPW4VoidMethod(iws[0],"getNpts");
			int[] npts2=(int[])jutils.runPW4VoidMethod(iws[1],"getNpts");
			if(npts1.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals1,yvals1,npts1,0);
				xvals1=(float[][])temp[0]; yvals1=(float[][])temp[1]; npts1=(int[])temp[2];
			}
			if(npts2.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals2,yvals2,npts2,0);
				xvals2=(float[][])temp[0]; yvals2=(float[][])temp[1]; npts2=(int[])temp[2];
			}

			Object[] temp=getClosestPoints(xvals1,yvals1,xvals2,yvals2);
			float[] mindist=(float[])temp[0];
			int[] dest=(int[])temp[1];
			StringBuffer sb=new StringBuffer();
			int npairs=0;
			for(int i=0;i<xvals1.length;i++){
				sb.append(mindist[i]+"\t"+dest[i]+"\n");
				if(mindist[i]<=distthresh) npairs++;
			}
			sb.deleteCharAt(sb.length()-1);
			new TextWindow("Min_Distances","Dist\tDest_Index",sb.toString(),200,200);
			showPairs(xvals1,yvals1,mindist,distthresh,npairs);
			if(outdest){
				showDest(xvals2,yvals2,mindist,distthresh,dest);
			}
		}
		}
	}

	public void showDest(float[][] xvals2,float[][] yvals2,float[][][] zvals2,float[] mindist,float distthresh,int[] dest){
		int npairs=dest.length;
		float[][] xvals=new float[npairs][1];
		float[][] yvals=new float[npairs][1];
		float[][] zvals=new float[npairs][1];
		int[] newnpts=new int[npairs];
		int counter=0;
		for(int i=0;i<xvals2.length;i++){
			if(mindist[i]<=distthresh){
				int d=dest[i];
				xvals[counter][0]=xvals2[d][0]; yvals[counter][0]=yvals2[d][0]; zvals[counter][0]=zvals2[0][d][0]; newnpts[counter]=1;
				counter++;
			}
		}
		Traj3D t3D=new Traj3D("x","y","z",xvals,yvals,zvals,newnpts);
		int[] shapes=t3D.getShapes(); for(int i=0;i<shapes.length;i++) shapes[i]=1;
		new PlotWindow3D("Filtered 3D Destination Points",t3D).draw();
	}

	public void showDest(float[][] xvals2,float[][] yvals2,float[] mindist,float distthresh,int[] dest){
		int npairs=dest.length;
		float[][] xvals=new float[npairs][1];
		float[][] yvals=new float[npairs][1];
		int[] newnpts=new int[npairs];
		int counter=0;
		for(int i=0;i<xvals2.length;i++){
			if(mindist[i]<=distthresh){
				int d=dest[i];
				xvals[counter][0]=xvals2[d][0]; yvals[counter][0]=yvals2[d][0]; newnpts[counter]=1;
				counter++;
			}
		}
		Plot4 plot=new Plot4("x","y",xvals,yvals,newnpts);
		int[] shapes=plot.getShapes(); for(int i=0;i<shapes.length;i++) shapes[i]=1;
		new PlotWindow4("Filtered Destination Points",plot).draw();
	}

	public void showPairs(float[][] xvals1,float[][] yvals1,float[][][] zvals1,float[] mindist,float distthresh,int npairs){
		float[][] xvals=new float[npairs][1];
		float[][] yvals=new float[npairs][1];
		float[][] zvals=new float[npairs][1];
		int[] newnpts=new int[npairs];
		int counter=0;
		for(int i=0;i<xvals1.length;i++){
			if(mindist[i]<=distthresh){
				xvals[counter][0]=xvals1[i][0]; yvals[counter][0]=yvals1[i][0]; zvals[counter][0]=zvals1[0][i][0]; newnpts[counter]=1;
				counter++;
			}
		}
		Traj3D t3D=new Traj3D("x","y","z",xvals,yvals,zvals,newnpts);
		int[] shapes=t3D.getShapes(); for(int i=0;i<shapes.length;i++) shapes[i]=1;
		new PlotWindow3D("Filtered 3D Points",t3D).draw();
	}

	public void showPairs(float[][] xvals1,float[][] yvals1,float[] mindist,float distthresh,int npairs){
		float[][] xvals=new float[npairs][1];
		float[][] yvals=new float[npairs][1];
		int[] newnpts=new int[npairs];
		int counter=0;
		for(int i=0;i<xvals1.length;i++){
			if(mindist[i]<=distthresh){
				xvals[counter][0]=xvals1[i][0]; yvals[counter][0]=yvals1[i][0]; newnpts[counter]=1;
				counter++;
			}
		}
		Plot4 plot=new Plot4("x","y",xvals,yvals,newnpts);
		int[] shapes=plot.getShapes(); for(int i=0;i<shapes.length;i++) shapes[i]=1;
		new PlotWindow4("Filtered Points",plot).draw();
	}

	public Object[] getClosestPoints(float[][] xvals1,float[][] yvals1,float[][][] zvals1,float[][] xvals2,float[][] yvals2,float[][][] zvals2){
		//assume all points are singles
		float[] mindist=new float[xvals1.length];
		int[] dest=new int[xvals1.length];
		int npts=0;
		for(int i=0;i<xvals1.length;i++){
			mindist[i]=dist(xvals1[i][0],xvals2[0][0],yvals1[i][0],yvals2[0][0],zvals1[0][i][0],zvals2[0][0][0]);
			for(int j=1;j<xvals2.length;j++){
				float dist1=dist(xvals1[i][0],xvals2[j][0],yvals1[i][0],yvals2[j][0],zvals1[0][i][0],zvals2[0][j][0]);
				if(dist1<mindist[i]){
					mindist[i]=dist1;
					dest[i]=j;
				}
			}
		}
		return new Object[]{mindist,dest};
	}

	public Object[] getClosestPoints(float[][] xvals1,float[][] yvals1,float[][][] zvals1){
		//assume all points are singles
		float[] mindist=new float[xvals1.length];
		int[] dest=new int[xvals1.length];
		int npts=0;
		for(int i=0;i<xvals1.length;i++){
			if(i>0) mindist[i]=dist(xvals1[i][0],xvals1[i-1][0],yvals1[i][0],yvals1[i-1][0],zvals1[0][i][0],zvals1[0][i-1][0]);
			else mindist[i]=dist(xvals1[i][0],xvals1[i+1][0],yvals1[i][0],yvals1[i+1][0],zvals1[0][i][0],zvals1[0][i+1][0]);
			for(int j=0;j<xvals1.length;j++){
				if(j!=i){
					float dist1=dist(xvals1[i][0],xvals1[j][0],yvals1[i][0],yvals1[j][0],zvals1[0][i][0],zvals1[0][j][0]);
					if(dist1<mindist[i]){
						mindist[i]=dist1;
						dest[i]=j;
					}
				}
			}
		}
		return new Object[]{mindist,dest};
	}

	public Object[] getClosestPoints(float[][] xvals1,float[][] yvals1,float[][] xvals2,float[][] yvals2){
		//assume all points are singles
		float[] mindist=new float[xvals1.length];
		int[] dest=new int[xvals1.length];
		int npts=0;
		for(int i=0;i<xvals1.length;i++){
			mindist[i]=dist(xvals1[i][0],xvals2[0][0],yvals1[i][0],yvals2[0][0]);
			for(int j=1;j<xvals2.length;j++){
				float dist1=dist(xvals1[i][0],xvals2[j][0],yvals1[i][0],yvals2[j][0]);
				if(dist1<mindist[i]){
					mindist[i]=dist1;
					dest[i]=j;
				}
			}
		}
		return new Object[]{mindist,dest};
	}

	public Object[] getClosestPoints(float[][] xvals1,float[][] yvals1){
		//assume all points are singles
		float[] mindist=new float[xvals1.length];
		int[] dest=new int[xvals1.length];
		int npts=0;
		for(int i=0;i<xvals1.length;i++){
			if(i>0) mindist[i]=dist(xvals1[i][0],xvals1[i-1][0],yvals1[i][0],yvals1[i-1][0]);
			else mindist[i]=dist(xvals1[i][0],xvals1[i+1][0],yvals1[i][0],yvals1[i+1][0]);
			for(int j=0;j<xvals1.length;j++){
				if(j!=i){
					float dist1=dist(xvals1[i][0],xvals1[j][0],yvals1[i][0],yvals1[j][0]);
					if(dist1<mindist[i]){
						mindist[i]=dist1;
						dest[i]=j;
					}
				}
			}
		}
		return new Object[]{mindist,dest};
	}

	public Object[] separatePoints(float[][] xpts,float[][] ypts,float[][][] zpts,int[] npts,int sel){
		//here we take trajectories and turn them into points
		float[][] newxpts=new float[npts[sel]][1];
		float[][] newypts=new float[npts[sel]][1];
		float[][][] newzpts=new float[1][npts[sel]][1];
		int[] newnpts=new int[npts[sel]];
		for(int i=0;i<npts[sel];i++){
			newxpts[i][0]=xpts[sel][i];
			newypts[i][0]=ypts[sel][i];
			newzpts[0][i][0]=zpts[0][sel][i];
			newnpts[i]=1;
		}
		return new Object[]{newxpts,newypts,newzpts,newnpts};
	}

	public Object[] separatePoints(float[][] xpts,float[][] ypts,int[] npts,int sel){
		//here we take trajectories and turn them into points
		float[][] newxpts=new float[npts[sel]][1];
		float[][] newypts=new float[npts[sel]][1];
		int[] newnpts=new int[npts[sel]];
		for(int i=0;i<npts[sel];i++){
			newxpts[i][0]=xpts[sel][i];
			newypts[i][0]=ypts[sel][i];
			newnpts[i]=1;
		}
		return new Object[]{newxpts,newypts,newnpts};
	}

	public float dist(float x1,float x2,float y1,float y2,float z1,float z2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	}

	public float dist(float x1,float x2,float y1,float y2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}

}
