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

public class filter_point_triples_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImageWindow[] iws=jutils.selectPlotFamily(false,3,new String[]{"source","filter1","filter2"});
		if(iws==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Distance Threshold",10.0f,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float distthresh=(float)gd.getNextNumber();
		float[][] xvals1=(float[][])jutils.runPW4VoidMethod(iws[0],"getXValues");
		float[][] yvals1=(float[][])jutils.runPW4VoidMethod(iws[0],"getYValues");
		float[][] xvals2=(float[][])jutils.runPW4VoidMethod(iws[1],"getXValues");
		float[][] yvals2=(float[][])jutils.runPW4VoidMethod(iws[1],"getYValues");
		float[][] xvals3=(float[][])jutils.runPW4VoidMethod(iws[2],"getXValues");
		float[][] yvals3=(float[][])jutils.runPW4VoidMethod(iws[2],"getYValues");
		if(iws[0].getClass().getName().equals("jguis.PlotWindow3D")){
			int[] npts1=((int[][])jutils.runPW4VoidMethod(iws[0],"getNpts"))[0];
			int[] npts2=((int[][])jutils.runPW4VoidMethod(iws[1],"getNpts"))[0];
			int[] npts3=((int[][])jutils.runPW4VoidMethod(iws[2],"getNpts"))[0];
			float[][][] zvals1=(float[][][])jutils.runPW4VoidMethod(iws[0],"getZValues");
			float[][][] zvals2=(float[][][])jutils.runPW4VoidMethod(iws[1],"getZValues");
			float[][][] zvals3=(float[][][])jutils.runPW4VoidMethod(iws[2],"getZValues");
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
			if(npts3.length<=1){
				//need to convert to single points
				Object[] temp=separatePoints(xvals3,yvals3,zvals3,npts3,0);
				xvals3=(float[][])temp[0]; yvals3=(float[][])temp[1]; zvals3=(float[][][])temp[2]; npts3=(int[])temp[3];
			}
			float[] mindist=new float[npts1.length];
			int[] dest=new int[npts1.length];
			float[] mindist2=new float[npts1.length];
			int[] dest2=new int[npts1.length];
			StringBuffer sb=new StringBuffer();
			int npts=0;
			//search through the first list of points and find the closest point in the other lists and its index
			for(int i=0;i<npts1.length;i++){
				//search through the second list of points and find the closest one
				mindist[i]=dist(xvals1[i][0],xvals2[0][0],yvals1[i][0],yvals2[0][0],zvals1[0][i][0],zvals2[0][0][0]);
				for(int j=1;j<npts2.length;j++){
					float dist1=dist(xvals1[i][0],xvals2[j][0],yvals1[i][0],yvals2[j][0],zvals1[0][i][0],zvals2[0][j][0]);
					if(dist1<mindist[i]){
						mindist[i]=dist1;
						dest[i]=j;
					}
				}
				//search through the third list of points and find the closest one
				mindist2[i]=dist(xvals1[i][0],xvals3[0][0],yvals1[i][0],yvals3[0][0],zvals1[0][i][0],zvals3[0][0][0]);
				for(int j=1;j<npts3.length;j++){
					float dist1=dist(xvals1[i][0],xvals3[j][0],yvals1[i][0],yvals3[j][0],zvals1[0][i][0],zvals3[0][j][0]);
					if(dist1<mindist2[i]){
						mindist2[i]=dist1;
						dest2[i]=j;
					}
				}
				if(mindist[i]<=distthresh && mindist2[i]<=distthresh){
					float dist3=dist(xvals2[dest[i]][0],xvals3[dest2[i]][0],yvals2[dest[i]][0],yvals3[dest2[i]][0],zvals2[0][dest[i]][0],zvals3[0][dest2[i]][0]);
					float[] triplecoords={xvals1[i][0],yvals1[i][0],zvals1[0][i][0],xvals2[dest[i]][0],yvals2[dest[i]][0],zvals2[0][dest[i]][0],xvals3[dest2[i]][0],yvals3[dest2[i]][0],zvals3[0][dest2[i]][0]};
					sb.append(""+i+"\t"+mindist[i]+"\t"+dest[i]+"\t"+mindist2[i]+"\t"+dest2[i]+"\t"+dist3+"\t"+table_tools.print_float_array(triplecoords)+"\n");
					npts++;
				}
			}
			sb.deleteCharAt(sb.length()-1);
			new TextWindow("Min_Distances (under thresh)","Index\tDist1\tDest1_Index\tDist2\tDest2_Index\tDist3\tx1\ty1\tz1\tx2\ty2\tz2\tx3\ty3\tz3",sb.toString(),200,200);
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
			float[] mindist=new float[npts1.length];
			StringBuffer sb=new StringBuffer();
			int npts=0;
			for(int i=0;i<npts1.length;i++){
				mindist[i]=dist(xvals1[i][0],xvals2[0][0],yvals1[i][0],yvals2[0][0]);
				for(int j=1;j<npts2.length;j++){
					float dist1=dist(xvals1[i][0],xvals2[j][0],yvals1[i][0],yvals2[j][0]);
					if(dist1<mindist[i]) mindist[i]=dist1;
				}
				sb.append(mindist[i]+"\n");
				if(mindist[i]<=distthresh) npts++;
			}
			sb.deleteCharAt(sb.length()-1);
			new TextWindow("Min_Distances","Dist",sb.toString(),200,200);
			float[][] xvals=new float[npts][1];
			float[][] yvals=new float[npts][1];
			int[] newnpts=new int[npts];
			int counter=0;
			for(int i=0;i<npts1.length;i++){
				if(mindist[i]<=distthresh){
					xvals[counter][0]=xvals1[i][0]; yvals[counter][0]=yvals1[i][0]; newnpts[counter]=1;
					counter++;
				}
			}
			Plot4 plot=new Plot4("x","y",xvals,yvals,newnpts);
			int[] shapes=plot.getShapes(); for(int i=0;i<shapes.length;i++) shapes[i]=1;
			new PlotWindow4("Filtered Points",plot).draw();
		}
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
