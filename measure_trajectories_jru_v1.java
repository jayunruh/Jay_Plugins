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
import jalgs.*;

public class measure_trajectories_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we measure trajectories on a stack
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Image","Trajectory"});
		if(imps==null) return;
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Statistic",jstatistics.stats,"Avg");
		gd.addNumericField("Measure_Radius",5.0,5,15,null);
		gd.addNumericField("Measure_Radius_Z (slices, optional)",5.0,5,15,null);
		gd.addNumericField("Z_Ratio",1.0,5,15,null);
		gd.addCheckbox("Measure_Circ?",false);
		gd.addNumericField("Circ_Radius",1.0,5,15,null);
		gd.addChoice("Circ_Stat",jstatistics.stats,"Min");
		gd.addCheckbox("Output_Objects?",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		float rad=(float)gd.getNextNumber();
		float zrad=(float)gd.getNextNumber();
		float zratio=(float)gd.getNextNumber();
		boolean circmeas=gd.getNextBoolean();
		float circrad=(float)gd.getNextNumber();
		String circstat=jstatistics.stats[gd.getNextChoiceIndex()];
		boolean outobj=gd.getNextBoolean();
		ImageWindow iw=imps[1].getWindow();
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		String[] starts=(String[])jutils.runPW4VoidMethod(iw,"getAnnotations");
		boolean threed=jutils.is3DPlot(iw);
		int[] npts=null;
		if(threed) npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
		else npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
		float[][] zvals=null;
		if(threed) zvals=((float[][][])jutils.runPW4VoidMethod(iw,"getZValues"))[0];

		ImageStack stack=imps[0].getStack();
		int width=imps[0].getWidth(); int height=imps[0].getHeight();
		int channels=imps[0].getNChannels();
		int currchan=imps[0].getC();
		int slices=imps[0].getNSlices();
		int frames=imps[0].getNFrames();
		if(threed && slices==1){slices=frames; frames=1;}
		//int frames=stack.getSize(); //assume single channel for now
		//IJ.log(""+npts.length);
		int totlength=0;
		for(int i=0;i<npts.length;i++) totlength+=npts[i];
		float[][] stats=new float[totlength][5+channels];
		if(threed) stats=new float[totlength][6+channels];
		if(circmeas){
			int temp=stats[0].length;
			stats=new float[totlength][temp+channels];
		}
		float[][] objstack=null;
		if(outobj){
			if(threed) objstack=new float[frames*slices][width*height];
			else objstack=new float[frames][width*height];
		}
		int counter=0;
		float circzrad=(zrad/rad)*(rad+circrad);
		//stats are id,start,frame,x,y,(z),stat
		for(int i=0;i<npts.length;i++){
			int start=0;
			if(starts!=null){
				start=(int)Float.parseFloat(starts[i]);
			}
			for(int j=start;j<(start+npts[i]);j++){
				int tpos=j;
				if(frames==1) tpos=0;
				float[][] circ=new float[channels][];
				float[][] circ2=new float[channels][];
				if(threed){
					for(int k=0;k<channels;k++){
						Object[] frame=jutils.get3DZSeries(stack,k,tpos,frames,slices,channels);
						circ[k]=getSphereVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],zvals[i][j-start]/zratio,rad,zrad);
						if(circmeas) circ2[k]=getSphereVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],zvals[i][j-start]/zratio,rad,rad+circrad,circzrad);
					}
					if(outobj){
						Object[] tempobjstack=algutils.get3DZSeries(objstack,tpos,0,frames,slices,1);
						setSphereObj(tempobjstack,(float)(i+1),width,height,xvals[i][j-start],yvals[i][j-start],zvals[i][j-start]/zratio,rad,zrad);
					}
				} else { 
					for(int k=0;k<channels;k++){
						Object frame=jutils.get3DSlice(stack,tpos,0,k,frames,slices,channels);
						circ[k]=getCircleVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],rad);
						if(circmeas) circ2[k]=getCircleVals(frame,width,height,xvals[i][j-start],yvals[i][j-start],rad,rad+circrad);
					}
					if(outobj) setCircleObj((float[])objstack[tpos],(float)(i+1),width,height,xvals[i][j-start],yvals[i][j-start],rad);
				}
				stats[counter][0]=i;
				stats[counter][1]=start;
				stats[counter][2]=j-start;
				stats[counter][3]=xvals[i][j-start];
				stats[counter][4]=yvals[i][j-start];
				if(threed){
					stats[counter][5]=zvals[i][j-start];
					for(int k=0;k<channels;k++){
						stats[counter][6+k]=jstatistics.getstatistic(stat,circ[k],null);
					}
				} else {
					for(int k=0;k<channels;k++){
						stats[counter][5+k]=jstatistics.getstatistic(stat,circ[k],null);
					}
				}
				if(circmeas){
					int off=stats[0].length-channels;
					for(int k=0;k<channels;k++){
						stats[counter][off+k]=jstatistics.getstatistic(circstat,circ2[k],null);
					}
				}
				counter++;
			}
		}
		//now output the stats
		String[] statlabels=new String[stats[0].length];
		statlabels[0]="id"; statlabels[1]="start"; statlabels[2]="frame"; statlabels[3]="x"; statlabels[4]="y";
		int off=5;
		if(threed){statlabels[off]="z"; off++;}
		for(int i=0;i<channels;i++){statlabels[off]=stat+(i+1); off++;}
		if(circmeas){
			for(int i=0;i<channels;i++){statlabels[off]="circ_"+circstat+(i+1); off++;}
		}
		table_tools.create_table("Traj Stats",stats,statlabels);
		//and the object image if required
		if(outobj){
			new ImagePlus("Objects",jutils.array2stack(objstack,width,height)).show();
		}
	}

	public void setCircleObj(float[] image,float objval,int width,int height,float xc,float yc,float rad){
		setCircleObj(image,objval,width,height,xc,yc,0.0f,rad);
	}

	public void setCircleObj(float[] image,float objval,int width,int height,float xc,float yc,float rad1,float rad2){
		int size=(int)(2.0f*rad2+1.0f);
		float[] square=new float[size*size];
		int startx=(int)(xc-rad2);
		int starty=(int)(yc-rad2);
		for(int i=starty;i<(starty+size);i++){
			float y=yc-(float)i;
			for(int j=startx;j<(startx+size);j++){
				float x=xc-(float)j;
				float dist2=(float)Math.sqrt(x*x+y*y);
				if(dist2<=rad2 && dist2>=rad1 && i>=0 && i<height && j>=0 && j<width){
					if(image[j+i*width]<=0.0f) image[j+i*width]=objval;
				}
			}
		}
	}

	public float[] getCircleVals(Object image,int width,int height,float xc,float yc,float rad){
		return getCircleVals(image,width,height,xc,yc,0.0f,rad);	
	}

	public float[] getCircleVals(Object image,int width,int height,float xc,float yc,float rad1,float rad2){
		int size=(int)(2.0f*rad2+1.0f);
		float[] square=new float[size*size];
		int startx=(int)(xc-rad2);
		int starty=(int)(yc-rad2);
		int counter=0;
		for(int i=starty;i<(starty+size);i++){
			float y=yc-(float)i;
			for(int j=startx;j<(startx+size);j++){
				float x=xc-(float)j;
				float dist2=(float)Math.sqrt(x*x+y*y);
				if(dist2<=rad2 && dist2>=rad1 && i>=0 && i<height && j>=0 && j<width){
					int pos=j+i*width;
					Object temp=image;
					if(temp instanceof float[]) square[counter]=((float[])temp)[pos];
					if(temp instanceof short[]) square[counter]=(float)(((short[])temp)[pos]&0xffff);
					if(temp instanceof byte[]) square[counter]=(float)(((byte[])temp)[pos]&0xff);
					counter++;
				}
			}
		}
		return (float[])algutils.get_subarray(square,0,counter);		
	}

	public void setSphereObj(Object[] stack,float objval,int width,int height,float xc,float yc,float zc,float rad,float zrad){
		setSphereObj(stack,objval,width,height,xc,yc,zc,0.0f,rad,zrad);
	}

	public void setSphereObj(Object[] stack,float objval,int width,int height,float xc,float yc,float zc,float rad1,float rad2,float zrad2){
		//this puts the object number in the sphere if there is not already something there
		if(rad2==0.0f && zrad2==0.0f){ //0 radius, place a point
			int zpos=(int)(zc+0.5f); if(zpos<0) zpos=0; if(zpos>=stack.length) zpos=stack.length-1;
			int xpos=(int)(xc+0.5f); if(xpos<0) xpos=0; if(xpos>=width) zpos=width-1;
			int ypos=(int)(yc+0.5f); if(ypos<0) ypos=0; if(ypos>=height) ypos=height-1;
			if(((float[])stack[zpos])[xpos+ypos*width]<=0.0f) ((float[])stack[zpos])[xpos+ypos*width]=objval;
			return;
		}
		int size=(int)(2.0f*rad2+1.0f);
		int zsize=(int)(2.0f*zrad2+1.0f);
		int slices=stack.length;
		float[] square=new float[size*size*zsize];
		int startx=(int)(xc-rad2);
		int starty=(int)(yc-rad2);
		int startz=(int)(zc-zrad2);
		float zratio=zrad2/rad2;
		for(int i=starty;i<(starty+size);i++){
			float y=yc-(float)i;
			for(int j=startx;j<(startx+size);j++){
				float x=xc-(float)j;
				for(int k=startz;k<(startz+zsize);k++){
					float z=zc-(float)k;
					float dist2=(float)Math.sqrt(x*x+y*y+z*z/(zratio*zratio));
					if(dist2<=rad2 && dist2>=rad1 && i>=0 && i<height && j>=0 && j<width && k>=0 && k<slices){
						if(((float[])stack[k])[j+i*width]<=0.0f) ((float[])stack[k])[j+i*width]=objval;
					}
				}
			}
		}
	}

	//gets the intensity values in a sphere of xy radius rad and z radius zrad
	public float[] getSphereVals(Object[] image,int width,int height,float xc,float yc,float zc,float rad,float zrad){
		return getSphereVals(image,width,height,xc,yc,zc,0.0f,rad,zrad);	
	}

	//gets the intensity values in a sphere or shell between rad1 and rad2
	public float[] getSphereVals(Object[] image,int width,int height,float xc,float yc,float zc,float rad1,float rad2,float zrad2){
		if(rad2==0.0f && zrad2==0.0f){ //0 radius, return a value
			float temp=interpolation.interp3D(image,width,height,xc,yc,zc);
			return new float[]{temp};
		}
		int size=(int)(2.0f*rad2+1.0f);
		int zsize=(int)(2.0f*zrad2+1.0f);
		int slices=image.length;
		float[] square=new float[size*size*zsize];
		int startx=(int)(xc-rad2);
		int starty=(int)(yc-rad2);
		int startz=(int)(zc-zrad2);
		float zratio=zrad2/rad2;
		int counter=0;
		for(int i=starty;i<(starty+size);i++){
			float y=yc-(float)i;
			for(int j=startx;j<(startx+size);j++){
				float x=xc-(float)j;
				for(int k=startz;k<(startz+zsize);k++){
					float z=zc-(float)k;
					float dist2=(float)Math.sqrt(x*x+y*y+z*z/(zratio*zratio));
					if(dist2<=rad2 && dist2>=rad1 && i>=0 && i<height && j>=0 && j<width && k>=0 && k<slices){
						Object temp=image[k];
						int pos=j+i*width;
						if(temp instanceof float[]) square[counter]=((float[])temp)[pos];
						if(temp instanceof short[]) square[counter]=(float)(((short[])temp)[pos]&0xffff);
						if(temp instanceof byte[]) square[counter]=(float)(((byte[])temp)[pos]&0xff);
						counter++;
					}
				}
			}
		}
		return (float[])algutils.get_subarray(square,0,counter);		
	}

	public PolygonRoi traj2roi(float[] xvals,float[] yvals,int npts){
		int[] xvals2=new int[npts];
		int[] yvals2=new int[npts];
		for(int i=0;i<npts;i++){
			xvals2[i]=(int)xvals[i];
			yvals2[i]=(int)yvals[i];
		}
		return new PolygonRoi(xvals2,yvals2,npts,Roi.POLYLINE);
	}

	public int[][] traj2int(float[] xvals,float[] yvals,int npts){
		int[] xvals2=new int[npts];
		int[] yvals2=new int[npts];
		for(int i=0;i<npts;i++){
			xvals2[i]=(int)xvals[i];
			yvals2[i]=(int)yvals[i];
		}
		return new int[][]{xvals2,yvals2};
	}

}
