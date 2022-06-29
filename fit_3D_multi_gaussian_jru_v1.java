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
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;
import ij.plugin.frame.*;
import ij.text.*;

public class fit_3D_multi_gaussian_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	int xpts,ypts,zpts,ngaus;
	float zratio;
	gausfunc gf;

	//this plugin fits a 3D stack to multiple 3D gaussians
	//a 3D profile plot is provided for interactive visualization
	//positions are initialized from a set of rois
	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		xpts=imp.getWidth(); ypts=imp.getHeight();
		ImageStack stack=imp.getStack();
		zpts=stack.getSize();
		zratio=(float)jutils.get_zratio(imp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addNumericField("XY_Stdev",0.95,5,15,null);
		gd.addNumericField("Z_Stdev",1.2,5,15,null);
		gd.addCheckbox("Calibrate?",true);
		gd.addCheckbox("Show 3D Plot",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		zratio=(float)gd.getNextNumber();
		float startstdev=(float)gd.getNextNumber();
		float startzstdev=(float)gd.getNextNumber();
		boolean cal=gd.getNextBoolean();
		boolean showplot=gd.getNextBoolean();
		float psize=(float)jutils.get_psize(imp);
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			IJ.error("need point rois");
		}
		Roi[] rois=rman.getRoisAsArray();
		//Object[] stack2=jutils.stack2array(stack);
		int currt=imp.getT()-1; int currc=imp.getC()-1;
		int nframes=imp.getNFrames(); int nchans=imp.getNChannels(); zpts=imp.getNSlices();
		if(zpts==1){zpts=nframes; nframes=1;}
		Object[] stack2=jutils.get3DZSeries(stack,currc,currt,nframes,zpts,nchans);
		float[][] fstack=algutils.get_region2(stack2,0,0,xpts,ypts,xpts,ypts);
		ngaus=rois.length;
		//IJ.log(""+ngaus)
		double[] params=new double[ngaus*6+1];
		double[][] constraints=new double[2][ngaus*6+1];
		params[0]=get3DMinMax(fstack)[0]; //this is the baseline
		float max=get3DMinMax(fstack)[1];
		constraints[0][0]=-max; constraints[1][0]=max;
		for(int i=0;i<rois.length;i++){
			Rectangle r=rois[i].getBounds();
			int xc=r.x; int yc=r.y;
			float[] zprofile=getAvgProfile(fstack,xc,yc);
			//float zc=zratio*centerOfMass(zprofile);
			float zc=zratio*maxPos(zprofile);
			float[][] tempstack=algutils.get_region(fstack,xc,yc,3,3,xpts,ypts);
			float amp=get3DMinMax(tempstack)[1]-(float)params[0];
			params[i*6+1]=amp;
			params[i*6+2]=xc; params[i*6+3]=yc; params[i*6+4]=zc;
			params[i*6+5]=startstdev; params[i*6+6]=startzstdev*zratio;
			constraints[0][i*6+1]=0.2*params[i*6+1]; constraints[1][i*6+1]=5.0*params[i*6+1]; 
			constraints[0][i*6+2]=params[i*6+2]-2; constraints[1][i*6+2]=params[i*6+2]+2;
			constraints[0][i*6+3]=params[i*6+3]-2; constraints[1][i*6+3]=params[i*6+3]+2;
			constraints[0][i*6+4]=params[i*6+4]-2*zratio; constraints[1][i*6+4]=params[i*6+4]+2*zratio;
			constraints[0][i*6+5]=0.25; constraints[1][i*6+5]=3.0*startstdev; 
			constraints[0][i*6+6]=0.25*zratio; constraints[1][i*6+6]=3.0*zratio*startzstdev; 
		}
		double[] origparams=params.clone();
		String[] paramsnames=new String[params.length];
		paramsnames[0]="baseline";
		for(int i=0;i<ngaus;i++){
			paramsnames[6*i+1]="Amp"+(i+1);
			paramsnames[6*i+2]="xc"+(i+1);
			paramsnames[6*i+3]="yc"+(i+1);
			paramsnames[6*i+4]="zc"+(i+1);
			paramsnames[6*i+5]="xystdev"+(i+1);
			paramsnames[6*i+6]="zstdev"+(i+1);
		}
		//for(int i=0;i<params.length;i++) IJ.log(""+params[i]);
		gf=new gausfunc();
		int[] fixes=new int[params.length];
		if(zpts==1){
			for(int i=0;i<ngaus;i++){
				fixes[6*i+4]=1;
				fixes[6*i+6]=1;
			}
		}
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0.0001,10,0.1);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,constraints,convertTo1D(fstack),null,stats,true);
		float[][] fitstack=convertTo2D(fit,xpts*ypts,zpts);
		
		ImageStack fitstack2=jutils.array2stack(fitstack,xpts,ypts);
		new ImagePlus("fit",fitstack2).show();
		float[][][] zprofiles=getZProfiles(origparams,fstack,fitstack);
		new PlotWindow4("Z Profiles","z","Intensity",zprofiles[0],zprofiles[1],null).draw();
		float[] maxstack=algutils.get_stack_proj_stat("Avg",fstack,xpts,ypts,zpts,null);
		float[][] maxstack2=convertTo2D(maxstack,xpts,ypts);
		float[] maxfit=algutils.get_stack_proj_stat("Avg",fitstack,xpts,ypts,zpts,null);
		float[][] maxfit2=convertTo2D(maxfit,xpts,ypts);
		new PlotWindow3D("Avg Projections","x","y","Max",new float[][][]{maxstack2,maxfit2},0,null).draw();

		double[] calparams=params.clone();
		for(int i=0;i<ngaus;i++){
			calparams[6*i+2]*=(double)psize;
			calparams[6*i+3]*=(double)psize;
			calparams[6*i+4]*=(double)psize;
			calparams[6*i+5]*=(double)psize;
			calparams[6*i+6]*=(double)psize;
		}

		StringBuffer sb=new StringBuffer();
		sb.append(imp.getTitle());
		IJ.log("Chi Squared = "+(float)stats[1]);
		sb.append("\t"+(float)stats[1]);
		IJ.log("Iterations = "+(int)stats[0]);
		sb.append("\t"+(int)stats[0]);
		for(int i=0;i<params.length;i++){
			IJ.log(paramsnames[i]+" : "+params[i]+" : Fixed : "+(fixes[i]==1));
			if(cal) sb.append("\t"+(float)calparams[i]);
			else sb.append("\t"+(float)params[i]);
		}
		TextWindow outtable=jutils.selectTable("MultiGaus Fits");
		if(outtable==null){
			outtable=FitDialog_v2.make_outtable("MultiGaus Fits",paramsnames);
		} else {
			FitDialog_v2.adapt_outtable(outtable,paramsnames);
		}
		outtable.append(sb.toString());
		if(showplot){
			float[] xvals=new float[ngaus];
			float[] yvals=new float[ngaus];
			float[] zvals=new float[ngaus];
			for(int i=0;i<ngaus;i++){
				xvals[i]=(float)calparams[6*i+2];
				yvals[i]=(float)calparams[6*i+3];
				zvals[i]=(float)calparams[6*i+4];
			}
			Traj3D t3D=new Traj3D("x","y","z",xvals,yvals,zvals);
			new PlotWindow3D("MultiGaus Plot",t3D).draw();
		}
	}

	public float[][][] getZProfiles(double[] params,float[][] fstack,float[][] fit){
		float[][] yvals=new float[2*ngaus][];
		float[][] xvals=new float[2*ngaus][];
		for(int i=0;i<ngaus;i++){
			int xc=(int)params[6*i+2]; int yc=(int)params[6*i+3];
			yvals[2*i]=getAvgProfile(fstack,xc,yc);
			yvals[2*i+1]=getAvgProfile(fit,xc,yc);
			xvals[2*i]=getXVals(zratio,zpts);
			xvals[2*i+1]=getXVals(zratio,zpts);
		}
		return new float[][][]{xvals,yvals};
	}

	public float[] getXVals(float xinc,int npts){
		float[] temp=new float[npts];
		for(int i=0;i<npts;i++) temp[i]=xinc*(float)i;
		return temp;
	}

	public float[] getAvgProfile(float[][] fstack,int xc,int yc){
		float[] profile=new float[zpts];
		float[][] tempstack=algutils.get_region(fstack,xc,yc,3,3,xpts,ypts);
		for(int j=0;j<zpts;j++) profile[j]=jstatistics.getstatistic("Avg",tempstack[j],null);
		return profile;
	}

	public float[] convertTo1D(float[][] data){
		float[] temp=new float[data.length*data[0].length];
		for(int i=0;i<data.length;i++){
			System.arraycopy(data[i],0,temp,i*data[0].length,data[0].length);
		}
		return temp;
	}

	public float[][] convertTo2D(float[] data,int sublength,int slices){
		float[][] temp=new float[slices][sublength];
		for(int i=0;i<slices;i++){
			System.arraycopy(data,i*sublength,temp[i],0,sublength);
		}
		return temp;
	}

	public float centerOfMass(float[] dist){
		float sum=0.0f;
		float xsum=0.0f;
		for(int i=0;i<dist.length;i++){
			sum+=dist[i];
			xsum+=dist[i]*(float)i;
		}
		return xsum/sum;
	}

	public float maxPos(float[] dist){
		float max=dist[0];
		int maxpos=0;
		for(int i=0;i<dist.length;i++){
			if(dist[i]>max){
				max=dist[i];
				maxpos=i;
			}
		}
		return (float)maxpos;
	}	

	public float[] get3DMinMax(float[][] box3D){
		float min=jstatistics.getstatistic("Min",box3D[0],null);
		float max=jstatistics.getstatistic("Max",box3D[0],null);
		for(int i=1;i<box3D.length;i++){
			float tempmin=jstatistics.getstatistic("Min",box3D[i],null);
			if(tempmin<min) min=tempmin;
			float tempmax=jstatistics.getstatistic("Max",box3D[i],null);
			if(tempmax>max) max=tempmax;
		}
		return new float[]{min,max};
	}

	public double[] fitfunc(double[] params){
		//params are baseline,amp1,xc1,yc1,zc1,xystdev1,zstdev1,amp2,...
		//6 per gaussian
		float[][] func=new float[zpts][xpts*ypts];
		for(int i=0;i<ngaus;i++){
			double zc=params[i*6+4]/zratio;
			double stdevz=params[i*6+6]/zratio;
			gf.draw_3D_func(func,params[i*6+2],params[i*6+3],zc,xpts,ypts,params[i*6+5],stdevz,(float)params[i*6+1]);
		}
		//now convert to double and add the baseline
		double[] func2=new double[zpts*xpts*ypts];
		int imgsize=xpts*ypts;
		for(int i=0;i<zpts;i++){
			for(int j=0;j<imgsize;j++) func2[j+i*imgsize]=(double)func[i][j]+params[0];
		}
		return func2;
	}

	public void showresults(String results){
		IJ.log(results);
	}

}
