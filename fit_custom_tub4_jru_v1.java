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

public class fit_custom_tub4_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	int xpts,ypts,zpts,profwidth,ngaus;
	float zratio;
	gausfunc gf;

	//this plugin fits elongated 3D spots with a custom strategy
	//first we create a 2D xz profile summing over the spots
	//then we fit each spot at its center z position to a lateral gaussian to get width
	//positions are initialized from a set of rois
	//latest version fits only the current frame and doesn't fit z if there is only one slice
	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		xpts=imp.getWidth(); ypts=imp.getHeight();
		ImageStack stack=imp.getStack();
		//zpts=stack.getSize();
		zratio=(float)jutils.get_zratio(imp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addNumericField("X_Stdev",0.95,5,15,null);
		gd.addNumericField("Z_Stdev",1.2,5,15,null);
		gd.addCheckbox("Calibrate?",true);
		gd.addCheckbox("Hide Realigned",false);
		gd.addNumericField("XZ_Width (pix)",10,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		zratio=(float)gd.getNextNumber();
		float startstdev=(float)gd.getNextNumber();
		float startzstdev=(float)gd.getNextNumber();
		boolean cal=gd.getNextBoolean();
		boolean hidereal=gd.getNextBoolean();
		int xzwidth=(int)gd.getNextNumber();
		float psize=(float)jutils.get_psize(imp);
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			IJ.error("need point rois");
		}
		Roi[] rois=rman.getRoisAsArray();
		int slices=imp.getNSlices();
		zpts=slices;
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int currframe=imp.getT();
		int currchan=imp.getC();
		Object[] stack2=jutils.get3DZSeries(stack,currchan-1,currframe-1,frames,slices,channels);
		//Object[] stack2=jutils.stack2array(stack);
		//duplicate the stack to get a float array
		float[][] fstack=algutils.get_region2(stack2,0,0,xpts,ypts,xpts,ypts);
		ngaus=rois.length; //this should always be 2

		//start by getting the xz profile
		int xc1=rois[0].getBounds().x;
		int yc1=rois[0].getBounds().y;
		int xc2=rois[1].getBounds().x;
		int yc2=rois[1].getBounds().y;
		//extend the profile the length of the spacing on either end
		int xinc=xc2-xc1;
		int x1=xc1-xinc;
		int yinc=yc2-yc1;
		int y1=yc1-yinc;
		int x2=xc2+xinc;
		int y2=yc2+yinc;
		float dist=(float)Math.sqrt(xinc*xinc+yinc*yinc);
		//IJ.log(""+dist);
		int newc1=(int)(dist+0.5f);
		int newc2=(int)(2.0f*dist+0.5f);
		float[] zprof1=new float[zpts];
		float[] zprof2=new float[zpts];
		float[] zreg1=new float[3*zpts];
		float[] zreg2=new float[3*zpts];
		//int xzwidth=10;
		float[] temp=profiler.get2DLineProfile(new float[]{x1,y1,x2,y2},fstack[0],1,0,xzwidth,xpts,ypts);
		for(int i=0;i<temp.length;i++) temp[i]*=(float)xzwidth;
		zprof1[0]=temp[newc1-1]; zprof1[0]+=temp[newc1]; zprof1[0]+=temp[newc1+1];
		zprof2[0]=temp[newc2-1]; zprof2[0]+=temp[newc2]; zprof2[0]+=temp[newc2+1];
		zreg1[0]=temp[newc1-1]; zreg1[1]=temp[newc1]; zreg1[2]=temp[newc1+1];
		zreg2[0]=temp[newc2-1]; zreg2[1]=temp[newc2]; zreg2[2]+=temp[newc2+1];
		profwidth=temp.length;
		float[] xzprof=new float[profwidth*zpts];
		System.arraycopy(temp,0,xzprof,0,profwidth);
		for(int i=1;i<zpts;i++){
			temp=profiler.get2DLineProfile(new float[]{x1,y1,x2,y2},fstack[i],1,0,xzwidth,xpts,ypts);
			for(int j=0;j<temp.length;j++) temp[j]*=(float)xzwidth;
			zprof1[i]=temp[newc1-1]; zprof1[i]+=temp[newc1]; zprof1[i]+=temp[newc1+1];
			zprof2[i]=temp[newc2-1]; zprof2[i]+=temp[newc2]; zprof2[i]+=temp[newc2+1];
			zreg1[3*i]=temp[newc1-1]; zreg1[3*i+1]=temp[newc1]; zreg1[3*i+2]=temp[newc1+1];
			zreg2[3*i]=temp[newc2-1]; zreg2[3*i+1]=temp[newc2]; zreg2[3*i+2]+=temp[newc2+1];
			System.arraycopy(temp,0,xzprof,i*profwidth,profwidth);
		}
		//show the profile
		new ImagePlus("XZ Profile",new FloatProcessor(profwidth,zpts,xzprof,null)).show();
		//now get the z maxima
		float zc1=zratio*maxPos(zprof1);
		float zc2=zratio*maxPos(zprof2);
		IJ.log(""+zc1+" , "+zc2);
		double[] params=new double[11];
		double[][] constraints=new double[2][11];
		params[0]=jstatistics.getstatistic("Min",xzprof,null); //baseline
		float max=jstatistics.getstatistic("Max",xzprof,null);
		constraints[0][0]=-max; constraints[1][0]=max;
		float amp1=jstatistics.getstatistic("Max",zreg1,null)-(float)params[0];
		params[1]=amp1;
		params[2]=newc1;
		params[3]=zc1;
		params[4]=startstdev; params[5]=startzstdev*zratio;
		constraints[0][1]=0.2*params[1]; constraints[1][1]=5.0*params[1];
		constraints[0][2]=params[2]-2; constraints[1][2]=params[2]+2;
		constraints[0][3]=params[3]-2*zratio; constraints[1][3]=params[3]+2*zratio;
		constraints[0][4]=0.25; constraints[1][4]=3.0*startstdev; 
		constraints[0][5]=0.25*zratio; constraints[1][5]=3.0*zratio*startzstdev; 
		float amp2=jstatistics.getstatistic("Max",zreg2,null)-(float)params[0];
		params[6]=amp2;
		params[7]=newc2;
		params[8]=zc2;
		params[9]=startstdev; params[10]=startzstdev*zratio;
		constraints[0][6]=0.2*params[6]; constraints[1][6]=5.0*params[6];
		constraints[0][7]=params[7]-2; constraints[1][7]=params[7]+2;
		constraints[0][8]=params[8]-2*zratio; constraints[1][8]=params[8]+2*zratio;
		constraints[0][9]=0.25; constraints[1][9]=3.0*startstdev; 
		constraints[0][10]=0.25*zratio; constraints[1][10]=3.0*zratio*startzstdev; 
		gf=new gausfunc();
		int[] fixes=new int[params.length];
		if(zpts==1){
			fixes[3]=1;
			fixes[5]=1;
			fixes[8]=1;
			fixes[10]=1;
		}
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0.0001,10,0.1);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,constraints,xzprof,null,stats,true);
		new ImagePlus("fit",new FloatProcessor(profwidth,zpts,fit,null)).show();
		//change the xc param back to cartesian coordinates
		//params are 0baseline,1amp1,2xc1,3zc1,4xystdev1,5zstdev1,6amp2,...
		float dx=(float)xinc/dist; float dy=(float)yinc/dist;
		double fxc1=(double)x1+(double)dx*params[2];
		double fyc1=(double)y1+(double)dy*params[2];
		double fxc2=(double)x1+(double)dx*params[7];
		double fyc2=(double)y1+(double)dy*params[7];

		//now calculate the profiles along the "y" axis
		//the perpendicular direction is -dy, dx
		double hl=6.0; //the half length of the profile
		double px11=fxc1-hl*dy;
		double py11=fyc1+hl*dx;
		double px12=fxc1+hl*dy;
		double py12=fyc1-hl*dx;
		double px21=fxc2-hl*dy;
		double py21=fyc2+hl*dx;
		double px22=fxc2+hl*dy;
		double py22=fyc2-hl*dx;

		//interpolate the plane for zc1
		float[] zcplane1=null;
		double tempz=params[3]/(double)zratio;
		if(tempz<=0.0) zcplane1=fstack[0];
		else if(tempz>=(double)(zpts-1)) zcplane1=fstack[zpts-1];
		else zcplane1=interpolation.interpz(fstack[(int)tempz],fstack[(int)tempz+1],xpts,ypts,(float)(tempz-(int)tempz));
		float[] p1=profiler.get2DLineProfile(new float[]{(float)px11,(float)py11,(float)px12,(float)py12},zcplane1,1,0,5,xpts,ypts);
		float[] zcplane2=null;
		tempz=params[8]/(double)zratio;
		if(tempz<=0.0) zcplane2=fstack[0];
		else if(tempz>=(double)(zpts-1)) zcplane2=fstack[zpts-1];
		else zcplane2=interpolation.interpz(fstack[(int)tempz],fstack[(int)tempz+1],xpts,ypts,(float)(tempz-(int)tempz));
		float[] p2=profiler.get2DLineProfile(new float[]{(float)px21,(float)py21,(float)px22,(float)py22},zcplane2,1,0,5,xpts,ypts);
		//now fit the y profiles
		ijpluginrunner ipr=new ijpluginrunner("fit_gaussian_jru_v1");
		float[] xvals1=new float[p1.length]; for(int i=0;i<p1.length;i++) xvals1[i]=(float)i;
		double[] params1=(double[])ipr.runPluginMethod("guessParams",new Object[]{xvals1,p1,p1.length});
		//IJ.log(""+params1[0]);
		constraints=(double[][])ipr.runPluginMethod("getConstraints",new Object[]{xvals1,params1});
		//IJ.log(""+constraints[0][0]);
		double[] stats1=new double[2];
		fixes=new int[4];
		float[] fit1=(float[])ipr.runPluginMethod("runFit",new Object[]{xvals1,p1,params1,stats1,constraints,fixes});
		//new PlotWindow4("Profile 1b","pixel","intensity",xvals1,p1).draw();
		new PlotWindow4("Profile 1","pixel","intensity",new float[][]{xvals1,xvals1},new float[][]{p1,fit1},null).draw();
		float[] xvals2=new float[p2.length]; for(int i=0;i<p2.length;i++) xvals2[i]=(float)i;
		double[] params2=(double[])ipr.runPluginMethod("guessParams",new Object[]{xvals2,p2,p2.length});
		constraints=(double[][])ipr.runPluginMethod("getConstraints",new Object[]{xvals2,params2});
		double[] stats2=new double[2];
		//params2 are baseline,xc,stdev,amp
		float[] fit2=(float[])ipr.runPluginMethod("runFit",new Object[]{xvals2,p2,params2,stats2,constraints,fixes});
		new PlotWindow4("Profile 2","pixel","intensity",new float[][]{xvals2,xvals2},new float[][]{p2,fit2},null).draw();

		//now output a realigned 2D profile along the spindle ("x") axis
		int realignwidth=15;
		fxc1=px11+params1[1]*dy;
		fyc1=py11-params1[1]*dx;
		fxc2=px21+params2[1]*dy;
		fyc2=py21-params2[1]*dx;
		float[] tempxvals={(float)(fxc1-(fxc2-fxc1)),(float)(fxc2+(fxc2-fxc1))};
		float[] tempyvals={(float)(fyc1-(fyc2-fyc1)),(float)(fyc2+(fyc2-fyc1))};
		float[] tempzvals={(float)(params[3]-(params[8]-params[3])),(float)(params[8]+(params[8]-params[3]))};
		int realignlength=profiler.get3DPolygonLength(tempxvals,tempyvals,tempzvals,false);
		float[] realigned=profiler.get3DStraightened(fstack,xpts,ypts,tempxvals,tempyvals,tempzvals,false,realignwidth,0,zratio);
		String rename=imp.getTitle();
		int dotpos=rename.lastIndexOf(".");
		if(dotpos>=0) rename=rename.substring(0,dotpos);
		rename=rename+"_realigned.tif";
		if(!hidereal) new ImagePlus(rename,new FloatProcessor(realignwidth,realignlength,realigned,null)).show();
		Traj3D t3D=new Traj3D("x","y","z",tempxvals,tempyvals,tempzvals);
		new PlotWindow3D("MultiGaus Plot",t3D).draw();

		double[] calparams=params.clone();
		double perpstdev1=params1[2];
		double perpstdev2=params2[2];
		//note that the z parameters are already in xy pixel units
		if(cal){
			calparams[2]*=(double)psize;
			calparams[3]*=(double)psize;
			calparams[4]*=(double)psize;
			calparams[5]*=(double)psize;
			calparams[7]*=(double)psize;
			calparams[8]*=(double)psize;
			calparams[9]*=(double)psize;
			calparams[10]*=(double)psize;
			fxc1*=(double)psize; fyc1*=(double)psize;
			fxc2*=(double)psize; fyc2*=(double)psize;
			perpstdev1*=(double)psize;
			perpstdev2*=(double)psize;
		}
		//we already integrated the ampliudes over the perpendicular axis, need to integrate over the parallel axis and the z axis
		double int1=params[1]*params[4]*params[5]*2.0*Math.PI/zratio;
		double int2=params[6]*params[9]*params[10]*2.0*Math.PI/zratio;

		//params are 0baseline,1amp1,2xc1,3zc1,4xystdev1,5zstdev1,6amp2,...
		String[] paramsnames={"base","amp1","xc1","yc1","zc1","parstdev1","zstdev1","perpstdev1","amp2","xc2","yc2","zc2","parstdev2","zstdev2","perpstdev2","int1","int2"};
		TextWindow outtable=jutils.selectTable("Tub4 Fits");
		if(outtable==null) outtable=FitDialog_v2.make_outtable("Tub4 Fits",paramsnames);
		StringBuffer sb=new StringBuffer();
		sb.append(imp.getTitle()+"\t"+(float)stats[1]+"\t"+(int)stats[0]+"\t"+calparams[0]+"\t");
		sb.append(""+(float)calparams[1]+"\t"+(float)fxc1+"\t"+(float)fyc1+"\t"+(float)calparams[3]+"\t"+(float)calparams[4]+"\t"+(float)calparams[5]+"\t"+(float)perpstdev1+"\t");
		sb.append(""+(float)calparams[6]+"\t"+(float)fxc2+"\t"+(float)fyc2+"\t"+(float)calparams[8]+"\t"+(float)calparams[9]+"\t"+(float)calparams[10]+"\t"+(float)perpstdev2+"\t"+(float)int1+"\t"+int2);
		outtable.append(sb.toString());
	}

	public float[][][] getZProfiles(double[] params,float[][] fstack,float[][] fit){
		//this is never used
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
		//params are baseline,amp1,xc1,zc1,xystdev1,zstdev1,amp2,...
		//5 per gaussian
		float[] func=new float[zpts*profwidth];
		double zc=params[3]/zratio;
		double stdevz=params[5]/zratio;
		gf.draw_2D_func(func,params[2],zc,profwidth,zpts,params[4],stdevz,(float)params[1]);
		zc=params[8]/zratio;
		stdevz=params[10]/zratio;
		gf.draw_2D_func(func,params[7],zc,profwidth,zpts,params[9],stdevz,(float)params[6]);
		//now convert to double and add the baseline
		double[] func2=new double[zpts*profwidth];
		for(int i=0;i<zpts*profwidth;i++){
			func2[i]=(double)func[i]+params[0];
		}
		return func2;
	}

	public void showresults(String results){
		IJ.log(results);
	}

}
