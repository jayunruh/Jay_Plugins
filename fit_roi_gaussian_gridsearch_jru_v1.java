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
import jalgs.*;
import jalgs.jfit.*;
import ij.text.*;
import jguis.*;

public class fit_roi_gaussian_gridsearch_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	int xpts,ypts;
	float[] data;
	gausfunc gf;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int nchannels=imp.getNChannels();
		Rectangle r=imp.getRoi().getBounds();
		Color tempc=imp.getRoi().getStrokeColor();
		String rcolor="yellow";
		if(tempc!=null) rcolor=jutils.get_closest_color_name(tempc);
		float[][] charrays=new float[2][];
		xpts=r.width; ypts=r.height;
		int centerx=r.x+r.width/2;
		int centery=r.y+r.height/2;
		if(r.width<2 && r.height<2){
			//assume we have a point roi
			xpts=10; r.width=10;
			ypts=10; r.height=10;
			r.x=centerx-xpts/2;
			r.y=centery-ypts/2;
		}
		boolean findmax=true;
		GenericDialog gd2=new GenericDialog("Options");
		gd2.addCheckbox("Find_Z_Max",findmax);
		gd2.addNumericField("X_Size (pixels)",xpts,0);
		gd2.addNumericField("Y_Size (pixels)",ypts,0);
		gd2.addCheckbox("Single_Channel",true);
		gd2.addCheckbox("Output_Roi_Color",false);
		gd2.showDialog(); if(gd2.wasCanceled()){return;}
		findmax=gd2.getNextBoolean();
		xpts=(int)gd2.getNextNumber(); r.width=xpts;
		ypts=(int)gd2.getNextNumber(); r.height=ypts;
		boolean singchan=gd2.getNextBoolean();
		boolean outcolor=gd2.getNextBoolean();
		r.x=centerx-xpts/2;
		r.y=centery-ypts/2;
		int currz=imp.getSlice();
		int currt=imp.getFrame();
		int currc=imp.getChannel();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int maxslice=0;
		if(nchannels>1 && !singchan){
			Object[] temp=jutils.get3DCSeries(imp.getStack(),currz-1,currt-1,frames,slices,nchannels);
			int ch1=0; int ch2=1;
			if(nchannels>2){
				float[] avgs={jstatistics.getstatistic("Avg",temp[0],null),jstatistics.getstatistic("Avg",temp[1],null),jstatistics.getstatistic("Avg",temp[2],null)};
				//the max intensity channel is always the transmitted light channel
				if(avgs[1]>avgs[0] && avgs[1]>avgs[2]){ch2=2;}
				if(avgs[0]>avgs[1] && avgs[0]>avgs[2]){ch1=1; ch2=2;}
			}
			if(findmax){
				Object[] temp1=jutils.get3DZSeries(imp.getStack(),ch1,currt-1,frames,slices,nchannels);
				Object[] temp2=jutils.get3DZSeries(imp.getStack(),ch2,currt-1,frames,slices,nchannels);
				float maxint=jstatistics.getstatistic("Avg",temp1[0],width,height,r,null);
				maxint+=jstatistics.getstatistic("Avg",temp2[0],width,height,r,null);
				maxslice=0;
				for(int i=1;i<slices;i++){
					float tempint=jstatistics.getstatistic("Avg",temp1[i],width,height,r,null);
					tempint+=jstatistics.getstatistic("Avg",temp2[i],width,height,r,null);
					if(tempint>maxint){
						maxint=tempint;
						maxslice=i;
					}
				}
				temp[ch1]=temp1[maxslice];
				temp[ch2]=temp2[maxslice];
			}
			charrays[0]=algutils.convert_arr_float(algutils.get_region(temp[ch1],centerx,centery,xpts,ypts,width,height));
			charrays[1]=algutils.convert_arr_float(algutils.get_region(temp[ch2],centerx,centery,xpts,ypts,width,height));
			data=new float[xpts*ypts];
			for(int i=0;i<data.length;i++) data[i]=charrays[0][i]+charrays[1][i];
		} else {
			Object temp=imp.getProcessor().getPixels();
			if(findmax){
				Object[] temp1=jutils.get3DZSeries(imp.getStack(),currc-1,currt-1,frames,slices,nchannels);
				float maxint=jstatistics.getstatistic("Avg",temp1[0],width,height,r,null);
				maxslice=0;
				for(int i=1;i<slices;i++){
					float tempint=jstatistics.getstatistic("Avg",temp1[i],width,height,r,null);
					if(tempint>maxint){
						maxint=tempint;
						maxslice=i;
					}
				}
				temp=temp1[maxslice];
			}
			data=algutils.convert_arr_float(algutils.get_region(temp,r.x+r.width/2,r.y+r.height/2,r.width,r.height,width,height));
		}
		if(maxslice==0) maxslice=currz;
		float[] xymax=getxymax(data,xpts,ypts);
		int xmax=(int)xymax[1]; int ymax=(int)xymax[2]; float max=xymax[0];
		//ymax=ypts-ymax;
		GenericDialog gd=new GenericDialog("Options");
		double xcmin=(double)xmax-2.0;
		gd.addNumericField("xc_min",xcmin,5,10,null);
		double xcmax=(double)xmax+2.0;
		gd.addNumericField("xc_max",xcmax,5,10,null);
		double ycmin=(double)ymax-2.0;
		gd.addNumericField("yc_min",ycmin,5,10,null);
		double ycmax=(double)ymax+2.0;
		gd.addNumericField("yc_max",ycmax,5,10,null);
		double stdevmin=0.5;
		gd.addNumericField("stdev_min",stdevmin,5,10,null);
		double stdevmax=(double)xpts*0.25;
		gd.addNumericField("stdev_max",stdevmax,5,10,null);
		gd.addCheckbox("Table_Output",false);
		gd.addCheckbox("Calibrate",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		xcmin=gd.getNextNumber();
		xcmax=gd.getNextNumber();
		ycmin=gd.getNextNumber();
		ycmax=gd.getNextNumber();
		stdevmin=gd.getNextNumber();
		stdevmax=gd.getNextNumber();
		boolean tableout=gd.getNextBoolean();
		boolean calibrate=gd.getNextBoolean();
		gridfit fitclass=new gridfit(this);
		fitclass.output=false;
		double[] params=new double[3];
		double[][] constraints=new double[3][3];
		constraints[0][0]=stdevmin; constraints[1][0]=stdevmax; constraints[2][0]=0.1;
		constraints[0][1]=xcmin; constraints[1][1]=xcmax; constraints[2][1]=0.25;
		constraints[0][2]=ycmin; constraints[1][2]=ycmax; constraints[2][2]=0.25;
		double[] stats=new double[2];
		int[] fixes={0,0,0};
		gf=new gausfunc();
		float[] fit=fitclass.fitdata(params,fixes,constraints,data,null,stats);
		double[] ampoffset=(new linleastsquares()).get_amp_offset(get_gaussian(params),data,true);
		double[][] champoffset=new double[2][];
		if(nchannels>2 && !singchan){
			champoffset[0]=(new linleastsquares()).get_amp_offset(get_gaussian(params),charrays[0],true);
			champoffset[1]=(new linleastsquares()).get_amp_offset(get_gaussian(params),charrays[1],true);
		}
		float xoff=(float)r.x;
		float yoff=(float)r.y;
		if(calibrate){
			double psize=jutils.get_psize(imp);
			params[0]*=psize;
			params[1]*=psize;
			params[2]*=psize;
			xoff*=(float)psize;
			yoff*=(float)psize;
		}
		if(tableout){
			TextWindow tw=jutils.selectTable("Gaussian_Output");
			if(nchannels>1 && !singchan){
				if(outcolor){
					if(tw==null){
						tw=new TextWindow("Gaussian_Output","baseline\tamplitude\tstdev\txc\tyc\tzslice\tb1\ta1\tb2\ta2\troicolor","",400,200);
					}
					tw.append(""+(float)ampoffset[1]+"\t"+(float)ampoffset[0]+"\t"+(float)params[0]+"\t"+((float)params[1]+xoff)+"\t"+((float)params[2]+yoff)+"\t"+(float)(maxslice)+"\t"+(float)champoffset[0][1]+"\t"+(float)champoffset[0][0]+"\t"+(float)champoffset[1][1]+"\t"+(float)champoffset[1][0]+"\t"+rcolor+"\n");
				} else {
					if(tw==null){
						tw=new TextWindow("Gaussian_Output","baseline\tamplitude\tstdev\txc\tyc\tzslice\tb1\ta1\tb2\ta2","",400,200);
					}
					tw.append(""+(float)ampoffset[1]+"\t"+(float)ampoffset[0]+"\t"+(float)params[0]+"\t"+((float)params[1]+xoff)+"\t"+((float)params[2]+yoff)+"\t"+(float)(maxslice)+"\t"+(float)champoffset[0][1]+"\t"+(float)champoffset[0][0]+"\t"+(float)champoffset[1][1]+"\t"+(float)champoffset[1][0]+"\n");
				}
			} else {
				if(outcolor){
					if(tw==null){
						tw=new TextWindow("Gaussian_Output","baseline\tamplitude\tstdev\txc\tyc\tzslice\troicolor","",400,200);
					}
					tw.append(""+(float)ampoffset[1]+"\t"+(float)ampoffset[0]+"\t"+(float)params[0]+"\t"+((float)params[1]+xoff)+"\t"+((float)params[2]+yoff)+"\t"+(float)(maxslice)+"\t"+rcolor+"\n");
				} else {
					if(tw==null){
						tw=new TextWindow("Gaussian_Output","baseline\tamplitude\tstdev\txc\tyc\tzslice\tc2","",400,200);
					}
					tw.append(""+(float)ampoffset[1]+"\t"+(float)ampoffset[0]+"\t"+(float)params[0]+"\t"+((float)params[1]+xoff)+"\t"+((float)params[2]+yoff)+"\t"+(float)(maxslice)+"\t"+(float)stats[1]+"\n");
				}
			}
		} else {
			IJ.log("baseline = "+(float)ampoffset[1]);
			IJ.log("amplitude = "+(float)ampoffset[0]);
			IJ.log("stdev = "+(float)params[0]);
			IJ.log("x center = "+(float)params[1]);
			IJ.log("y center = "+(float)params[2]);
			if(nchannels>1 && !singchan){
				IJ.log("baseline1 = "+(float)champoffset[0][1]);
				IJ.log("amp1 = "+(float)champoffset[0][0]);
				IJ.log("baseline2 = "+(float)champoffset[1][1]);
				IJ.log("amp2 = "+(float)champoffset[1][0]);
			}
		}
		//IJ.log(""+r.width+" , "+r.height+" , "+data.length);
		new ImagePlus("fit",new FloatProcessor(r.width,r.height,fit,null)).show();		
	}

	public float[] getxymax(float[] image,int width,int height){
		float max=image[0];
		int xmax=0;
		int ymax=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if(image[j+i*width]>max){
					max=image[j+i*width];
					xmax=j;
					ymax=i;
				}
			}
		}
		return new float[]{max,(float)xmax,(float)ymax};
	}

	public double[] fitfunc(double[] params)
	{
		//the params list is stdev,xc,yc
		double[] func=get_gaussian(params);
		double[] ampoffset=(new linleastsquares()).get_amp_offset(func,data,true);
		for(int i=0;i<(xpts*ypts);i++){func[i]*=ampoffset[0]; func[i]+=ampoffset[1];}
		return func;
	}

	public double[] get_gaussian(double[] params){
		float[] func= gf.get_2D_func2(-params[1],xpts,1.0,-params[2],ypts,1.0,params[0]);
		return algutils.convert_arr_double(func);
		/*double[] func=new double[xpts*ypts];
		for(int i=0;i<ypts;i++){
			double ydiff=((double)i-params[2]);
			for(int j=0;j<xpts;j++){
				double xdiff=((double)j-params[1]);
				func[j+i*xpts]= Math.exp(((-0.5)*(xdiff*xdiff))/(params[0]*params[0]));
				func[j+i*xpts]*= Math.exp(((-0.5)*(ydiff*ydiff))/(params[0]*params[0]));
			}
		}
		return func;*/
	}

	public void showresults(String results){
		IJ.log(results);
	}


}
