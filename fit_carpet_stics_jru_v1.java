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
import ij.plugin.frame.*;
import jalgs.jfit.*;
import ij.text.*;

public class fit_carpet_stics_jru_v1 implements PlugIn {
	gausfunc gf;
	linleastsquares lls;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		float[] pixels=(float[])imp.getProcessor().getPixels();
		GenericDialog gd=new GenericDialog("Options");
		boolean skipfirst=true;
		gd.addCheckbox("Skip_First_G(0)",skipfirst);
		boolean skiprest=true;
		gd.addCheckbox("Skip_Other_G(0)s",skiprest);
		int npts=width/2;
		gd.addNumericField("Pts_to_fit_in_x",npts,0);
		boolean fixstdev=false;
		gd.addCheckbox("Fix_Stdev",fixstdev);
		float stdev=4.0f;
		gd.addNumericField("Stdev(pixels)",stdev,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		skipfirst=gd.getNextBoolean();
		skiprest=gd.getNextBoolean();
		npts=(int)gd.getNextNumber();
		fixstdev=gd.getNextBoolean();
		stdev=(float)gd.getNextNumber();
		gf=new gausfunc();
		lls=new linleastsquares();
		float[] fit=new float[width*height];
		System.arraycopy(pixels,0,fit,0,width*height);
		int startr=0; if(skipfirst){startr=1;}
		float[] data=new float[npts-startr];
		for(int i=startr;i<npts;i++) data[i-startr]=pixels[width/2+i];
		float[][] temp=null;
		if(fixstdev){
			temp=fit_data(data,stdev,stdev,(double)startr);
		} else {
			temp=fit_data(data,1.5,0.5*(double)npts,(double)startr);
		}
		for(int i=startr;i<npts;i++) {
			fit[width/2+i]=temp[1][i-startr];
			fit[width/2-i]=temp[1][i-startr];
		}
		TextWindow tw=new TextWindow("Fit Results","Stdev\tAmp\tBaseline\tc2","",400,200);
		tw.append(""+temp[0][0]+"\t"+temp[0][1]+"\t"+temp[0][2]+"\t"+temp[0][3]+"\n");
		for(int j=1;j<height;j++){
			startr=0; if(skiprest){startr=1;}
			data=new float[npts-startr];
			for(int i=startr;i<npts;i++) data[i-startr]=pixels[j*width+width/2+i];
			if(fixstdev){
				temp=fit_data(data,stdev,stdev,(double)startr);
			} else {
				temp=fit_data(data,1.5,0.5*(double)npts,(double)startr);
			}
			for(int i=startr;i<npts;i++) {
				fit[j*width+width/2+i]=temp[1][i-startr];
				fit[j*width+width/2-i]=temp[1][i-startr];
			}
			tw.append(""+temp[0][0]+"\t"+temp[0][1]+"\t"+temp[0][2]+"\t"+temp[0][3]);
		}
		new ImagePlus("Fit",new FloatProcessor(width,height,fit,null)).show();
	}

	public float[][] fit_data(float[] data,double minstdev,double maxstdev,double startr){
		float[][] outvals=new float[2][];
		outvals[0]=new float[4];
		double minc2=0.0;
		int npts=data.length;
		if(maxstdev>minstdev){
			for(double stdev=minstdev;stdev<maxstdev;stdev+=0.02){
				float[] func=gf.get_func(startr,data.length,1.0,stdev);
				float[] coef=lls.get_amp_offset(func,data,true);
				float[] fit=get_fit(func,coef);
				double c2=this_c2(data,fit);
				if(c2<minc2 || stdev==minstdev){
					minc2=c2;
					outvals[0][0]=(float)stdev; outvals[0][1]=coef[0]; outvals[0][2]=coef[1]; outvals[0][3]=(float)c2;
					outvals[1]=fit;
				}
			}
		} else {
			double stdev=minstdev;
			float[] func=gf.get_func(startr,data.length,1.0,stdev);
			float[] coef=lls.get_amp_offset(func,data,true);
			float[] fit=get_fit(func,coef);
			double c2=this_c2(data,fit);
			outvals[0][0]=(float)stdev; outvals[0][1]=coef[0]; outvals[0][2]=coef[1]; outvals[0][3]=(float)c2;
			outvals[1]=fit;
		}
		return outvals;
	}

	public float[] get_fit(float[] func,float[] coef){
		float[] temp=new float[func.length];
		for(int i=0;i<func.length;i++) temp[i]=coef[0]*func[i]+coef[1];
		return temp;
	}

	public double this_c2(float[] data,float[] fit){
		double c2=0.0;
		double npts=(double)data.length;
		for(int i=0;i<data.length;i++) c2+=(double)((data[i]-fit[i])*(data[i]-fit[i]))/npts;
		return c2;
	}

}
