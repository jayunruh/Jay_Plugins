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

public class carpet_detrend_linear_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Detrend Options");
		int segments=1;
		gd.addNumericField("Segments",segments,0);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		segments=(int)gd.getNextNumber();
		ImagePlus imp = WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int segheight=height/segments;
		float[] result=new float[width*segments*segheight];
		ImageProcessor ip=(ImageProcessor)imp.getProcessor();
		if(ip instanceof FloatProcessor){
			float[] pixels=(float[])ip.getPixels();
			for(int i=0;i<width;i++){
				float avg=0.0f;
				for(int k=0;k<segments;k++){
					float[] temp=new float[segheight];
					float[] params=new float[3];
					for(int j=0;j<segheight;j++){
						temp[j]=pixels[j*width+i+k*segheight*width];
						avg+=pixels[j*width+i+k*segheight*width];
					}
					gettrend(temp,params);
					for(int j=0;j<segheight;j++){
						result[j*width+i+k*segheight*width]=pixels[j*width+i+k*segheight*width]-params[0]*(float)j-params[1];
					}
				}
				avg/=(float)(segments*segheight);
				for(int j=0;j<(segments*segheight);j++){
					result[j*width+i]+=avg;
				}
				IJ.showProgress(i,width);
			}
		}
		if(ip instanceof ShortProcessor){
			short[] pixels=(short[])ip.getPixels();
			for(int i=0;i<width;i++){
				float avg=0.0f;
				for(int k=0;k<segments;k++){
					float[] temp=new float[segheight];
					float[] params=new float[3];
					for(int j=0;j<segheight;j++){
						temp[j]=pixels[j*width+i+k*segheight*width]&0xffff;
						avg+=temp[j];
					}
					gettrend(temp,params);
					for(int j=0;j<segheight;j++){
						float temp5=pixels[j*width+i+k*segheight*width]&0xffff;
						result[j*width+i+k*segheight*width]=temp5-params[0]*(float)j-params[1];
					}
				}
				avg/=(float)(segments*segheight);
				for(int j=0;j<(segments*segheight);j++){
					result[j*width+i]+=avg;
				}
				IJ.showProgress(i,width);
			}
		}
		FloatProcessor fp2=new FloatProcessor(width,segments*segheight,result,null);
		ImagePlus imp2=new ImagePlus("Detrended",fp2);
		imp2.show();
	}

	public void gettrend(float[] data, float[] params){
		//assume x values start from zero and increment by 1
		int length=data.length;
		float xsum=0.0f;
		float ysum=0.0f;
		float xsqsum=0.0f;
		float xysum=0.0f;
		for(int i=0;i<length;i++){
			xsum+=(float)i;
			ysum+=data[i];
			xsqsum+=(float)(i*i);
			xysum+=data[i]*(float)i;
		}
		float dumflt=1.0f/(xsqsum*(float)length-xsum*xsum);
		params[1]=dumflt*(xsqsum*ysum-xsum*xysum);
		params[0]=dumflt*(xysum*(float)length-xsum*ysum);
		params[2]=ysum/(float)length;
		return;
	}
}
