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

public class stack_detrend_linear_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean addavg;
		GenericDialog gd=new GenericDialog("Detrend Options");
		int segments=1;
		gd.addNumericField("Segments",segments,0);
		addavg=true;
		gd.addCheckbox("Maintain Temporal Average?",addavg);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		segments=(int)gd.getNextNumber();
		addavg=gd.getNextBoolean();
		ImagePlus imp = WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		int segheight=slices/segments;
		ImageStack resstack=new ImageStack(width,height);
		for(int i=0;i<(segments*segheight);i++){
			float[] pixels=new float[width*height];
			Object temp=stack.getPixels(i+1);
			if(temp instanceof float[]){
				for(int j=0;j<width*height;j++){
					pixels[j]=((float[])temp)[j];
				}
			}else{
				for(int j=0;j<width*height;j++){
					pixels[j]=(float)(((short[])temp)[j]&0xffff);
				}
			}
			resstack.addSlice("",(Object)pixels);
		}
		for(int i=0;i<width*height;i++){
			float avg=0.0f;
			for(int k=0;k<segments;k++){
				float[] temp=new float[segheight];
				for(int j=0;j<segheight;j++){
					temp[j]=((float[])resstack.getPixels(j+1+k*segheight))[i];
				}
				float[] params=new float[3];
				gettrend(temp,params);
				avg+=params[2];
				for(int j=0;j<segheight;j++){
					float[] pixels=(float[])resstack.getPixels(j+1+k*segheight);
					pixels[i]-=(params[0]*(float)j+params[1]);
				}
			}
			avg/=(float)segments;
			if(addavg){
				for(int j=0;j<(segments*segheight);j++){
					float[] pixels=(float[])resstack.getPixels(j+1);
					pixels[i]+=avg;
				}
			}
			IJ.showProgress(i,width*height);
		}
		ImagePlus imp2=new ImagePlus("Detrended",resstack);
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
