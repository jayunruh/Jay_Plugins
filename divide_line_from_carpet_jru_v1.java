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

public class divide_line_from_carpet_jru_v1 implements PlugIn {

	public void run(String arg) {
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(int i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){
				titles[i]=imp.getTitle();
			} else {
				titles[i]="NA";
			}
		}
		GenericDialog gd = new GenericDialog("Options");
		gd.addChoice("Plot Window",titles,titles[0]);
		gd.addChoice("Carpet",titles,titles[0]);
		gd.addCheckbox("Horizontal_Line",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int index1 = gd.getNextChoiceIndex();
		int index2 = gd.getNextChoiceIndex();
		boolean horizontal=gd.getNextBoolean();
		ImagePlus imp = WindowManager.getImage(wList[index1]);
		ImagePlus imp2 = WindowManager.getImage(wList[index2]);

		ImageWindow iw=imp.getWindow();
		float[] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0];
		Object pixels=imp2.getProcessor().getPixels();
		int width=imp2.getWidth();
		int height=imp2.getHeight();
		if(horizontal){
			int length=yvals.length;
			if(pixels instanceof float[]){
				float[] pixels2=(float[])pixels;
				float[] newpixels=new float[pixels2.length];
				for(int i=0;i<height;i++){
					for(int j=0;j<length;j++){
						newpixels[j+i*width]=pixels2[j+i*width]/yvals[j];
					}
				}
				(new ImagePlus("Divided",new FloatProcessor(width,height,newpixels,null))).show();
			} else if (pixels instanceof short[]){
				short[] pixels2=(short[])pixels;
				short[] newpixels=new short[pixels2.length];
				for(int i=0;i<height;i++){
					for(int j=0;j<length;j++){
						float temp=(float)(pixels2[j+i*width]&0xffff)/yvals[j];
						int temp2=(int)temp; if(temp2<0) temp2=0; if(temp2>65535) temp2=65535;
						newpixels[j+i*width]=(short)(temp2);
					}
				}
				(new ImagePlus("Divided",new ShortProcessor(width,height,newpixels,null))).show();
			} else {
				byte[] pixels2=(byte[])pixels;
				byte[] newpixels=new byte[pixels2.length];
				for(int i=0;i<height;i++){
					for(int j=0;j<length;j++){
						float temp=(float)(pixels2[j+i*width]&0xff)/yvals[j];
						int temp2=(int)temp; if(temp2<0) temp2=0; if(temp2>65535) temp2=65535;
						newpixels[j+i*width]=(byte)(temp2);
					}
				}
				(new ImagePlus("Divided",new ByteProcessor(width,height,newpixels))).show();
			}
		} else {
			int length=yvals.length;
			if(pixels instanceof float[]){
				float[] pixels2=(float[])pixels;
				float[] newpixels=new float[pixels2.length];
				for(int i=0;i<length;i++){
					for(int j=0;j<width;j++){
						newpixels[j+i*width]=pixels2[j+i*width]/yvals[i];
					}
				}
				(new ImagePlus("Divided",new FloatProcessor(width,length,newpixels,null))).show();
			} else if (pixels instanceof short[]){
				short[] pixels2=(short[])pixels;
				short[] newpixels=new short[pixels2.length];
				for(int i=0;i<length;i++){
					for(int j=0;j<width;j++){
						float temp=(float)(pixels2[j+i*width]&0xffff)/yvals[i];
						int temp2=(int)temp; if(temp2<0) temp2=0; if(temp2>65535) temp2=65535;
						newpixels[j+i*width]=(short)(temp2);
					}
				}
				(new ImagePlus("Divided",new ShortProcessor(width,length,newpixels,null))).show();
			} else {
				byte[] pixels2=(byte[])pixels;
				byte[] newpixels=new byte[pixels2.length];
				for(int i=0;i<length;i++){
					for(int j=0;j<width;j++){
						float temp=(float)(pixels2[j+i*width]&0xff)/yvals[i];
						int temp2=(int)temp; if(temp2<0) temp2=0; if(temp2>65535) temp2=65535;
						newpixels[j+i*width]=(byte)(temp2);
					}
				}
				(new ImagePlus("Divided",new ByteProcessor(width,length,newpixels))).show();
			}
		}
	}

}
