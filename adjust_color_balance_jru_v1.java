/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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

public class adjust_color_balance_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Red_Multiplier",1.077f,5,15,null);
		gd.addNumericField("Green_Multiplier",1.057f,5,15,null);
		gd.addNumericField("Blue_Multiplier",1.27f,5,15,null);
		gd.addNumericField("Red_Offset",0.0f,5,15,null);
		gd.addNumericField("Green_Offset",0.0f,5,15,null);
		gd.addNumericField("Blue_Offset",0.0f,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float rmult=(float)gd.getNextNumber();
		float gmult=(float)gd.getNextNumber();
		float bmult=(float)gd.getNextNumber();
		float roff=(float)gd.getNextNumber();
		float goff=(float)gd.getNextNumber();
		float boff=(float)gd.getNextNumber();
		if(imp.getProcessor() instanceof ColorProcessor){
			ColorProcessor cp=(ColorProcessor)imp.getProcessor();
			int[] pix=(int[])cp.getPixels();
			for(int i=0;i<pix.length;i++){
				int[] rgb=jutils.intval2rgb(pix[i]);
				rgb[0]=(int)(roff+(float)rgb[0]*rmult); 
				if(rgb[0]>255) rgb[0]=255;
				if(rgb[0]<0) rgb[0]=0;
				rgb[1]=(int)(goff+(float)rgb[1]*gmult);
				if(rgb[1]>255) rgb[1]=255;
				if(rgb[1]<0) rgb[1]=0;
				rgb[2]=(int)(boff+(float)rgb[2]*bmult);
				if(rgb[2]>255) rgb[2]=255;
				if(rgb[2]<0) rgb[2]=0;
				pix[i]=jutils.rgb2intval(rgb[0],rgb[1],rgb[2]);
			}
		} else {
			//if it's not a color processor, needs to have exactly 3 channels
			ImageStack stack=imp.getStack();
			for(int j=0;j<stack.getSize();j++){
				if(imp.getProcessor() instanceof ByteProcessor){
					byte[] pix=(byte[])stack.getPixels(j+1);
					for(int i=0;i<pix.length;i++){
						int temp=(int)(pix[i]&0xff);
						if(j%3==0) temp=(int)(roff+(float)temp*rmult);
						if(j%3==1) temp=(int)(goff+(float)temp*gmult);
						if(j%3==2) temp=(int)(boff+(float)temp*bmult);
						if(temp>255) temp=255;
						if(temp<0) temp=0;
						pix[i]=(byte)temp;
					}
				}
				if(imp.getProcessor() instanceof ShortProcessor){
					short[] pix=(short[])stack.getPixels(j+1);
					for(int i=0;i<pix.length;i++){
						int temp=(int)(pix[i]&0xffff);
						if(j%3==0) temp=(int)(roff+(float)temp*rmult);
						if(j%3==1) temp=(int)(goff+(float)temp*gmult);
						if(j%3==2) temp=(int)(boff+(float)temp*bmult);
						if(temp>65535) temp=65535;
						if(temp<0) temp=0;
						pix[i]=(short)temp;
					}
				}
			}
		}
		imp.updateAndDraw();
	}

}
