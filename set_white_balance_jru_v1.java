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
import jalgs.*;

public class set_white_balance_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		int offset=-5;
		int whiteval=240;
		if(imp.getProcessor() instanceof ShortProcessor){
			offset=-10;
			whiteval=60000;
		}
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Offset (new 0 value)",offset,0);
		gd.addNumericField("White Value",whiteval,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		offset=(int)gd.getNextNumber();
		whiteval=(int)gd.getNextNumber();
		Roi roi=imp.getRoi();
		Polygon poly=roi.getPolygon();
		boolean[] mask=jstatistics.poly2mask(poly,width,height);
		float rmean=0.0f; float gmean=0.0f; float bmean=0.0f;
		int currslice=0; int currframe=0; int currchan=0;

		//start by measuring the intensity in the red, green, and blue channels
		if(imp.getProcessor() instanceof ColorProcessor){
			byte[][] rgbvals=jutils.intval2rgb((int[])imp.getProcessor().getPixels());
			rmean=jstatistics.getstatistic("Avg",rgbvals[0],width,height,mask,null);
			gmean=jstatistics.getstatistic("Avg",rgbvals[1],width,height,mask,null);
			bmean=jstatistics.getstatistic("Avg",rgbvals[2],width,height,mask,null);
		} else {
			currslice=imp.getZ()-1;
			currframe=imp.getT()-1;
			currchan=imp.getC()-1;
			Object[] cstack=jutils.get3DCSeries(imp.getStack(),currslice,currframe,imp.getNFrames(),imp.getNSlices(),imp.getNChannels());
			float[] spectrum=jstatistics.getspectrum("Avg",cstack,width,height,mask,null);
			rmean=spectrum[0]; gmean=spectrum[1]; bmean=spectrum[2];
		}

		//now calculate the offsets and multipliers
		float rmult=(float)whiteval/(rmean-(float)offset);
		float gmult=(float)whiteval/(gmean-(float)offset);
		float bmult=(float)whiteval/(bmean-(float)offset);
		float roff=0.0f-(float)offset-rmult;
		float goff=0.0f-(float)offset-gmult;
		float boff=0.0f-(float)offset-bmult;
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
		if(imp.getProcessor() instanceof ShortProcessor){
			for(int i=0;i<3;i++){
				imp.setPositionWithoutUpdate(i+1,currslice+1,currframe+1);
				imp.setDisplayRange(0,65535);
			}
			imp.setPosition(currchan+1,currslice+1,currframe+1);
			imp.updateChannelAndDraw();
		}
		imp.updateAndDraw();
	}

}
