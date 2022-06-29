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
import java.awt.event.*;
import java.awt.image.*;
import ij.plugin.*;
import ij.measure.*;
import java.io.*;
import jguis.*;
import jalgs.*;

public class histogram_2D_phasor_jru_v2 implements PlugIn {

	public void run(String arg) {

		int i;
		int[] wList = WindowManager.getIDList();
		String[] titles = new String[wList.length];
		for(i=0;i<wList.length;i++){
			ImagePlus imp = WindowManager.getImage(wList[i]);
			if(imp!=null){titles[i]=imp.getTitle();}
			else{titles[i]="";}
		}
		GenericDialog gd = new GenericDialog("Choose Images");
		boolean calcgs=true;
		gd.addCheckbox("Calculate G and S?",calcgs);
		int harmonic=1;
		gd.addNumericField("Harmonic (for G S Calc)",harmonic,0);
		gd.addChoice("FLIM Stack (for G S Calc)",titles,titles[0]);
		gd.addChoice("G image",titles,titles[0]);
		gd.addChoice("S image",titles,titles[0]);
		gd.addChoice("display image",titles,titles[0]);
		gd.addNumericField("Phase_Shift (deg)",0.0f,5,15,null);
		gd.addNumericField("Mod_Scale (fraction)",1.0f,5,15,null);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		int[] index = new int[4];
		calcgs=gd.getNextBoolean();
		harmonic=(int)gd.getNextNumber();
		for(i=0;i<4;i++){
			index[i]=gd.getNextChoiceIndex();
		}
		float pshift=(float)gd.getNextNumber();
		float mscale=(float)gd.getNextNumber();
		ImagePlus impstack=WindowManager.getImage(wList[index[0]]);
		ImagePlus imp1 = WindowManager.getImage(wList[index[1]]);
		ImagePlus imp2 = WindowManager.getImage(wList[index[2]]);
		ImagePlus imp3 = WindowManager.getImage(wList[index[3]]);

		if(calcgs){
			ImagePlus[] imps=stack2gsi(impstack.getStack(),harmonic,pshift,mscale);
			imp1=imps[0];
			imp2=imps[1];
			imp3=imps[2];
		} else {
			impstack=null;
		}

		final Hist2DWindow_v2 cw=new Hist2DWindow_v2();
		//cw.init(imp1,imp2,null,impstack,imp3,5);
		float[] xpix=(float[])imp1.getProcessor().convertToFloat().getPixels();
		float[] ypix=(float[])imp2.getProcessor().convertToFloat().getPixels();
		float[] dpix=(float[])imp3.getProcessor().convertToFloat().getPixels();
		cw.xlab="G";
		cw.ylab="S";
		cw.init(imp3.getWidth(),imp3.getHeight(),xpix,ypix,dpix,impstack);
		Hist2DWindow_v2.launch_frame(cw);
	}

	ImagePlus[] stack2gsi(ImageStack stack,int harmonic,float pshift,float mscale){
		int slices=stack.getSize();
		int width=stack.getWidth();
		int height=stack.getHeight();
		float[] intensity=new float[width*height];
		float[] G=new float[width*height];
		float[] S=new float[width*height];
		float[] cosvals=new float[slices];
		float[] sinvals=new float[slices];
		float pshiftfrac=pshift/360.0f;
		if(slices==4){
			//should correct for phase shift here somehow in the future
			cosvals=new float[]{mscale,-mscale,-mscale,mscale};
			sinvals=new float[]{mscale,mscale,-mscale,-mscale};
		} else {
			for(int i=0;i<slices;i++){
				cosvals[i]=mscale*(float)Math.cos(2.0*Math.PI*((double)(i*harmonic)/(double)slices-pshiftfrac));
				sinvals[i]=mscale*(float)Math.sin(2.0*Math.PI*((double)(i*harmonic)/(double)slices-pshiftfrac));
			}
		}
		if(stack.getProcessor(1) instanceof FloatProcessor){
			for(int i=0;i<(width*height);i++){
				for(int j=0;j<slices;j++){
					float temp=((float[])stack.getPixels(j+1))[i];
					intensity[i]+=temp;
					G[i]+=temp*cosvals[j];
					S[i]+=temp*sinvals[j];
				}
				if(intensity[i]>0.0f){
					G[i]/=intensity[i];
					S[i]/=intensity[i];
				}
			}
		}
		if(stack.getProcessor(1) instanceof ShortProcessor){
			for(int i=0;i<(width*height);i++){
				for(int j=0;j<slices;j++){
					float temp=(float)(((short[])stack.getPixels(j+1))[i]&0xffff);
					intensity[i]+=temp;
					G[i]+=temp*cosvals[j];
					S[i]+=temp*sinvals[j];
				}
				if(intensity[i]>0.0f){
					G[i]/=intensity[i];
					S[i]/=intensity[i];
				}
			}
		}
		if(stack.getProcessor(1) instanceof ByteProcessor){
				for(int i=0;i<(width*height);i++){
					for(int j=0;j<slices;j++){
						float temp=(float)(((byte[])stack.getPixels(j+1))[i]&0xff);
						intensity[i]+=temp;
						G[i]+=temp*cosvals[j];
						S[i]+=temp*sinvals[j];
					}
					if(intensity[i]>0.0f){
						G[i]/=intensity[i];
						S[i]/=intensity[i];
					}
				}
			}
		ImagePlus[] imps=new ImagePlus[3];
		imps[0]=new ImagePlus("G",new FloatProcessor(width,height,G,null));
		imps[1]=new ImagePlus("S",new FloatProcessor(width,height,S,null));
		imps[2]=new ImagePlus("Intensity",new FloatProcessor(width,height,intensity,null));
		return imps;
	}
}
