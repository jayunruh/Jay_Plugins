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
import jalgs.jfft.*;

public class carpet_stics_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int height=imp.getHeight();
		int width=imp.getWidth();
		GenericDialog gd = new GenericDialog("Options");
		int slength=height/2;
		gd.addNumericField("STICS length",slength,0);
		boolean extrapg0=true;
		gd.addCheckbox("Extrapolate G(0)?",extrapg0);
		boolean allstics=false;
		gd.addCheckbox("Use 2D FFT?",allstics);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		slength=(int)gd.getNextNumber();
		extrapg0=gd.getNextBoolean();
		allstics=gd.getNextBoolean();
		FloatProcessor fp=(FloatProcessor)imp.getProcessor();
		float[] pixels=(float[])fp.getPixels();
		if(!allstics){
			float[] stics=new float[width*slength];
			crosscorr ccclass=new crosscorr(width);
			for(int k=0;k<slength;k++){
				for(int i=0;i<(height-k);i++){
					float[] temp1=new float[width];
					float[] temp2=new float[width];
					for(int j=0;j<width;j++){
						temp1[j]=pixels[i*width+j];
						temp2[j]=pixels[(i+k)*width+j];
					}
					float[] temp3=ccclass.docrosscorr(temp1,temp2);
					for(int j=0;j<width;j++){
						int shift=j-width/2;
						if(shift<0){shift+=width;}
						stics[k*width+j]+=temp3[shift]/(float)(height-k);
					}
				}
				IJ.showProgress(k,slength);
			}
			if(extrapg0){stics[width/2]=(stics[width/2-1]+stics[width/2+1])/2.0f;}
			FloatProcessor fp2=new FloatProcessor(width,slength,stics,null);
			ImagePlus imp2=new ImagePlus("STICS",fp2);
			imp2.copyScale(imp);
			imp2.show();
		} else {
			autocorr2D acclass=new autocorr2D(width,height);
			float[] stics=acclass.doautocorr2D(pixels,true,false);
			float[] newstics=new float[width*(height/2)];
			System.arraycopy(stics,0,newstics,0,width*(height/2));
			stics=newstics;
			if(extrapg0){stics[width/2]=(stics[width/2-1]+stics[width/2+1])/2.0f;}
			FloatProcessor fp2=new FloatProcessor(width,height/2,stics,null);
			ImagePlus imp2=new ImagePlus("STICS",fp2);
			imp2.copyScale(imp);
			imp2.show();
		}

	}
}
