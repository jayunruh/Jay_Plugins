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
import java.util.*;
import jalgs.jsim.*;

public class sim_FLIM_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin simulates photon counting or analog lifetime data
		GenericDialog gd=new GenericDialog("Options");
		int channels=32;
		gd.addNumericField("Channels",channels,0);
		int width=256;
		gd.addNumericField("Image Width",width,0);
		int height=256;
		gd.addNumericField("Image Height",height,0);
		float tau1=4.0f;
		gd.addNumericField("tau1 (channels)",tau1,5,10,null);
		float amp1=100.0f;
		gd.addNumericField("amp1 (channels)",amp1,5,10,null);
		float tau2=2.0f;
		gd.addNumericField("tau2 (channels)",tau2,5,10,null);
		float amp2=0.0f;
		gd.addNumericField("amp2 (channels)",amp2,5,10,null);
		float tau3=1.0f;
		gd.addNumericField("tau3 (channels)",tau3,5,10,null);
		float amp3=0.0f;
		gd.addNumericField("amp3 (channels)",amp3,5,10,null);
		float irfstdev=0.0f;
		gd.addNumericField("irfstdev (channels)",irfstdev,5,10,null);
		float irfoffset=0.0f;
		gd.addNumericField("irfoffset (channels)",irfoffset,5,10,null);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		channels=(int)gd.getNextNumber();
		width=(int)gd.getNextNumber();
		height=(int)gd.getNextNumber();
		tau1=(float)gd.getNextNumber();
		amp1=(float)gd.getNextNumber();
		tau2=(float)gd.getNextNumber();
		amp2=(float)gd.getNextNumber();
		tau3=(float)gd.getNextNumber();
		amp3=(float)gd.getNextNumber();
		irfstdev=(float)gd.getNextNumber();
		irfoffset=(float)gd.getNextNumber();
		float[] profile=time_profile(tau1,amp1,tau2,amp2,tau3,amp3,channels,irfstdev,irfoffset);
		ImageStack result_stack=new ImageStack(width,height);
		rngs rngclass=new rngs();
		for(int j=0;j<channels;j++){
			float[] temp=new float[width*height];
			for(int i=0;i<width*height;i++){
				temp[i]=(float)rngclass.poidev((double)profile[j]);
				//temp[i]=profile[j];
			}
			result_stack.addSlice("",(Object)temp);
			IJ.showProgress(j,channels);
		}
		(new ImagePlus("Lifetime Sim",result_stack)).show();
	}

	private float[] time_profile(float tau1,float amp1,float tau2,float amp2,float tau3,float amp3,int channels,float irfstdev,float irfoffset){
		float[] irf=new float[channels];
		float[] result=new float[channels];
		if(irfstdev==0.0f){irf[(int)irfoffset]=1.0f;}
		else{
			float amp=1.0f/(irfstdev*(float)Math.sqrt(2.0*Math.PI));
			for(int i=-channels;i<0;i++){
				irf[i+channels]+=amp*(float)Math.exp((((float)i-irfoffset)*((float)i-irfoffset))/(-2.0f*irfstdev*irfstdev));
			}
			for(int i=0;i<channels;i++){
				irf[i]+=amp*(float)Math.exp((((float)i-irfoffset)*((float)i-irfoffset))/(-2.0f*irfstdev*irfstdev));
			}
			for(int i=channels;i<2*channels;i++){
				irf[i-channels]+=amp*(float)Math.exp((((float)i-irfoffset)*((float)i-irfoffset))/(-2.0f*irfstdev*irfstdev));
			}
			float sum=0.0f;
			for(int i=0;i<channels;i++){
				sum+=irf[i];
			}
			for(int i=0;i<channels;i++){
				irf[i]/=sum;
			}
		}
		for(int j=0;j<channels;j++){
			for(int i=0;i<10*channels;i++){
				int chanval=(j+i)%channels;
				result[chanval]+=irf[j]*amp1*(float)Math.exp(-((float)i)/tau1);
				if(amp2>0.0f){result[chanval]+=irf[j]*amp2*(float)Math.exp(-((float)i)/tau2);}
				if(amp3>0.0f){result[chanval]+=irf[j]*amp3*(float)Math.exp(-((float)i)/tau3);}
			}
		}
		return result;
		//return irf;
	}

}
