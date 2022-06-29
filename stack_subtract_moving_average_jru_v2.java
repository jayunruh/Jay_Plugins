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

public class stack_subtract_moving_average_jru_v2 implements PlugIn {
	int slices,channels;
	ImageStack stack;

	public void run(String arg) {
		boolean outmovavg,addspatialconst,addtemporalconst;
		GenericDialog gd = new GenericDialog("subtract moving average");
		gd.addCheckbox("Output (Moving) Average?",true);
		gd.addCheckbox("Subtract Static Average?",false);
		String[] sitems={"No Const Add","Maintain Spatial Average","Maintain Temporal Average"};
		gd.addChoice("Avg Correction",sitems,sitems[0]);
		gd.addNumericField("period (must be odd)",5,0);
		gd.addCheckbox("Output_Sub_Stack",true);
		gd.showDialog();
		if(gd.wasCanceled()) return;
		outmovavg = gd.getNextBoolean();
		boolean isstatic=gd.getNextBoolean();
		int index=gd.getNextChoiceIndex();
		addspatialconst=false;
		addtemporalconst=false;
		if(index==1){addspatialconst=true;}
		if(index==2){addtemporalconst=true;}
		int period = (int)gd.getNextNumber();
		boolean outsub=gd.getNextBoolean();
		ImagePlus imp = WindowManager.getCurrentImage();
		Object[] outputs=exec(imp,outmovavg,isstatic,addspatialconst,addtemporalconst,period,outsub);
		for(int i=0;i<outputs.length;i++){
			((ImagePlus)outputs[i]).show();
		}
	}

	public Object[] exec(ImagePlus imp,boolean outmovavg,boolean isstatic,boolean addspatialconst,boolean addtemporalconst,int period,boolean outsub){
		int width = imp.getWidth();
		int height = imp.getHeight();
		channels=imp.getNChannels();
		stack = imp.getStack();
		slices=stack.getSize();
		slices/=channels;
		boolean isfloat=(getpix(0,0) instanceof float[]);
		ImageStack result_stack = new ImageStack(width,height);
		ImageStack movavg_stack = new ImageStack(width,height);
		int periodshift = (period-1)/2;
		if(isstatic){periodshift=0;}
		int startslice = periodshift;
		int endslice = slices-startslice-1;
		float[][] oldcurravg = null;
		float[][] totavg=null;
		if(!isstatic){
			oldcurravg = new float[channels][width*height];
			totavg=new float[channels][width*height];
		}
		float[][] newcurravg=new float[channels][width*height];;
		for(int i=startslice;i<=endslice;i++){
			if(!isstatic){
				newcurravg=new float[channels][width*height];
			}
			for(int l=0;l<channels;l++){
				//calculate the moving average
				if(!isstatic){
					if(i==startslice){
						//at the beginning calculate it directly
						for(int j=i-periodshift;j<=i+periodshift;j++){
							Object tempimg=getpix(l,j);
							for(int k=0;k<width*height;k++){
								float temp=0.0f;
								if(isfloat){temp=((float[])tempimg)[k];}
								else {temp=(float)(((short[])tempimg)[k]&0xffff);}
								newcurravg[l][k]+=temp/(float)period;
							}
						}
						if(addtemporalconst){
							for(int k=0;k<width*height;k++){
								totavg[l][k]=(newcurravg[l][k]*(float)period)/(float)slices;
							}
						}
					} else{
						//later, just subtract the first image from the average and add a new image
						Object tempsubimg=getpix(l,i-periodshift-1);
						Object tempaddimg=getpix(l,i+periodshift);
						for(int k=0;k<width*height;k++){
							float temp1,temp2;
							if(isfloat){temp1=((float[])tempsubimg)[k]; temp2=((float[])tempaddimg)[k];}
							else {temp1=(float)(((short[])tempsubimg)[k]&0xffff); temp2=(float)(((short[])tempaddimg)[k]&0xffff);}
							newcurravg[l][k]=oldcurravg[l][k]-temp1/(float)period+temp2/(float)period;
							if(addtemporalconst){
								totavg[l][k]+=temp2/(float)slices;
							}
						}
					}
				} else {
					if(i==startslice && l==0){
						for(int j=0;j<channels;j++){
							for(int k=0;k<slices;k++){
								Object tempimg=getpix(j,k);
								for(int m=0;m<width*height;m++){
									float temp=0.0f;
									if(isfloat){temp=((float[])tempimg)[m];}
									else {temp=(float)(((short[])tempimg)[m]&0xffff);}
									newcurravg[j][m]+=temp/(float)slices;
								}
							}
							if(outmovavg){
								movavg_stack.addSlice(null,newcurravg[j]);
							}
						}
						totavg=newcurravg;
					}
				}
				//now subtract it from the current image
				float[] currimage=new float[width*height];
				//IJ.log(""+newcurravg[l][0]);
				Object oldcurrimage=getpix(l,i);
				for(int k=0;k<width*height;k++){
					float temp;
					if(isfloat){temp=((float[])oldcurrimage)[k];}
					else{temp=(float)(((short[])oldcurrimage)[k]&0xffff);}
					currimage[k]=temp-newcurravg[l][k];
				}
				if(addspatialconst){
					float avgint=jstatistics.getstatistic("Avg",oldcurrimage,null);
					for(int k=0;k<width*height;k++){currimage[k]+=avgint;}
				}
				if(outsub) result_stack.addSlice(null,currimage);
				if(outmovavg && !isstatic){
					movavg_stack.addSlice(null,newcurravg[l]);
				}
				if(!isstatic){
					for(int k=0;k<width*height;k++){oldcurravg[l][k]=newcurravg[l][k];}
				}
			}
			IJ.showProgress(i-startslice,endslice-startslice);
		}
		int nouts=1; if(outsub && outmovavg) nouts=2;
		Object[] output=new Object[nouts];
		int outcounter=0;
		if(outsub){
			int outslices=result_stack.getSize();
			outslices/=channels;
			if(addtemporalconst){
				for(int i=0;i<outslices;i++){
					for(int j=0;j<channels;j++){
						float[] pixels=(float[])result_stack.getPixels(j+i*channels+1);
						for(int k=0;k<width*height;k++){
							pixels[k]+=totavg[j][k];
						}
					}
				}
			}
			//output the results
			ImagePlus imp2 = new ImagePlus("Subtracted Image",result_stack);
			imp2.copyScale(imp);
			imp2.setOpenAsHyperStack(true);
			imp2.setDimensions(channels,1,outslices);
			if(channels>1){
				output[outcounter]=new CompositeImage(imp2,CompositeImage.COMPOSITE); outcounter++;
			} else {
				output[outcounter]=imp2; outcounter++;
			}
		}
		if(outmovavg){
			int outslices=movavg_stack.getSize();
			outslices/=channels;
			ImagePlus imp3 = new ImagePlus("(Moving) Average",movavg_stack);
			imp3.copyScale(imp);
			imp3.setOpenAsHyperStack(true);
			imp3.setDimensions(channels,1,outslices);
			if(channels>1){
				output[outcounter]=new CompositeImage(imp3,CompositeImage.COMPOSITE);
			} else {
				output[outcounter]=imp3;
			}
		}
		return output;
	}

	public Object getpix(int channel,int slice){
		return stack.getPixels(channel+slice*channels+1);
	}

}
