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

public class stack_subtract_moving_statistic_jru_v1 implements PlugIn {

	public void run(String arg) {
		boolean outmovavg,addspatialconst,addtemporalconst;
		GenericDialog gd = new GenericDialog("subtract moving statistic");
		gd.addCheckbox("Output Moving Statistic?",true);
		String[] sitems={"No Const Add","Maintain Spatial Avg","Maintain Temporal Avg"};
		gd.addChoice("Avg Correction",sitems,sitems[0]);
		String[] statistics={"Avg","Median","Max","Min","Mode","ConditionalAvg"};
		gd.addChoice("Statistic to subtract",statistics,statistics[0]);
		int period=5;
		gd.addNumericField("period (must be odd)",period,0);
		boolean ends=false;
		gd.addCheckbox("include ends",ends);
		boolean divide=false;
		gd.addCheckbox("divide (instead of subtract)",divide);
		gd.showDialog();
		if(gd.wasCanceled()) return;
		outmovavg = gd.getNextBoolean();
		int index=gd.getNextChoiceIndex();
		addspatialconst=false;
		addtemporalconst=false;
		if(index==1){addspatialconst=true;}
		if(index==2){addtemporalconst=true;}
		int statindex=gd.getNextChoiceIndex();
		period = (int)gd.getNextNumber();
		ends=gd.getNextBoolean();
		divide=gd.getNextBoolean();
		ImagePlus imp = WindowManager.getCurrentImage();
		int width = imp.getWidth();
		int height = imp.getHeight();
		ImageStack stack = imp.getStack();
		ImageStack result_stack = new ImageStack(width,height);
		ImageStack movavg_stack = new ImageStack(width,height);
		int slices = stack.getSize();
		int periodshift = (period-1)/2;
		int startslice = periodshift;
		int endslice = slices-startslice-1;

		float[] extras=null;
		if(statindex==4){
			extras=getmodeextras();
			if(extras==null){return;}
		}
		if(statindex==5){
			extras=getcondavgextras();
			if(extras==null){return;}
		}

		//start by calculating the averages for the correction
		float[] corr_avg=new float[1];
		Rectangle nullrect=null;
		float[] nullfltarray=null;
		if(addspatialconst){
			corr_avg=new float[slices];
			for(int i=0;i<slices;i++){
				corr_avg[i]=jstatistics.getstatistic("Avg",stack.getPixels(i+1),0,0,nullrect,nullfltarray);
			}
		}
		if(addtemporalconst){
			corr_avg=new float[width*height];
			for(int i=0;i<width*height;i++){
				corr_avg[i]=getzstatistic(stack,0,slices-1,i,"Avg",null);
			}
		}

		//initialize the first statistic
		float[] newcurrstat=new float[width*height];
		for(int i=0;i<width*height;i++){
			newcurrstat[i]=getzstatistic(stack,0,startslice+periodshift,i,statistics[statindex],extras);
		}
		if(outmovavg){
			float[] temp5=new float[width*height];
			System.arraycopy(newcurrstat,0,temp5,0,width*height);
			movavg_stack.addSlice("",temp5);
		}

		if(ends){
			//for the ends, use the first period frames for the statistic
			for(int j=0;j<startslice;j++){
				float[] temp=new float[width*height];
				Object temp2=stack.getPixels(j+1);
				if(temp2 instanceof float[]){
					for(int i=0;i<width*height;i++){
						if(!divide) temp[i]=((float[])temp2)[i]-newcurrstat[i];
						else temp[i]=((float[])temp2)[i]/newcurrstat[i];
					}
				} else {
					for(int i=0;i<width*height;i++){
						float temp3=((short[])temp2)[i]&0xffff;
						if(!divide) temp[i]=temp3-newcurrstat[i];
						else temp[i]=temp3/newcurrstat[i];
					}
				}
				//add the average intensity back in if called for
				if(addspatialconst){
					if(statindex==0){
						for(int k=0;k<width*height;k++){
							if(!divide) temp[k]+=corr_avg[j];
							else temp[k]*=corr_avg[j];
						}
					} else {
						float temp1=jstatistics.getstatistic("Avg",(Object)temp,0,0,nullrect,nullfltarray);
						for(int k=0;k<width*height;k++){
							if(!divide) temp[k]+=corr_avg[j]-temp1;
							else temp[k]*=corr_avg[j]-temp1;
						}
					}
				}
				result_stack.addSlice("",(Object)temp);
			}
		}

		for(int i=startslice;i<=endslice;i++){
			//get the moving average image
			if(i!=startslice){
				if(statindex==0 && stack.getProcessor(1) instanceof FloatProcessor){
					float[] pixels1=(float[])stack.getPixels(i-periodshift+1);
					float[] pixels2=(float[])stack.getPixels(i+periodshift+1);
					//just subtract the first image from the average and add a new image
					for(int k=0;k<width*height;k++){newcurrstat[k]-=pixels1[k]/((float)period)-pixels2[k]/((float)period);}
				} else {
					for(int k=0;k<width*height;k++){
						newcurrstat[k]=getzstatistic(stack,i-periodshift,i+periodshift,k,statistics[statindex],extras);
					}
				}
			}
			//now subtract it from the current image
			float[] currimage = new float[width*height];
			if(stack.getProcessor(1) instanceof FloatProcessor){
				for(int k=0;k<width*height;k++){currimage[k]=((float[])stack.getPixels(i+1))[k]-newcurrstat[k];}
			} else {
				for(int k=0;k<width*height;k++){
					float temp=((short[])stack.getPixels(i+1))[k]&0xffff;
					if(!divide) currimage[k]=temp-newcurrstat[k];
					else currimage[k]=temp/newcurrstat[k];
				}
			}
			//add the average intensity back in if called for
			if(addspatialconst){
				if(statindex==0){
					for(int k=0;k<width*height;k++){
						if(!divide) currimage[k]+=corr_avg[i];
						else currimage[k]*=corr_avg[i];
					}
				} else {
					float temp=jstatistics.getstatistic("Avg",(Object)currimage,0,0,nullrect,nullfltarray);
					for(int k=0;k<width*height;k++){
						if(!divide) currimage[k]+=corr_avg[i]-temp;
						else currimage[k]*=corr_avg[i]-temp;
					}
				}
			}
			result_stack.addSlice("",(Object)currimage);
			if(outmovavg && i!=startslice){
				float[] temp5=new float[width*height];
				System.arraycopy(newcurrstat,0,temp5,0,width*height);
				movavg_stack.addSlice(null,temp5);
			}
			IJ.showProgress(i-startslice,endslice-startslice);
		}

		if(ends){
			//for the ends, use the first period frames for the statistic
			for(int j=endslice+1;j<slices;j++){
				float[] temp=new float[width*height];
				Object temp2=stack.getPixels(j+1);
				if(temp2 instanceof float[]){
					for(int i=0;i<width*height;i++){
						if(!divide) temp[i]=((float[])temp2)[i]-newcurrstat[i];
						else temp[i]=((float[])temp2)[i]/newcurrstat[i];
					}
				} else {
					for(int i=0;i<width*height;i++){
						float temp3=((short[])temp2)[i]&0xffff;
						if(!divide) temp[i]=temp3-newcurrstat[i];
						else temp[i]=temp3/newcurrstat[i];
					}
				}
				if(addspatialconst){
					if(statindex==0){
						for(int k=0;k<width*height;k++){
							if(!divide) temp[k]+=corr_avg[j];
							else temp[k]*=corr_avg[j];
						}
					} else {
						float temp1=jstatistics.getstatistic("Avg",(Object)temp,0,0,nullrect,nullfltarray);
						for(int k=0;k<width*height;k++){
							if(!divide) temp[k]+=corr_avg[j]-temp1;
							else temp[k]*=corr_avg[j]-temp1;
						}
					}
				}
				result_stack.addSlice("",temp);
			}
		}

		if(addtemporalconst){
			int outslices=result_stack.getSize();
			if(statindex==0){
				for(int i=0;i<outslices;i++){
					float[] pixels=(float[])result_stack.getPixels(i+1);
					for(int j=0;j<width*height;j++){
						if(!divide) pixels[j]+=corr_avg[j];
						else pixels[j]*=corr_avg[j];
					}
				}
			} else {
				float[] curravg=new float[width*height];
				for(int i=0;i<width*height;i++){
					curravg[i]=getzstatistic(result_stack,0,outslices-1,i,"Avg",null);
				}
				for(int i=0;i<outslices;i++){
					float[] pixels=(float[])result_stack.getPixels(i+1);
					for(int j=0;j<width*height;j++){
						if(!divide) pixels[j]+=corr_avg[j]-curravg[i];
						else pixels[j]*=corr_avg[j]-curravg[i];
					}
				}
			}
		}
		//output the results
		ImagePlus imp2 = new ImagePlus("Subtracted Image",result_stack);
		imp2.show();
		if(outmovavg){
			ImagePlus imp3 = new ImagePlus("Moving Statistic",movavg_stack);
			imp3.show();
		}
	}

	private float[] getmodeextras(){
		float[] extras=new float[3];
		GenericDialog gd=new GenericDialog("Mode Options");
		extras[0]=0.0f;
		gd.addNumericField("Histogram Bins",(int)extras[0],0);
		extras[1]=0.0f;
		gd.addNumericField("Histogram Start",extras[1],5,10,null);
		extras[2]=100.0f;
		gd.addNumericField("Histogram End",extras[2],5,10,null);
		gd.showDialog();  if(gd.wasCanceled()){return null;}
		extras[0]=(float)gd.getNextNumber();
		extras[1]=(float)gd.getNextNumber();
		extras[2]=(float)gd.getNextNumber();
		return extras;
	}

	private float[] getcondavgextras(){
		float[] extras=new float[2];
		GenericDialog gd=new GenericDialog("ConditionalAvg Options");
		extras[0]=100.0f;
		gd.addNumericField("Upper Limit",extras[0],5,10,null);
		extras[1]=0.0f;
		gd.addNumericField("Lower Limit",extras[1],5,10,null);
		gd.showDialog();  if(gd.wasCanceled()){return null;}
		extras[0]=(float)gd.getNextNumber();
		extras[1]=(float)gd.getNextNumber();
		return extras;
	}

	private float getzstatistic(ImageStack stack,int zstart,int zend,int x,int y,String stat,float[] extras){
		int pixel=(y*stack.getWidth())+x;
		return getzstatistic(stack,zstart,zend,pixel,stat,extras);
	}

	private float getzstatistic(ImageStack stack,int zstart,int zend,int pixel,String stat,float[] extras){
		Rectangle nullrect=null;
		if(stack.getProcessor(1) instanceof FloatProcessor){
			float[] temp=new float[zend-zstart+1];
			for(int j=zstart;j<=zend;j++){temp[j-zstart]=((float[])stack.getPixels(j+1))[pixel];}
			return jstatistics.getstatistic(stat,(Object)temp,0,0,nullrect,extras);
		} else {
			short[] temp=new short[zend-zstart+1];
			for(int j=zstart;j<=zend;j++){temp[j-zstart]=((short[])stack.getPixels(j+1))[pixel];}
			return jstatistics.getstatistic(stat,(Object)temp,0,0,nullrect,extras);
		}
	}
}
