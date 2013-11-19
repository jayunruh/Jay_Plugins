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
import java.io.*;

public class stack_statistics_jru_v1 implements PlugIn {
		//this plugin calculates statistics for each pixel in a stack
		//if frames to avg is less than the total number of frames, a movie of the statistic is created
		//this can be used for N and B analysis
		int avgframes,index;
		boolean avgall,batchmode,incthird,incfourth;

	public void run(String arg) {
		init_options();
		GenericDialog gd = new GenericDialog("Options");
		gd.addNumericField("Frames to Avg",avgframes,0);
		gd.addCheckbox("Avg all frames",avgall);
		String[] vartypes={"Variance","Covariance","Y Spatial Covar","X Spatial Covar","Temporal Covar"};
		gd.addChoice("Variance Type",vartypes,vartypes[index]);
		gd.addCheckbox("Batch mode",batchmode);
		gd.addCheckbox("Include third moment",incthird);
		gd.addCheckbox("Include fourth moment",incfourth);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		avgframes=(int)gd.getNextNumber();
		avgall=gd.getNextBoolean();
		index=gd.getNextChoiceIndex();
		batchmode=gd.getNextBoolean();
		incthird=gd.getNextBoolean();
		incfourth=gd.getNextBoolean();
		if(index==1){batchmode=false;}
		if(index!=0){incthird=false; incfourth=false;}
		
		int[] wList=WindowManager.getIDList();
		int images=1;
		if(batchmode){images=wList.length;}
		for(int im=0;im<images;im++){
			//get the current image and its info
			ImagePlus imp;
			if(!batchmode){imp = WindowManager.getCurrentImage();}
			else{imp=WindowManager.getImage(wList[im]);}
			String otitle=imp.getTitle();
			int namelen=otitle.lastIndexOf(".");
			String trunctitle;
			if(namelen>0){
				trunctitle=otitle.substring(0,namelen);
			} else{
				trunctitle=otitle;
			}
			int height = imp.getHeight();
			int width = imp.getWidth();
			ImageStack stack = imp.getStack();
			int slices = stack.getSize();
			int result_slices=(int)((float)slices/(float)avgframes);
			if(avgall){result_slices=1; avgframes=slices;}
			ImageProcessor ip=(ImageProcessor)imp.getProcessor();
			
			//now that we have the data, calculate the average and variance
			if(index==0){
				ImageStack avg_stack=new ImageStack(width,height); ImageStack var_stack=new ImageStack(width,height);
				ImageStack third_stack=new ImageStack(width,height); ImageStack fourth_stack=new ImageStack(width,height);
				for(int i=0;i<result_slices;i++){
					float[] avg=new float[width*height];
					float[] var=new float[width*height];
					float[] third=new float[width*height];
					float[] fourth=new float[width*height];
					for(int j=0;j<avgframes;j++){
						if(ip instanceof FloatProcessor){
							float[] temp=(float[])stack.getPixels(j+i*avgframes+1);
							for(int k=0;k<(width*height);k++){
								avg[k]+=temp[k]/(float)avgframes;
								var[k]+=(temp[k]*temp[k])/(float)avgframes;
								if(incthird || incfourth){third[k]+=(temp[k]*temp[k]*temp[k])/(float)avgframes;}
								if(incfourth){fourth[k]+=(temp[k]*temp[k]*temp[k]*temp[k])/(float)avgframes;}
							}
						}
						if(ip instanceof ShortProcessor){
							short[] temp=(short[])stack.getPixels(j+i*avgframes+1);
							for(int k=0;k<(width*height);k++){
								float temp2=temp[k]&0xffff;
								avg[k]+=temp2/(float)avgframes;
								var[k]+=(temp2*temp2)/(float)avgframes;
								if(incthird || incfourth){third[k]+=(temp2*temp2*temp2)/(float)avgframes;}
								if(incfourth){fourth[k]+=(temp2*temp2*temp2*temp2)/(float)avgframes;}
							}
						}
					}
					for(int j=0;j<(width*height);j++){
						float n=(float)avgframes;
						if(incfourth){
							fourth[j]=-(-3.0f*var[j]*var[j]+4.0f*avg[j]*third[j]+n*(6.0f*avg[j]*avg[j]*avg[j]*avg[j]-12.0f*avg[j]*avg[j]*var[j]+3.0f*var[j]*var[j]+4.0f*avg[j]*
							third[j]-fourth[j])-fourth[j])*((n*n)/(-6.0f+11.0f*n-6.0f*n*n+n*n*n));
						} 
						if(incthird){
							third[j]+=(2.0f*avg[j]*avg[j]*avg[j]-3.0f*var[j]*avg[j]);
							third[j]*=n*n/(n*n-3.0f*n+2.0f);
						}
						var[j]-=avg[j]*avg[j];
						var[j]*=n/(n-1.0f);
					}
					IJ.showProgress(i,result_slices);
					avg_stack.addSlice("",(Object)avg);
					var_stack.addSlice("",(Object)var);
					if(incthird){third_stack.addSlice("",(Object)third);}
					if(incfourth){fourth_stack.addSlice("",(Object)fourth);}
				}
				ImagePlus imp2=new ImagePlus(trunctitle+"_avg",avg_stack);
				imp2.show();
				ImagePlus imp3=new ImagePlus(trunctitle+"_var",var_stack);
				imp3.show();
				if(incthird){
					ImagePlus imp4=new ImagePlus(trunctitle+"_third",third_stack);
					imp4.show();
				}
				if(incfourth){
					new ImagePlus(trunctitle+"_fourth",fourth_stack).show();
				}
			}
			if(index==1){
				ImageStack avg_stack=new ImageStack(width,height); ImageStack var_stack=new ImageStack(width,height);
				String[] titles = new String[wList.length];
				for(int i=0;i<wList.length;i++){
					ImagePlus imp10 = WindowManager.getImage(wList[i]);
					if(imp10!=null){titles[i]=imp10.getTitle();}
					else{titles[i]="";}
				}
				GenericDialog gd1 = new GenericDialog("Options");
				gd1.addChoice("Stack 1",titles,titles[0]);
				gd1.addChoice("Stack 2",titles,titles[0]);
				gd1.showDialog();
				if(gd1.wasCanceled()){return;}
				int index1 = gd1.getNextChoiceIndex();
				int index2 = gd1.getNextChoiceIndex();
				ImagePlus imp1 = WindowManager.getImage(wList[index1]);
				ImagePlus imp2 = WindowManager.getImage(wList[index2]);
				ImageStack stack1 = imp1.getStack();
				ImageStack stack2=imp2.getStack();
				width=imp1.getWidth();
				height=imp1.getHeight();
				slices=stack1.getSize();
				result_slices=(int)((float)slices/(float)avgframes);
				if(avgall){result_slices=1; avgframes=slices;}
				ip=(ImageProcessor)imp1.getProcessor();
				for(int i=0;i<result_slices;i++){
					float[] avg1=new float[width*height];
					float[] avg2=new float[width*height];
					float[] var=new float[width*height];
					for(int j=0;j<avgframes;j++){
						if(ip instanceof FloatProcessor){
							float[] temp1=(float[])stack1.getPixels(j+i*avgframes+1);
							float[] temp2=(float[])stack2.getPixels(j+i*avgframes+1);
							for(int k=0;k<(width*height);k++){
								avg1[k]+=temp1[k]/(float)avgframes;
								avg2[k]+=temp2[k]/(float)avgframes;
								var[k]+=(temp1[k]*temp2[k])/(float)avgframes;
							}
						}
						if(ip instanceof ShortProcessor){
							short[] temp1=(short[])stack1.getPixels(j+i*avgframes+1);
							short[] temp2=(short[])stack2.getPixels(j+i*avgframes+1);
							for(int k=0;k<(width*height);k++){
								float temp3=temp1[k]&0xffff;
								float temp4=temp2[k]&0xffff;
								avg1[k]+=temp3/(float)avgframes;
								avg2[k]+=temp4/(float)avgframes;
								var[k]+=(temp3*temp4)/(float)avgframes;
							}
						}
					}
					for(int j=0;j<(width*height);j++){
						var[j]-=avg1[j]*avg2[j];
						avg1[j]=(float)Math.sqrt((double)(avg1[j]*avg2[j]));
					}
					IJ.showProgress(i,result_slices);
					avg_stack.addSlice("",(Object)avg1);
					var_stack.addSlice("",(Object)var);
				}
				ImagePlus imp3=new ImagePlus("covariance",var_stack);
				imp3.show();
				ImagePlus imp4=new ImagePlus("average",avg_stack);
				imp4.show();
			}
			
			if(index==2){
				ImageStack avg_stack=new ImageStack(width,height-1); ImageStack var_stack=new ImageStack(width,height-1);
				for(int i=0;i<result_slices;i++){
					float[] avg=new float[width*(height-1)];
					float[] var=new float[width*(height-1)];
					for(int j=0;j<avgframes;j++){
						if(ip instanceof FloatProcessor){
							float[] temp=(float[])stack.getPixels(j+i*avgframes+1);
							for(int k=0;k<(width*(height-1));k++){
								avg[k]+=(temp[k]+temp[k+width])/(float)(avgframes*2);
								var[k]+=(temp[k]*temp[k+width])/(float)avgframes;
							}
						}
						if(ip instanceof ShortProcessor){
							short[] temp=(short[])stack.getPixels(j+i*avgframes+1);
							for(int k=0;k<(width*(height-1));k++){
								float temp2=temp[k]&0xffff;
								float temp3=temp[k+width]&0xffff;
								avg[k]+=(temp2+temp3)/(float)(2*avgframes);
								var[k]+=(temp2*temp3)/(float)avgframes;
							}
						}
					}
					for(int j=0;j<(width*(height-1));j++){
						var[j]-=avg[j]*avg[j];
					}
					IJ.showProgress(i,result_slices);
					avg_stack.addSlice("",(Object)avg);
					var_stack.addSlice("",(Object)var);
				}
				ImagePlus imp2=new ImagePlus(trunctitle+"_g1yvar",var_stack);
				imp2.show();
				ImagePlus imp3=new ImagePlus(trunctitle+"_avg",avg_stack);
				imp3.show();
			}
			
			if(index==3){
				ImageStack avg_stack=new ImageStack(width-1,height); ImageStack var_stack=new ImageStack(width-1,height);
				for(int i=0;i<result_slices;i++){
					float[] avg=new float[(width-1)*height];
					float[] var=new float[(width-1)*height];
					for(int j=0;j<avgframes;j++){
						if(ip instanceof FloatProcessor){
							float[] temp=(float[])stack.getPixels(j+i*avgframes+1);
							for(int k=0;k<height;k++){
								for(int l=0;l<(width-1);l++){
									avg[k]+=(temp[k*width+l]+temp[k*width+l+1])/(float)(avgframes*2);
									var[k]+=(temp[k*width+l]*temp[k*width+l+1])/(float)avgframes;
								}
							}
						}
						if(ip instanceof ShortProcessor){
							short[] temp=(short[])stack.getPixels(j+i*avgframes+1);
							for(int k=0;k<height;k++){
								for(int l=0;l<(width-1);l++){
									float temp2=temp[k*width+l]&0xffff;
									float temp3=temp[k*width+l+1]&0xffff;
									avg[k]+=(temp2+temp3)/(float)(2*avgframes);
									var[k]+=(temp2*temp3)/(float)avgframes;
								}
							}
						}
					}
					for(int j=0;j<((width-1)*height);j++){
						var[j]-=avg[j]*avg[j];
					}
					IJ.showProgress(i,result_slices);
					avg_stack.addSlice("",(Object)avg);
					var_stack.addSlice("",(Object)var);
				}
				ImagePlus imp2=new ImagePlus(trunctitle+"_g1xvar",var_stack);
				imp2.show();
				ImagePlus imp3=new ImagePlus(trunctitle+"_avg",avg_stack);
				imp3.show();
			}
			
			if(index==4){
				ImageStack avg_stack=new ImageStack(width,height); ImageStack var_stack=new ImageStack(width,height);
				for(int i=0;i<result_slices;i++){
					float[] avg=new float[width*height];
					float[] var=new float[width*height];
					for(int j=0;j<(avgframes-1);j++){
						if(ip instanceof FloatProcessor){
							float[] temp=(float[])stack.getPixels(j+i*avgframes+1);
							float[] temp2=(float[])stack.getPixels(j+i*avgframes+2);
							for(int k=0;k<(width*height);k++){
								avg[k]+=(temp[k]+temp2[k])/(float)(2*avgframes);
								var[k]+=(temp[k]*temp2[k])/(float)avgframes;
							}
						}
						if(ip instanceof ShortProcessor){
							short[] temp=(short[])stack.getPixels(j+i*avgframes+1);
							short[] temp2=(short[])stack.getPixels(j+i*avgframes+2);
							for(int k=0;k<(width*height);k++){
								float temp3=temp[k]&0xffff;
								float temp4=temp2[k]&0xffff;
								avg[k]+=(temp3+temp4)/(float)(2*avgframes);
								var[k]+=(temp3*temp4)/(float)avgframes;
							}
						}
					}
					for(int j=0;j<(width*height);j++){
						var[j]-=avg[j]*avg[j];
					}
					IJ.showProgress(i,result_slices);
					avg_stack.addSlice("",(Object)avg);
					var_stack.addSlice("",(Object)var);
				}
				ImagePlus imp2=new ImagePlus(trunctitle+"_g1var",var_stack);
				imp2.show();
				ImagePlus imp3=new ImagePlus(trunctitle+"_avg",avg_stack);
				imp3.show();
			}
		}
		set_options();

	}
	
	void init_options(){
		String dir=System.getProperty("user.home");
		try{
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"stack_statistics_jru_v1.jrn");
			BufferedReader d=new BufferedReader(new FileReader(b));
			avgframes=Integer.parseInt(d.readLine());
			avgall = (Integer.parseInt(d.readLine())==0) ? false : true;
			index=Integer.parseInt(d.readLine());
			batchmode = (Integer.parseInt(d.readLine())==0) ? false : true;
			incthird = (Integer.parseInt(d.readLine())==0) ? false : true;
			incfourth = (Integer.parseInt(d.readLine())==0) ? false : true;
			d.close();
		}
		catch(IOException e){
			index=0;
			avgframes=50;
			avgall=false;
			index=0;
			batchmode=false;
			incthird=false;
			incfourth=false;
			set_options();
		}
		return;
	}
	
	void set_options(){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){a.mkdir();}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"stack_statistics_jru_v1.jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			d.write(""+avgframes+"\n");
			d.write(""+(avgall ? 1:0)+"\n");
			d.write(""+index+"\n");
			d.write(""+(batchmode ? 1:0)+"\n");
			d.write(""+(incthird ? 1:0)+"\n");
			d.write(""+(incfourth ? 1:0)+"\n");
			d.close();
		}
		catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}

}
