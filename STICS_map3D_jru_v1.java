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
import jalgs.jseg.*;
import jguis.*;

public class STICS_map3D_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		int subsize=32;
		gd.addNumericField("Subregion_Size_(pixels)?",subsize,0);
		int stepsize=32;
		gd.addNumericField("Step_Size?",stepsize,0);
		int subsizez=8;
		gd.addNumericField("Z_Subregion_Size_(pixels)?",subsizez,0);
		int stepsizez=8;
		gd.addNumericField("Z_Step_Size?",stepsizez,0);
		int shift=3;
		gd.addNumericField("STICS_temporal_Shift?",shift,0);
		float xoffset=0.0f;
		gd.addNumericField("X_Offset",xoffset,5,15,null);
		float yoffset=0.0f;
		gd.addNumericField("Y_Offset",yoffset,5,15,null);
		float zoffset=0.0f;
		gd.addNumericField("Z_Offset",zoffset,5,15,null);
		float multiplier=8.0f;
		gd.addNumericField("Velocity_Multiplier",multiplier,5,15,null);
		boolean norm=false;
		gd.addCheckbox("Normalize_Vector_lengths?",norm);
		boolean centered=true;
		gd.addCheckbox("Center_Vectors?",centered);
		float magthresh=0.5f;
		gd.addNumericField("Magnitude_Threshhold?",magthresh,5,15,null);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		subsize=(int)gd.getNextNumber();
		stepsize=(int)gd.getNextNumber();
		subsizez=(int)gd.getNextNumber();
		stepsizez=(int)gd.getNextNumber();
		shift=(int)gd.getNextNumber();
		xoffset=(float)gd.getNextNumber();
		yoffset=(float)gd.getNextNumber();
		zoffset=(float)gd.getNextNumber();
		multiplier=(float)gd.getNextNumber();
		norm=gd.getNextBoolean();
		centered=gd.getNextBoolean();
		magthresh=(float)gd.getNextNumber();

		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int xregions=1+(int)(((float)width-(float)subsize)/(float)stepsize);
		int newwidth=xregions*subsize;
		int yregions=1+(int)(((float)height-(float)subsize)/(float)stepsize);
		int newheight=yregions*subsize;
		int zregions=1+(int)(((float)slices-(float)subsizez)/(float)stepsizez);
		int newslices=zregions*subsizez;
		int zratio=(int)((float)subsize/(float)subsizez);
		float[][][] results=new float[3][zregions][xregions*yregions];
		//ImageStack ccoutput=new ImageStack(newwidth,newheight);
		//for(int i=0;i<newslices;i++) ccoutput.addSlice("",new float[newwidth*newheight]);
		crosscorr3D cc3D=new crosscorr3D(subsize,subsize,subsizez);
		int counter=0;
		for(int n=0;n<zregions;n++){
			int zpos=n*stepsizez;
			for(int m=0;m<yregions;m++){
				int ypos=m*stepsize;
				for(int l=0;l<xregions;l++){
					int xpos=l*stepsize;
					for(int k=0;k<(frames-shift);k++){
						//copy the subarray for correlation
						Object[] subarray1=get_subregion(stack,xpos,ypos,zpos,subsize,subsizez,slices,k);
						Object[] subarray2=get_subregion(stack,xpos,ypos,zpos,subsize,subsizez,slices,k+shift);
						subarray1=cc3D.docrosscorr3D(subarray2,subarray1,true,true);
						//set_subregion(ccoutput,subarray1,l*subsize,m*subsize,n*subsizez,subsize,subsizez);
						//find the 3D maximum of the cross correlation
						int maxj=subsize/2; int maxi=subsize/2; int maxz=subsizez/2; float max=((float[])subarray1[subsizez/2])[subsize/2+subsize*subsize/2];
						for(int z=0;z<subsizez;z++){
							int scalez=z*zratio;
							for(int i=0;i<subsize;i++){
								for(int j=0;j<subsize;j++){
									float dist=(float)Math.sqrt((i-subsize/2)*(i-subsize/2)+(j-subsize/2)*(j-subsize/2)+(scalez-subsize/2)*(scalez-subsize/2));
									if(dist<((float)subsize/4.0f)){
										if(((float[])subarray1[z])[j+i*subsize]>max){
											max=((float[])subarray1[z])[j+i*subsize]; maxj=j; maxi=i; maxz=z;
										}
									}
								}
							}
						}
						if(maxi<=0){maxi=1;}
						if(maxi>=(subsize-1)){maxi=subsize-2;}
						if(maxj<=0){maxj=1;}
						if(maxj>=(subsize-1)){maxj=subsize-2;}
						if(maxz<=0){maxz=1;}
						if(maxz>=(subsizez-1)){maxz=subsizez-2;}
						float[] tempmax=interpolation.get_local_max3D(subarray1,maxj,maxi,maxz,subsize,subsize);
						float xsum=tempmax[0]; float ysum=tempmax[1]; float zsum=tempmax[2];
						xsum+=xoffset-0.5f*(float)subsize;
						ysum+=yoffset-0.5f*(float)subsize;
						zsum+=zoffset-0.5f*(float)subsizez;
						results[0][n][l+m*xregions]=xsum;
						results[1][n][l+m*xregions]=ysum;
						results[2][n][l+m*xregions]=zsum;
					}
					counter++;
					IJ.showProgress(counter,xregions*yregions*zregions);
				}
			}
		}
		ImageStack retstack=new ImageStack(xregions,yregions);
		for(int i=0;i<zregions;i++){
			for(int j=0;j<3;j++){
				retstack.addSlice("",results[j][i]);
			}
		}
		ImagePlus imp2=new ImagePlus("Velocities",retstack);
		imp2.setOpenAsHyperStack(true);
		imp2.setDimensions(3,zregions,1);
		new CompositeImage(imp2,CompositeImage.COLOR).show();
		//new ImagePlus("STICS",ccoutput).show();
	}

	public Object[] get_subregion(ImageStack stack,int xpos,int ypos,int zpos,int subsize,int subsizez,int slices,int frame){
		int width=stack.getWidth();
		Object[] subarray=new Object[subsizez];
		for(int z=zpos;z<(zpos+subsizez);z++){
			float[] temp1=new float[subsize*subsize];
			if(stack.getPixels(1) instanceof float[]){
				for(int i=ypos;i<(ypos+subsize);i++){
					for(int j=xpos;j<(xpos+subsize);j++){
						temp1[(i-ypos)*subsize+(j-xpos)]=((float[])stack.getPixels(frame*slices+z+1))[i*width+j];
					}
				}
			} else {
				if(stack.getPixels(1) instanceof short[]){
					for(int i=ypos;i<(ypos+subsize);i++){
						for(int j=xpos;j<(xpos+subsize);j++){
							temp1[(i-ypos)*subsize+(j-xpos)]=((short[])stack.getPixels(frame*slices+z+1))[i*width+j]&0xffff;
						}
					}
				} else {
					for(int i=ypos;i<(ypos+subsize);i++){
						for(int j=xpos;j<(xpos+subsize);j++){
							temp1[(i-ypos)*subsize+(j-xpos)]=((byte[])stack.getPixels(frame*slices+z+1))[i*width+j]&0xff;
						}
					}
				}
			}
			subarray[z-zpos]=temp1;
		}
		return subarray;
	}

	public void set_subregion(ImageStack stack,Object[] subarray,int xpos,int ypos,int zpos,int subsize,int subsizez){
		int width=stack.getWidth();
		for(int z=zpos;z<(zpos+subsizez);z++){
			float[] temp1=(float[])subarray[z-zpos];
			float[] temp2=(float[])stack.getPixels(z+1);
			for(int i=ypos;i<(ypos+subsize);i++){
				for(int j=xpos;j<(xpos+subsize);j++){
					temp2[i*width+j]=temp1[(i-ypos)*subsize+(j-xpos)];
				}
			}
		}
	}

}
