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
import java.io.*;
import ij.plugin.*;
import ij.measure.*;

public class bin_image_jru_v1 implements PlugIn {
	int index,binby,binx,biny,binc,binz;
	boolean cumulative;

	public void run(String arg) {
		init_options();
		GenericDialog gd1=new GenericDialog("Binning Options");
		String[] bintypes={"Temporal_Bin","Spatial Bin","Z Bin","Channel Bin"};
		gd1.addChoice("Bin_Type",bintypes,bintypes[index]);
		gd1.addCheckbox("Cumulative?",cumulative);
		gd1.showDialog(); if(gd1.wasCanceled()){return;}
		index=gd1.getNextChoiceIndex();
		cumulative=gd1.getNextBoolean();
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageStack is=imp.getStack();
		int width=imp.getWidth();
		int height=imp.getHeight();
		int slices=imp.getNSlices();
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		ImageStack result_stack=new ImageStack(width,height);
		boolean isfloat=(is.getPixels(1) instanceof float[]);
		boolean isbyte=(is.getPixels(1) instanceof byte[]);
		Calibration cal=imp.getCalibration().copy();
		ImagePlus imp2=null;
		GenericDialog gd=null;
		switch(index){
			case 0:
				if(frames==1){IJ.showMessage("#frames = 1; try z bin"); break;}
				gd = new GenericDialog("Options");
				gd.addNumericField("bin_by_t?",binby,0);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				binby = (int)gd.getNextNumber();
				if(binby>frames) binby=frames;
				int binframes=(int)(frames/binby);
				cal.frameInterval*=(double)binby;
				for(int i=0;i<binframes;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							float[] result=new float[width*height];
							for(int l=0;l<binby;l++){
								Object frame=get3DSlice(is,i*binby+l,j,k,frames,slices,channels);
								for(int m=0;m<width*height;m++){
									if(isfloat){result[m]+=((float[])frame)[m];}
									else if(isbyte){result[m]+=(float)(((byte[])frame)[m]&0xff);}
									else{result[m]+=(float)(((short[])frame)[m]&0xffff);}
								}
							}
							if(!cumulative){
								for(int m=0;m<width*height;m++){
									result[m]/=(float)binby;
								}
							}
							result_stack.addSlice("",result);
						}
					}
				}
				imp2 = new ImagePlus("Binned Stack",result_stack);
				imp2.setOpenAsHyperStack(true);
				imp2.setDimensions(channels,slices,binframes);
				imp2.setCalibration(cal);
				if(channels>1 && imp.isComposite()){
					CompositeImage ci=new CompositeImage(imp2,((CompositeImage)imp).getMode());
					LUT[] lut=((CompositeImage)imp).getLuts();
					ci.setLuts(lut);
					ci.resetDisplayRanges();
					ci.show();
				} else {
					imp2.show();
				}
				break;
			case 1:
				gd = new GenericDialog("Options");
				gd.addNumericField("bin_by_x?",binx,0);
				gd.addNumericField("bin_by_y?",biny,0);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				binx = (int)gd.getNextNumber();
				biny = (int)gd.getNextNumber();
				int binwidth=(int)(width/binx);
				int binheight=(int)(height/biny);
				cal.pixelWidth*=(double)binx;
				cal.pixelHeight*=(double)biny;
				result_stack=new ImageStack(binwidth,binheight);
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							Object frame=get3DSlice(is,i,j,k,frames,slices,channels);
							float[] result=spatial_bin(frame,width,height,binx,biny);
							result_stack.addSlice("",result);
						}
					}
				}
				imp2 = new ImagePlus("Binned Stack",result_stack);
				imp2.setCalibration(cal);
				imp2.setOpenAsHyperStack(true);
				imp2.setDimensions(channels,slices,frames);
				if(channels>1 && imp.isComposite()){
					CompositeImage ci=new CompositeImage(imp2,((CompositeImage)imp).getMode());
					LUT[] lut=((CompositeImage)imp).getLuts();
					ci.setLuts(lut);
					ci.resetDisplayRanges();
					ci.show();
				} else {
					imp2.show();
				}
				break;
			case 2:
				if(slices==1){IJ.showMessage("#slices = 1; try temporal bin"); break;}
				gd = new GenericDialog("Options");
				gd.addNumericField("bin_by_z?",binz,0);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				binz = (int)gd.getNextNumber();
				if(binz>slices) binz=slices;
				int binslices=(int)(slices/binz);
				cal.pixelDepth*=(double)binz;
				for(int i=0;i<frames;i++){
					for(int j=0;j<binslices;j++){
						for(int k=0;k<channels;k++){
							float[] result=new float[width*height];
							for(int l=0;l<binz;l++){
								Object frame=get3DSlice(is,i,l+binz*j,k,frames,slices,channels);
								for(int m=0;m<width*height;m++){
									if(isfloat){result[m]+=((float[])frame)[m];}
									else if(isbyte){result[m]+=(float)(((byte[])frame)[m]&0xff);}
									else{result[m]+=(float)(((short[])frame)[m]&0xffff);}
								}
							}
							if(!cumulative){
								for(int m=0;m<width*height;m++){
									result[m]/=(float)binz;
								}
							}
							result_stack.addSlice("",result);
						}
					}
				}
				imp2 = new ImagePlus("Binned Stack",result_stack);
				imp2.setOpenAsHyperStack(true);
				imp2.setDimensions(channels,binslices,frames);
				imp2.setCalibration(cal);
				if(channels>1 && imp.isComposite()){
					CompositeImage ci=new CompositeImage(imp2,((CompositeImage)imp).getMode());
					LUT[] lut=((CompositeImage)imp).getLuts();
					ci.setLuts(lut);
					ci.resetDisplayRanges();
					ci.show();
				} else {
					imp2.show();
				}
				break;
			case 3:
				gd = new GenericDialog("Options");
				gd.addNumericField("bin_by_channel?",binc,0);
				gd.showDialog(); if(gd.wasCanceled()){return;}
				binc = (int)gd.getNextNumber();
				if(binc>channels) binc=channels;
				int binchannels=(int)(channels/binc);
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<binchannels;k++){
							float[] result=new float[width*height];
							for(int l=0;l<binc;l++){
								Object frame=get3DSlice(is,i,j,k*binc+l,frames,slices,channels);
								for(int m=0;m<width*height;m++){
									if(isfloat){result[m]+=((float[])frame)[m];}
									else if(isbyte){result[m]+=(float)(((byte[])frame)[m]&0xff);}
									else{result[m]+=(float)(((short[])frame)[m]&0xffff);}
								}
							}
							if(!cumulative){
								for(int m=0;m<width*height;m++){
									result[m]/=(float)binc;
								}
							}
							result_stack.addSlice("",result);
						}
					}
				}
				imp2 = new ImagePlus("Binned Stack",result_stack);
				imp2.setOpenAsHyperStack(true);
				imp2.setDimensions(binchannels,slices,frames);
				imp2.setCalibration(cal);
				if(binchannels>1){
					new CompositeImage(imp2,((CompositeImage)imp).getMode()).show();
				} else {
					imp2.show();
				}
				break;
			default:
				break;
		}
		set_options();
	}

	Object get3DSlice(ImageStack is,int frame,int slice,int channel,int frames,int slices,int channels){
		return is.getPixels(1+channel+slice*channels+frame*channels*slices);
	}

	public float[] spatial_bin(Object image,int width,int height,int binx,int biny){
		if(image instanceof float[]) return spatial_bin((float[])image,width,height,binx,biny);
		if(image instanceof short[]) return spatial_bin((short[])image,width,height,binx,biny);
		if(image instanceof byte[]) return spatial_bin((byte[])image,width,height,binx,biny);
		return null;
	}

	float[] spatial_bin(float[] image,int width,int height,int binx,int biny){
		int binwidth=(int)(width/binx);
		int binheight=(int)(height/biny);
		float[] result=new float[binwidth*binheight];
		for(int i=0;i<binheight;i++){
			int ystart=i*biny;
			for(int j=0;j<binwidth;j++){
				int xstart=j*binx;
				for(int k=0;k<biny;k++){
					for(int l=0;l<binx;l++){
						result[i*binwidth+j]+=image[(ystart+k)*width+(xstart+l)];
					}
				}
			}
		}
		if(!cumulative){
			for(int i=0;i<binwidth*binheight;i++){
				result[i]/=(float)(binx*biny);
			}
		}
		return result;
	}

	float[] spatial_bin(short[] image,int width,int height,int binx,int biny){
		int binwidth=(int)(width/binx);
		int binheight=(int)(height/biny);
		float[] result=new float[binwidth*binheight];
		for(int i=0;i<binheight;i++){
			int ystart=i*biny;
			for(int j=0;j<binwidth;j++){
				int xstart=j*binx;
				for(int k=0;k<biny;k++){
					for(int l=0;l<binx;l++){
						result[i*binwidth+j]+=(float)(image[(ystart+k)*width+(xstart+l)]&0xffff);
					}
				}
			}
		}
		if(!cumulative){
			for(int i=0;i<binwidth*binheight;i++){
				result[i]/=(float)(binx*biny);
			}
		}
		return result;
	}

	float[] spatial_bin(byte[] image,int width,int height,int binx,int biny){
		int binwidth=(int)(width/binx);
		int binheight=(int)(height/biny);
		float[] result=new float[binwidth*binheight];
		for(int i=0;i<binheight;i++){
			int ystart=i*biny;
			for(int j=0;j<binwidth;j++){
				int xstart=j*binx;
				for(int k=0;k<biny;k++){
					for(int l=0;l<binx;l++){
						result[i*binwidth+j]+=(float)(image[(ystart+k)*width+(xstart+l)]&0xff);
					}
				}
			}
		}
		if(!cumulative){
			for(int i=0;i<binwidth*binheight;i++){
				result[i]/=(float)(binx*biny);
			}
		}
		return result;
	}
	
	void init_options(){
		String dir=System.getProperty("user.home");
		try{
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"bin_image_jru_v1.jrn");
			BufferedReader d=new BufferedReader(new FileReader(b));
			index=Integer.parseInt(d.readLine());
			binby=Integer.parseInt(d.readLine());
			binx=Integer.parseInt(d.readLine());
			biny=Integer.parseInt(d.readLine());
			binc=Integer.parseInt(d.readLine());
			binz=Integer.parseInt(d.readLine());
			if((Integer.parseInt(d.readLine()))==0){cumulative=false;}
			else{cumulative=true;}
			d.close();
		}
		catch(IOException e){
			index=0;
			binby=2;
			binx=2;
			biny=2;
			binc=2;
			binz=2;
			set_options();
		}
		catch(NumberFormatException e){
			index=0;
			binby=2;
			binx=2;
			biny=2;
			binc=2;
			binz=2;
			set_options();
		}
		return;
	}
	
	void set_options(){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){a.mkdir();}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"bin_image_jru_v1.jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			d.write(""+index+"\n");
			d.write(""+binby+"\n");
			d.write(""+binx+"\n");
			d.write(""+biny+"\n");
			d.write(""+binc+"\n");
			d.write(""+binz+"\n");
			if(!cumulative){d.write("0\n");}
			else{d.write("1\n");}
			d.close();
		}
		catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}

}
