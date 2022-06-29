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
import jguis.*;
import ij.measure.*;

public class avg_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd= new GenericDialog("Options");
		boolean vert=true;
		gd.addCheckbox("Vertical?",vert);
		String[] proftype={"T profile","Z profile","C profile"};
		gd.addChoice("Profile_Type",proftype,proftype[0]);
		gd.addChoice("Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		vert=gd.getNextBoolean();
		int profindex=gd.getNextChoiceIndex();
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		ImagePlus imp = WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int frames=imp.getNFrames();
		int slices=imp.getNSlices();
		int channels=imp.getNChannels();
		int thissize=frames; if(profindex==1){thissize=slices;} if(profindex==2){thissize=channels;}
		Rectangle r=null;
		Roi roi=imp.getRoi();
		if(roi!=null){r=roi.getBounds();}
		else{r=new Rectangle(0,0,width,height);}
		int length=r.width;
		if(vert) length=r.height;
		float[][][][] data=null;
		float[][] xvals=new float[channels][];
		for(int i=0;i<channels;i++){
			xvals[i]=get_xvals(r,imp.getCalibration(),vert);
		}
		if(profindex==0){
			data=new float[1][slices][channels][length*thissize];
			for(int i=0;i<frames;i++){
				for(int j=0;j<slices;j++){
					for(int k=0;k<channels;k++){
						float[] temp=get_profile(stack.getPixels(i*channels*slices+j*channels+k+1),width,height,vert,r,stat);
						System.arraycopy(temp,0,data[0][j][k],i*length,length);
					}
				}
			}
		} else{
			if(profindex==1){
				data=new float[frames][1][channels][length*thissize];
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							float[] temp=get_profile(stack.getPixels(i*channels*slices+j*channels+k+1),width,height,vert,r,stat);
							System.arraycopy(temp,0,data[i][0][k],j*length,length);
						}
					}
				}
			} else {
				data=new float[frames][slices][1][length*thissize];
				for(int i=0;i<frames;i++){
					for(int j=0;j<slices;j++){
						for(int k=0;k<channels;k++){
							float[] temp=get_profile(stack.getPixels(i*channels*slices+j*channels+k+1),width,height,vert,r,stat);
							System.arraycopy(temp,0,data[i][j][0],k*length,length);
						}
					}
				}
			}
		}
		if(frames==1 && slices==1){
			new PlotWindow4(stat+" Profile","distance","intensity",xvals,data[0][0],null).draw();
		} else {
			ImageStack retstack=new ImageStack(length,thissize);
			for(int i=0;i<data.length;i++){
				for(int j=0;j<data[0].length;j++){
					for(int k=0;k<data[0][0].length;k++){
						retstack.addSlice("",data[i][j][k]);
					}
				}
			}
			ImagePlus tempimp=new ImagePlus(stat+" Profiles",retstack);
			tempimp.copyScale(imp);
			tempimp.setOpenAsHyperStack(true);
			tempimp.setDimensions(data[0][0].length,data[0].length,data.length);
			if(data[0][0].length>1 && imp.isComposite()){
				CompositeImage ci=new CompositeImage(tempimp,((CompositeImage)imp).getMode());
				ci.copyLuts(imp);
				ci.show();
			} else {
				tempimp.show();
			}
		}
	}

	private float[] get_xvals(Rectangle r,Calibration cal,boolean vert){
		if(vert){
			float[] temp=new float[r.height];
			for(int i=0;i<r.height;i++){
				temp[i]=(float)cal.pixelHeight*(float)i;
			}
			return temp;
		} else {
			float[] temp=new float[r.width];
			for(int i=0;i<r.width;i++){
				temp[i]=(float)cal.pixelWidth*(float)i;
			}
			return temp;
		}
	}

	private float[] get_profile(Object pixels,int width,int height,boolean vert,Rectangle r,String stat){
		if(pixels instanceof float[]){
			return get_profile((float[])pixels,width,height,vert,r,stat);
		} else {
			if(pixels instanceof short[]){
				return get_profile((short[])pixels,width,height,vert,r,stat);
			} else {
				return get_profile((byte[])pixels,width,height,vert,r,stat);
			}
		}
	}

	private float[] get_profile(short[] pixels,int width,int height,boolean vert,Rectangle r,String stat){
		float[] data=null;
		if(vert){
			data=new float[r.height];
			for(int i=r.y;i<(r.y+r.height);i++){
				short[] temp=new short[r.width];
				for(int j=r.x;j<(r.x+r.width);j++){
					temp[j-r.x]=pixels[i*width+j];
				}
				data[i-r.y]=jstatistics.getstatistic(stat,temp,null);
			}
		}
		else {
			data=new float[r.width];
			for(int i=r.x;i<(r.x+r.width);i++){
				short[] temp=new short[r.height];
				for(int j=r.y;j<(r.y+r.height);j++){
					temp[j-r.y]=pixels[j*width+i];
				}
				data[i-r.x]=jstatistics.getstatistic(stat,temp,null);
			}
		}
		return data;
	}

	private float[] get_profile(byte[] pixels,int width,int height,boolean vert,Rectangle r,String stat){
		float[] data=null;
		jstatistics jstats=new jstatistics();
		if(vert){
			data=new float[r.height];
			for(int i=r.y;i<(r.y+r.height);i++){
				byte[] temp=new byte[r.width];
				for(int j=r.x;j<(r.x+r.width);j++){
					temp[j-r.x]=pixels[i*width+j];
				}
				data[i-r.y]=jstatistics.getstatistic(stat,temp,null);
			}
		}
		else {
			data=new float[r.width];
			for(int i=r.x;i<(r.x+r.width);i++){
				byte[] temp=new byte[r.height];
				for(int j=r.y;j<(r.y+r.height);j++){
					temp[j-r.y]=pixels[j*width+i];
				}
				data[i-r.x]=jstatistics.getstatistic(stat,temp,null);
			}
		}
		return data;
	}

	private float[] get_profile(float[] pixels,int width,int height,boolean vert,Rectangle r,String stat){
		float[] data=null;
		jstatistics jstats=new jstatistics();
		if(vert){
			data=new float[r.height];
			for(int i=r.y;i<(r.y+r.height);i++){
				float[] temp=new float[r.width];
				for(int j=r.x;j<(r.x+r.width);j++){
					temp[j-r.x]=pixels[i*width+j];
				}
				data[i-r.y]=jstatistics.getstatistic(stat,temp,null);
			}
		} else {
			data=new float[r.width];
			for(int i=r.x;i<(r.x+r.width);i++){
				float[] temp=new float[r.height];
				for(int j=r.y;j<(r.y+r.height);j++){
					temp[j-r.y]=pixels[j*width+i];
				}
				data[i-r.x]=jstatistics.getstatistic(stat,temp,null);
			}
		}
		return data;
	}

}
