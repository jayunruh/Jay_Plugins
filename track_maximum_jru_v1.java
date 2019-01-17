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
import ij.plugin.frame.*;
import jguis.*;

public class track_maximum_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("interpolate_maximum?",false);
		String[] centeroptions={"use_upper_left","use_center","use_initial"};
		gd.addChoice("Traj_Origin",centeroptions,centeroptions[0]);
		gd.addNumericField("Edge_Buffer",10,0);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean interp=gd.getNextBoolean();
		int origindex=gd.getNextChoiceIndex();
		int edgebuf=(int)gd.getNextNumber();
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth();
		int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		Roi roi=imp.getRoi();
		Rectangle r=null;
		if(roi!=null) r=roi.getBounds();
		float[][] traj=new float[2][slices];
		for(int i=0;i<slices;i++){
			float[] image=(float[])stack.getPixels(i+1);
			float max=image[edgebuf+edgebuf*width];
			int maxj=edgebuf;
			int maxk=edgebuf;
			if(r==null){
				for(int j=0;j<height;j++){
					for(int k=0;k<width;k++){
						if(inbounds(width,height,edgebuf,k,j) && image[k+j*width]>max){
							max=image[k+j*width];
							maxk=k;
							maxj=j;
						}
					}
				}
			} else {
				max=image[r.x+r.y*width];
				maxj=r.x; maxk=r.y;
				for(int j=r.y;j<(r.y+r.height);j++){
					for(int k=r.x;k<(r.x+r.width);k++){
						if(inbounds(width,height,edgebuf,k,j) && image[k+j*width]>max){
							max=image[k+j*width];
							maxk=k;
							maxj=j;
						}
					}
				}
				r.translate(maxk-r.x-r.width/2,maxj-r.y-r.height/2);
			}
			if(interp){
				float[] interpmax=interp_max(image,width,height,maxk,maxj);
				traj[0][i]=interpmax[0];
				traj[1][i]=interpmax[1];
			} else {
				traj[0][i]=(float)maxk;
				traj[1][i]=(float)maxj;
			}
		}
		float subx=0.0f; float suby=0.0f;
		if(origindex==1){
			subx=(float)(width/2); suby=(float)(height/2);
		}
		if(origindex==2){
			subx=traj[0][0]; suby=traj[1][0];
		}
		if(slices>1){
			for(int i=0;i<slices;i++){
				traj[0][i]-=subx;
				traj[1][i]-=suby;
			}
			new PlotWindow4("Max Track","x","y",traj[0],traj[1]).draw();
		} else {
			IJ.log("maxx = "+(traj[0][0]-subx)+" , maxy = "+(traj[1][0]-suby));
		}
	}

	public boolean inbounds(int width,int height,int buf,int x,int y){
		if(x<buf) return false;
		if(y<buf) return false;
		if(x>=(width-buf)) return false;
		if(y>=(height-buf)) return false;
		return true;
	}

	public float[] interp_max(float[] image,int width,int height,int maxx,int maxy){
		if(maxx<=0){maxx=1;}
		if(maxx>=(width-1)){maxx=width-2;}
		if(maxy<=0){maxy=1;}
		if(maxy>=(height-1)){maxy=height-2;}
		float sum=0.0f; float xsum=0.0f; float ysum=0.0f;
		for(int j=(maxy-1);j<=(maxy+1);j++){
			for(int k=(maxx-1);k<=(maxx+1);k++){
				sum+=image[k+j*width];
				xsum+=((float)k)*image[k+j*width];
				ysum+=((float)j)*image[k+j*width];
			}
		}
		return new float[]{xsum/sum,ysum/sum};
	}

}
