/*******************************************************************************
 * Copyright (c) 2019 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import ij.plugin.frame.RoiManager;
import jalgs.*;
import jalgs.jseg.*;
import jalgs.jfit.*;

public class manual_transform_3D_images_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin takes two 3D images and a set of alternating Rois in the roi manager
		//it creates the best scaled rotation transform from the second image Roi points to the first and transforms the image
		//these steps are individually described by coord transformation jru v1 and custom image transformation jru v1
		//rois must be ordered dest, source, dest, source, ...
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Reference Image (dest)","Image To Transform (src)"});
		if(imps==null || imps.length<2) return;
		int destwidth=imps[0].getWidth(); int destheight=imps[0].getHeight(); int destslices=imps[0].getNSlices();
		int srcwidth=imps[1].getWidth(); int srcheight=imps[1].getHeight();
		int srcframes=imps[1].getNFrames();
		int srcslices=imps[1].getNSlices();
		int srcchans=imps[1].getNChannels();
		ImageStack srcstack=imps[1].getStack();
		float destzratio=(float)jutils.get_zratio(imps[0]);
		float srczratio=(float)jutils.get_zratio(imps[1]);
		//start by getting the rois
		RoiManager rman=RoiManager.getInstance();
		Roi[] rois=rman.getRoisAsArray();
		float[][] centroids1=new float[3][rois.length/2];
		float[][] centroids2=new float[3][rois.length/2];
		for(int i=0;i<rois.length/2;i++){
			float[] temp=getRoiCentroid(rois[2*i]);
			centroids1[0][i]=temp[0]; centroids1[1][i]=temp[1];
			centroids1[2][i]=rois[2*i].getZPosition()-1; centroids1[2][i]*=destzratio;
			temp=getRoiCentroid(rois[2*i+1]);
			centroids2[0][i]=temp[0]; centroids2[1][i]=temp[1];
			centroids2[2][i]=rois[2*i+1].getZPosition()-1; centroids2[2][i]*=srczratio;
		}
		//now get the transformation (order?)
		Object[] trans=(new jreg()).scaled_rotation_fiducials(centroids1,centroids2);
		float[] t=(float[])trans[1];
		float[][] r=(float[][])trans[0];
		float s=((Float)trans[2]).floatValue();
		table_tools.create_table("Rotation Matrix",r,null);
		table_tools.create_table("Translation Matrix",new float[][]{t},new String[]{"x","y","z"});
		table_tools.create_table("Centroids",new float[][]{(float[])trans[3],(float[])trans[4]},new String[]{"x","y","z"});
		IJ.log("scaling = "+s);
		//now transform the alignment coordinates
		//first add scaling
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				r[i][j]*=s;
			}
		}
		//shift the coordinates to the centroid
		shift_points(centroids2,(float[])trans[4],-1.0f);
		for(int i=0;i<centroids2[0].length;i++){
			float[][] rt=(new jreg()).mattrans(r);
			float[] temp=matmult(new float[]{centroids2[0][i],centroids2[1][i],centroids2[2][i]},rt);
			centroids2[0][i]=temp[0]; centroids2[1][i]=temp[1]; centroids2[2][i]=temp[2];
		}
		//shift to the reference centroid
		shift_points(centroids2,(float[])trans[3],1.0f);
		//table_tools.create_table("Transformed Coords",(new jreg()).mattrans(centroids2),new String[]{"x","y","z"});
		Traj3D traj3D=new Traj3D("x","y","z",centroids2[0],centroids2[1],centroids2[2]);
		traj3D.addPoints(centroids1[0],centroids1[1],centroids1[2],true);
		new PlotWindow3D("Transformed vs Ref Coordinates",traj3D).draw();
		float[] cdest=(float[])trans[3];
		float[] csource=(float[])trans[4];
		//now perform the transformation
		float[][] transformed=new float[destslices*srcframes*srcchans][];
		for(int l=0;l<srcframes;l++){
			for(int m=0;m<srcchans;m++){
				IJ.showStatus("transforming frame "+(l+1)+", chan "+(m+1));
				float[][] t3D=new float[destslices][destheight*destwidth];
				Object[] srcpix=jutils.get3DZSeries(srcstack,m,l,srcframes,srcslices,srcchans);
				for(int i=0;i<destslices;i++){
					for(int j=0;j<destheight;j++){
						for(int k=0;k<destwidth;k++){
							float[][] coords={{(float)k},{(float)j},{destzratio*(float)i}};
							shift_points(coords,csource,-1.0f);
							float[] temp=matmult(new float[]{coords[0][0],coords[1][0],coords[2][0]},r);
							coords[0][0]=temp[0]; coords[1][0]=temp[1]; coords[2][0]=temp[2];
							shift_points(coords,cdest,1.0f);
							coords[2][0]/=srczratio;
							t3D[i][k+j*destwidth]=interpolation.interp3D(srcpix,srcwidth,srcheight,coords[0][0],coords[1][0],coords[2][0]);
						}
					}
					IJ.showProgress(i,destslices);
				}
				for(int i=0;i<destslices;i++) transformed[l*destslices*srcchans+i*srcchans+m]=t3D[i];
			}
		}
		ImagePlus imp=new ImagePlus("transformed",jutils.array2stack(transformed,destwidth,destheight));
		imp.setOpenAsHyperStack(true);
		imp.setDimensions(srcchans,destslices,srcframes);
		new CompositeImage(imp,CompositeImage.COLOR).show();
	}

	public float[] getRoiCentroid(Roi roi){
		if(roi instanceof PointRoi){
			Rectangle r=roi.getBounds();
			float x=(float)r.x;
			float y=(float)r.y;
			return new float[]{x,y};
		} else {
			float[] centroid=measure_object.centroid(roi.getPolygon());
			return centroid;
		}
	}

	public void shift_points(float[][] points,float[] centroid,float sign){
		for(int i=0;i<points[0].length;i++){
			for(int j=0;j<3;j++) points[j][i]+=sign*centroid[j];
		}
	}

	public float[] matmult(float[] vec,float[][] mat){
		float[] newvec=new float[vec.length];
		for(int i=0;i<vec.length;i++){
			for(int j=0;j<vec.length;j++){
				newvec[i]+=vec[j]*mat[i][j];
			}
		}
		return newvec;
	}

}
