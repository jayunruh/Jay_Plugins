/*******************************************************************************
 * Copyright (c) 2017 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;

public class get_crypt_xz_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		float zratio=(float)jutils.get_zratio(imp);
		float zoff=50f;
		int zsize=250;
		float expansion=1.5f;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addNumericField("Z_Offset (pixel units)",zoff,5,15,null);
		gd.addNumericField("Z_Size (pixel units)",zsize,0);
		gd.addNumericField("XY_Size (fraction of diameter)",expansion,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		zratio=(float)gd.getNextNumber();
		zoff=(float)gd.getNextNumber();
		zsize=(int)gd.getNextNumber();
		expansion=(float)gd.getNextNumber();
		ImageStack stack=imp.getStack();
		int nchans=imp.getNChannels();
		int nslices=imp.getNSlices();
		Object[][] stack2=new Object[nchans][nslices];
		Object[] temp=jutils.stack2array(stack);
		for(int i=0;i<nslices;i++){
			for(int j=0;j<nchans;j++){
				stack2[j][i]=temp[j+i*nchans];
			}
		}
		RoiManager rman=RoiManager.getInstance();
		if(rman==null) {
			IJ.error("need crypt neck ellipses and vertex points in roi manager");
			return;
		}
		//each pair of items in the roi manager will be an ellipse around the crypt neck and then a point at the crypt vertex
		//need to transform each crypt so it is oriented vertically
		//then crop it (use max ellipse diameter) and create the radial profile in all channels
		Roi[] rois=rman.getRoisAsArray();
		int ncrypts=rois.length/2;
		for(int i=0;i<ncrypts;i++){
			IJ.log("analyzing crypt "+(i+1));
			//for each crypt get the centroid of the ellipse and the vertex and the longest dimension of the ellipse
			Roi vertexroi=rois[2*i+1];
			Rectangle r=vertexroi.getBounds();
			Roi neckroi=rois[2*i];
			int vertexz=vertexroi.getZPosition(); if(vertexz==0) vertexz=vertexroi.getPosition();
			float[] vertex={r.x,r.y,zratio*(float)(vertexz-1)};
			IJ.log("vertex pos = \t"+table_tools.print_float_array(vertex));
			double[] params=((EllipseRoi)neckroi).getParams();
			int neckz=neckroi.getZPosition(); if(neckz==0) neckz=neckroi.getPosition();
			float[] neck={0.5f*(float)(params[0]+params[2]),0.5f*(float)(params[1]+params[3]),zratio*(float)(neckz-1)};
			IJ.log("neck center = \t"+table_tools.print_float_array(neck));
			float maxd=(float)Math.sqrt((params[2]-params[0])*(params[2]-params[0])+(params[3]-params[1])*(params[3]-params[1]));
			//for the transformation, need to rotate about the cross-product between the crypt vector and the z axis
			//rotation angle is given by the dot product between the crypt vector and the z axis
			float[] cryptvec={(neck[0]-vertex[0]),(neck[1]-vertex[1]),(neck[2]-vertex[2])};
			cryptvec=measure_object.norm_vector(cryptvec);
			float[] zvec={0.0f,0.0f,1.0f};
			float angle=measure_object.get_inner_angle(zvec,cryptvec); //order?
			IJ.log("angle = \t"+angle);
			//now get the cross product
			float[] crossprod=measure_object.crossProd(zvec,cryptvec);
			crossprod=measure_object.norm_vector(crossprod);
			//IJ.log("rot vector = \t"+table_tools.print_float_array(crossprod));
			float[][] rotmat=measure_object.getRotationMatrix(crossprod,angle);
			IJ.log("rot matrix = \t"+table_tools.print_float_array(rotmat));
			int rotsize=(int)(expansion*maxd);
			int rsize=(int)(0.5f*(float)rotsize-0.5f);
			float[][] xzstack=new float[nchans][];
			for(int j=0;j<nchans;j++){
				Object[] chanstack=stack2[j];
				float[][] rotated=profiler.getRotated3DImage(chanstack,width,height,vertex,zratio,rotmat,(int)(expansion*maxd),zoff,zsize);
				//new ImagePlus("rotated",jutils.array2stack(rotated,rotsize,rotsize)).show();
				xzstack[j]=getxzProfile(rotated,rotsize,rotsize,rsize);
			}
			jutils.create_hyperstack("crypt "+i+" xzprofile",jutils.array2stack(xzstack,rsize*2-1,zsize),imp,1,1,nchans).show();
		}
	}

	//here we do the circular average through an entire stack
	public float[] getxzProfile(float[][] stack,int width,int height,int rsize){
		int newsize=2*rsize-1;
		float[] xzprof=new float[newsize*stack.length];
		for(int i=0;i<stack.length;i++){
			float[] circavg=interpolation.circavg(stack[i],width,height,rsize,width/2,height/2,true);
			System.arraycopy(circavg,0,xzprof,i*newsize,newsize);
		}
		return xzprof;
	}

}
