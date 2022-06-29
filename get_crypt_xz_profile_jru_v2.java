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

public class get_crypt_xz_profile_jru_v2 implements PlugIn {
	//this version uses three points per crypt instead of two

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		int width=imp.getWidth(); int height=imp.getHeight();
		float zratio=(float)jutils.get_zratio(imp);
		float zoff=50f;
		int zsize=400;
		float expansion=1.5f;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addNumericField("Z_Offset (pixel units)",zoff,5,15,null);
		gd.addNumericField("Z_Size (pixel units)",zsize,0);
		gd.addNumericField("XY_Size (fraction of diameter)",expansion,5,15,null);
		gd.addNumericField("N_Crypt_Points (including vertex)",3,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		zratio=(float)gd.getNextNumber();
		zoff=(float)gd.getNextNumber();
		zsize=(int)gd.getNextNumber();
		expansion=(float)gd.getNextNumber();
		int npts=(int)gd.getNextNumber();
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
		//each npts set of items in the roi manager will be ellipses moving down the crypt neck and then a point at the crypt vertex
		//need to transform each crypt (segment by segment so it is oriented vertically
		//then crop it (use max ellipse diameter) and create the radial profile in all channels
		Roi[] rois=rman.getRoisAsArray();
		int ncrypts=rois.length/npts;
		for(int i=0;i<ncrypts;i++){
			IJ.log("analyzing crypt "+(i+1));
			//for each crypt get the centroid of the ellipse and the vertex and the longest dimension of the ellipse
			Roi vertexroi=rois[npts*i+(npts-1)];
			Rectangle r=vertexroi.getBounds();
			Roi[] neckrois=new Roi[npts-1];
			for(int j=0;j<(npts-1);j++) neckrois[j]=rois[npts*i+j];
			//Roi neckroi=rois[npts*i];
			int vertexz=vertexroi.getZPosition(); if(vertexz==0) vertexz=vertexroi.getPosition();
			float[] vertex={r.x,r.y,zratio*(float)(vertexz-1)};
			IJ.log("vertex pos = \t"+table_tools.print_float_array(vertex));
			//double[] params=((EllipseRoi)neckroi).getParams();
			//now get the neck centroids
			float[][] neck=new float[npts][];
			float maxd=0; //the max neck diameter
			for(int j=0;j<(npts-1);j++){
				double[] params=((EllipseRoi)neckrois[j]).getParams();
				int neckz=neckrois[j].getZPosition(); if(neckz==0) neckz=neckrois[j].getPosition();
				neck[j]=new float[]{0.5f*(float)(params[0]+params[2]),0.5f*(float)(params[1]+params[3]),zratio*(float)(neckz-1)};
				IJ.log("neck center "+(j+1)+" = \t"+table_tools.print_float_array(neck[j]));
				float tempd=(float)Math.sqrt((params[2]-params[0])*(params[2]-params[0])+(params[3]-params[1])*(params[3]-params[1]));
				if(tempd>maxd) maxd=tempd;
			}
			neck[npts-1]=vertex;

			//for the transformation, need to rotate about the cross-product between the crypt vector and the z axis
			//rotation angle is given by the dot product between the crypt vector and the z axis
			Object[] xzstacks=new Object[npts-1];
			int[] sizes=new int[npts-1];
			int rotsize=(int)(expansion*maxd);
			int rsize=(int)(0.5f*(float)rotsize-0.5f);
			int newwidth=rsize*2-1;
			int sizeleft=zsize;
			//realign all segments but the first and the last
			for(int j=1;j<(npts-2);j++){
			//for(int j=0;j<(npts-1);j++){
				float[] cryptvec={(neck[j][0]-neck[j+1][0]),(neck[j][1]-neck[j+1][1]),(neck[j][2]-neck[j+1][2])};
				float tempdist=(float)Math.sqrt(cryptvec[0]*cryptvec[0]+cryptvec[1]*cryptvec[1]+cryptvec[2]*cryptvec[2]);
				cryptvec=measure_object.norm_vector(cryptvec);
				float[][] xzstack=new float[nchans][];
				for(int k=0;k<nchans;k++){
					float[][] rotated=profiler.getRotated3DImage(stack2[k],width,height,neck[j+1],zratio,cryptvec,(int)(expansion*maxd),0.0f,(int)tempdist);
					//new ImagePlus("crypt "+i+" rotated "+j,jutils.array2stack(rotated,rotsize,rotsize)).show();
					xzstack[k]=getxzProfile(rotated,rotsize,rotsize,rsize);
				}
				sizeleft-=(int)tempdist;
				sizes[j]=(int)tempdist;
				xzstacks[j]=algutils.clone_multidim_array(xzstack);
				//new ImagePlus("crypt "+i+" profile "+j,jutils.array2stack((float[][])xzstacks[j],newwidth,(int)tempdist)).show();
			}

			//now realign the last segment along with the z offset
			float[] cryptvec={(neck[npts-2][0]-neck[npts-1][0]),(neck[npts-2][1]-neck[npts-1][1]),(neck[npts-2][2]-neck[npts-1][2])};
			float tempdist=(float)Math.sqrt(cryptvec[0]*cryptvec[0]+cryptvec[1]*cryptvec[1]+cryptvec[2]*cryptvec[2]);
			tempdist+=zoff;
			cryptvec=measure_object.norm_vector(cryptvec);
			float[][] xzstack=new float[nchans][];
			for(int k=0;k<nchans;k++){
				float[][] rotated=profiler.getRotated3DImage(stack2[k],width,height,neck[npts-1],zratio,cryptvec,(int)(expansion*maxd),zoff,(int)tempdist);
				//new ImagePlus("crypt "+i+" profile "+(npts-2),jutils.array2stack(rotated,rotsize,rotsize)).show();
				xzstack[k]=getxzProfile(rotated,rotsize,rotsize,rsize);
			}
			sizeleft-=(int)tempdist;
			sizes[npts-2]=(int)tempdist;
			xzstacks[npts-2]=algutils.clone_multidim_array(xzstack);
			//new ImagePlus("crypt "+i+" profile "+(npts-2),jutils.array2stack((float[][])xzstacks[npts-2],newwidth,(int)tempdist)).show();

			//now realign the first segment along with the leftover zdist
			cryptvec=new float[]{(neck[0][0]-neck[1][0]),(neck[0][1]-neck[1][1]),(neck[0][2]-neck[1][2])};
			cryptvec=measure_object.norm_vector(cryptvec);
			xzstack=new float[nchans][];
			for(int k=0;k<nchans;k++){
				float[][] rotated=profiler.getRotated3DImage(stack2[k],width,height,neck[1],zratio,cryptvec,(int)(expansion*maxd),0.0f,sizeleft);
				//new ImagePlus("crypt "+i+" profile 0",jutils.array2stack(rotated,rotsize,rotsize)).show();
				xzstack[k]=getxzProfile(rotated,rotsize,rotsize,rsize);
			}
			sizes[0]=sizeleft;
			xzstacks[0]=algutils.clone_multidim_array(xzstack);
			//new ImagePlus("crypt "+i+" profile 0",jutils.array2stack((float[][])xzstacks[0],newwidth,(int)tempdist)).show();

			//now we need to concatenate all of the profiles together
			int totsize=0;
			for(int j=0;j<sizes.length;j++){totsize+=sizes[j]; IJ.log(""+sizes[j]);}
			xzstack=new float[nchans][newwidth*totsize];
			int pos=0;
			for(int j=(npts-2);j>=0;j--){
				float[][] temp2=(float[][])xzstacks[j];
				for(int k=0;k<nchans;k++){
					System.arraycopy(temp2[k],0,xzstack[k],pos,sizes[j]*newwidth);
				}
				pos+=sizes[j]*newwidth;
			}
			jutils.create_hyperstack("crypt "+i+" xzprofile",jutils.array2stack(xzstack,newwidth,totsize),imp,1,1,nchans).show();
		}
	}

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
