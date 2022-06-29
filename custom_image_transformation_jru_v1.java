/*******************************************************************************
 * Copyright (c) 2016 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import jguis.*;
import jalgs.*;
import jalgs.jseg.*;
import jalgs.jfit.*;
import java.util.*;
import ij.text.*;

public class custom_image_transformation_jru_v1 implements PlugIn {

	public void run(String arg) {
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Reference Image (dest)","Image To Transform (src)"});
		if(imps==null || imps.length<2) return;
		TextWindow[] tw=jutils.selectTables(false,2,new String[]{"Rot_matrix (rotate dest into src)","Centroids (src then dest)"});
		if(tw==null || tw.length<2) return;
		int destwidth=imps[0].getWidth(); int destheight=imps[0].getHeight(); int destslices=imps[0].getNSlices();
		int srcwidth=imps[1].getWidth(); int srcheight=imps[1].getHeight();
		int srcframes=imps[1].getNFrames();
		int srcslices=imps[1].getNSlices();
		int srcchans=imps[1].getNChannels();
		ImageStack srcstack=imps[1].getStack();
		float destzratio=6.06f;
		float srczratio=6.06f;
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Dest_Z_Ratio",destzratio,5,15,null);
		gd.addNumericField("Src_Z_Ratio",srczratio,5,15,null);
		gd.addNumericField("Scaling (dest/src)",1.0,5,15,null);
		gd.addNumericField("Pixel_Size (1 if centroid is pixel units)",1.0,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		destzratio=(float)gd.getNextNumber();
		srczratio=(float)gd.getNextNumber();
		float s=(float)gd.getNextNumber();
		float psize=(float)gd.getNextNumber();
		List<List<String>> rlisttable=table_tools.table2listtable(tw[0].getTextPanel());
		List<List<String>> clisttable=table_tools.table2listtable(tw[1].getTextPanel());
		//double[][] rd={{0.999046,-0.02942476,0.03226281},{0.02566373,0.9934486,0.11135952},{-0.035328194,-0.11042526,0.9932563}};
		//IJ.log(table_tools.print_double_array(rd));
		//float[][] r=algutils.convert_arr_float(rd);
		float[][] rt=table_tools.get_matrix(rlisttable);
		for(int i=0;i<rt.length;i++){
			for(int j=0;j<rt[i].length;j++) rt[i][j]*=s;
		}
		//IJ.log(table_tools.print_float_array(rt));
		float[][] r=(new jreg()).mattrans(rt);
		//float[] csource={1306.101889f,583.3888889f,129.96633f};
		//float[] cdest={1589.6304f,1075.0339f,130.63972f};
		float[] cdest=table_tools.get_row_array(clisttable,0);
		//IJ.log(table_tools.print_float_array(cdest));
		float[] csource=table_tools.get_row_array(clisttable,1);
		for(int i=0;i<csource.length;i++) csource[i]/=psize;
		for(int i=0;i<cdest.length;i++) cdest[i]/=psize;
		//IJ.log(table_tools.print_float_array(csource));
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
							float[] temp=matmult(new float[]{coords[0][0],coords[1][0],coords[2][0]},rt);
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
