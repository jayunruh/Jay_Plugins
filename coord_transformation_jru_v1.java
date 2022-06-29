/*******************************************************************************
 * Copyright (c) 2015 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.text.*;
import java.util.*;
import jguis.*;
import jalgs.jfit.*;
import jalgs.*;

public class coord_transformation_jru_v1 implements PlugIn {
	//this plugin aligns two 3 point data sets in 3D space
	//the aligned data set can have more coordinates which will be similarly aligned and scaled
	//the translation, rotation, and scaling are output as well as transformed dataset
	//the tables must have 3 columns: x, y, and z

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Affine_Transform?",false);
		gd.addCheckbox("Allow_Scaling",false);
		gd.addCheckbox("Plot_Transformed",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean affine=gd.getNextBoolean();
		boolean allows=gd.getNextBoolean();
		boolean plot=gd.getNextBoolean();
		TextWindow[] tw=jutils.selectTables(false,2,new String[]{"Ref_Table","Align_Table"});
		if(tw==null || tw.length<2) return;
		List<List<String>> rlisttable=table_tools.table2listtable(tw[0].getTextPanel());
		List<List<String>> alisttable=table_tools.table2listtable(tw[1].getTextPanel());
		int reflength=rlisttable.size();
		float[][] coordsr=new float[3][reflength];
		float[][] coordsa=new float[3][reflength];
		float[][] coordsa2=new float[3][alisttable.size()];

		for(int i=0;i<reflength;i++){ //over the data sets
			for(int j=0;j<3;j++){ //over the dimensions
				coordsr[j][i]=table_tools.get_number(rlisttable,i,j);
				coordsa[j][i]=table_tools.get_number(alisttable,i,j);
			}
		}
		for(int i=0;i<alisttable.size();i++){ //over the data sets
			for(int j=0;j<3;j++){ //over the dimensions
				coordsa2[j][i]=table_tools.get_number(alisttable,i,j);
			}
		}
		if(affine){
			Object[] trans=(new jreg()).affine_transformation_fiducials(coordsr,coordsa);
			table_tools.create_table("Centroids",new float[][]{(float[])trans[1],(float[])trans[2]},new String[]{"x","y","z"});
			float[][] r=algutils.convert_arr_float((double[][])trans[0]);
			table_tools.create_table("Transformation Matrix",r,null);
			//shift the coordinates to the centroid
			shift_points(coordsa2,(float[])trans[2],-1.0f);
			for(int i=0;i<coordsa2[0].length;i++){
				float[][] rt=(new jreg()).mattrans(r);
				float[] temp=matmult(new float[]{coordsa2[0][i],coordsa2[1][i],coordsa2[2][i]},rt);
				coordsa2[0][i]=temp[0]; coordsa2[1][i]=temp[1]; coordsa2[2][i]=temp[2];
				//coordsa2[0][i]+=t[0]; //coordsa2[1][i]+=t[1]; //coordsa2[2][i]+=t[2];
			}
			//shift to the reference centroid
			shift_points(coordsa2,(float[])trans[1],1.0f);
			//subtract the translation
			//shift_points(coordsa2,t,-1.0f);
			table_tools.create_table("Transformed Coords",(new jreg()).mattrans(coordsa2),new String[]{"x","y","z"});
			if(plot){
				Traj3D t3D=new Traj3D("x","y","z",coordsa2[0],coordsa2[1],coordsa2[2]);
				t3D.addPoints(coordsr[0],coordsr[1],coordsr[2],true);
				new PlotWindow3D("Transformed vs Ref Coordinates",t3D).draw();
			}
		} else {
			Object[] trans=(new jreg()).scaled_rotation_fiducials(coordsr,coordsa);
			//Object[] trans=(new jreg()).rigid_body_fiducials(coordsr,coordsa);
			//output the transformation data
			float[] t=(float[])trans[1];
			float[][] r=(float[][])trans[0];
			float s=((Float)trans[2]).floatValue(); //factor of 2?
			//float s=1.0f;
			table_tools.create_table("Rotation Matrix",r,null);
			table_tools.create_table("Translation Matrix",new float[][]{t},new String[]{"x","y","z"});
			table_tools.create_table("Centroids",new float[][]{(float[])trans[3],(float[])trans[4]},new String[]{"x","y","z"});
			IJ.log("scaling = "+s);
			//now transform the alignment coordinates
			//first reverse the rotation matrix and add scaling
			if(allows){
				for(int i=0;i<3;i++){
					for(int j=0;j<3;j++){
						r[i][j]*=s;
					}
				}
			}
			//shift the coordinates to the centroid
			shift_points(coordsa2,(float[])trans[4],-1.0f);
			for(int i=0;i<coordsa2[0].length;i++){
				float[][] rt=(new jreg()).mattrans(r);
				float[] temp=matmult(new float[]{coordsa2[0][i],coordsa2[1][i],coordsa2[2][i]},rt);
				coordsa2[0][i]=temp[0]; coordsa2[1][i]=temp[1]; coordsa2[2][i]=temp[2];
				//coordsa2[0][i]+=t[0]; //coordsa2[1][i]+=t[1]; //coordsa2[2][i]+=t[2];
			}
			//shift to the reference centroid
			shift_points(coordsa2,(float[])trans[3],1.0f);
			//subtract the translation
			//shift_points(coordsa2,t,-1.0f);
			table_tools.create_table("Transformed Coords",(new jreg()).mattrans(coordsa2),new String[]{"x","y","z"});
			if(plot){
				Traj3D t3D=new Traj3D("x","y","z",coordsa2[0],coordsa2[1],coordsa2[2]);
				t3D.addPoints(coordsr[0],coordsr[1],coordsr[2],true);
				new PlotWindow3D("Transformed vs Ref Coordinates",t3D).draw();
			}
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
