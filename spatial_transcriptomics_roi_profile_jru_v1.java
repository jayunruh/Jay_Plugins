/*******************************************************************************
 * Copyright (c) 2020 Jay Unruh, Stowers Institute for Medical Research.
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
import ij.plugin.frame.RoiManager;
import ij.text.*;
import java.util.*;
import jguis.*;
import jalgs.*;

public class spatial_transcriptomics_roi_profile_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we are going to get a gene profile from slide seq data
		//need 3 inputs
		//an image containing the barcode x gene count matrix
		//a table of gene vs. tissue count
		//a table of x y positions
		//use an previously generated image with x and y offset to find the beads contained within our roi
		//need to some those gene profiles for those beads
		TextWindow[] tw=jutils.selectTables(false,2,new String[]{"tissue_count_table","xy_position_table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels1=table_tools.getcollabels(tp);
		List<List<String>> tctable=table_tools.table2listtable(tp);
		int namecol=1;
		if(tctable.get(0).size()<3) namecol=0;
		boolean hasquotes=false;
		if(tctable.get(0).get(0).indexOf("\"")>=0) hasquotes=true;
		String[] genenames=new String[tctable.size()];
		for(int i=0;i<genenames.length;i++){
			genenames[i]=tctable.get(i).get(namecol);
			//trim out quotes if they are there
			//genenames[i]=genenames[i].replace("\"","");
			//genenames[i]=genenames[i].replaceAll("^\"|\"$", "");
			if(hasquotes) genenames[i]=genenames[i].substring(1,genenames[i].length()-1);
		}
		TextPanel tp2=tw[1].getTextPanel();
		String[] col_labels2=table_tools.getcollabels(tp2);
		List<List<String>> postable=table_tools.table2listtable(tp2);
		float[] xpos=table_tools.get_column_array(postable,1);
		float[] ypos=table_tools.get_column_array(postable,2);
		ImagePlus[] imps=jutils.selectImages(false,2,new String[]{"Roi_Image","Gene_Expression_Matrix"});
		if(imps==null) return;
		int nbeads=xpos.length;
		int ngenes=genenames.length;
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			IJ.error("add Rois to the roi manager");
		}
		Roi[] rois=rman.getRoisAsArray();
		int imwidth=imps[0].getWidth();
		int imheight=imps[0].getHeight();
		short[] gematrix=(short[])imps[1].getProcessor().getPixels();
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X_Offset",2000.0f,5,15,null);
		gd.addNumericField("Y_Offset",1500.0f,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float xmin=(float)gd.getNextNumber();
		float ymin=(float)gd.getNextNumber();

		//the first gene expression profile is always the whole image
		int[][] geprofile=new int[rois.length+1][ngenes];
		for(int i=0;i<ypos.length;i++){
			//if the point is in the image, put it in the first column
			if(xpos[i]>=xmin && xpos[i]<(imwidth+xmin) && ypos[i]>=ymin && ypos[i]<(imheight+ymin)){
				int[] subarr=getImageColumn(gematrix,nbeads,ngenes,i);
				addToArray(geprofile[0],subarr);
				for(int j=0;j<rois.length;j++){
					if(rois[j].contains((int)(xpos[i]-xmin),(int)(ypos[i]-ymin))){
						addToArray(geprofile[j+1],subarr);
					}
				}
			}
		}

		//now write the output to a table
		String[] outlabels=new String[rois.length+2];
		outlabels[0]="Gene";
		outlabels[1]="Whole_Image";
		for(int i=2;i<(rois.length+2);i++) outlabels[i]="roi"+(i-1);
		String tablevals=printGEProfiles(geprofile,genenames);
		new TextWindow("GE_Profiles",table_tools.print_string_array(outlabels),tablevals,400,200);
	}

	public String printGEProfiles(int[][] image,String[] genelist){
		StringBuffer retvals=new StringBuffer();
		retvals.append(genelist[0]);
		retvals.append("\t"+table_tools.print_int_array(getImageColumn(image,0)));
		for(int i=1;i<image[0].length;i++){
			retvals.append("\n"+genelist[i]);
			retvals.append("\t"+table_tools.print_int_array(getImageColumn(image,i)));
		}
		return retvals.toString();
	}

	public int[] getImageColumn(int[][] image,int col){
		int[] retarr=new int[image.length];
		for(int i=0;i<image.length;i++) retarr[i]=image[i][col];
		return retarr;
	}

	public void addToArray(int[] arr,int[] addvals){
		for(int i=0;i<addvals.length;i++) arr[i]+=addvals[i];
		return;
	}

	public int[] getImageColumn(short[] image,int width,int height,int col){
		int[] colvals=new int[height];
		for(int i=0;i<height;i++){
			colvals[i]=image[i*width+col]&0xffff;
		}
		return colvals;
	}

}
