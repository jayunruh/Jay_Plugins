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
import ij.text.*;
import java.util.*;
import jguis.*;
import jalgs.*;

public class spatial_transcriptomics_image_maker_jru_v1 implements PlugIn {

	public void run(String arg) {
		//here we are going to generate images for slide seq data
		//need 3 inputs
		//an image containing the barcode x gene count matrix
		//a table of gene vs. tissue count
		//a table of x y positions
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
		ImagePlus[] imps=jutils.selectImages(false,1,new String[]{"Gene_Expression_Matrix"});
		if(imps==null) return;
		int nbeads=xpos.length;
		int ngenes=genenames.length;
		short[] gematrix=(short[])imps[0].getProcessor().getPixels();
		float xmin=jstatistics.getstatistic("Min",xpos,null);
		float xmax=jstatistics.getstatistic("Max",xpos,null);
		float ymin=jstatistics.getstatistic("Min",ypos,null);
		float ymax=jstatistics.getstatistic("Max",ypos,null);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("X_Min",xmin,5,15,null);
		gd.addNumericField("X_Max",xmax,5,15,null);
		gd.addNumericField("Y_Min",ymin,5,15,null);
		gd.addNumericField("Y_Max",ymax,5,15,null);
		gd.addNumericField("Spot_Radius",5.0,5,15,null);
		gd.addNumericField("Number_of_genes",3,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		xmin=(float)gd.getNextNumber();
		xmax=(float)gd.getNextNumber();
		ymin=(float)gd.getNextNumber();
		ymax=(float)gd.getNextNumber();
		float rad=(float)gd.getNextNumber();
		int genecounter=(int)gd.getNextNumber();
		int[] geneindices=new int[genecounter];
		GenericDialog gd2=new GenericDialog("Select Gene");
		for(int i=0;i<genecounter;i++) gd2.addChoice("Gene"+(i+1),genenames,genenames[i]);
		gd2.showDialog(); if(gd2.wasCanceled()) return;
		for(int i=0;i<genecounter;i++) geneindices[i]=gd2.getNextChoiceIndex();
		float[] limits=new float[]{xmin,xmax,ymin,ymax};
		int width=(int)(limits[1]-limits[0]);
		int height=(int)(limits[3]-limits[2]);
		ImageStack stack=new ImageStack(width,height);
		for(int i=0;i<genecounter;i++){
			short[] selcounts=new short[nbeads];
			System.arraycopy(gematrix,geneindices[i]*nbeads,selcounts,0,nbeads);
			short[] image=drawImage(xpos,ypos,selcounts,rad,limits);
			String gname=genenames[geneindices[i]];
			stack.addSlice(gname,image);
		}
		jutils.create_hyperstack("Slide_Seq Images",stack,1,1,genecounter,true,null).show();
	}

	public short[] drawImage(float[] xpos,float[] ypos,short[] counts,float spotrad,float[] limits){
		int width=(int)(limits[1]-limits[0]);
		float xmin=limits[0];
		int height=(int)(limits[3]-limits[2]);
		float ymin=limits[2];
		short[] image=new short[width*height];
		for(int i=0;i<counts.length;i++){
			if(counts[i]!=(short)0){
				float xc=xpos[i]-xmin;
				float yc=ypos[i]-ymin;
				drawCircle(image,width,height,xc,yc,counts[i],spotrad);
			}
		}
		return image;
	}

	public void drawCircle(short[] image,int width,int height,float xc,float yc,short value,float rad){
		int xmin=(int)(xc-rad);
		int xmax=xmin+(int)(2.0f*rad+1.0f);
		int ymin=(int)(yc-rad);
		int ymax=ymin+(int)(2.0f*rad+1.0f);
		float rad2=rad*rad;
		for(int i=ymin;i<=ymax;i++){
			float ypos=0.0f;
			if(i>=0 && i<height) ypos=(float)i;
			else continue;
			for(int j=xmin;j<=xmax;j++){
				float xpos=0.0f;
				if(j>=0 && j<width) xpos=(float)j;
				else continue;
				float dist2=(xpos-xc)*(xpos-xc)+(ypos-yc)*(ypos-yc);
				if(dist2<=rad2) image[j+i*width]=value;
			}
		}
		return;
	}

}
