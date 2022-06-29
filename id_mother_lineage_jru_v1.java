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

public class id_mother_lineage_jru_v1 implements PlugIn {

	public void run(String arg) {
		float derthresh=5000.0f; //minimum intensity derivative for mother cell
		float velthresh=2.0f; //minimum velocity threshold for mother cell
		float distthresh=15.0f; //maximum distance for mother cell
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Table"});
		if(tw==null || tw.length<1) return;
		TextPanel tp=tw[0].getTextPanel();
		String[] col_labels=table_tools.getcollabels(tp);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		//need to have id, x,y,frame,start, and intensity
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("id_column",col_labels,col_labels[0]);
		gd.addChoice("start_column",col_labels,col_labels[1]);
		gd.addChoice("frame_column",col_labels,col_labels[2]);
		gd.addChoice("x_column",col_labels,col_labels[3]);
		gd.addChoice("y_column",col_labels,col_labels[4]);
		gd.addChoice("intensity_column",col_labels,col_labels[5]);
		gd.addNumericField("Intensity_Der_Thresh",derthresh,5,15,null);
		gd.addNumericField("Velocity_Thresh",velthresh,5,15,null);
		gd.addNumericField("Max_Distance",distthresh,5,15,null);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int idcol=gd.getNextChoiceIndex();
		int startcol=gd.getNextChoiceIndex();
		int framecol=gd.getNextChoiceIndex();
		int xcol=gd.getNextChoiceIndex();
		int ycol=gd.getNextChoiceIndex();
		int intcol=gd.getNextChoiceIndex();
		derthresh=(float)gd.getNextNumber();
		velthresh=(float)gd.getNextNumber();
		distthresh=(float)gd.getNextNumber();
		List<List<List<String>>> celllist=table_tools.get_cells_listtable(listtable,idcol,1);
		String[] mothers=new String[celllist.size()];
		mothers[0]="NaN";
		List<List<String>> cell1=celllist.get(0);
		for(int i=0;i<cell1.size();i++) cell1.get(i).add(mothers[0]);
		for(int i=1;i<celllist.size();i++){
			//for each new cell, search for its mother
			//if none is found, set it to "NaN"
			//if multiple are found, use the one with the closest distance
			List<List<String>> qcell=celllist.get(i);
			float queryframe=table_tools.get_number(qcell,0,startcol)+table_tools.get_number(qcell,0,framecol);
			float qx=table_tools.get_number(qcell,0,xcol);
			float qy=table_tools.get_number(qcell,0,ycol);
			int bestmother=-1;
			float bestdist=Float.MAX_VALUE;
			for(int j=0;j<celllist.size();j++){
				if(j==i) continue;
				List<List<String>> bcell=celllist.get(j);
				//see if the cell is present at the same time
				float startframe=table_tools.get_number(bcell,0,startcol)+table_tools.get_number(bcell,0,framecol);
				float endframe=table_tools.get_number(bcell,bcell.size()-1,startcol)+table_tools.get_number(bcell,bcell.size()-1,framecol);
				if(queryframe>startframe && queryframe<=endframe){
					//find the equivalent frame
					int eqframe=(int)(queryframe-startframe);
					float dist=get_dist(qx,qy,table_tools.get_number(bcell,eqframe,xcol),table_tools.get_number(bcell,eqframe,ycol));
					IJ.log(""+i+" , "+j+" =\t"+dist);
					if(dist<=distthresh){
						float vel=get_dist(table_tools.get_number(bcell,eqframe-1,xcol),table_tools.get_number(bcell,eqframe-1,ycol),table_tools.get_number(bcell,eqframe,xcol),table_tools.get_number(bcell,eqframe,ycol));
						float intder=(float)Math.abs(table_tools.get_number(bcell,eqframe,intcol)-table_tools.get_number(bcell,eqframe-1,intcol));
						if(vel>=velthresh && intder>=derthresh){
							if(dist<bestdist){
								bestmother=j;
								bestdist=dist;
							}
						}
					}
				}
			}
			if(bestmother==-1) mothers[i]="NaN";
			else mothers[i]=celllist.get(bestmother).get(0).get(idcol);
			for(int j=0;j<qcell.size();j++) qcell.get(j).add(mothers[i]);
		}
		//now rearrange back into a single table
		List<List<String>> outtable=new ArrayList<List<String>>();
		for(int i=0;i<celllist.size();i++){
			List<List<String>> selcell=celllist.get(i);
			for(int j=0;j<selcell.size();j++) outtable.add(selcell.get(j));
		}
		String[] newcollabels=new String[col_labels.length+1];
		for(int i=0;i<col_labels.length;i++) newcollabels[i]=col_labels[i];
		newcollabels[col_labels.length]="mother_id";
		table_tools.create_table("Mother ID table",outtable,newcollabels);
	}

	private static float get_dist(float x1,float y1,float x2,float y2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}

}
