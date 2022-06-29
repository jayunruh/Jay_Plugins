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
import ij.text.*;
import java.util.*;
import jguis.*;
import jalgs.*;
import jalgs.jsim.*;
import ij.io.*;
import java.io.*;

public class custom_yeast_nuclear_table_jru_v2 implements PlugIn {
	//this plugin custom analyzes Jaspersen lab nuclear volume measurement data
	//this version has three tables per image: Volumes, greenQuantif, and redQuantif
	//file names are M or Q_NDExpWellXXX_PointXXXX_SeqXXXX...
	//in the final averaged plate, the first row and the 7th row are wt, the others should be expressed as ratios to wt if cell count is high enough?
	//propagate errors through the ratios as always

	public void run(String arg) {
		//start by getting the parameters
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Min_Nuc_Vol (pix)",250,0);
		gd.addNumericField("Max_Nuc_Vol (pix)",2000,0);
		gd.addChoice("Red_Center_Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd.addCheckbox("Filter_On_Red",false);
		gd.addChoice("Red_Shift_Statistic",jstatistics.stats,jstatistics.stats[0]);
		gd.addNumericField("Min_Red_Multiplier",1.0,5,15,null);
		gd.addNumericField("Max_Red_Multiplier",1.0,5,15,null);
		gd.addNumericField("Count_Threshold (downstream analysis)",20,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float minvol=(float)gd.getNextNumber();
		float maxvol=(float)gd.getNextNumber();
		String stat=jstatistics.stats[gd.getNextChoiceIndex()];
		boolean filterred=gd.getNextBoolean();
		String stat2=jstatistics.stats[gd.getNextChoiceIndex()];
		float minrmult=(float)gd.getNextNumber();
		float maxrmult=(float)gd.getNextNumber();
		int ctthresh=(int)gd.getNextNumber();

		//now get the files and combine them into a table
		OpenDialog od = new OpenDialog("Open Image Sequence...", arg);
        		String directory = od.getDirectory();
		String fname=od.getFileName();
		if(fname==null || fname.length()==0){return;}
		int period=fname.lastIndexOf('.');
		String extension=fname.substring(period+1);
		String[] mask={extension,"Volumes"};
		String[] flist=(new jdataio()).get_sorted_string_list(directory,mask,null);
		List<List<String>> combined=new ArrayList<List<String>>();
		String[] combined_labels={"Image","Object","vol","volpix","surf","surfpix","greensum","greenavg","redsum","redavg","well"};
		for(int i=0;i<flist.length;i++){
			String prefix=flist[i].substring(2,flist[i].length()-16);
			int wellpos=flist[i].indexOf("Well");
			String well=flist[i].substring(wellpos+4,wellpos+7);
			IJ.log(prefix+" , "+well);
			String gqname="Q_"+prefix+".nd2_greenQuantif.xls";
			String rqname="Q_"+prefix+".nd2_redQuantif.xls";
			List<List<String>> vtable=table_tools.getTableFromFile(new File(directory+flist[i]),"\t",false);
			List<List<String>> gtable=table_tools.getTableFromFile(new File(directory+gqname),"\t",false);
			List<List<String>> rtable=table_tools.getTableFromFile(new File(directory+rqname),"\t",false);
			for(int j=1;j<vtable.size();j++){
				List<String> row=new ArrayList<String>();
				row.add(prefix+".nd2"); row.add(""+j); 
				row.add(vtable.get(j).get(3)); row.add(vtable.get(j).get(4)); row.add(vtable.get(j).get(5)); row.add(vtable.get(j).get(6)); 
				row.add(gtable.get(j).get(4)); row.add(gtable.get(j).get(3)); row.add(rtable.get(j).get(4)); row.add(rtable.get(j).get(3));
				row.add(well);
				combined.add(row);
			}
		}
		//sort the combined table by well
		table_tools.sort_listtable(combined,10);
		//now optionally output the combined tables for testing
		table_tools.create_table("Combined Table",combined,combined_labels);
		//probably should save this table somewhere
		//now we need to filter the results
		//filter on volume, pus1 intensity
		List<List<String>> filtered=new ArrayList<List<String>>();
		String[] filtered_labels={"Image","Object","vol","volpix","surf","surfpix","greensum","greenavg","redsum","redavg","well","green/surf","red/vol"};
		//first filter on volume
		for(int i=0;i<combined.size();i++){
			float vol=table_tools.get_number(combined,i,3);
			if(!Float.isNaN(vol) && vol>=(float)minvol && vol<=(float)maxvol){
				List<String> row=combined.get(i);
				float green=table_tools.get_number(combined,i,6);
				float red=table_tools.get_number(combined,i,8);
				float surf=table_tools.get_number(combined,i,4);
				float vol2=table_tools.get_number(combined,i,2);
				row.add(""+(green/surf)); row.add(""+(red/vol2));
				filtered.add(row);
				
			}
		}
		table_tools.create_table("Volume Filtered",filtered,filtered_labels);
		//now get the well statistics for pus1
		List<List<String>> wellstats=table_tools.get_cell_stat_list(filtered,10,stat,true);
		List<List<String>> wellstats2=table_tools.get_cell_stat_list(filtered,10,stat2,true);
		List<List<String>> avgstats=wellstats;
		if(!stat.equals("Avg")); avgstats=table_tools.get_cell_stat_list(filtered,10,"Avg",true);
		String[] avg_labels={"Image","Object","vol","volpix","surf","surfpix","greensum","greenavg","redsum","redavg","well","green/surf","red/vol","count"};
		table_tools.create_table("Volume Filtered_Avg",avgstats,avg_labels);
		//make a filter set for each well based on pus1 values
		List<List<String>> filtered2=new ArrayList<List<String>>();
		float[] cellstatval=table_tools.get_column_array(wellstats,12);
		float[] cellstat2val=table_tools.get_column_array(wellstats2,12);
		String[] cellnames=table_tools.get_listtable_column(wellstats,10);
		if(!filterred){
			filtered2=filtered;
		} else {
			float[][] limits=new float[cellnames.length][2];
			for(int i=0;i<cellnames.length;i++){
				limits[i][0]=cellstatval[i]-cellstat2val[i]*minrmult;
				limits[i][1]=cellstatval[i]+cellstat2val[i]*maxrmult;
			}
			//now filter on a cell by cell basis based on the pus1 statistics
			String oldwell=filtered.get(0).get(10);
			float[] oldfilter=getWellFilter(oldwell,cellnames,limits);
			for(int i=0;i<filtered.size();i++){
				String newwell=filtered.get(i).get(10);
				if(!newwell.equals(oldwell)){
					oldwell=newwell;
					oldfilter=getWellFilter(oldwell,cellnames,limits);
				}
				float filterval=table_tools.get_number(filtered,i,12);
				if(filterval>=oldfilter[0] && filterval<=oldfilter[1]) filtered2.add(filtered.get(i));
			}
			//now optionally output the combined tables for testing
			table_tools.create_table("Vol&Intensity Filtered",filtered2,filtered_labels);
		}
		//finally get the avg, sterr, and count per well of the newly filtered set
		List<List<String>> filteredcellstats1=table_tools.get_cell_stat_list(filtered2,10,"Avg",true);
		List<List<String>> filteredcellstats2=table_tools.get_cell_stat_list(filtered2,10,"StErr",true);
		String[] stat_labels={"Image","Object","vol","volpix","surf","surfpix","greensum","greenavg","redsum","redavg","well","green/surf","red/vol","count","row","col","volsem","volpixsem","surfsem","surfpixsem","greensumsem","greenavgsem","redsumsem","redavgsem","green/surf_sem","red/vol_sem"};
		for(int i=0;i<filteredcellstats1.size();i++){
			String well=filteredcellstats1.get(i).get(10);
			int[] rowcol=well2RowCol(well);
			filteredcellstats1.get(i).add(""+rowcol[0]);
			filteredcellstats1.get(i).add(""+rowcol[1]);
		}
		//append the sem values to the avg table
		String[] tempcol=table_tools.get_listtable_column(filteredcellstats2,2);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,3);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,4);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,5);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,6);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,7);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,8);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,9);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,11);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		tempcol=table_tools.get_listtable_column(filteredcellstats2,12);
		table_tools.add_listtable_column(filteredcellstats1,tempcol,1000);
		table_tools.create_table("Final Filtered Well_Avg",filteredcellstats1,stat_labels);

		//now create another table where values are expressed as ratios to the wt (both red and green)
		//wt is in columns 1 and 7
		//there are at most 16 labels (1 wt and 5 mutants each)
		//start by filling the mutant values tables
		int nlabels=16;
		int halflabels=nlabels/2;
		float[][] gvals=new float[nlabels][6];
		float[][] gsem=new float[nlabels][6];
		float[][] rvals=new float[nlabels][6];
		float[][] rsem=new float[nlabels][6];
		int[][] counts=new int[nlabels][6];
		for(int i=0;i<filteredcellstats1.size();i++){
			int row=table_tools.get_integer(filteredcellstats1,i,14);
			int col=table_tools.get_integer(filteredcellstats1,i,15);
			int label=row-1;
			if(col>6) label+=halflabels;
			int mut=col-1;
			if(col>6) mut-=6;
			//IJ.log(""+row+" , "+col+" , "+label+" , "+mut);
			counts[label][mut]=table_tools.get_integer(filteredcellstats1,i,13);
			gvals[label][mut]=table_tools.get_number(filteredcellstats1,i,11);
			rvals[label][mut]=table_tools.get_number(filteredcellstats1,i,12);
			gsem[label][mut]=table_tools.get_number(filteredcellstats1,i,24);
			rsem[label][mut]=table_tools.get_number(filteredcellstats1,i,25);
		}
		//now calculate the ratios and their sem values if possible
		float[][] gratios=new float[nlabels][6];
		float[][] grsem=new float[nlabels][6];
		float[][] rratios=new float[nlabels][6];
		float[][] rrsem=new float[nlabels][6];
		for(int i=0;i<nlabels;i++){
			if(counts[i][0]>=ctthresh){
				gratios[i][0]=1.0f; grsem[i][0]=gsem[i][0]/gvals[i][0];
				rratios[i][0]=1.0f; rrsem[i][0]=rsem[i][0]/rvals[i][0];
				for(int j=1;j<6;j++){
					if(counts[i][j]>=ctthresh){
						gratios[i][j]=gvals[i][j]/gvals[i][0];
						grsem[i][j]=propRatioErrs(gvals[i][j],gvals[i][0],gsem[i][j],gsem[i][0]);
						rratios[i][j]=rvals[i][j]/rvals[i][0];
						rrsem[i][j]=propRatioErrs(rvals[i][j],rvals[i][0],rsem[i][j],rsem[i][0]);
					}
				}
			}
		}
		//finally output a new table listing the ratios
		List<List<String>> rlist=new ArrayList<List<String>>();
		String[] ratio_labels=new String[27]; ratio_labels[0]="Strain_id"; ratio_labels[1]="row"; ratio_labels[2]="col";
		ratio_labels[3]="wt_green"; ratio_labels[9]="wt_greensem"; ratio_labels[15]="wt_red"; ratio_labels[21]="wt_redsem";
		for(int i=1;i<6;i++) ratio_labels[i+3]="mut"+i+"_green";
		for(int i=1;i<6;i++) ratio_labels[i+9]="mut"+i+"_greensem";
		for(int i=1;i<6;i++) ratio_labels[i+15]="mut"+i+"_red";
		for(int i=1;i<6;i++) ratio_labels[i+21]="mut"+i+"_redsem";
		for(int i=0;i<nlabels;i++){
			List<String> row=new ArrayList<String>();
			row.add(""+i);
			int col=1;
			if(i>=halflabels) col=7;
			int rtemp=i+1;
			if(i>=halflabels) rtemp-=halflabels;
			row.add(""+rtemp);
			row.add(""+col);
			for(int j=0;j<6;j++) row.add(""+gratios[i][j]);
			for(int j=0;j<6;j++) row.add(""+grsem[i][j]);
			for(int j=0;j<6;j++) row.add(""+rratios[i][j]);
			for(int j=0;j<6;j++) row.add(""+rrsem[i][j]);
			rlist.add(row);
		}
		table_tools.create_table("Ratio Well_Avg",rlist,ratio_labels);
	}

	public float propRatioErrs(float num,float den,float numsem,float densem){
		rngs random=new rngs();
		float[] temp=new float[100];
		for(int i=0;i<100;i++) temp[i]=(float)random.gasdev(num,numsem)/(float)random.gasdev(den,densem);
		return jstatistics.getstatistic("StDev",temp,null);
	}

	public int[] well2RowCol(String well){
		String lett=well.substring(0,1);
		String col1=well.substring(1,well.length());
		String ulett=lett.toUpperCase();
		int row=(int)ulett.charAt(0)-64;
		int col=Integer.parseInt(col1);
		return new int[]{row,col};
	}

	public float[] getWellFilter(String well,String[] wellnames,float[][] limits){
		for(int i=0;i<wellnames.length;i++){
			if(well.equals(wellnames[i])){
				return limits[i];
			}
		}
		return null;
	}

}
