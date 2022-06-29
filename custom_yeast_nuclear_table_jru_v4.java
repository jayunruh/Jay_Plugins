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

public class custom_yeast_nuclear_table_jru_v4 implements PlugIn {
	//this plugin custom analyzes Jaspersen lab nuclear volume measurement data
	//this version starts with a single (image avg) table with 0plate,1row (letter), 2col, 3image, 4cell_number, 5pus1, 6nucleus, 7nucleoplasm (center), 8cyto (circ), 9nuc/cyto, 10nuc/center, 11nuc/pus1
	//actually now have 0plate, 1image, 2cell_number, 3pus1, 4nucleus, 5nucleoplasm (center), 6cyto (circ), 7nuc/cyto, 8nuc/center, 9nuc/pus1
	//the image contains "WellA01_" at the beginning
	//need to average all images together (when n>count threshold) and get sem values (where appropriate)
	//then measure mutant ratios and propagate errors for those
	String[] letters={"A","B","C","D","E","F","G","H"};

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1,new String[]{"Measurements_Table"});
		if(tw==null) return;
		//start by getting the parameters
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Count_Threshold (downstream analysis)",5,0);
		gd.showDialog(); if(gd.wasCanceled()) return;
		int ctthresh=(int)gd.getNextNumber();

		List<List<String>> temp=table_tools.table2listtable(tw[0].getTextPanel());
		//start by adding a unique id to end of each table row (image average)--will be the 12th row
		//also massage the table to the previous format
		for(int i=0;i<temp.size();i++){
			List<String> row=temp.get(i);
			String imname=row.get(1);
			String rowlet=imname.substring(4,5);
			String coltext=imname.substring(5,7);
			row.add(1,rowlet);
			row.add(2,coltext);
			String id=row.get(0)+"_"+row.get(1)+row.get(2);
			row.add(id);
		}
		//sort by row? may already be sorted
		table_tools.sort_listtable(temp,12);
		//now get the filtered cell stats for each well id
		List<String> cellnames=table_tools.get_cell_list(temp,12);
		//this is a list of filter column, filter type (<,=,>), and filter value
		float[][] filters={{4.0f,2.0f,(float)ctthresh}};
		List<List<String>> filtered=table_tools.get_filtered_listtable(temp,filters);
		List<List<String>> stattable=table_tools.get_cell_stat_list(filtered,12,"Avg",true); //this adds an image count column at the end (13)
		List<List<String>> errtable=table_tools.get_cell_stat_list(filtered,12,"StErr",false);
		//output a combined avg, count, and error table
		String[] combined_labels={"id","cell_count","pus1","nuc","center","cyto","nuc/cyto","nuc/center","nuc/pus1","pus1err","nucerr","centererr","cytoerr","nuc/cytoerr","nuc/centererr","nuc/pus1err"};
		List<List<String>> combined=new ArrayList<List<String>>();
		for(int i=0;i<stattable.size();i++){
			List<String> statrow=stattable.get(i);
			List<String> errrow=errtable.get(i);
			List<String> row=new ArrayList<String>();
			float cellct=table_tools.get_number(stattable,i,4);
			float imgct=table_tools.get_number(stattable,i,13);
			cellct*=imgct;
			row.add(statrow.get(12)); row.add(""+cellct); 
			for(int j=0;j<7;j++) row.add(statrow.get(j+5));
			for(int j=0;j<7;j++) row.add(errrow.get(j+5));
			combined.add(row);
		}
		List<List<String>> missing=table_tools.get_missing_cells(combined,0,cellnames," ");
		combined.addAll(missing);
		
		table_tools.sort_listtable(combined,0); //sort by wellid
		table_tools.create_table("Well_Summary",combined,combined_labels);
		//now measure the ratio to wt and output it
		//wt is in columns 1 and 7
		//first need to find out how many plates we have
		String oldplate="";
		int nplates=0;
		for(int i=0;i<combined.size();i++){
			String id=combined.get(i).get(0);
			//String plate=id.substring(0,id.indexOf("_"));
			String[] split=table_tools.split(id,"_");
			if(!split[0].equals(oldplate)){
				nplates++;
				oldplate=split[0];
			}
		}
		//there are 16 strains per plate (8 per column and 2 columns) and 5 mutants each (6 columns total each strain)
		//we have 7 measurements plus keep sem values for nuc/pus1 ratio
		int nstrains=16;
		float[][][][] vals=new float[nplates][nstrains][6][8]; //measurement order is 0pus1,1nuc,2center,3cyto,4nuc/cyto,5nuc/center,6nuc/pus1,7nuc/pus1err
		float[][][] counts=new float[nplates][nstrains][6];
		String[] platenames=new String[nplates];
		platenames[0]=table_tools.split(combined.get(0).get(0),"_")[0];
		int platecounter=0;
		for(int i=0;i<combined.size();i++){
			String id=combined.get(i).get(0);
			String[] split=table_tools.split(id,"_");
			if(!split[0].equals(platenames[platecounter])){
				platecounter++;
				platenames[platecounter]=split[0];
			}
			int[] rowcol=well2RowCol(split[1]);
			int strain=rowcol[0]-1;
			if(rowcol[1]>6) strain+=nstrains/2;
			int mut=rowcol[1]-1;
			if(rowcol[1]>6) mut-=6;
			//IJ.log(platenames[platecounter]+" , "+id+" , "+rowcol[0]+" , "+rowcol[1]+" , "+strain+" , "+mut);
			counts[platecounter][strain][mut]=table_tools.get_number(combined,i,1);
			for(int j=0;j<7;j++) vals[platecounter][strain][mut][j]=table_tools.get_number(combined,i,2+j);
			vals[platecounter][strain][mut][7]=table_tools.get_number(combined,i,15);
		}
		//now that we have all of our measurements together, calculate the ratios
		List<List<String>> rlist=new ArrayList<List<String>>();
		String[] ratio_labels=new String[55]; ratio_labels[0]="Plate"; ratio_labels[1]="Row"; ratio_labels[2]="Wt_col";
		int pos=3;
		for(int i=0;i<5;i++) ratio_labels[i+pos]="mut"+(i+1)+"/wt";
		pos+=5;
		for(int i=0;i<5;i++) ratio_labels[i+pos]="mut"+(i+1)+"/wt_sem";
		pos+=5;
		ratio_labels[pos]="wt_nuc/center";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_nuc/center";
		pos+=6;
		ratio_labels[pos]="wt_nuc/pus1";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_nuc/pus1";
		pos+=6;
		ratio_labels[pos]="wt_nuc/cyto";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_nuc/cyto";
		pos+=6;
		ratio_labels[pos]="wt_pus1";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_pus1";
		pos+=6;
		ratio_labels[pos]="wt_nuc";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_nuc";
		pos+=6;
		ratio_labels[pos]="wt_center";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_center";
		pos+=6;
		ratio_labels[pos]="wt_cyto";
		for(int i=1;i<6;i++) ratio_labels[i+pos]="mut"+i+"_cyto";
		for(int i=0;i<nplates;i++){
			for(int j=0;j<nstrains;j++){
				List<String> rvals=new ArrayList<String>();
				rvals.add(platenames[i]);
				int rownum=j; 
				if(j>=nstrains/2){rvals.add(letters[rownum-nstrains/2]); rvals.add("7");}
				else{rvals.add(letters[rownum]); rvals.add("1");}
				float wtcount=counts[i][j][0];
				float[] ratios=new float[5];
				float[] rerrs=new float[5];
				if(wtcount>(float)ctthresh) for(int k=0;k<5;k++){ratios[k]=vals[i][j][k+1][6]/vals[i][j][0][6]; rerrs[k]=propRatioErrs(vals[i][j][k+1][6],vals[i][j][0][6],vals[i][j][k+1][7],vals[i][j][0][7]);}
				for(int k=0;k<5;k++) rvals.add(""+ratios[k]);
				for(int k=0;k<5;k++) rvals.add(""+rerrs[k]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][5]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][6]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][4]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][0]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][1]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][2]);
				for(int k=0;k<6;k++) rvals.add(""+vals[i][j][k][3]);
				rlist.add(rvals);
			}
		}
		table_tools.sort_listtable(rlist,0,true); //sort by plate number
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
