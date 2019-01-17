/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.text.TextPanel;
import ij.text.TextWindow;
import jalgs.jstatistics;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class table_tools{
	// here we have static utility methods for working with ImageJ table panels
	// and multidimensional string arrays
	public static String[] split_string_tab(String line){
		return split(line,"\t");
	}

	public static String[] split(String line,String delim){
		String temp;
		if(line.endsWith(delim)){
			temp=line.substring(0,line.length()-delim.length());
		}else{
			temp=line.substring(0);
		}
		return temp.split(delim);
	}
	
	public static String[] split(String line,String delim,boolean noconsec){
		String temp;
		if(line.endsWith(delim)){
			temp=line.substring(0,line.length()-delim.length());
		}else{
			temp=line.substring(0);
		}
		String[] temp2=temp.split(delim);
		if(!noconsec) return temp2;
		List<String> nonnulls=new ArrayList<String>();
		for(int i=0;i<temp2.length;i++){
			if(temp2[i].length()>0) nonnulls.add(temp2[i]);
		}
		return list2stringarray(nonnulls);
	}
	
	public static int[] wellName2RowCol(String wellname) {
		//turns a wellname (in A1 or A01 format) to row and column numbers (base 1)
		int slen=wellname.length();
		if(slen>3) return new int[] {-1,-1};
		int col=Integer.parseInt(wellname.substring(1,slen));
		String temp=wellname.toUpperCase();
		int row=1+(int)temp.charAt(0)-(int)('A');
		return new int[] {row,col};
	}

	public static boolean is_number(String snumber){
		try{
			Float.parseFloat(snumber);
		}catch(NumberFormatException e){
			return false;
		}
		return true;
	}

	public static float get_number(List<List<String>> table,int row,int column){
		if(row<table.size()){
			List<String> r=table.get(row);
			if(column<r.size()){
				String temp=r.get(column);
				if(is_number(temp)){
					return Float.parseFloat(temp);
				}
			}
		}
		return Float.NaN;
	}

	public static double get_double(List<List<String>> table,int row,int column){
		if(row<table.size()){
			List<String> r=table.get(row);
			if(column<r.size()){
				String temp=r.get(column);
				if(is_number(temp)){
					return Double.parseDouble(temp);
				}
			}
		}
		return Double.NaN;
	}

	public static int get_integer(List<List<String>> table,int row,int column){
		if(row<table.size()){
			List<String> r=table.get(row);
			if(column<r.size()){
				String temp=r.get(column);
				if(is_number(temp)){
					return Integer.parseInt(temp);
				}
			}
		}
		return 0;
	}

	public static String[] getcollabels(TextPanel tp){
		String headings=tp.getColumnHeadings();
		String[] col_labels=split_string_tab(headings);
		return col_labels;
	}

	public static String createcollabels(int ncols){
		return createcollabels(ncols,"col");
	}
	
	public static String[] make_labels_unique(String[] labels){
		String[] new_labels=new String[labels.length];
		for(int i=0;i<labels.length;i++){
			new_labels[i]=(labels[i].trim()).replace(' ','_');
			if(new_labels[i].equals("")) new_labels[i]="_";
			new_labels[i]=new_labels[i].replace('*','_');
			new_labels[i]=new_labels[i].replace('.','_');
			new_labels[i]=new_labels[i].replace('-','_');
			new_labels[i]=new_labels[i].replace('/','_');
			new_labels[i]=new_labels[i].replace('#','_');
			new_labels[i]=new_labels[i].replace('[','_');
			new_labels[i]=new_labels[i].replace(']','_');
			new_labels[i]=new_labels[i].replace('\"','_');
			//IJ.log(new_col_labels[i]);
		}
		return new_labels;
	}

	public static String createcollabels(int ncols,String prefix){
		return createcollabels(ncols,prefix,0);
	}
	
	public static String createcollabels(int ncols,String prefix,int delim){
		StringBuffer sb=new StringBuffer();
		sb.append(prefix+"1");
		String sep=get_delim(delim);
		for(int i=1;i<ncols;i++){
			sb.append(sep+prefix+(i+1));
		}
		return sb.toString();
	}

	public static String[][] tp2array(TextPanel tp){
		int nlines=tp.getLineCount();
		String[][] retarray=new String[nlines][];
		for(int i=0;i<nlines;i++){
			String line=tp.getLine(i);
			retarray[i]=split_string_tab(line);
		}
		return retarray;
	}

	public static String[][] tp2array(String[] lines){
		int nlines=lines.length;
		String[][] retarray=new String[nlines][];
		for(int i=0;i<nlines;i++){
			String line=lines[i];
			retarray[i]=split_string_tab(line);
		}
		return retarray;
	}
	
	public static List<List<String>> getTableFromFile(File infile,String delim,boolean noconsec){
		try{
			BufferedReader b=new BufferedReader(new FileReader(infile));
			List<String> lines=new ArrayList<String>();
			String temp=b.readLine();
			while(temp!=null && temp.length()>1){
				lines.add(temp);
				temp=b.readLine();
			}
			List<List<String>> listtable=table_tools.table2listtable(lines,delim,noconsec);
			lines=null;
			b.close();
			return listtable;
		}
		catch(IOException e){
			IJ.error("File Reader",e.getMessage());
			return null;
		}
	}
	
	public static boolean writeTableToFile(String path,String[] collabels,List<List<String>> listtable,int delim){
		String temp=null;
		if(collabels!=null) temp=print_string_array(collabels,delim);
		if(temp==null) temp=createcollabels(listtable.get(0).size(),"Col",delim);
		return writeTableToFile(path,temp,print_listtable(listtable,delim));
	}
	
	public static boolean writeTableToFile(String path,String collabels,String listtable){
		try{
			BufferedWriter b=new BufferedWriter(new FileWriter(new File(path)));
			if(collabels!=null) b.write(collabels+"\n");
			b.write(listtable);
			b.close();
			return true;
		} catch(IOException e){
			IJ.error("error writing file",e.getMessage());
			return false;
		}
	}

	public static List<List<String>> table2listtable(TextPanel tp){
		int nlines=tp.getLineCount();
		List<List<String>> retvals=new ArrayList<List<String>>();
		int longest=0;
		for(int i=0;i<nlines;i++){
			String line=tp.getLine(i);
			String[] temp2=split_string_tab(line);
			List<String> temp=stringarray2list(temp2);
			if(i==0){
				longest=temp.size();
			}
			for(int j=temp.size();j<longest;j++){
				temp.add("");
			}
			retvals.add(temp);
		}
		return retvals;
	}

	public static List<List<String>> table2listtable(String data,String delim){
		String[] lines=split(data,"\n");
		return table2listtable(lines,delim);
	}

	public static List<List<String>> table2listtable(String[] lines,String delim){
		int nlines=lines.length;
		// eliminate empty lines at the end
		while(lines[nlines-1]==null||lines[nlines-1].equals("")){nlines--;}
		return table2listtable(lines,delim,nlines);
	}

	public static List<List<String>> table2listtable(String[] lines,String delim,int nlines){
		List<List<String>> retvals=new ArrayList<List<String>>();
		int longest=0;
		for(int i=0;i<nlines;i++){
			String line=lines[i];
			String[] temp2=split(line,delim);
			List<String> temp=stringarray2list(temp2);
			if(i==0){
				longest=temp.size();
			}
			for(int j=temp.size();j<longest;j++){
				temp.add("");
			}
			retvals.add(temp);
		}
		return retvals;
	}
	
	public static List<List<String>> table2listtable(String[] lines,String delim,boolean noconsec){
		int nlines=lines.length;
		// eliminate empty lines at the end
		while(lines[nlines-1]==null||lines[nlines-1].equals("")){nlines--;}
		return table2listtable(lines,delim,nlines,noconsec);
	}
	
	public static List<List<String>> table2listtable(String[] lines,String delim,int nlines,boolean noconsec){
		List<List<String>> retvals=new ArrayList<List<String>>();
		int longest=0;
		for(int i=0;i<nlines;i++){
			String line=lines[i];
			String[] temp2=split(line,delim,noconsec);
			List<String> temp=stringarray2list(temp2);
			if(i==0){
				longest=temp.size();
			}
			for(int j=temp.size();j<longest;j++){
				temp.add("");
			}
			retvals.add(temp);
		}
		return retvals;
	}

	public static List<List<String>> table2listtable(List<String> lines,String delim){
		int nlines=lines.size();
		// eliminate empty lines at the end
		while(lines.get(nlines-1)==null||lines.get(nlines-1).equals("")){nlines--;}
		return table2listtable(lines,delim,nlines);
	}

	public static List<List<String>> table2listtable(List<String> lines,String delim,int nlines){
		List<List<String>> retvals=new ArrayList<List<String>>();
		int longest=0;
		for(int i=0;i<nlines;i++){
			String line=lines.get(i);
			String[] temp2=split(line,delim);
			List<String> temp=stringarray2list(temp2);
			if(i==0){
				longest=temp.size();
			}
			for(int j=temp.size();j<longest;j++){
				temp.add("");
			}
			retvals.add(temp);
		}
		return retvals;
	}
	
	public static List<List<String>> table2listtable(List<String> lines,String delim,boolean noconsec){
		int nlines=lines.size();
		// eliminate empty lines at the end
		while(lines.get(nlines-1)==null||lines.get(nlines-1).equals("")){nlines--;}
		return table2listtable(lines,delim,nlines,noconsec);
	}

	public static List<List<String>> table2listtable(List<String> lines,String delim,int nlines,boolean noconsec){
		List<List<String>> retvals=new ArrayList<List<String>>();
		int longest=0;
		for(int i=0;i<nlines;i++){
			String line=lines.get(i);
			String[] temp2=split(line,delim,noconsec);
			List<String> temp=stringarray2list(temp2);
			if(i==0){
				longest=temp.size();
			}
			for(int j=temp.size();j<longest;j++){
				temp.add("");
			}
			retvals.add(temp);
		}
		return retvals;
	}

	public static List<List<String>> table2listtable(String[] lines){
		int nlines=lines.length;
		List<List<String>> retvals=new ArrayList<List<String>>();
		int longest=0;
		for(int i=0;i<nlines;i++){
			String line=lines[i];
			String[] temp2=split_string_tab(line);
			List<String> temp=stringarray2list(temp2);
			if(i==0){
				longest=temp.size();
			}
			for(int j=temp.size();j<longest;j++){
				temp.add("");
			}
			retvals.add(temp);
		}
		return retvals;
	}
	
	public static List<List<String>> table2listtable(float[][] tabdata){
		List<List<String>> temp=new ArrayList<List<String>>();
		for(int i=0;i<tabdata.length;i++) {
			List<String> temp2=new ArrayList<String>();
			for(int j=0;j<tabdata[i].length;j++) {
				temp2.add(""+tabdata[i][j]);
			}
			temp.add(temp2);
		}
		return temp;
	}

	public static List<String> stringarray2list(String[] arr){
		List<String> temp=new ArrayList<String>();
		for(int i=0;i<arr.length;i++){
			temp.add(arr[i]);
		}
		return temp;
	}
	
	public static String[] list2stringarray(List<String> arr){
		String[] temp=new String[arr.size()];
		for(int i=0;i<arr.size();i++){
			temp[i]=arr.get(i);
		}
		return temp;
	}

	public static void sort_listtable(List<List<String>> list,int sortcolumn){
		final int tempcolumn=sortcolumn;
		Collections.sort(list,new Comparator<List<String>>(){
			public int compare(List<String> o1,List<String> o2){
				return o1.get(tempcolumn).compareTo(o2.get(tempcolumn));
				//String temp1=o1.get(tempcolumn).trim().toUpperCase();
				//String temp2=o2.get(tempcolumn).trim().toUpperCase();
				//return temp1.compareTo(temp2);
			}
		});
	}
	
	public static void sort_listtable(List<List<String>> list,int sortcolumn,boolean numeric){
		if(!numeric){sort_listtable(list,sortcolumn); return;}
		final int tempcolumn=sortcolumn;
		Collections.sort(list,new Comparator<List<String>>(){
			public int compare(List<String> o1,List<String> o2){
				float temp1=Float.parseFloat(o1.get(tempcolumn));
				float temp2=Float.parseFloat(o2.get(tempcolumn));
				return Float.compare(temp1,temp2);
			}
		});
	}
	
	public static void reverse_listtable(List<List<String>> list){
		Collections.reverse(list);
	}

	public static boolean filter_listtable(List<List<String>> list,float[][] filters){
		// here filters have three values--filter column,filter type, and filter
		// value
		List<List<String>> temp=new ArrayList<List<String>>();
		for(int i=0;i<filters.length;i++){
			int column=(int)filters[i][0];
			int type=(int)filters[i][1];
			for(int j=0;j<list.size();j++){
				List<String> row=list.get(j);
				String value=row.get(column);
				float fval;
				try{
					fval=Float.parseFloat(value);
				}catch(NumberFormatException e){
					return false;
				}
				boolean filtered=false;
				if(type==0){
					filtered=(fval>=filters[i][2]);
				}
				if(type==1){
					filtered=(fval!=filters[i][2]);
				}
				if(type==2){
					filtered=(fval<=filters[i][2]);
				}
				if(!filtered)
					temp.add(row);
			}
		}
		list=temp;
		return true;
	}

	public static List<List<String>> get_filtered_listtable(List<List<String>> list,float[][] filters){
		// here filters have three values--filter column,filter type, and filter
		// value
		List<List<String>> temp=new ArrayList<List<String>>();
		for(int i=0;i<list.size();i++){
			List<String> row=list.get(i);
			boolean filtered=false;
			for(int j=0;j<filters.length;j++){
				int column=(int)filters[j][0];
				int type=(int)filters[j][1];
				String value=row.get(column);
				float fval=0.0f;
				try{
					fval=Float.parseFloat(value);
				}catch(NumberFormatException e){
					filtered=true; break;
				}
				if(!filtered){
					if(type==0) filtered=(fval>=filters[j][2]);
					if(type==1) filtered=(fval!=filters[j][2]);
					if(type==2) filtered=(fval<=filters[j][2]);
				}
				if(filtered) break;
			}
			if(!filtered) temp.add(row);
		}
		return temp;
	}

	public static List<List<String>> get_cell_stat_list(List<List<String>> list,int cellcolumn,String stat){
		return get_cell_stat_list(list,cellcolumn,stat,false);
	}
	
	public static List<List<String>> get_cell_stat_list(List<List<String>> list,int cellcolumn,String stat,boolean addct){
		return get_cell_stat_list(list,cellcolumn,stat,addct,null);
	}
	
	public static List<List<String>> get_cell_stat_list(List<List<String>> list,int cellcolumn,String stat,boolean addct,float[] options){
		List<String> celllist=get_cell_list(list,cellcolumn);
		List<List<String>> rettable=new ArrayList<List<String>>();
		for(int i=0;i<celllist.size();i++){
			List<String> cellstat=get_cell_stat(list,cellcolumn,celllist.get(i),stat,addct,options);
			rettable.add(cellstat);
		}
		return rettable;
	}

	public static List<String> get_cell_list(List<List<String>> list,int cellcolumn){
		List<String> celllist=new ArrayList<String>();
		String currcell=list.get(0).get(cellcolumn);
		celllist.add(currcell);
		for(int i=1;i<list.size();i++){
			String temp=list.get(i).get(cellcolumn);
			if(!temp.equals(currcell)){
				celllist.add(temp);
				currcell=temp;
			}
		}
		return celllist;
	}
	
	public static List<List<String>> get_missing_cells(List<List<String>> listtable,int cellcolumn,List<String> celllist,String emptystring){
		//the listtable must be sorted for this
		List<List<String>> extralist=new ArrayList<List<String>>();
		int ncols=listtable.get(0).size();
		for(int i=0;i<celllist.size();i++){
			String cellname=celllist.get(i);
			int foundindex=find_sorted_listtable_string(listtable,cellcolumn,cellname);
			if(foundindex<0){
				List<String> temp=new ArrayList<String>();
				for(int j=0;j<ncols;j++){
					if(j==cellcolumn) temp.add(cellname);
					else temp.add(emptystring);
				}
				extralist.add(temp);
				//stattable.add(temp);
			}
		}
		return extralist;
	}

	public static List<List<String>> get_cell_listtable(List<List<String>> list,int cellcolumn,String cellid){
		// first find the starting row
		boolean found=false;
		int startindex=0;
		for(int i=0;i<list.size();i++){
			String cellcolid=list.get(i).get(cellcolumn);
			if(cellcolid.equals(cellid)){
				startindex=i;
				found=true;
				break;
			}
		}
		if(!found){
			return null;
		}
		// now find the ending row
		int endindex=startindex;
		found=false;
		for(int i=startindex+1;i<list.size();i++){
			String cellcolid=list.get(i).get(cellcolumn);
			if(!cellcolid.equals(cellid)){
				endindex=i-1;
				found=true;
				break;
			}
		}
		if(!found){
			endindex=list.size()-1;
		}
		// now return the sublist
		return list.subList(startindex,endindex+1);
	}

	public static List<List<List<String>>> get_cells_listtable(List<List<String>> table,int cellcolumn,int minsize){
		// here nest each "cell" or contiguous region of the cell column into a
		// sublist
		String currid=table.get(0).get(cellcolumn);
		List<List<List<String>>> table2=new ArrayList<List<List<String>>>();
		List<List<String>> temptable=new ArrayList<List<String>>();
		temptable.add(table.get(0));
		for(int i=1;i<table.size();i++){
			if(!table.get(i).get(cellcolumn).equals(currid)){
				currid=table.get(i).get(cellcolumn);
				if(temptable.size()>minsize){
					table2.add(temptable);
				}
				temptable=new ArrayList<List<String>>();
			}
			temptable.add(table.get(i));
		}
		if(temptable.size()>minsize){
			table2.add(temptable);
		}
		return table2;
	}
	
	public static List<List<String>> cells2listtable(List<List<List<String>>> cellstable){
		List<List<String>> temptable=new ArrayList<List<String>>();
		for(int i=0;i<cellstable.size();i++){
			List<List<String>> celltable=cellstable.get(i);
			for(int j=0;j<celltable.size();j++){
				temptable.add(celltable.get(j));
			}
		}
		return temptable;
	}

	public static List<String> get_cell_stat(List<List<String>> list,int cellcolumn,String cellid,String stat){
		return get_cell_stat(list,cellcolumn,cellid,stat,false);
	}
	
	public static List<String> get_cell_stat(List<List<String>> list,int cellcolumn,String cellid,String stat,boolean addct){
		return get_cell_stat(list,cellcolumn,cellid,stat,false,null);
	}
	
	public static List<String> get_cell_stat(List<List<String>> list,int cellcolumn,String cellid,String stat,boolean addct,float[] options){
		List<List<String>> celltable=get_cell_listtable(list,cellcolumn,cellid);
		List<String> stats=table_tools.get_table_stat(celltable,stat,options);
		stats.set(cellcolumn,cellid);
		if(addct) stats.add(""+celltable.size());
		return stats;
	}

	public static List<String> get_table_stat(List<List<String>> list,String stat){
		return get_table_stat(list,stat,null);
	}
	
	public static List<String> get_table_stat(List<List<String>> list,String stat,float[] options){
		int ncolumns=list.get(0).size();
		List<String> stattable=new ArrayList<String>();
		for(int i=0;i<ncolumns;i++){
			if(is_number(list.get(0).get(i))){
				float[] data=new float[list.size()];
				for(int j=0;j<list.size();j++){
					String temp=list.get(j).get(i);
					if(is_number(temp)) data[j]=Float.parseFloat(list.get(j).get(i));
					else data[j]=Float.NaN;
				}
				float[] tempoptions=null;
				if(options!=null) tempoptions=options.clone();
				float temp=jstatistics.getstatistic(stat,data,tempoptions);
				stattable.add(""+temp);
			}else{
				stattable.add(""+list.get(0).get(i));
			}
		}
		return stattable;
	}
	
	public static float[] get_row_array(List<List<String>> list,int row){
		List<String> temp=list.get(row);
		float[] temp2=new float[temp.size()];
		for(int i=0;i<temp.size();i++){
			String temp3=temp.get(i);
			if(is_number(temp3)){
				temp2[i]=Float.parseFloat(temp3);
			} else {
				temp2[i]=Float.NaN;
			}
		}
		return temp2;
	}

	public static float[] get_column_array(List<List<String>> list,int col){
		float[] temp=new float[list.size()];
		for(int i=0;i<temp.length;i++){
			if(col<list.get(i).size()){
				String temp2=list.get(i).get(col);
				if(is_number(temp2)){
					temp[i]=Float.parseFloat(temp2);
				}else{
					temp[i]=Float.NaN;
				}
			}else{
				temp[i]=Float.NaN;
			}
		}
		return temp;
	}
	
	public static float[][] get_matrix(List<List<String>> list){
		float[][] temp=new float[list.size()][];
		for(int i=0;i<list.size();i++) temp[i]=get_column_array(list,i);
		return temp;
	}
	
	public static String[] get_listtable_column(List<List<String>> list,int col){
		String[] temp=new String[list.size()];
		for(int i=0;i<temp.length;i++){
			if(col<list.get(i).size()){
				temp[i]=list.get(i).get(col);
			}else{
				temp[i]="";
			}
		}
		return temp;
	}

	public static void add_listtable_column(List<List<String>> list,String[] column,int pos){
		for(int i=0;i<list.size();i++){
			List<String> row=list.get(i);
			int pos2=pos;
			if(pos2>row.size())
				pos2=row.size();
			row.add(pos2,column[i]);
		}
	}
	
	public static void add_listtable_column(List<List<String>> list,float[] column,int pos){
		for(int i=0;i<list.size();i++){
			List<String> row=list.get(i);
			int pos2=pos;
			if(pos2>row.size())
				pos2=row.size();
			row.add(pos2,""+column[i]);
		}
	}

	public static void delete_listtable_column(List<List<String>> list,int delcolumn){
		for(int i=0;i<list.size();i++){
			list.get(i).remove(delcolumn);
		}
	}

	public static void delete_listtable_row(List<List<String>> list,int delrow){
		list.remove(delrow);
	}

	public static String delete_col_label(String headings,int delcol){
		String[] labels=split_string_tab(headings);
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<delcol;i++){
			sb.append(labels[i]+"\t");
		}
		for(int i=delcol+1;i<labels.length;i++){
			sb.append(labels[i]+"\t");
		}
		return sb.substring(0,sb.length()-1);
	}

	/*****************
	 * prints a 2D arraylist with tab, comma, or space delimiting
	 * @param list
	 * @param delim: 0=tab,1=comma, 2=space
	 * @return
	 */
	public static String print_listtable(List<List<String>> list,int delim){
		// delims are tab,comma,space
		StringBuffer retvals=new StringBuffer();
		for(int i=0;i<list.size();i++){
			if(i!=0){
				retvals.append("\n");
			}
			retvals.append(print_string_array(list.get(i),delim));
		}
		return retvals.toString();
	}

	public static String print_listtable(List<List<String>> list){
		StringBuffer retvals=new StringBuffer();
		for(int i=0;i<list.size();i++){
			if(i!=0){
				retvals.append("\n");
			}
			retvals.append(print_string_array(list.get(i)));
		}
		return retvals.toString();
	}
	
	public static void replace_table(TextPanel tp,List<List<String>> list,String[] collabels){
		//set the new column labels and clear the table
		tp.setColumnHeadings(print_string_array(collabels));
		tp.append(print_listtable(list));
	}
	
	public static void replace_table(TextPanel tp,List<List<String>> list,String collabels){
		//set the new column labels and clear the table
		tp.setColumnHeadings(collabels);
		tp.append(print_listtable(list));
	}
	
	public static void change_table_labels(TextPanel tp,String[] newlabels){
		change_table_labels(tp,print_string_array(newlabels));
	}
	
	public static void change_table_labels(TextPanel tp,String newlabels){
		if(!newlabels.equals(tp.getColumnHeadings())){
			List<List<String>> listtable=table2listtable(tp);
			String[] headings2=split_string_tab(newlabels);
			String[] oldheadings=getcollabels(tp);
			if(headings2.length>oldheadings.length){
				List<String> row1=listtable.get(0);
				int start=row1.size();
				for(int i=start;i<headings2.length;i++) row1.add("");
			}
			replace_table(tp,listtable,newlabels);
		}
	}

	public static void create_table(String title,List<List<String>> list,String[] collabels){
		int ncolumns=0;
		if(list.size()>0) ncolumns=list.get(0).size();
		if(ncolumns==0 && collabels!=null) ncolumns=collabels.length; //blank table
		String[] labels=collabels;
		if(collabels==null){
			if(ncolumns==0) return; //empty file
			labels=new String[ncolumns];
			for(int i=0;i<ncolumns;i++)
				labels[i]="col"+(i+1);
		}
		String title2=title;
		if(title==null){
			title2="Table";
		}
		new TextWindow(title2,print_string_array(labels),print_listtable(list),400,200);
	}
	
	public static void create_table(String title,List<List<String>> list){
		//here the column labels are in the first row
		//warning: the first row is removed by this subroutine
		String[] collabels=list2stringarray(list.get(0));
		list.remove(0);
		create_table(title,list,collabels);
		int ncolumns=list.get(0).size();
	}

	public static void create_table(String title,float[][] list,String[] collabels){
		int ncolumns=list[0].length;
		String[] labels=collabels;
		if(collabels==null){
			labels=new String[ncolumns];
			for(int i=0;i<ncolumns;i++)
				labels[i]="col"+(i+1);
		}
		String title2=title;
		if(title==null){
			title2="Table";
		}
		new TextWindow(title2,print_string_array(labels),print_float_array(list),400,200);
	}
	
	public static void create_table(String title,double[][] list,String[] collabels){
		int ncolumns=list[0].length;
		String[] labels=collabels;
		if(collabels==null){
			labels=new String[ncolumns];
			for(int i=0;i<ncolumns;i++)
				labels[i]="col"+(i+1);
		}
		String title2=title;
		if(title==null){
			title2="Table";
		}
		new TextWindow(title2,print_string_array(labels),print_double_array(list),400,200);
	}
	
	public static void create_table(String title,int[][] list,String[] collabels){
		int ncolumns=list[0].length;
		String[] labels=collabels;
		if(collabels==null){
			labels=new String[ncolumns];
			for(int i=0;i<ncolumns;i++)
				labels[i]="col"+(i+1);
		}
		String title2=title;
		if(title==null){
			title2="Table";
		}
		new TextWindow(title2,print_string_array(labels),print_int_array(list),400,200);
	}
	
	public static void create_table(String title,String[][] list,String[] collabels){
		int ncolumns=list[0].length;
		String[] labels=collabels;
		if(collabels==null){
			labels=new String[ncolumns];
			for(int i=0;i<ncolumns;i++)
				labels[i]="col"+(i+1);
		}
		String title2=title;
		if(title==null){
			title2="Table";
		}
		new TextWindow(title2,print_string_array(labels),print_string_array(list),400,200);
	}

	public static String print_string_array(String[] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(data[0]);
		for(int i=1;i<data.length;i++){
			retvals.append("\t"+data[i]);
		}
		return retvals.toString();
	}
	
	public static String print_string_array(String[][] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(print_string_array(data[0]));
		for(int i=1;i<data.length;i++){
			retvals.append("\n"+print_string_array(data[i]));
		}
		return retvals.toString();
	}

	public static String print_string_array(List<String> data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(data.get(0));
		for(int i=1;i<data.size();i++){
			retvals.append("\t"+data.get(i));
		}
		return retvals.toString();
	}

	/*****************
	 * prints a 1D arraylist with tab, comma, or space delimiting
	 * @param data
	 * @param delim: 0=tab,1=comma, 2=space
	 * @return
	 */
	public static String print_string_array(List<String> data,int delim){
		// delims are tab,comma,space
		StringBuffer retvals=new StringBuffer();
		retvals.append(data.get(0));
		String sep=get_delim(delim);
		for(int i=1;i<data.size();i++){
			retvals.append(sep+data.get(i));
		}
		return retvals.toString();
	}
	
	public static String[] delims= {"tab","comma","space","newline"};
	
	public static String get_delim(int delim){
		String sep="\t";
		if(delim==1) sep=",";
		if(delim==2) sep=" ";
		if(delim==3) sep="\n";
		return sep;
	}

	/*****************
	 * prints a 1D arraylist with tab, comma, or space delimiting
	 * @param data
	 * @param delim: 0=tab,1=comma, 2=space, 3=newline
	 * @return
	 */
	public static String print_string_array(String[] data,int delim){
		// delims are tab,comma,space
		StringBuffer retvals=new StringBuffer();
		retvals.append(data[0]);
		String sep=get_delim(delim);
		for(int i=1;i<data.length;i++){
			retvals.append(sep+data[i]);
		}
		return retvals.toString();
	}

	public static String print_int_array(int[][] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(print_int_array(data[0]));
		for(int i=1;i<data.length;i++){
			retvals.append("\n"+print_int_array(data[i]));
		}
		return retvals.toString();
	}

	public static String print_int_array(int[] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(""+data[0]);
		for(int i=1;i<data.length;i++){
			retvals.append("\t"+data[i]);
		}
		return retvals.toString();
	}

	public static String print_float_array(float[][] data,boolean addnumbers){
		StringBuffer retvals=new StringBuffer();
		if(addnumbers)
			retvals.append("1\t");
		retvals.append(print_float_array(data[0]));
		for(int i=1;i<data.length;i++){
			retvals.append("\n");
			if(addnumbers)
				retvals.append(""+(i+1)+"\t");
			retvals.append(print_float_array(data[i]));
		}
		return retvals.toString();
	}

	public static String print_float_array(float[][] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(print_float_array(data[0]));
		for(int i=1;i<data.length;i++){
			retvals.append("\n"+print_float_array(data[i]));
		}
		return retvals.toString();
	}
	
	public static String print_float_array(float[][] data,int delim){
		StringBuffer retvals=new StringBuffer();
		retvals.append(print_float_array(data[0],delim));
		for(int i=1;i<data.length;i++){
			retvals.append("\n"+print_float_array(data[i],delim));
		}
		return retvals.toString();
	}

	public static String print_float_array(float[] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(""+data[0]);
		for(int i=1;i<data.length;i++){
			retvals.append("\t"+data[i]);
		}
		return retvals.toString();
	}
	
	public static String print_float_array(float[] data,int delim){
		StringBuffer retvals=new StringBuffer();
		retvals.append(""+data[0]);
		String sep=get_delim(delim);
		for(int i=1;i<data.length;i++){
			retvals.append(sep+data[i]);
		}
		return retvals.toString();
	}

	public static String print_double_array(double[][] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(print_double_array(data[0]));
		for(int i=1;i<data.length;i++){
			retvals.append("\n"+print_double_array(data[i]));
		}
		return retvals.toString();
	}

	public static String print_double_array(double[] data){
		StringBuffer retvals=new StringBuffer();
		retvals.append(""+(float)data[0]);
		for(int i=1;i<data.length;i++){
			retvals.append("\t"+(float)data[i]);
		}
		return retvals.toString();
	}

	public static int find_string_array_index(String[] starr,String target){
		for(int i=0;i<starr.length;i++){
			if(starr[i].equals(target)){
				return i;
			}
		}
		return -1;
	}

	public static int find_sorted_listtable_string(List<List<String>> list,int column,String target){
		// note that the listtable must be sorted for this to work
		List<String> key=new ArrayList<String>();
		int size=list.get(0).size();
		for(int i=0;i<size;i++){
			key.add("");
		}
		key.set(column,target);
		final int tempcolumn=column;
		int index=Collections.binarySearch(list,key,new Comparator<List<String>>(){
			public int compare(List<String> o1,List<String> o2){
				return o1.get(tempcolumn).compareTo(o2.get(tempcolumn));
				//String temp1=o1.get(tempcolumn).trim().toUpperCase();
				//String temp2=o2.get(tempcolumn).trim().toUpperCase();
				//return temp1.compareTo(temp2);
			}
		});
		return index;
	}
}
