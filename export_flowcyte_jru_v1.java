/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
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
import java.util.*;
import ij.io.*;
import java.io.*;
import jalgs.*;
import jguis.*;
import ij.text.*;

public class export_flowcyte_jru_v1 implements PlugIn {
	private static String labeldelim="|\u0024";
	private static String delim2="|";

	public void run(String arg) {
		TextWindow[] tw=jutils.selectTables(false,1);
		if(tw==null) return;
		String title=tw[0].getTitle();
		if(!title.endsWith(".fcs")) title+=".fcs";
		SaveDialog sd=new SaveDialog("Save FCS File",title,".fcs");
		String dir=sd.getDirectory();
		String fname=sd.getFileName();
		TextPanel tp=tw[0].getTextPanel();
		IJ.log(""+fname);
		List<List<String>> listtable=table_tools.table2listtable(tp);
		String[] labels=table_tools.getcollabels(tp);
		//labels=table_tools.make_labels_unique(labels);
		(new export_flowcyte()).write_table(listtable,labels,dir,fname);
		/*int nch=labels.length;
		int npts=listtable.size();
		float[][] data=new float[npts][];
		for(int i=0;i<npts;i++){
			data[i]=table_tools.get_row_array(listtable,i);
		}
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(dir+fname));
			//now write the header
			int offset=0;
			String header=get_FCS_header(nch,npts);
			os.write(header.getBytes()); offset+=header.length();
			//now the TEXT header
			String texthead=get_TEXT_header(nch,npts,fname);
			os.write(texthead.getBytes()); offset+=texthead.length();
			//now the parameter header
			int[] ranges=get_ranges(data,nch);
			String parambuf=get_PARAM_buffer(labels,ranges);
			os.write(parambuf.getBytes()); offset+=parambuf.length();
			//pad with spaces
			for(int i=offset;i<20048;i++) os.write((byte)32);
			//now the table data
			(new jdataio()).writeintelfloatarray(os,data);
			os.close();
		} catch(IOException e){
			IJ.log(e.getMessage());
		}*/
	}

	private int[] get_ranges(float[][] data,int nch){
		int[] ranges=new int[nch];
		for(int i=0;i<data.length;i++){
			for(int j=0;j<data[i].length;j++){
				if(!Float.isInfinite(data[i][j])){
					if(!Float.isNaN(data[i][j])){
						if((int)data[i][j]>ranges[j]) ranges[j]=(int)data[i][j];
					}
				}
			}
		}
		for(int i=0;i<nch;i++) ranges[i]+=1;
		return ranges;				
	}

	private String get_FCS_header(int nch,int npts){
		StringBuffer sb=new StringBuffer();
		sb.append("FCS3.0    ");
		sb.append(pad_integer(58));
		sb.append(pad_integer(20047));
		sb.append(pad_integer(20048));
		int totpts=nch*npts*4;
		sb.append(pad_integer(20048+totpts-1));
		sb.append(pad_integer(0));
		sb.append(pad_integer(0));
		return sb.toString();
	}

	private String get_TEXT_header(int nch,int npts,String filename){
		StringBuffer sb=new StringBuffer();
		sb.append(labeldelim+"TOT"+delim2+npts);
		sb.append(labeldelim+"PAR"+delim2+nch);
		sb.append(labeldelim+"BYTEORD"+delim2+"1,2,3,4");
		sb.append(labeldelim+"DATATYPE"+delim2+"F");
		sb.append(labeldelim+"MODE"+delim2+"L");
		sb.append(labeldelim+"BEGINDATA"+delim2+"20048");
		int dataend=nch*npts*8+20048;
		sb.append(labeldelim+"ENDDATA"+delim2+dataend);
		sb.append(labeldelim+"BEGINANALYSIS"+delim2+"0");
		sb.append(labeldelim+"ENDANALYSIS"+delim2+"0");
		sb.append(labeldelim+"BEGINSTEXT"+delim2+"58");
		sb.append(labeldelim+"ENDSTEXT"+delim2+"20048");
		sb.append(labeldelim+"COM"+delim2+"Jay Unruh Plugins");
		sb.append(labeldelim+"FIL"+delim2+filename);
		sb.append(labeldelim+"NEXTDATA"+delim2+"0");
		Calendar calendar = new GregorianCalendar();
		Date today = new Date();
		calendar.setTime(today);
		String date=""+(calendar.get(Calendar.MONTH)+1)+"-"+calendar.get(Calendar.DAY_OF_MONTH)+"-"+calendar.get(Calendar.YEAR);
		sb.append(labeldelim+"DATE"+delim2+date);
		sb.append(labeldelim+"CYT"+delim2+"ImageJ_Table");
		//sb.append(labeldelim+"WELL ID"+delim2+wellID);
		//sb.append(labeldelim+"PLATE ID"+delim2+plateID);
		return sb.toString();
	}

	private String get_PARAM_buffer(String[] labels,int[] ranges){
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<labels.length;i++){
			sb.append(labeldelim+"P"+(i+1)+"S"+delim2+labels[i]);
			//sb.append(labeldelim+"P"+(i+1)+"N"+delim2+"P"+(i+1));
			sb.append(labeldelim+"P"+(i+1)+"N"+delim2+labels[i]);
			sb.append(labeldelim+"P"+(i+1)+"E"+delim2+"0,0");
			sb.append(labeldelim+"P"+(i+1)+"G"+delim2+"1");
			sb.append(labeldelim+"P"+(i+1)+"R"+delim2+ranges[i]);
			sb.append(labeldelim+"P"+(i+1)+"B"+delim2+"32");
		}
		sb.append(delim2);
		return sb.toString();
	}

	private String pad_integer(int val){
		if(val<0) return null;
		if(val<10) return "       "+val;
		if(val<100) return "      "+val;
		if(val<1000) return "     "+val;
		if(val<10000) return "    "+val;
		if(val<100000) return "   "+val;
		if(val<1000000) return "  "+val;
		if(val<10000000) return " "+val;
		if(val<100000000) return ""+val;
		return null;
	}

}
