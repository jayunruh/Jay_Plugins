package jguis;

import ij.IJ;
import jalgs.jdataio;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.List;

public class export_flowcyte{
	private static String labeldelim="|\u0024";
	private static String delim2="|";
	public int textheadlength=20048;
	
	public boolean write_table(List<List<String>> table,String[] labels,String dir,String fname){
		int npts=table.size();
		float[][] data=new float[npts][];
		for(int i=0;i<npts;i++){
			data[i]=table_tools.get_row_array(table,i);
		}
		return write_table(data,labels,dir,fname);
	}
	
	public boolean write_table(float[][] data,String[] labels,String dir,String fname){
		return write_table(data,labels,dir+fname);
	}
	
	public boolean write_table(float[][] data,String[] labels,String path) {
		String[][] meta=null;
		return write_table(data,labels,path,meta);
	}
	
	public boolean write_table(float[][] data,String[] labels,String path,String[][] meta){
		int nch=labels.length;
		int npts=data.length;
		/*if(meta!=null) {
			for(int i=0;i<meta.length;i++) {
				textheadlength+=3+meta[i][0].length()+meta[i][1].length();
			}
		}*/
		String fname=(new File(path)).getName();
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(path));
			//now write the header
			int offset=0;
			String header=get_FCS_header(nch,npts);
			os.write(header.getBytes()); offset+=header.length();
			//now the TEXT header
			String texthead=get_TEXT_header(nch,npts,fname,meta);
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
			return false;
		}
		return true;
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
		sb.append(pad_integer(textheadlength-1));
		sb.append(pad_integer(textheadlength));
		int totpts=nch*npts*4;
		sb.append(pad_integer(textheadlength+totpts-1));
		sb.append(pad_integer(0));
		sb.append(pad_integer(0));
		return sb.toString();
	}

	private String get_TEXT_header(int nch,int npts,String filename,String[][] meta){
		StringBuffer sb=new StringBuffer();
		sb.append(labeldelim+"TOT"+delim2+npts);
		sb.append(labeldelim+"PAR"+delim2+nch);
		sb.append(labeldelim+"BYTEORD"+delim2+"1,2,3,4");
		sb.append(labeldelim+"DATATYPE"+delim2+"F");
		sb.append(labeldelim+"MODE"+delim2+"L");
		sb.append(labeldelim+"BEGINDATA"+delim2+textheadlength);
		int dataend=nch*npts*8+textheadlength;
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
		if(meta!=null) {
    		for(int i=0;i<meta[0].length;i++) {
    			sb.append(labeldelim+meta[0][i]+delim2+meta[1][i]);
    		}
		}
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
