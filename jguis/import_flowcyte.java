package jguis;

import java.io.InputStream;

import ij.IJ;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;

import jalgs.jdataio;

public class import_flowcyte{
	
	/**************
	 * this version just gets the header info
	 * @param path
	 * @return
	 */
	public Object[] getChNames(String path) {
		try{
			jdataio jdio=new jdataio();
			//first read the header
			InputStream instream=new BufferedInputStream(new FileInputStream(path));
			String label=jdio.readstring(instream,6);
			jdio.skipstreambytes(instream,4);
			String toff=jdio.readstring(instream,8); int textoff=Integer.parseInt(toff.trim()); //start of TEXT segment
			String teoff=jdio.readstring(instream,8); int texteoff=Integer.parseInt(teoff.trim()); //end of TEXT segment
			//IJ.log(""+textoff+" , "+texteoff+" , "+dataoff+" , "+dataeoff+" , "+analoff+" , "+analeoff);
			instream.close();
			instream=new BufferedInputStream(new FileInputStream(path));
			jdio.skipstreambytes(instream,textoff);
			//now read all parameters and values
			int textlength=texteoff-textoff;
			String text=jdio.readstring(instream,textlength);
			jdio.skipstreambytes(instream,1);
			text=text.replace('\\','/');
			text=text.replace('|','/');
			text=text.replace("\u000C","/");
			text=text.replace("\u001E","/");
			text=text.replace("*","/");
			String[] params=text.split("/\\u0024");
			//note that sometimes $ (x24) occurs before the keyword and sometimes it doesn't
			//String[] params=text.split("/\u0c24");
			//IJ.log(text);
			String[] labels=new String[params.length-1];
			String[] values=new String[params.length-1];
			boolean[] flags=new boolean[params.length-1];
			int numextra=0;
			for(int i=1;i<params.length;i++){
				String[] temp=params[i].split("/");
				labels[i-1]=temp[0]; values[i-1]=temp[1];
				//if temp is longer than 2 values, need to add those values
				if(temp.length>2) {
					flags[i-1]=true;
					numextra+=(temp.length-2)/2;
				}
				//if(showmetadata) IJ.log(labels[i-1]+" , "+values[i-1]);
			}
			//now add on the extras if they exist
			if(numextra>0) {
    			String[] newlabels=new String[labels.length+numextra];
    			String[] newvalues=new String[labels.length+numextra];
    			int counter=flags.length;
    			for(int i=0;i<flags.length;i++) {
    				newlabels[i]=labels[i];
    				newvalues[i]=values[i];
    				if(flags[i]) {
    					String[] temp=params[i+1].split("/");
    					int temp1=(temp.length-2)/2;
    					for(int j=0;j<temp1;j++) {
    						newlabels[counter]=temp[2*j+2];
    						newvalues[counter]=temp[2*j+3];
    						counter++;
    					}
    				}
    			}
    			labels=newlabels;
    			values=newvalues;
			}
			String[][] meta={labels,values};
			//now get the number of data points and number of channels
			int npts=get_label_number(labels,values,"TOT");
			int nch=get_label_number(labels,values,"PAR");
			String byteord=get_label_value(labels,values,"BYTEORD");
			//IJ.log(byteord);
			boolean motorola=true;
			if(byteord.startsWith("1,")) motorola=false;
			String dtype=get_label_value(labels,values,"DATATYPE");
			boolean isints=dtype.startsWith("I");
			//now go through each channel getting its name, bits, and range
			String[] chnames=new String[nch];
			int[] chbits=new int[nch];
			int[] range=new int[nch];
			for(int i=0;i<nch;i++){
				chnames[i]=get_label_value(labels,values,"P"+(i+1)+"N");
				if(chnames[i]==null) chnames[i]="P"+(i+1)+"N";
				chbits[i]=get_label_number(labels,values,"P"+(i+1)+"B");
				if(chbits[i]<0) chbits[i]=32;
				range[i]=get_label_number(labels,values,"P"+(i+1)+"R");
			}
			instream.close();
			return new Object[] {chnames,meta};
		} catch(IOException e){
			IJ.log(e.getMessage());
			return null;
		}
	}

	public Object[] getFCSFile(String path){
		try{
			jdataio jdio=new jdataio();
			//first read the header
			InputStream instream=new BufferedInputStream(new FileInputStream(path));
			String label=jdio.readstring(instream,6);
			jdio.skipstreambytes(instream,4);
			String toff=jdio.readstring(instream,8); int textoff=Integer.parseInt(toff.trim()); //start of TEXT segment
			String teoff=jdio.readstring(instream,8); int texteoff=Integer.parseInt(teoff.trim()); //end of TEXT segment
			String doff=jdio.readstring(instream,8); int dataoff=Integer.parseInt(doff.trim()); //start of DATA segment
			String deoff=jdio.readstring(instream,8); int dataeoff=Integer.parseInt(deoff.trim()); //end of DATA segment
			String aoff=jdio.readstring(instream,8); int analoff=Integer.parseInt(aoff.trim()); //start of ANALYSIS segment
			String aeoff=jdio.readstring(instream,8); int analeoff=Integer.parseInt(aeoff.trim()); //end of ANALYSIS segment
			//IJ.log(""+textoff+" , "+texteoff+" , "+dataoff+" , "+dataeoff+" , "+analoff+" , "+analeoff);
			instream.close();
			instream=new BufferedInputStream(new FileInputStream(path));
			jdio.skipstreambytes(instream,textoff);
			//now read all parameters and values
			int textlength=texteoff-textoff;
			String text=jdio.readstring(instream,textlength);
			jdio.skipstreambytes(instream,1);
			text=text.replace('\\','/');
			text=text.replace('|','/');
			text=text.replace("\u000C","/");
			text=text.replace("\u001E","/");
			text=text.replace("*","/");
			String[] params=text.split("/\\u0024");
			//note that sometimes $ (x24) occurs before the keyword and sometimes it doesn't
			//String[] params=text.split("/\u0c24");
			//IJ.log(text);
			String[] labels=new String[params.length-1];
			String[] values=new String[params.length-1];
			boolean[] flags=new boolean[params.length-1];
			int numextra=0;
			for(int i=1;i<params.length;i++){
				String[] temp=params[i].split("/");
				labels[i-1]=temp[0]; values[i-1]=temp[1];
				//if temp is longer than 2 values, need to add those values
				if(temp.length>2) {
					flags[i-1]=true;
					numextra+=(temp.length-2)/2;
				}
				//if(showmetadata) IJ.log(labels[i-1]+" , "+values[i-1]);
			}
			//now add on the extras if they exist
			if(numextra>0) {
    			String[] newlabels=new String[labels.length+numextra];
    			String[] newvalues=new String[labels.length+numextra];
    			int counter=flags.length;
    			for(int i=0;i<flags.length;i++) {
    				newlabels[i]=labels[i];
    				newvalues[i]=values[i];
    				if(flags[i]) {
    					String[] temp=params[i+1].split("/");
    					int temp1=(temp.length-2)/2;
    					for(int j=0;j<temp1;j++) {
    						newlabels[counter]=temp[2*j+2];
    						newvalues[counter]=temp[2*j+3];
    						counter++;
    					}
    				}
    			}
    			labels=newlabels;
    			values=newvalues;
			}
			String[][] meta={labels,values};
			//now get the number of data points and number of channels
			int npts=get_label_number(labels,values,"TOT");
			int nch=get_label_number(labels,values,"PAR");
			String byteord=get_label_value(labels,values,"BYTEORD");
			//IJ.log(byteord);
			boolean motorola=true;
			if(byteord.startsWith("1,")) motorola=false;
			String dtype=get_label_value(labels,values,"DATATYPE");
			boolean isints=dtype.startsWith("I");
			//now go through each channel getting its name, bits, and range
			String[] chnames=new String[nch];
			int[] chbits=new int[nch];
			int[] range=new int[nch];
			for(int i=0;i<nch;i++){
				chnames[i]=get_label_value(labels,values,"P"+(i+1)+"N");
				if(chnames[i]==null) chnames[i]="P"+(i+1)+"N";
				chbits[i]=get_label_number(labels,values,"P"+(i+1)+"B");
				if(chbits[i]<0) chbits[i]=32;
				range[i]=get_label_number(labels,values,"P"+(i+1)+"R");
			}
			//String headings=table_tools.print_string_array(chnames);
			//TextWindow tw=new TextWindow(name,headings,"",400,200);
			float[][] temp=new float[npts][nch];
			instream.close();
			if(npts<=0){
				return new Object[]{chnames,null,meta};
			}
			instream=new BufferedInputStream(new FileInputStream(path));
			jdio.skipstreambytes(instream,dataoff);
			for(int i=0;i<npts;i++){ //need to allow for integer data types
				for(int j=0;j<nch;j++){
					if(chbits[j]==16){
						if(!motorola) temp[i][j]=(float)jdio.readintelshort(instream);
						else temp[i][j]=(float)jdio.readmotorolashort(instream);
					} else {
						if(chbits[j]==32){
							if(isints) {
								if(!motorola) temp[i][j]=(float)jdio.readintelint(instream);
								else temp[i][j]=(float)jdio.readmotorolaint(instream);
							} else {
								if(!motorola) temp[i][j]=(float)jdio.readintelfloat(instream);
								else temp[i][j]=(float)jdio.readmotorolafloat(instream);
							}
						} else {
							if(chbits[j]==8){
								temp[i][j]=(float)jdio.readintelbyte(instream);
							}
						}
					}
				}
			}
			//String data=table_tools.print_float_array(temp);
			//tw.append(data);
			instream.close();
			return new Object[]{chnames,temp,meta};
		} catch(IOException e){
			IJ.log(e.getMessage());
			return null;
		}
	}

	private String get_label_value(String[] labels,String[] values,String match){
		for(int i=0;i<labels.length;i++){
			if(labels[i].equals(match)){
				return values[i];
			}
		}
		return null;
	}

	private int get_label_number(String[] labels,String[] values,String match){
		String value=get_label_value(labels,values,match);
		if(value!=null){
			try{
				return Integer.parseInt(value.trim());
			}catch(NumberFormatException e){
				return -1;
			}
		}
		return -1;
	}
	
}
