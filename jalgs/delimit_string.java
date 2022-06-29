/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.util.ArrayList;
import java.util.List;

public class delimit_string{
	public char delimiter;

	// this class contains utility classes to delimit strings

	public delimit_string(char delimiter1){
		delimiter=delimiter1;
	}

	public int count_lines(String input){
		int temp=-1;
		int counter=0;
		do{
			temp=input.indexOf('\n',temp+1);
			counter++;
		}while(temp>=0);
		if(input.charAt(input.length()-1)=='\n'){
			return counter-1;
		}else{
			return counter;
		}
	}

	public String[] getrows(String input){
		String[] output=input.split("\n");
		if(output[output.length-1]==null){
			String[] output2=new String[output.length-1];
			for(int i=0;i<(output.length-1);i++){
				output2[i]=output[i].substring(0);
			}
			return output2;
		}else{
			return output;
		}
	}

	public String[] split_string(String line){
		if(line.endsWith(""+delimiter)){
			return line.substring(0,line.length()-1).split(""+delimiter);
		}else{
			return line.split(""+delimiter);
		}
	}

	public int getnumcolumns(String line){
		int temp=-1;
		int counter=0;
		do{
			temp=line.indexOf(delimiter,temp+1);
			counter++;
		}while(temp>=0);
		if(line.charAt(line.length()-1)==delimiter){
			return counter-1;
		}else{
			return counter;
		}
	}

	public String[] delim2string(String line,int numcolumns){
		String[] templine=line.split(""+delimiter);
		if(templine.length>numcolumns){
			String[] templine2=new String[numcolumns];
			for(int i=0;i<numcolumns;i++){
				templine2[i]=templine[i].substring(0);
			}
			return templine2;
		}else{
			if(templine.length==numcolumns){
				return templine;
			}else{
				String[] templine3=new String[numcolumns];
				for(int i=0;i<templine.length;i++){
					templine3[i]=templine[i].substring(0);
				}
				return templine3;
			}
		}
	}

	public List<String> delim2list(String line,int numcolumns){
		String[] temp=delim2string(line,numcolumns);
		List<String> list=new ArrayList<String>();
		for(int i=0;i<temp.length;i++)
			list.add(temp[i]);
		return list;
	}

	public float[] delim2float(String line,int numcolumns){
		String[] delim_string=delim2string(line,numcolumns);
		float[] data=new float[numcolumns];
		for(int i=0;i<delim_string.length;i++){
			try{
				data[i]=Float.parseFloat(delim_string[i]);
			}catch(Exception e){
				data[i]=Float.NaN;
			}
		}
		return data;
	}

	public int[] delim2int(String line,int numcolumns){
		String[] delim_string=delim2string(line,numcolumns);
		int[] data=new int[numcolumns];
		for(int i=0;i<delim_string.length;i++){
			try{
				data[i]=Integer.parseInt(delim_string[i]);
			}catch(NumberFormatException e){
				data[i]=0;
			}
		}
		return data;
	}
}
