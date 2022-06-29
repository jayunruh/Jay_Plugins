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
import jguis.*;
import jalgs.*;
import java.io.*;
import java.util.*;

public class edit_KISS_codes_jru_v1 implements PlugIn {

	public void run(String arg) {
		//this plugin creates and edits KISS color codes
		//first find out which codes we want to edit
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Mouse?",false);
		gd.addChoice("Action",new String[]{"New","Edit","Delete"},"New");
		gd.showDialog(); if(gd.wasCanceled()) return;
		boolean mouse=gd.getNextBoolean();
		int action=gd.getNextChoiceIndex();
		Object[] codes=KissPanel.getCustomCodes();
		String[] mouse_names=null;
		String[] human_names=null;
		Object[] mouse_codes=null;
		Object[] human_codes=null;
		int[] human_indices=null;
		int[] mouse_indices=null;
		if(codes==null){
			human_names=new String[]{"Human_Default"};
			human_codes=new Object[]{KissPanel.humancode};
			mouse_names=new String[]{"Mouse_Default"};
			mouse_codes=new Object[]{KissPanel.mousecode};
			human_indices=new int[]{-1};
			mouse_indices=new int[]{-1};
		} else {
			String[] names=(String[])codes[0];
			boolean[] ismouse=(boolean[])codes[1];
			int nmouse=0;
			for(int i=0;i<ismouse.length;i++) if(ismouse[i]) nmouse++;
			//IJ.log("object length = "+codes.length);
			//IJ.log("ncodes = "+ismouse.length);
			//IJ.log("nmouse = "+nmouse);
			mouse_names=new String[nmouse+1];
			human_names=new String[ismouse.length-nmouse+1];
			mouse_codes=new Object[nmouse+1];
			human_codes=new Object[ismouse.length-nmouse+1];
			human_indices=new int[ismouse.length-nmouse+1]; human_indices[0]=-1;
			mouse_indices=new int[nmouse+1]; mouse_indices[0]=-1;
			mouse_names[0]="Mouse_Default";
			human_names[0]="Human_Default";
			mouse_codes[0]=KissPanel.mousecode;
			human_codes[0]=KissPanel.humancode;
			mouse_indices[0]=-1;
			human_indices[0]=-1;
			int counter1=1; int counter2=1;
			for(int i=0;i<ismouse.length;i++){
				if(ismouse[i]){
					mouse_names[counter1]=names[i];
					mouse_codes[counter1]=codes[i+2]; 
					mouse_indices[counter1]=i;
					counter1++;
				} else {
					human_names[counter2]=names[i];
					human_codes[counter2]=codes[i+2]; 
					human_indices[counter2]=i;
					counter2++;
				}
			}
		}
		if(mouse){
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Select Code",mouse_names,mouse_names[0]);
			gd2.showDialog(); if(gd2.wasCanceled()) return;
			int editcode=gd2.getNextChoiceIndex();
			if(editcode==0 && action==2) return; //can't delete the default codes, exit
			if(editcode==0 && action==1) action=0; //can't edit the default codes, add them
			if(action==2){
				int delcode=mouse_indices[editcode];
				setCustomCodes(codes,delcode);
				return;
			}
			int[][] selcode=algutils.clone_multidim_array((int[][])mouse_codes[editcode]);
			Object[][] tabledata=new Object[22][7];
			String tempname=mouse_names[editcode];
			if(editcode==0) tempname="Mouse-1";
			else tempname+="-1";
			tabledata[0][0]="Name"; tabledata[0][1]=tempname;
			for(int i=0;i<21;i++){
				tabledata[i+1][0]=""+(i+1);
				for(int j=0;j<5;j++){
					tabledata[i+1][j+1]=new Boolean(selcode[i][j]==1);
				}
				tabledata[i+1][6]=new Integer(selcode[i][7]);
			}
			tabledata[20][0]="X";
			tabledata[21][0]="Y";
			Object[][] out=TableDialog2.showDialog(null,null,"Code_Editor",new String[]{"Chrom","A","B","C","D","E","Size"},tabledata,null);
			if(out==null) return;
			tempname=(String)out[0][1];
			for(int i=0;i<21;i++){
				for(int j=0;j<5;j++){
					selcode[i][j]=((Boolean)out[i+1][j+1]).booleanValue()?1:0;
				}
				selcode[i][7]=((Double)out[i+1][6]).intValue();
			}
			if(action==1){ //edit the code
				int edcode2=mouse_indices[editcode];
				setCustomCodes(codes,edcode2,selcode);
				return;
			} else {
				setCustomCodes(codes,tempname,true,selcode);
				return;
			}
		} else {
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Select Code",human_names,human_names[0]);
			gd2.showDialog(); if(gd2.wasCanceled()) return;
			int editcode=gd2.getNextChoiceIndex();
			if(editcode==0 && action==2) return; //can't delete the default codes, exit
			if(editcode==0 && action==1) action=0; //can't edit the default codes, add them
			if(action==2){
				int delcode=human_indices[editcode];
				setCustomCodes(codes,delcode);
				return;
			}
			int[][] selcode=algutils.clone_multidim_array((int[][])human_codes[editcode]);
			Object[][] tabledata=new Object[25][7];
			String tempname=human_names[editcode];
			if(editcode==0) tempname="Human-1";
			else tempname+="-1";
			tabledata[0][0]="Name"; tabledata[0][1]=tempname;
			for(int i=0;i<24;i++){
				tabledata[i+1][0]=""+(i+1);
				for(int j=0;j<5;j++){
					tabledata[i+1][j+1]=new Boolean(selcode[i][j]==1);
				}
				tabledata[i+1][6]=new Integer(selcode[i][7]);
			}
			tabledata[23][0]="X";
			tabledata[24][0]="Y";
			Object[][] out=TableDialog2.showDialog(null,null,"Code_Editor",new String[]{"Chrom","A","B","C","D","E","Size"},tabledata,null);
			if(out==null) return;
			tempname=(String)out[0][1];
			for(int i=0;i<24;i++){
				for(int j=0;j<5;j++){
					selcode[i][j]=((Boolean)out[i+1][j+1]).booleanValue()?1:0;
				}
				selcode[i][7]=((Double)out[i+1][6]).intValue();
			}
			if(action==1){ //edit the code
				int edcode2=human_indices[editcode];
				setCustomCodes(codes,edcode2,selcode);
				return;
			} else {
				setCustomCodes(codes,tempname,false,selcode);
				return;
			}
		}
	}

	public void setCustomCodes(Object[] data,int edindex,Object newdata){
		//here we edit a code
		data[edindex+2]=newdata;
		KissPanel.setCustomCodes(data);
	}

	public void setCustomCodes(Object[] data,String newname,boolean newmouse,Object newdata){
		//here we add a new code
		if(data==null){
			data=new Object[3];
			data[0]=new String[]{newname};
			data[1]=new boolean[]{newmouse};
			data[2]=newdata;
			KissPanel.setCustomCodes(data);
			return;
		}
		if(newdata==null){KissPanel.setCustomCodes(data); return;}
		String[] names=(String[])data[0];
		boolean[] mouse=(boolean[])data[1];
		String[] newnames=new String[names.length+1];
		boolean[] newmouses=new boolean[names.length+1];
		Object[] newdatas=new Object[names.length+3];
		for(int i=0;i<names.length;i++){
			newnames[i]=names[i];
			newmouses[i]=mouse[i];
			newdatas[i+2]=data[i+2];
		}
		newnames[names.length]=newname;
		newmouses[names.length]=newmouse;
		newdatas[0]=newnames;
		newdatas[1]=newmouses;
		newdatas[names.length+2]=newdata;
		KissPanel.setCustomCodes(newdatas);
	}

	public void setCustomCodes(Object[] data,int delindex){
		//here we delete an index
		String[] names=(String[])data[0];
		boolean[] mouse=(boolean[])data[1];
		String[] newnames=new String[names.length-1];
		boolean[] newmouses=new boolean[names.length-1];
		Object[] newdatas=new Object[names.length+1];
		if(data.length>3){
			for(int i=0;i<delindex;i++){
				newnames[i]=names[i];
				newmouses[i]=mouse[i];
				newdatas[i+2]=data[i+2];
			}
			for(int i=delindex+1;i<names.length;i++){
				newnames[i-1]=names[i];
				newmouses[i-1]=mouse[i];
				newdatas[i+1]=data[i+2];
			}
		}
		newdatas[0]=newnames;
		newdatas[1]=newmouses;
		KissPanel.setCustomCodes(newdatas);
	}


}
