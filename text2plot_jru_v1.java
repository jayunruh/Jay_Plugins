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
import java.awt.*;
import ij.plugin.*;
import ij.text.*;
import ij.plugin.frame.*;
import jguis.*;
import jalgs.*;
import ij.io.*;
import java.io.*;
import java.awt.datatransfer.*;

public class text2plot_jru_v1 implements PlugIn {

	public void run(String arg) {
		Frame[] niframes=WindowManager.getNonImageWindows();
		String[] titles=new String[niframes.length+1];
		for(int i=0;i<niframes.length;i++){
			titles[i]=niframes[i].getTitle();
		}
		titles[niframes.length]="Clipboard";
		GenericDialog gd=new GenericDialog("Windows");
		boolean importfile=false;
		gd.addCheckbox("Import from file?",importfile);
		gd.addChoice("Windows",titles,titles[0]);
		boolean hasxvals=false;
		gd.addCheckbox("X Vals Column?",hasxvals);
		boolean multix=false;
		gd.addCheckbox("Multi_X_Columns?",multix);
		boolean skipendzeros=false;
		gd.addCheckbox("Skip_end_zeros?",skipendzeros);
		String[] delimiters={"Tab","Comma","Space"};
		gd.addChoice("Delimiter",delimiters,delimiters[0]);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		importfile=gd.getNextBoolean();
		int index=gd.getNextChoiceIndex();
		hasxvals=gd.getNextBoolean();
		multix=gd.getNextBoolean();
		skipendzeros=gd.getNextBoolean();
		int delimindex=gd.getNextChoiceIndex();
		if(multix) hasxvals=true;
		String textdata="";
		if(importfile){
			OpenDialog od = new OpenDialog("Open File","",".txt");
			String directory=od.getDirectory();
			String name=od.getFileName();
			if(name==null){return;}
			try{
				File infile=new File(directory+name);
				BufferedReader b=new BufferedReader(new FileReader(infile));
				textdata=(new jdataio()).readstringfile(b);
				b.close();
			} catch(IOException e){
				return;
			}
		} else {
			if(index==niframes.length){
				//here we get the data from the clipboard
				Transferable t=Toolkit.getDefaultToolkit().getSystemClipboard().getContents(null);
				try {
 					if (t != null && t.isDataFlavorSupported(DataFlavor.stringFlavor)) {
						textdata = (String)t.getTransferData(DataFlavor.stringFlavor);
					}
				} catch (UnsupportedFlavorException e) {}
				catch (IOException e) {}
				if(textdata.equals("")){
					IJ.error("Error copying from clipboard.");
					return;
				}	
			} else {
				if(niframes[index] instanceof Editor){
					Editor tw=(Editor)niframes[index];
					textdata=tw.getText();
				} else {
					if(niframes[index] instanceof TextWindow){
						TextWindow tw=(TextWindow)niframes[index];
						textdata=tw.getTextPanel().getText();
					} else {
						IJ.showMessage("Not a valid text window");
						return;
					}
				}
			}
		}
		if(textdata==null){IJ.showMessage("Error in Obtaining String"); return;}
		if(textdata.indexOf("\r")>=0){
			textdata=textdata.replace('\r','\n');
		}
		char[] delims={'\t',',',' '};
		delimit_string ds=new delimit_string(delims[delimindex]);
		String[] rows=ds.getrows(textdata);
		int lines=rows.length;
		int columns=ds.getnumcolumns(rows[0]);
		int ycolumns=columns;
		if(hasxvals){
			if(multix){ycolumns/=2;}
			else{ycolumns--;}
		}
		if(multix){
			float[][] ydata=new float[ycolumns][lines];
			float[][] xdata=new float[ycolumns][lines];
			for(int i=0;i<lines;i++){
				float[] temp=ds.delim2float(rows[i],columns);
				for(int j=0;j<ycolumns;j++){
					ydata[j][i]=temp[2*j+1];
					xdata[j][i]=temp[2*j];
				}
			}
			int[] npts=new int[ycolumns];
			for(int i=0;i<ycolumns;i++){npts[i]=lines;}
			if(skipendzeros){
				for(int i=0;i<ycolumns;i++){
					int counter=lines-1;
					while((xdata[i][counter]==0.0f || Float.isNaN(xdata[i][counter])) && counter>0){
						xdata[i][counter]=0.0f;
						ydata[i][counter]=0.0f;
						npts[i]--;
						counter--;
					}
				}
			}
			(new PlotWindow4("Text Plot","x","y",xdata,ydata,npts)).draw();
		} else {
			float[][] tempydata=new float[ycolumns][lines];
			float[] tempxdata=new float[lines];
			float[][] xdata=null;
			float[][] ydata=null;
			int startcolumn=0;
			if(hasxvals) startcolumn=1;
			for(int i=0;i<lines;i++){
				float[] temp=ds.delim2float(rows[i],columns);
				if(hasxvals){tempxdata[i]=temp[0];}
				else{tempxdata[i]=(float)(i+1);}
				for(int j=0;j<ycolumns;j++){
					tempydata[j][i]=temp[j+startcolumn];
				}
			}
			int[] npts=new int[ycolumns];
			npts[0]=lines;
			if(skipendzeros){
				int maxpts=0;
				for(int i=0;i<ycolumns;i++){
					int counter=lines-1;
					npts[i]=lines;
					while((tempydata[i][counter]==0.0f || Float.isNaN(tempydata[i][counter])) && counter>0){
						npts[i]--;
						counter--;
					}
					if(npts[i]>maxpts) maxpts=npts[i];
					IJ.log(""+npts[i]);
				}
				ydata=new float[ycolumns][maxpts];
				xdata=new float[ycolumns][maxpts];
				for(int i=0;i<ycolumns;i++){
					//npts[i]=npts[0];
					System.arraycopy(tempxdata,0,xdata[i],0,npts[i]);
					System.arraycopy(tempydata[i],0,ydata[i],0,npts[i]);
				}
			} else {
				ydata=tempydata;
				xdata=new float[ycolumns][];
				for(int i=0;i<ycolumns;i++){
					npts[i]=npts[0];
					xdata[i]=tempxdata.clone();
				}
			}
			(new PlotWindow4("Text Plot","x","y",xdata,ydata,npts)).draw();
		}
	}
			
}
