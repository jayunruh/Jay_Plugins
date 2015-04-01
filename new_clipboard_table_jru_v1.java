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
import ij.util.*;
import ij.text.*;
import java.io.*;
import java.awt.datatransfer.*;

public class new_clipboard_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Contains_Column_Titles?",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean titles=gd.getNextBoolean();
		String textdata="";
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
		if(textdata.indexOf("\r")>=0){
			textdata=textdata.replace('\r','\n');
		}
		if(titles){
			int first=textdata.indexOf("\n");
			String headings=textdata.substring(0,first);
			String rest=textdata.substring(first+1);
			new TextWindow("Clipboard",headings,rest,400,200);
		} else {
			String first=textdata.substring(0,textdata.indexOf("\n"));
			String[] temp=first.split("\t");
			StringBuffer sb=new StringBuffer();
			for(int i=0;i<temp.length;i++){
				sb.append("Col"+(i+1)+"\t");
			}
			new TextWindow("Clipboard",sb.toString(),textdata,400,200);
		}
	}
}
