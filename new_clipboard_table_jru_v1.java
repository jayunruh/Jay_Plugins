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
import java.awt.Toolkit;
import java.util.*;
import ij.plugin.*;
import ij.util.*;
import ij.text.*;
import java.io.*;
import java.awt.datatransfer.*;
import jguis.*;

public class new_clipboard_table_jru_v1 implements PlugIn {

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addCheckbox("Contains_Column_Titles?",true);
		String[] delims={"tab","comma","space"};
		gd.addChoice("Delimiter",delims,delims[0]);
		gd.addCheckbox("Treat_Consecutive_Delims_As_One",false);
		gd.addStringField("Title","Clipboard");
		gd.showDialog(); if(gd.wasCanceled()){return;}
		boolean titles=gd.getNextBoolean();
		int delimindex=gd.getNextChoiceIndex();
		boolean noconsec=gd.getNextBoolean();
		String title=gd.getNextString();
		String delim="\t";
		if(delimindex==1) delim=",";
		if(delimindex==2) delim=" ";
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
		String[] lines=table_tools.split(textdata,"\n");
		List<List<String>> listtable=table_tools.table2listtable(lines,delim,noconsec);
		String[] headings=null;
		if(titles){
			headings=table_tools.list2stringarray(listtable.get(0));
			listtable.remove(0);
		}
		table_tools.create_table(title,listtable,headings);
	}
}
