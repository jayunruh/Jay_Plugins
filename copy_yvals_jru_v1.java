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
import ij.plugin.PlugIn;
import jguis.*;
import java.io.*;
import java.awt.datatransfer.*;

public class copy_yvals_jru_v1 implements PlugIn, ClipboardOwner {

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		Clipboard systemClipboard = null;
		try {systemClipboard = IJ.getInstance().getToolkit().getSystemClipboard();}
		catch (Exception e) {systemClipboard = null; }
		if (systemClipboard==null)
			{IJ.error("Unable to copy to Clipboard."); return;}
		IJ.showStatus("Copying plot values...");
		StringBuffer sb = new StringBuffer();
		int length=yvals[0].length;
		int nseries=yvals.length;
		for (int i=0; i<length; i++) {
			for (int j=0; j<nseries; j++) {
				sb.append(""+yvals[j][i]);
				if(j<(nseries-1)){sb.append("\t");}
			}
			sb.append("\n");
		}
		String text = sb.toString();
		StringSelection contents = new StringSelection(text);
		systemClipboard.setContents(contents, this);
		IJ.showStatus(text.length() + " characters copied to Clipboard");
	}

	public void lostOwnership(Clipboard clipboard, Transferable contents) {}

}
