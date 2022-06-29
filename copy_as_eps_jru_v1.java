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
import ij.plugin.frame.*;
import jguis.*;
import java.awt.datatransfer.*;
import java.io.*;

public class copy_as_eps_jru_v1 implements PlugIn, Transferable {
	Object p4;
	DataFlavor epsflavor;

	public void run(String arg) {
		ImageWindow iw=WindowManager.getCurrentWindow();
		p4=jutils.runReflectionMethod(iw,"getPlot",null,null);
		if(p4!=null){
			try{
				epsflavor=new DataFlavor("image/eps","Encapsulated Postscript File");
				((SystemFlavorMap)SystemFlavorMap.getDefaultFlavorMap()).addUnencodedNativeForFlavor(epsflavor,"PSFILE");
				Toolkit.getDefaultToolkit().getSystemClipboard().setContents(this,null);
				IJ.showStatus("Plot Copied");
			}catch(Throwable e){
				IJ.log(e.toString());
			}
		}
	}

	public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException{
		if(isDataFlavorSupported(flavor)){
			byte[] bindata=(byte[])jutils.runReflectionMethod(p4,"getPSBinary",null,null);
			return new ByteArrayInputStream(bindata);
		}
		throw new UnsupportedFlavorException(flavor);
	}

	public DataFlavor[] getTransferDataFlavors(){
		return new DataFlavor[]{epsflavor};
	}

	public boolean isDataFlavorSupported(DataFlavor dataFlavor){
		if (dataFlavor.match(epsflavor)) return true;
		return false;
	}

}
