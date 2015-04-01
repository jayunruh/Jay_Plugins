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
import java.io.*;
import ij.io.*;

public class import_plot_jru_v1 implements PlugIn {

	public void run(String arg) {
		OpenDialog od=new OpenDialog("Open Plot Object File",arg);
		String dir=od.getDirectory();
		String name=od.getFileName();
		if(name==null){return;}
		String path=dir+File.separator+name;
		if(Plot4.is_this(path)){
			IJ.showStatus("Importing PlotWindow");
			Plot4 plot=new Plot4(path);
			new PlotWindow4(name,plot).draw();
		} else {
			if(Plot3D.is_this(path)){
				IJ.showStatus("Importing Plot3D");
				Plot3D plot=new Plot3D(path);
				new PlotWindow3D(name,plot).draw();
			} else {
				if(Traj3D.is_this(path)){
					IJ.showStatus("Importing Traj3D");
					Traj3D plot=new Traj3D(path);
					new PlotWindow3D(name,plot).draw();
				} else {
					if(PlotHist.is_this(path)){
						IJ.showStatus("Importing PlotHist");
						PlotHist plot=new PlotHist(path);
						new PlotWindowHist(name,plot).draw();
					} else {
						if(Plot2DHist.is_this(path)){
							IJ.showStatus("Importing Plot2DHist");
							Plot2DHist plot=new Plot2DHist(path);
							new PlotWindow2DHist(name,plot).draw();
						} else {
							 if(PlotColumn.is_this(path)){
								IJ.showStatus("Importing PlotColumn");
								PlotColumn plot=new PlotColumn(path);
								new PlotWindowColumn(name,plot).draw();
							} else {
								IJ.log("unsupported or corrupted pw2 file");
							}
						}
					}
				}
			}
		}
	}

}
