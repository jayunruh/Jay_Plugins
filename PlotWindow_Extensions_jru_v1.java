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
import ij.macro.Functions;
import ij.macro.MacroExtension;
import ij.macro.ExtensionDescriptor;
import java.awt.Color;
import java.awt.Frame;
import ij.plugin.*;
import jalgs.*;
import jguis.*;
import ij.text.*;
import java.util.*;


public class PlotWindow_Extensions_jru_v1 implements PlugIn, MacroExtension {
	private ExtensionDescriptor[] extensions;

	public void run(String arg) {
		if (!IJ.macroRunning()) {
			IJ.error("Cannot install extensions from outside a macro!");
			return;
		}

		extensions=new ExtensionDescriptor[]{
			ExtensionDescriptor.newDescriptor("getNSeries",this,MacroExtension.ARG_OUTPUT+MacroExtension.ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("selectSeries",this,MacroExtension.ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("scalePlotCoords",this,ARG_NUMBER,ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("setLimits",this,ARG_NUMBER,ARG_NUMBER,ARG_NUMBER,ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getLimits",this,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("setLogAxes",this,ARG_NUMBER,ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("setXLabel",this,ARG_STRING),
			ExtensionDescriptor.newDescriptor("setYLabel",this,ARG_STRING),
			ExtensionDescriptor.newDescriptor("setGridWhiteness",this,ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("deleteSelected",this),
			ExtensionDescriptor.newDescriptor("autoscaleX",this),
			ExtensionDescriptor.newDescriptor("autoscaleY",this),
			ExtensionDescriptor.newDescriptor("autoscaleZ",this),
			ExtensionDescriptor.newDescriptor("scaleROI",this),
			ExtensionDescriptor.newDescriptor("setMagnification",this,ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("setMagRatio",this,ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getSelNpts",this,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getSelIndexXYVals",this,ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getSelIndexXYZVals",this,ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getSelStat",this,ARG_STRING,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getSelStat2",this,ARG_STRING,ARG_NUMBER,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("addXYSeries",this,ARG_ARRAY,ARG_ARRAY),
			ExtensionDescriptor.newDescriptor("updateSelSeries",this,ARG_ARRAY,ARG_ARRAY),
			ExtensionDescriptor.newDescriptor("plot2List",this),
			ExtensionDescriptor.newDescriptor("createPlot",this,ARG_STRING,ARG_STRING,ARG_ARRAY,ARG_ARRAY),
			ExtensionDescriptor.newDescriptor("create3DTraj",this,ARG_ARRAY,ARG_ARRAY,ARG_ARRAY),
			ExtensionDescriptor.newDescriptor("convertToPW4",this),
			ExtensionDescriptor.newDescriptor("convertToPW",this),
			ExtensionDescriptor.newDescriptor("setBinSize",this,MacroExtension.ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getSelected",this,MacroExtension.ARG_OUTPUT+MacroExtension.ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("setCommand",this,ARG_STRING,ARG_STRING),
			ExtensionDescriptor.newDescriptor("getCommand",this,ARG_STRING,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("getXLabel",this,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("getYLabel",this,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("selectTable",this,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("tableExists",this,ARG_STRING,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("setTableValue",this,ARG_STRING,ARG_NUMBER,ARG_NUMBER,ARG_STRING),
			ExtensionDescriptor.newDescriptor("getTableValue",this,ARG_STRING,ARG_NUMBER,ARG_NUMBER,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("getTableHeadings",this,ARG_STRING,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("getTableSize",this,ARG_STRING,ARG_OUTPUT+ARG_NUMBER),
			ExtensionDescriptor.newDescriptor("getLastLog",this,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("selectImage",this,ARG_OUTPUT+ARG_STRING),
			ExtensionDescriptor.newDescriptor("testArray",this,ARG_ARRAY)
		};

		Functions.registerExtensions(this);
	}

	public ExtensionDescriptor[] getExtensionFunctions(){
		return extensions;
	}

	public String handleExtension(String name,Object[] args){
		if(name.equals("createPlot")){
			String xlabel=(String)args[0];
			String ylabel=(String)args[1];
			Object[] xvals=(Object[])args[2];
			Object[] yvals=(Object[])args[3];
			float[] yvals2=new float[yvals.length];
			for(int i=0;i<yvals.length;i++) yvals2[i]=((Double)yvals[i]).floatValue();
			if(xvals.length==yvals.length){
				float[] xvals2=new float[xvals.length];
				for(int i=0;i<xvals.length;i++) xvals2[i]=((Double)xvals[i]).floatValue();
				new PlotWindow4("Macro PlotWindow4",xlabel,ylabel,xvals2,yvals2).draw();
			} else {
				new PlotWindow4("Macro PlotWindow4",xlabel,ylabel,yvals2).draw();
			}
			return null;
		}
		if(name.equals("create3DTraj")){
			Object[] xvals=(Object[])args[0];
			Object[] yvals=(Object[])args[1];
			Object[] zvals=(Object[])args[2];
			float[] xvals2=new float[xvals.length];
			for(int i=0;i<xvals.length;i++) xvals2[i]=((Double)xvals[i]).floatValue();
			float[] yvals2=new float[yvals.length];
			for(int i=0;i<yvals.length;i++) yvals2[i]=((Double)yvals[i]).floatValue();
			float[] zvals2=new float[zvals.length];
			for(int i=0;i<zvals.length;i++) zvals2[i]=((Double)zvals[i]).floatValue();
			Traj3D t3d=new Traj3D("x","y","z",xvals2,yvals2,zvals2);
			new PlotWindow3D("Macro 3D Trajectory",t3d).draw();
			return null;
		}
		if(name.equals("selectImage")){
			ImagePlus[] imp=jutils.selectImages(false,1);
			if(imp!=null && imp.length>0){
				((String[])args[0])[0]=imp[0].getTitle();
			}
			return null;
		}
		if(name.equals("selectTable")){
			TextWindow[] tw=jutils.selectTables(false,1);
			if(tw!=null && tw.length>0){
				((String[])args[0])[0]=tw[0].getTitle();
			}
			return null;
		}
		if(name.equals("tableExists")){
			String title=(String)args[0];
			TextWindow tw=jutils.selectTable(title);
			if(tw!=null){
				((String[])args[1])[0]="true";
			} else {
				((String[])args[1])[0]="false";
			}
			return null;
		}
		if(name.equals("setTableValue")){
			String title=(String)args[0];
			int row=((Double)args[1]).intValue();
			int col=((Double)args[2]).intValue();
			TextWindow tw=jutils.selectTable(title);
			String val=(String)args[3];
			if(tw!=null){
				TextPanel tp=tw.getTextPanel();
				List<List<String>> listtable=table_tools.table2listtable(tp);
				List<String> row2=listtable.get(row);
				row2.set(col,val);
				table_tools.replace_table(tp,listtable,tp.getColumnHeadings());
			}
			return null;
		}
		if(name.equals("getTableValue")){
			//parameters are tabname,row,col,outvar
			String title=(String)args[0];
			int row=((Double)args[1]).intValue();
			int col=((Double)args[2]).intValue();
			TextWindow tw=jutils.selectTable(title);
			if(tw!=null){
				TextPanel tp=tw.getTextPanel();
				String line=tp.getLine(row);
				String[] temp=table_tools.split_string_tab(line);
				//((Double[])args[3])[0]=new Double(temp[col]);
				((String[])args[3])[0]=temp[col];
			}
			return null;
		}
		if(name.equals("getLastLog")){
			//this gets the last line in the log window
			String log=IJ.getLog();
			if(log!=null){
				String[] lines=log.split("\n");
				((String[])args[0])[0]=lines[lines.length-1];
			}
			return null;
		}
		if(name.equals("getTableHeadings")){
			String title=(String)args[0];
			TextWindow tw=jutils.selectTable(title);
			if(tw!=null){
				TextPanel tp=tw.getTextPanel();
				String headings=tp.getColumnHeadings();
				((String[])args[1])[0]=headings;
			}
			return null;
		}
		if(name.equals("getTableSize")){
			String title=(String)args[0];
			TextWindow tw=jutils.selectTable(title);
			if(tw!=null){
				TextPanel tp=tw.getTextPanel();
				int nlines=tp.getLineCount();
				((Double[])args[1])[0]=new Double(nlines);
			}
			return null;
		}
		if(name.equals("testArray")){
			for(int i=0;i<args.length;i++){
				if(args[i] instanceof String) IJ.log("args"+i+"="+(String)args[i]);
				if(args[i] instanceof Object[]){
					Object[] tempargs=(Object[])args[i];
					for(int j=0;j<tempargs.length;j++){ //these are either strings or doubles
						IJ.log("args"+i+","+j+"="+tempargs[j].toString());
						IJ.log("args"+i+","+j+" type "+tempargs[j].getClass().getName());
					}
				}
			}
			return null;
		}
		ImageWindow iw=WindowManager.getCurrentWindow();
		if(!jutils.isPlotFamily(iw)){
			IJ.error("Current Image is not Compatible");
			return null;
		}
		if(name.equals("convertToPW4")){
			jutils.pw2pw4((PlotWindow)iw).draw();
			return null;
		}
		if(name.equals("convertToPW")){
			jutils.pw42pw(iw);
			return null;
		}
		if(name.equals("plot2List")){
			//list the plot data in a table
			jutils.runPW4VoidMethod(iw,"showList");
			return null;
		}
		if(name.equals("getNSeries")){
			int nseries=((Integer)jutils.runPW4VoidMethod(iw,"getNSeries")).intValue();
			((Double[])args[0])[0]=new Double((double)nseries);
		}
		if(name.equals("selectSeries")){
			int series=((Double)args[0]).intValue();
			jutils.runReflectionMethod(iw,"selectSeries",new Object[]{series});
		}
		if(name.equals("setLimits")){
			float xmin=((Double)args[0]).floatValue();
			float xmax=((Double)args[1]).floatValue();
			float ymin=((Double)args[2]).floatValue();
			float ymax=((Double)args[3]).floatValue();
			float[] limits=(float[])jutils.runPW4VoidMethod(iw,"getLimits");
			limits[0]=xmin; limits[1]=xmax; limits[2]=ymin; limits[3]=ymax;
			jutils.runReflectionMethod(iw,"setLimits",new Object[]{limits});
		}
		if(name.equals("getLimits")){
			float[] limits=(float[])jutils.runPW4VoidMethod(iw,"getLimits");
			((Double[])args[0])[0]=new Double(limits[0]);
			((Double[])args[1])[0]=new Double(limits[1]);
			((Double[])args[2])[0]=new Double(limits[2]);
			((Double[])args[3])[0]=new Double(limits[3]);
		}
		if(name.equals("scalePlotCoords")){
			//need to make this equivalent to the "toScaled" macro function but with log plot support
			int x=((Double)args[0]).intValue();
			int y=((Double)args[1]).intValue();
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			float[] coords=(float[])jutils.runReflectionMethod(plot,"getPlotCoords",new Object[]{x,y});
			((Double[])args[2])[0]=new Double(coords[0]);
			((Double[])args[3])[0]=new Double(coords[1]);
		}
		if(name.equals("setLogAxes")){
			//0 means linear axis, 1 means log axis
			boolean logx=(((Double)args[0]).intValue()==1)?true:false;
			boolean logy=(((Double)args[1]).intValue()==1)?true:false;
			jutils.runReflectionMethod(iw,"setLogAxes",new Object[]{logx,logy});
		}
		if(name.equals("setXLabel")){
			String label=(String)args[0];
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"setxLabel",new Object[]{label});
		}
		if(name.equals("setYLabel")){
			String label=(String)args[0];
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"setyLabel",new Object[]{label});
		}
		if(name.equals("setGridWhiteness")){
			int whiteness=((Double)args[0]).intValue();
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"setGridWhiteness",new Object[]{whiteness});
		}
		if(name.equals("deleteSelected")){
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			int selindex=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
			jutils.runReflectionMethod(plot,"deleteSeries",new Object[]{selindex,false});
		}
		if(name.equals("autoscaleX")){
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"xautoscale",null);
		}
		if(name.equals("autoscaleY")){
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"yautoscale",null);
		}
		if(name.equals("autoscaleZ")){
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			if(plot instanceof Plot2DHist) jutils.runReflectionMethod(plot,"intautoscale",null);
			else jutils.runReflectionMethod(plot,"zautoscale",null);
		}
		if(name.equals("scaleROI")){
			jutils.runPW4VoidMethod(iw,"scaleroi");
		}
		if(name.equals("setMagnification")){
			float mag=((Double)args[0]).floatValue();
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"setmagnification",new Object[]{mag});
		}
		if(name.equals("setMagRatio")){
			float magratio=((Double)args[0]).floatValue();
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			jutils.runReflectionMethod(plot,"setmagratio",new Object[]{magratio});
		}
		if(name.equals("getSelNpts")){
			Object npts1=jutils.runPW4VoidMethod(iw,"getNpts");
			int[] npts=null;
			if(npts1 instanceof int[]) npts=(int[])npts1;
			else npts=((int[][])npts1)[0];
			int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
			if(sel<0) sel=0;
			int selnpts=npts[sel];
			((Double[])args[0])[0]=new Double(selnpts);
		}
		if(name.equals("getSelIndexXYVals")){
			float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
			if(sel<0) sel=0;
			int index=((Double)args[0]).intValue();
			int[] npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
			if(index<0 || index>=npts[sel]){
				((Double[])args[1])[0]=new Double(0.0);
				((Double[])args[2])[0]=new Double(0.0);
			} else {
				((Double[])args[1])[0]=new Double(xvals[sel][index]);
				((Double[])args[2])[0]=new Double(yvals[sel][index]);
			}
		}
		if(name.equals("getSelIndexXYZVals")){
			//assume this is a trajectory (not a surface)
			float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
			float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
			float[][] zvals=null;
			if(jutils.is3DPlot(iw)) zvals=((float[][][])jutils.runPW4VoidMethod(iw,"getZValues"))[0];
			int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
			if(sel<0) sel=0;
			int index=((Double)args[0]).intValue();
			int[] npts=null;
			if(jutils.is3DPlot(iw)) npts=((int[][])jutils.runPW4VoidMethod(iw,"getNpts"))[0];
			else npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
			((Double[])args[3])[0]=new Double(0.0);
			if(index<0 || index>=npts[sel]){
				((Double[])args[1])[0]=new Double(0.0);
				((Double[])args[2])[0]=new Double(0.0);
				//((Double[])args[3])[0]=new Double(0.0);
			} else {
				((Double[])args[1])[0]=new Double(xvals[sel][index]);
				((Double[])args[2])[0]=new Double(yvals[sel][index]);
				if(zvals!=null) ((Double[])args[3])[0]=new Double(zvals[sel][index]);
			}
		}
		if(name.equals("getSelStat")){
			float[][] yvals=null;
			int sel=0;
			int[] npts=null;
			if(jutils.isPlotHist(iw)){
				yvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
				npts=new int[]{yvals[0].length};
			} else {
				yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
				npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
				sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
				if(sel<0) sel=0;
			}
			String stat=(String)args[0];
			float[] temp=new float[npts[sel]];
			System.arraycopy(yvals[sel],0,temp,0,npts[sel]);
			float stat2=jstatistics.getstatistic(stat,temp,null);
			((Double[])args[1])[0]=new Double(stat2);
		}
		if(name.equals("getSelStat2")){
			float[][] yvals=null;
			int sel=0;
			int[] npts=null;
			if(jutils.isPlotHist(iw)){
				yvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
				npts=new int[]{yvals[0].length};
			} else {
				yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
				npts=(int[])jutils.runPW4VoidMethod(iw,"getNpts");
				sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
				if(sel<0) sel=0;
			}
			String stat=(String)args[0];
			float extra=((Double)args[1]).floatValue();
			float[] temp=new float[npts[sel]];
			System.arraycopy(yvals[sel],0,temp,0,npts[sel]);
			float stat2=jstatistics.getstatistic(stat,temp,new float[]{extra});
			((Double[])args[2])[0]=new Double(stat2);
		}
		if(name.equals("addXYSeries")){
			Object[] xvals=(Object[])args[0];
			Object[] yvals=(Object[])args[1];
			float[] yvals2=new float[yvals.length];
			for(int i=0;i<yvals.length;i++) yvals2[i]=((Double)yvals[i]).floatValue();
			if(xvals.length==yvals.length){
				float[] xvals2=new float[xvals.length];
				for(int i=0;i<xvals.length;i++) xvals2[i]=((Double)xvals[i]).floatValue();
				jutils.runReflectionMethod(iw,"addPoints",new Object[]{xvals2,yvals2,true});
			} else {
				jutils.runReflectionMethod(iw,"addPoints",new Object[]{yvals2,true});
			}
		}
		if(name.equals("updateSelSeries")){
			Object plot=jutils.runPW4VoidMethod(iw,"getPlot");
			int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected");
			if(sel<0) sel=0;
			Object[] xvals=(Object[])args[0];
			Object[] yvals=(Object[])args[1];
			float[] yvals2=new float[yvals.length];
			for(int i=0;i<yvals.length;i++) yvals2[i]=((Double)yvals[i]).floatValue();
			if(xvals.length==yvals.length){
				float[] xvals2=new float[xvals.length];
				for(int i=0;i<xvals.length;i++) xvals2[i]=((Double)xvals[i]).floatValue();
				jutils.runReflectionMethod(plot,"updateSeries",new Object[]{xvals2,yvals2,sel,true});
			} else {
				jutils.runReflectionMethod(plot,"updateSeries",new Object[]{yvals2,sel,true});
			}
		}
		if(name.equals("setBinSize")){
			float binsize=((Double)args[0]).floatValue();
			Object plot=jutils.runReflectionMethod(iw,"getPlot",null);
			if(iw.getClass().getName().equals("jguis.PlotWindow2DHist")) jutils.runReflectionMethod(plot,"setBinSize",new Object[]{(int)binsize});
			else jutils.runReflectionMethod(plot,"setBinSizeUnits",new Object[]{binsize});
		}
		if(name.equals("getSelected")){
			int sel=((Integer)jutils.runPW4VoidMethod(iw,"getSelected")).intValue();
			((Double[])args[0])[0]=new Double((double)sel);
		}
		if(name.equals("getXLabel")){
			String xlab=(String)jutils.runPW4VoidMethod(iw,"getxLabel");
			((String[])args[0])[0]=xlab;
		}
		if(name.equals("getYLabel")){
			String ylab=(String)jutils.runPW4VoidMethod(iw,"getyLabel");
			((String[])args[0])[0]=ylab;
		}
		jutils.runPW4VoidMethod(iw,"updatePlot");
		return null;
	}

}
