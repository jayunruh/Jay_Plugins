/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.plugin.frame.Recorder;
import ij.process.ColorProcessor;
import ij.text.TextWindow;
import ij.util.Tools;
import jalgs.jdataio;

import java.awt.Button;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.List;

/**
 * This class is an extended ImageWindow that displays line graphs. This class
 * was adapted from the PlotWindow class original to ImageJ Adaptations were
 * done by Jay Unruh (jru@stowers.org)
 */
public class PlotWindow4 extends ImageWindow implements ActionListener,ClipboardOwner{

	private Button list,save,copy,editbutton,selbutton,seldbutton;
	// private Label coordinates;
	protected static String defaultDirectory=null;
	public Plot4 p3;
	protected static ColorProcessor cp;

	public PlotWindow4(String title1){
		super(createImage(title1));
	}

	public PlotWindow4(String title1,String xLabel1,String yLabel1,float[][] xValues1,float[][] yValues1,Object npts1){
		super(createImage(title1));
		p3=new Plot4(xLabel1,yLabel1,xValues1,yValues1,npts1);
	}

	public PlotWindow4(String title1,String xLabel1,String yLabel1,float[][] yValues1,Object npts1){
		super(createImage(title1));
		p3=new Plot4(xLabel1,yLabel1,yValues1,npts1);
	}

	public PlotWindow4(String title1,String xLabel1,String yLabel1,float[] xValues1,float[] yValues1){
		super(createImage(title1));
		p3=new Plot4(xLabel1,yLabel1,xValues1,yValues1);
	}

	public PlotWindow4(String title1,String xLabel1,String yLabel1,double[] xValues1,double[] yValues1){
		this(title1,xLabel1,yLabel1,Tools.toFloat(xValues1),Tools.toFloat(yValues1));
	}

	public PlotWindow4(String title1,String xLabel1,String yLabel1,float[] yValues1){
		super(createImage(title1));
		p3=new Plot4(xLabel1,yLabel1,yValues1);
	}
	
	public PlotWindow4(String title1,String xLabel1,String yLabel1,List<List<float[]>> points){
		super(createImage(title1));
		int maxsize=0;
		int[] npts=new int[points.size()];
		for(int i=0;i<points.size();i++){
			npts[i]=points.get(i).size();
			if(npts[i]>maxsize) maxsize=npts[i];
		}
		float[][] xvals=new float[npts.length][maxsize];
		float[][] yvals=new float[npts.length][maxsize];
		for(int i=0;i<points.size();i++){
			List<float[]> series=points.get(i);
			for(int j=0;j<series.size();j++){
				float[] coords=series.get(j);
				xvals[i][j]=coords[0]; yvals[i][j]=coords[1];
			}
		}
		p3=new Plot4(xLabel1,yLabel1,xvals,yvals,npts);
	}

	public PlotWindow4(String title1,Plot4 plot){
		super(createImage(title1));
		p3=plot;
	}

	public PlotWindow4(ImagePlus imp,Plot4 plot){
		// turns the passed ImagePlus into a plot window
		// used for bioformats which wants to make the window for me
		super(imp);
		int width=Plot4.WIDTH+Plot4.LEFT_MARGIN+Plot4.RIGHT_MARGIN;
		int height=Plot4.HEIGHT+Plot4.TOP_MARGIN+Plot4.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		imp.setProcessor(cp);
		// cp=(ColorProcessor)this.imp.getProcessor();
		p3=plot;
	}

	public PlotWindow4 getInstance(){
		return this;
	}

	public void setmagnification(float newmag){
		p3.setmagnification(newmag);
		updatePlot();
	}

	public float getmagnification(){
		return p3.getmagnification();
	}

	public void setmagratio(float newmag){
		p3.setmagratio(newmag);
		updatePlot();
	}

	public float getmagratio(){
		return p3.getmagratio();
	}

	static ImagePlus createImage(String title1){
		int width=Plot4.WIDTH+Plot4.LEFT_MARGIN+Plot4.RIGHT_MARGIN;
		int height=Plot4.HEIGHT+Plot4.TOP_MARGIN+Plot4.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		return new ImagePlus(title1,cp);
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1){
		p3.setLimits(xMin1,xMax1,yMin1,yMax1);
		updatePlot();
	}

	public void setLimits(float[] limits){
		p3.setLimits(limits);
		updatePlot();
	}

	/** Sets the x-axis and y-axis range. */
	public void setLogAxes(boolean logx1,boolean logy1){
		p3.setLogAxes(logx1,logy1);
		updatePlot();
	}

	public void autoscale(){
		p3.autoscale();
		updatePlot();
	}

	public void xautoscale(){
		p3.xautoscale();
		updatePlot();
	}

	public void yautoscale(){
		p3.yautoscale();
		updatePlot();
	}

	public void scaleroi(){
		Rectangle rect=imp.getProcessor().getRoi();
		if(rect!=null){
			p3.scalerect(rect);
			imp.killRoi();
			updatePlot();
		}
	}

	public void updateSeries(float[] xValues1,float[] yValues1,int series,boolean rescale){
		p3.updateSeries(xValues1,yValues1,series,rescale);
		updatePlot();
	}

	public void updateSeries(float[] yValues1,int series,boolean rescale){
		p3.updateSeries(yValues1,series,rescale);
		updatePlot();
	}

	public void deleteSeries(int series,boolean rescale){
		p3.deleteSeries(series,rescale);
		updatePlot();
	}

	public void deleteMultiSeries(boolean[] fate,boolean rescale){
		p3.deleteMultiSeries(fate,rescale);
		updatePlot();
	}

	public void addPoints(float[] xValues1,float[] yValues1,boolean rescale){
		p3.addPoints(xValues1,yValues1,rescale);
		updatePlot();
	}

	public void addPoints(float[] xValues1,float[] yValues1,int shape){
		p3.addPoints(xValues1,yValues1,false);
		int[] shapes=p3.getShapes();
		shapes[shapes.length-1]=shape;
		updatePlot();
	}

	public void addPoints(double[] xValues1,double[] yValues1,int shape){
		addPoints(Tools.toFloat(xValues1),Tools.toFloat(yValues1),shape);
	}

	public void addPoints(float[] yValues1,boolean rescale){
		p3.addPoints(yValues1,rescale);
		updatePlot();
	}

	public void addErrors(float[][][] errors){
		p3.addErrors(errors);
		p3.setShowErrors(true);
		updatePlot();
	}

	public void addErrors(float[][] errors){
		p3.addErrors(errors);
		p3.setShowErrors(true);
		updatePlot();
	}

	public void addErrorBars(float[] errorBars){
		float[][] temp=new float[1][];
		temp[0]=errorBars;
		addErrors(temp);
	}
	
	public void addSeriesErrors(int series,float[][] errors){
		p3.addSeriesErrors(series,errors);
		p3.setShowErrors(true);
		updatePlot();
	}
	
	public void addSeriesErrors(int series,float[] errors){
		p3.addSeriesErrors(series,errors);
		p3.setShowErrors(true);
		updatePlot();
	}

	public void setShowErrors(boolean showerrors){
		p3.setShowErrors(showerrors);
	}

	public void updatePlot(){
		cp=p3.getProcessor();
		Calibration cal=p3.getCalibration();
		imp.setProcessor(null,cp);
		imp.setCalibration(cal);
		imp.updateAndDraw();
	}

	/** Displays the plot. */
	public void draw(){
		Panel buttons=new Panel();
		buttons.setLayout(new FlowLayout(FlowLayout.RIGHT));
		list=new Button(" List ");
		list.addActionListener(this);
		buttons.add(list);
		save=new Button("Save...");
		save.addActionListener(this);
		buttons.add(save);
		copy=new Button("Copy...");
		copy.addActionListener(this);
		buttons.add(copy);
		editbutton=new Button("Edit...");
		editbutton.addActionListener(this);
		buttons.add(editbutton);
		selbutton=new Button("Select+");
		selbutton.addActionListener(this);
		buttons.add(selbutton);
		seldbutton=new Button("Select-");
		seldbutton.addActionListener(this);
		buttons.add(seldbutton);
		// coordinates = new Label("X=12345678, Y=12345678");
		// coordinates.setFont(new Font("Monospaced", Font.PLAIN, 12));
		// buttons.add(coordinates);
		add(buttons);
		pack();
		// coordinates.setText("");
		updatePlot();
		IJ.register(this.getClass());
	}

	/**
	 * Updates the graph X and Y values when the mouse is moved. Overrides
	 * mouseMoved() in ImageWindow.
	 * 
	 * @see ij.gui.ImageWindow#mouseMoved
	 */
	public void mouseMoved(int x,int y){
		// super.mouseMoved(x, y);
		float[] temp=p3.getPlotCoords(x,y);
		if(temp!=null){
			IJ.showStatus("X="+temp[0]+", Y="+temp[1]);
		}
		// coordinates.setText("X="+temp[0]+", Y="+temp[1]);
	}

	public void showList(){
		/*StringBuffer sb=new StringBuffer();
		StringBuffer headings=new StringBuffer();
		int tempnseries=p3.getNSeries();
		float[][] tempxvals=p3.getXValues();
		float[][] tempyvals=p3.getYValues();
		for(int j=0;j<tempnseries;j++){
			headings.append("X"+(j+1)+"\t"+"Y"+(j+1));
			if(j<(tempnseries-1)){
				headings.append("\t");
			}
		}
		for(int i=0;i<p3.getmaxpts();i++){
			for(int j=0;j<tempnseries;j++){
				sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]);
				if(j<(tempnseries-1)){
					sb.append("\t");
				}
			}
			sb.append("\n");
		}
		new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);*/
		showList(true,true,true,true,false,false,true);
	}
	
	public void showList(boolean copyfirstx,boolean copyotherx,boolean copyfirsty,boolean copyothery,boolean copyerrs,boolean copybotherrs,boolean padnan){
		StringBuffer sb = new StringBuffer();
		StringBuffer headings=new StringBuffer();
		float[][] xvals=p3.getXValues();
		float[][] yvals=p3.getYValues();
		int[] npts=p3.getNpts();
		float[][][] errs=p3.getErrors();
		int length=yvals[0].length;
		int nseries=yvals.length;
		//first copy the column titles
		for (int j=0; j<nseries; j++) {
			if(j==0 && copyfirstx) headings.append("X"+(j+1)+"\t");
			if(j==0 && copyfirsty) headings.append("Y"+(j+1)+"\t");
			if(j>0 && copyotherx) headings.append("X"+(j+1)+"\t");
			if(j>0 && copyothery) headings.append("Y"+(j+1)+"\t");
			if(copyerrs && errs!=null){
				if(copybotherrs) headings.append("err1\terr2\t");
				else headings.append("err1\t");
			}
		}
		for (int i=0; i<length; i++) {
			for (int j=0; j<nseries; j++) {
				String xval=""+xvals[j][i];
				String yval=""+yvals[j][i];
				if(padnan){
					if(i>=npts[j]){xval="NaN"; yval="NaN";}
				}
				if(j==0 && copyfirstx) sb.append(xval+"\t");
				if(j==0 && copyfirsty) sb.append(yval+"\t");
				if(j>0 && copyotherx) sb.append(xval+"\t");
				if(j>0 && copyothery) sb.append(yval+"\t");
				if(copyerrs && errs!=null){
					String err1=""+errs[0][j][i];
					String err2=""+errs[1][j][i];
					if(padnan){
						if(i>=npts[j]){err1="NaN"; err2="NaN";}
					}
					if(copybotherrs) sb.append(err1+"\t"+err2+"\t");
					else sb.append(err1+"\t");
				}
				if(j==(nseries-1)){sb.deleteCharAt(sb.length()-1);}
			}
			sb.append("\n");
		}
		new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);
	}

	void saveAsText(){
		SaveDialog sd=new SaveDialog("Save as Text...",getTitle(),".txt");
		String name=sd.getFileName();
		String directory=sd.getDirectory();
		/*FileDialog fd=new FileDialog(this,"Save as Text...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();*/
		saveAsText(directory+name);
	}

	public void saveAsText(String fname){
		PrintWriter pw=null;
		try{
			FileOutputStream fos=new FileOutputStream(fname);
			BufferedOutputStream bos=new BufferedOutputStream(fos);
			pw=new PrintWriter(bos);
		}catch(IOException e){
			IJ.error(""+e);
			return;
		}
		IJ.wait(250); // give system time to redraw ImageJ window
		IJ.showStatus("Saving plot values...");
		int tempnseries=p3.getNSeries();
		float[][] tempxvals=p3.getXValues();
		float[][] tempyvals=p3.getYValues();
		String delim="\t";
		if(fname.endsWith(".csv")) delim=",";
		for(int i=0;i<p3.getmaxpts();i++){
			StringBuffer sb=new StringBuffer();
			for(int j=0;j<tempnseries;j++){
				sb.append(""+tempxvals[j][i]+delim+tempyvals[j][i]);
				if(j<(tempnseries-1)){
					sb.append(delim);
				}
			}
			pw.println(sb.toString());
		}
		pw.close();
	}

	void saveAsBinary(int saveseries,int typeflag){
		/*FileDialog fd=new FileDialog(this,"Save as Binary...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();*/
		SaveDialog sd=new SaveDialog("Save as Binary...",getTitle(),".bin");
		String name=sd.getFileName();
		String directory=sd.getDirectory();
		saveAsBinary(directory+name,saveseries,typeflag);
	}

	public void saveAsBinary(String fname,int saveseries,int typeflag){
		int[] tempnpts=p3.getNpts();
		float[][] tempyvals=p3.getYValues();
		try{
			OutputStream outstream=new BufferedOutputStream(new FileOutputStream(fname));
			for(int i=0;i<tempnpts[saveseries];i++){
				if(typeflag==0){
					int tmp=Float.floatToIntBits(tempyvals[saveseries][i]);
					byte[] dumbyte={(byte)tmp,(byte)(tmp>>8),(byte)(tmp>>16),(byte)(tmp>>24)};
					outstream.write(dumbyte[0]);
					outstream.write(dumbyte[1]);
					outstream.write(dumbyte[2]);
					outstream.write(dumbyte[3]);
				}else{
					if(typeflag==1){
						int tmp=(int)tempyvals[saveseries][i];
						byte[] dumbyte={(byte)tmp,(byte)(tmp>>8),(byte)(tmp>>16),(byte)(tmp>>24)};
						outstream.write(dumbyte[0]);
						outstream.write(dumbyte[1]);
						outstream.write(dumbyte[2]);
						outstream.write(dumbyte[3]);
					}else{
						if(typeflag==2){
							short tmp=(short)tempyvals[saveseries][i];
							byte[] dumbyte={(byte)tmp,(byte)(tmp>>8)};
							outstream.write(dumbyte[0]);
							outstream.write(dumbyte[1]);
						}
					}
				}
			}
			outstream.close();
		}catch(IOException e){
			IJ.error(""+e);
			return;
		}
	}

	void saveAsObject(){
		/*FileDialog fd=new FileDialog(this,"Save as Plot Object...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		String temptitle=getTitle();
		if(!temptitle.endsWith(".pw2")){
			temptitle+=".pw2";
		}
		fd.setFile(temptitle);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();*/
		SaveDialog sd=new SaveDialog("Save as Plot Object...",getTitle(),".pw2");
		String name=sd.getFileName();
		String directory=sd.getDirectory();
		if(name==null||name==""||directory==null||directory=="")
			return;
		if(!name.endsWith(".pw2")){
			if(name.endsWith(".pw")) name+="2";
			else name+=".pw2";
		}
		imp.setTitle(name);
		String dir2=directory.replace("\\","\\\\");
		if(Recorder.record && !IJ.isMacro()) Recorder.record("run","export plot jru v1", "save=["+dir2+name+"]");
		saveAsObject(directory+File.separator+name);
	}

	public void saveAsObject(String filename){
		p3.saveplot2file(filename);
	}

	public void saveAs(String filename,int outtype){
		if(outtype==0)
			saveAsObject(filename);
		if(outtype==1)
			saveAsText(filename);
		if(outtype==2)
			saveAsBinary(filename,0,0);
	}

	void copyToClipboard(){
		Clipboard systemClipboard=null;
		try{
			systemClipboard=getToolkit().getSystemClipboard();
		}catch(Exception e){
			systemClipboard=null;
		}
		if(systemClipboard==null){
			IJ.error("Unable to copy to Clipboard.");
			return;
		}
		IJ.showStatus("Copying plot values...");
		StringBuffer sb=new StringBuffer();
		int tempnseries=p3.getNSeries();
		float[][] tempxvals=p3.getXValues();
		float[][] tempyvals=p3.getYValues();
		for(int i=0;i<p3.getmaxpts();i++){
			for(int j=0;j<tempnseries;j++){
				sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]);
				if(j<(tempnseries-1)){
					sb.append("\t");
				}
			}
			sb.append("\n");
		}
		String text=sb.toString();
		StringSelection contents=new StringSelection(text);
		systemClipboard.setContents(contents,this);
		IJ.showStatus(text.length()+" characters copied to Clipboard");
	}
	
	/*************
	 * this is a more customizable version of copyToClipboard
	 * @param copytitles
	 * @param copyfirstx
	 * @param copyotherx
	 * @param copyfirsty
	 * @param copyothery
	 * @param copyerrs
	 * @param copybotherrs
	 * @param padnan
	 */
	void copyCustom(boolean copytitles,boolean copyfirstx,boolean copyotherx,boolean copyfirsty,boolean copyothery,boolean copyerrs,boolean copybotherrs,boolean padnan) {
		Clipboard systemClipboard = null;
		try {systemClipboard = IJ.getInstance().getToolkit().getSystemClipboard();}
		catch (Exception e) {systemClipboard = null; }
		if (systemClipboard==null)
			{IJ.error("Unable to copy to Clipboard."); return;}
		IJ.showStatus("Copying plot values...");
		StringBuffer sb = new StringBuffer();
		float[][] xvals=p3.getXValues();
		float[][] yvals=p3.getYValues();
		int[] npts=p3.getNpts();
		float[][][] errs=p3.getErrors();
		int length=yvals[0].length;
		int nseries=yvals.length;
		//first copy the column titles
		if(copytitles){
			for (int j=0; j<nseries; j++) {
				if(j==0 && copyfirstx) sb.append("x"+(j+1)+"\t");
				if(j==0 && copyfirsty) sb.append("y"+(j+1)+"\t");
				if(j>0 && copyotherx) sb.append("x"+(j+1)+"\t");
				if(j>0 && copyothery) sb.append("y"+(j+1)+"\t");
				if(copyerrs && errs!=null){
					if(copybotherrs) sb.append("err1\terr2\t");
					else sb.append("err1\t");
				}
			}
			sb.append("\n");
		}
		for (int i=0; i<length; i++) {
			for (int j=0; j<nseries; j++) {
				String xval=""+xvals[j][i];
				String yval=""+yvals[j][i];
				if(padnan){
					if(i>=npts[j]){xval="NaN"; yval="NaN";}
				}
				if(j==0 && copyfirstx) sb.append(xval+"\t");
				if(j==0 && copyfirsty) sb.append(yval+"\t");
				if(j>0 && copyotherx) sb.append(xval+"\t");
				if(j>0 && copyothery) sb.append(yval+"\t");
				if(copyerrs && errs!=null){
					String err1=""+errs[0][j][i];
					String err2=""+errs[1][j][i];
					if(padnan){
						if(i>=npts[j]){err1="NaN"; err2="NaN";}
					}
					if(copybotherrs) sb.append(err1+"\t"+err2+"\t");
					else sb.append(err1+"\t");
				}
				if(j==(nseries-1)){sb.deleteCharAt(sb.length()-1);}
			}
			sb.append("\n");
		}
		String text = sb.toString();
		StringSelection contents = new StringSelection(text);
		systemClipboard.setContents(contents, this);
		IJ.showStatus(text.length() + " characters copied to Clipboard");
	}

	void editPlot(){
		GenericDialog gd=new GenericDialog("Plot Options");
		float[] limits=p3.getLimits();
		gd.addNumericField("x min",limits[0],5,10,null);
		gd.addNumericField("x max",limits[1],5,10,null);
		gd.addNumericField("y min",limits[2],5,10,null);
		gd.addNumericField("y max",limits[3],5,10,null);
		boolean[] logs=p3.getLogAxes();
		gd.addCheckbox("Log x?",logs[0]);
		gd.addCheckbox("Log y?",logs[1]);
		gd.addStringField("x label",p3.getxLabel());
		gd.addStringField("y label",p3.getyLabel());
		gd.addNumericField("Grid Whiteness",p3.gridColor.getRed(),0);
		boolean delsel=false;
		gd.addCheckbox("Delete Selected",delsel);
		boolean ascalex=false;
		gd.addCheckbox("AutoScale x",ascalex);
		boolean ascaley=false;
		gd.addCheckbox("AutoScale y",ascalex);
		boolean doscaleroi=false;
		gd.addCheckbox("Scale Roi",doscaleroi);
		gd.addNumericField("Magnification",p3.getmagnification(),5,10,null);
		gd.addNumericField("Mag Ratio",p3.getmagratio(),5,10,null);
		gd.addNumericField("Selected Series",p3.getSelected(),0);
		gd.addCheckbox("Edit Selected Series",false);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		float[] newlims={(float)gd.getNextNumber(),(float)gd.getNextNumber(),(float)gd.getNextNumber(),(float)gd.getNextNumber()};
		if(limits[0]!=newlims[0] || limits[1]!=newlims[1] || limits[2]!=newlims[2] || limits[3]!=newlims[3]){
			p3.setLimits(newlims);
			if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.setLimits",""+newlims[0]+","+newlims[1]+","+newlims[2]+","+newlims[3]);
		}
		boolean[] newlogs={gd.getNextBoolean(),gd.getNextBoolean()};
		if(p3.getLogAxes()[0]!=newlogs[0] || p3.getLogAxes()[1]!=newlogs[1]){
			p3.setLogAxes(newlogs[0],newlogs[1]);
			if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.setLogAxes",""+newlogs[0]+","+newlogs[1]);
		}
		String[] newlabs={gd.getNextString(),gd.getNextString()};
		if(!newlabs[0].equals(p3.getxLabel()) || !newlabs[1].equals(p3.getyLabel())){
			p3.setxLabel(newlabs[0]);
			p3.setyLabel(newlabs[1]);
			if(Recorder.record && !IJ.isMacro()){
				Recorder.record("Ext.setXLabel",newlabs[0]);
				Recorder.record("Ext.setYLabel",newlabs[1]);
			}
		}
		int temp=(int)gd.getNextNumber();
		if(temp!=p3.gridColor.getRed()){
			p3.setGridWhiteness(temp);
			if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.setGridWhiteness",temp);
		}
		delsel=gd.getNextBoolean();
		ascalex=gd.getNextBoolean();
		ascaley=gd.getNextBoolean();
		doscaleroi=gd.getNextBoolean();
		p3.setmagnification((float)gd.getNextNumber());
		p3.setmagratio((float)gd.getNextNumber());
		p3.selectSeries((int)gd.getNextNumber());
		if(gd.getNextBoolean()){
			editSeries(p3.getSelected());
		}
		if(delsel){
			delsel=false;
			p3.deleteSeries(p3.getSelected(),false);
			Recorder.record("Ext.deleteSelected");
		}
		if(ascalex){
			ascalex=false;
			p3.xautoscale();
			if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.autoscaleX");
		}
		if(ascaley){
			ascaley=false;
			p3.yautoscale();
			if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.autoscaleY");
		}
		if(doscaleroi){
			doscaleroi=false;
			scaleroi();
			if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.scaleROI");
		}
		updatePlot();
	}

	private void editSeries(int series){
		if(series>=0){
			GenericDialog gd5=new GenericDialog("Edit Series "+series);
			int[] colors=p3.getColors();
			int[] shapes=p3.getShapes();
			int[] shapesizes=p3.getShapeSizes();
			String[] annotations=p3.getAnnotations();
			gd5.addChoice("Series Color",Plot4.color_names,Plot4.color_names[colors[series]%8]);
			gd5.addChoice("Series Shape",Plot4.shape_names,Plot4.shape_names[shapes[series]]);
			gd5.addNumericField("Shape Size",shapesizes[series],0);
			if(annotations!=null) gd5.addStringField("Annotation",annotations[series]);
			gd5.showDialog();
			if(gd5.wasCanceled()){
				return;
			}
			colors[series]=gd5.getNextChoiceIndex();
			shapes[series]=gd5.getNextChoiceIndex();
			shapesizes[series]=(int)gd5.getNextNumber();
			if(annotations!=null) annotations[series]=gd5.getNextString();
		}else{
			IJ.showMessage("You must select a series to edit");
		}
	}

	public void lostOwnership(Clipboard clipboard,Transferable contents){
	}

	public void actionPerformed(ActionEvent e){
		Object b=e.getSource();
		if(b==list){
			int mods=e.getModifiers();
			//byte[] mods2= {(byte)mods,(byte)(mods>>8),(byte)(mods>>16),(byte)(mods>>24)};
			//IJ.log("Modifiers = "+jdataio.printHexBytes(mods2));
			//if((mods&ActionEvent.CTRL_MASK)!=0 || (mods&ActionEvent.ALT_MASK)!=0 || (mods&ActionEvent.SHIFT_MASK)!=0) showList(true,true,true,true,true,false,true);
			//getModifiers should only be non-zero if a key was held during the click
			if(mods!=0) showList(true,true,true,true,true,false,true);
			else showList();
		}else{
			if(b==save){
				GenericDialog gd=new GenericDialog("Save Options");
				String[] savechoice={"Text","Binary","Plot Object"};
				gd.addChoice("File Type?",savechoice,savechoice[2]);
				int saveseries=0;
				gd.addNumericField("Save Series # (for binary)",saveseries,0);
				String[] binarytypechoice={"Float","Integer","Short"};
				gd.addChoice("Binary File Type?",binarytypechoice,binarytypechoice[0]);
				gd.showDialog();
				if(gd.wasCanceled()){
					return;
				}
				int choiceindex=gd.getNextChoiceIndex();
				if(choiceindex==0){
					saveAsText();
				}else{
					if(choiceindex==1){
						saveseries=(int)gd.getNextNumber();
						saveAsBinary(saveseries,gd.getNextChoiceIndex());
					}else{
						saveAsObject();
					}
				}
			}else{
				if(b==editbutton){
					editPlot();
				}else{
					if(b==copy){
						if(e.getModifiers()!=0) copyCustom(true,true,true,true,true,true,false,true);
						else copyToClipboard();
					}else{
						if(b==selbutton){
							p3.selectSeries(p3.getSelected()+1);
							String annot="";
							if(p3.getAnnotations()!=null) annot=":"+(p3.getAnnotations())[p3.getSelected()];
							IJ.showStatus("Series "+(p3.getSelected()+1)+" of "+p3.getNSeries()+" Selected"+annot);
							if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.selectSeries",p3.getSelected());
							updatePlot();
						}else{
							p3.selectSeries(p3.getSelected()-1);
							String annot="";
							if(p3.getAnnotations()!=null) annot=":"+(p3.getAnnotations())[p3.getSelected()];
							IJ.showStatus("Series "+(p3.getSelected()+1)+" of "+p3.getNSeries()+" Selected"+annot);
							if(Recorder.record && !IJ.isMacro()) Recorder.record("Ext.selectSeries",p3.getSelected());
							updatePlot();
						}
					}
				}
			}
		}
	}

	public void selectSeries(int series){
		p3.selectSeries(series);
		updatePlot();
	}

	public int getSelected(){
		return p3.getSelected();
	}

	public float[][] getXValues(){
		return p3.getXValues();
	}

	public float[] getXValues(int series){
		return p3.getXValues(series);
	}

	public float[][] getYValues(){
		return p3.getYValues();
	}

	public float[] getYValues(int series){
		return p3.getYValues(series);
	}

	public float[][][] getErrors(){
		return p3.getErrors();
	}

	public float[] getErrors(int series,boolean upper){
		return p3.getErrors(series,upper);
	}

	public boolean getShowErrors(){
		return p3.getShowErrors();
	}

	public String[] getAllLabels(){
		String[] temp={imp.getTitle(),p3.getxLabel(),p3.getyLabel()};
		return temp;
	}

	public String getPlotTitle(){
		return imp.getTitle();
	}

	public String getxLabel(){
		return p3.getxLabel();
	}

	public String getyLabel(){
		return p3.getyLabel();
	}

	public int[] getNpts(){
		return p3.getNpts();
	}

	public int getNSeries(){
		return p3.getNSeries();
	}

	public float[] getLimits(){
		return p3.getLimits();
	}

	public boolean[] getLogAxes(){
		return p3.getLogAxes();
	}

	public int[] getShapes(){
		return p3.getShapes();
	}

	public int[] getColors(){
		return p3.getColors();
	}
	
	public String[] getAnnotations(){
		return p3.getAnnotations();
	}

	public ColorProcessor getProcessor(){
		return p3.getProcessor();
	}

	public Plot4 getPlot(){
		return p3;
	}

	public Plot4 getPlotCopy(){
		return p3.getCopy();
	}

}
