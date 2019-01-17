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
import ij.process.ColorProcessor;
import ij.text.TextWindow;
import ij.util.Tools;

import java.awt.Button;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Frame;
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

/**
 * This class is a plot window that is independent of ImageJ's ImageWindow class
 * It uses a plot panel to display the plot This should provide compatibility
 * with ImageJ2 A couple of issues: 1. The ImageJ copy will no longer work 2.
 * The image will no longer be listed in any ImageJ window lists 3. The ImageJ
 * rename function is no longer useable and must now be listed from the drop
 * down 4. This window is no longer selectable by IJ macros not sure how it will
 * function with IJ2 macros in principle some of the features can be retained by
 * complex reflection routines Adaptations were done by Jay Unruh
 * (jru@stowers.org)
 */
public class PlotWindow5 extends Frame implements ActionListener,ClipboardOwner{

	private ImagePlus imp;
	private Button list,save,copy,editbutton,selbutton,seldbutton;
	// private Label coordinates;
	private static String defaultDirectory=null;
	public Plot4 p3;
	private static ColorProcessor cp;

	public PlotWindow5(String title1,String xLabel1,String yLabel1,float[][] xValues1,float[][] yValues1,Object npts1){
		super(title1);
		p3=new Plot4(xLabel1,yLabel1,xValues1,yValues1,npts1);
	}

	public PlotWindow5(String title1,String xLabel1,String yLabel1,float[][] yValues1,Object npts1){
		super(title1);
		p3=new Plot4(xLabel1,yLabel1,yValues1,npts1);
	}

	public PlotWindow5(String title1,String xLabel1,String yLabel1,float[] xValues1,float[] yValues1){
		super(title1);
		p3=new Plot4(xLabel1,yLabel1,xValues1,yValues1);
	}

	public PlotWindow5(String title1,String xLabel1,String yLabel1,double[] xValues1,double[] yValues1){
		this(title1,xLabel1,yLabel1,Tools.toFloat(xValues1),Tools.toFloat(yValues1));
	}

	public PlotWindow5(String title1,String xLabel1,String yLabel1,float[] yValues1){
		super(title1);
		p3=new Plot4(xLabel1,yLabel1,yValues1);
	}

	public PlotWindow5(String title1,Plot4 plot){
		super(title1);
		p3=plot;
	}

	public void initui(){
		setResizable(true);
	}

	public PlotWindow5 getInstance(){
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
		float[][] temp=new float[2][];
		temp[0]=errorBars;
		temp[1]=errorBars;
		addErrors(temp);
	}

	public void setShowErrors(boolean showerrors){
		p3.setShowErrors(showerrors);
	}

	public void updatePlot(){
		cp=p3.getProcessor();
		imp.setProcessor(null,cp);
		/*
		 * Calibration cal=imp.getCalibration(); float[] limits=p3.getLimits();
		 * cal
		 * .pixelWidth=(double)(limits[1]-limits[0])/(p3.getmagnification()*Plot4
		 * .WIDTH);
		 * cal.pixelHeight=(double)(limits[3]-limits[2])/(p3.getmagnification
		 * ()*Plot4.HEIGHT); //cal.xOrigin=-(double)xmin/cal.pixelWidth;
		 * //cal.yOrigin=(double)ymin/cal.pixelHeight+255.0;
		 * cal.setInvertY(true);
		 */
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

	void showList(){
		StringBuffer sb=new StringBuffer();
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
		new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);
	}

	void saveAsText(){
		FileDialog fd=new FileDialog(this,"Save as Text...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();
		PrintWriter pw=null;
		try{
			FileOutputStream fos=new FileOutputStream(directory+name);
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
		for(int i=0;i<p3.getmaxpts();i++){
			StringBuffer sb=new StringBuffer();
			for(int j=0;j<tempnseries;j++){
				sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]);
				if(j<(tempnseries-1)){
					sb.append("\t");
				}
			}
			pw.println(sb.toString());
		}
		pw.close();
	}

	void saveAsBinary(int saveseries,int typeflag){
		FileDialog fd=new FileDialog(this,"Save as Binary...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();
		int[] tempnpts=p3.getNpts();
		float[][] tempyvals=p3.getYValues();
		try{
			OutputStream outstream=new BufferedOutputStream(new FileOutputStream(directory+name));
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
		FileDialog fd=new FileDialog(this,"Save as Plot Object...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		String temptitle=getTitle();
		if(!temptitle.endsWith(".pw")){
			temptitle+=".pw";
		}
		fd.setFile(temptitle);
		/*
		 * fd.setFilenameFilter( new FilenameFilter(){ public boolean
		 * accept(File dir,String name){ if(name.endsWith(".pw")){ return true;
		 * } else { return false; } } } );
		 */
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();
		if(name==null||name==""||directory==null||directory=="")
			return;
		imp.setTitle(name);
		saveAsObject(directory+File.separator+name);
	}

	public void saveAsObject(String filename){
		p3.saveplot2file(filename);
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
		limits[0]=(float)gd.getNextNumber();
		limits[1]=(float)gd.getNextNumber();
		limits[2]=(float)gd.getNextNumber();
		limits[3]=(float)gd.getNextNumber();
		p3.setLimits(limits);
		logs[0]=gd.getNextBoolean();
		logs[1]=gd.getNextBoolean();
		p3.setLogAxes(logs[0],logs[1]);
		p3.setxLabel(gd.getNextString());
		p3.setyLabel(gd.getNextString());
		int temp=(int)gd.getNextNumber();
		if(temp>255)
			temp=255;
		if(temp<0)
			temp=0;
		p3.gridColor=new Color(temp,temp,temp);
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
		}
		if(ascalex){
			ascalex=false;
			p3.xautoscale();
		}
		if(ascaley){
			ascaley=false;
			p3.yautoscale();
		}
		if(doscaleroi){
			doscaleroi=false;
			scaleroi();
		}
		updatePlot();
	}

	private void editSeries(int series){
		if(series>=0){
			GenericDialog gd5=new GenericDialog("Edit Series "+series);
			int[] colors=p3.getColors();
			int[] shapes=p3.getShapes();
			gd5.addChoice("Series Color",Plot4.color_names,Plot4.color_names[colors[series]%8]);
			gd5.addChoice("Series Shape",Plot4.shape_names,Plot4.shape_names[shapes[series]]);
			gd5.showDialog();
			if(gd5.wasCanceled()){
				return;
			}
			colors[series]=gd5.getNextChoiceIndex();
			shapes[series]=gd5.getNextChoiceIndex();
		}else{
			IJ.showMessage("You must select a series to edit");
		}
	}

	public void lostOwnership(Clipboard clipboard,Transferable contents){
	}

	public void actionPerformed(ActionEvent e){
		Object b=e.getSource();
		if(b==list){
			showList();
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
						copyToClipboard();
					}else{
						if(b==selbutton){
							p3.selectSeries(p3.getSelected()+1);
							IJ.showStatus("Series "+(p3.getSelected()+1)+" of "+p3.getNSeries()+" Selected");
							updatePlot();
						}else{
							p3.selectSeries(p3.getSelected()-1);
							IJ.showStatus("Series "+(p3.getSelected()+1)+" of "+p3.getNSeries()+" Selected");
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
