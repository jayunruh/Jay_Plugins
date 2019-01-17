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
import ij.gui.Roi;
import ij.io.SaveDialog;
import ij.plugin.frame.Recorder;
import ij.process.ColorProcessor;
import ij.text.TextWindow;

import java.awt.Button;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Label;
import java.awt.Panel;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * This class is an extended ImageWindow that displays line graphs. This class
 * was adapted from the PlotWindow class original to ImageJ Adaptations were
 * done by Jay Unruh (jru@stowers.org)
 */
public class PlotWindowColumn extends ImageWindow implements ActionListener,ClipboardOwner,MouseMotionListener,MouseListener{

	private Button list,save,copy,editbutton;
	private Label coordinates;
	public boolean savehist,showerrors;
	protected static String defaultDirectory=null;
	public PlotColumn p3;
	protected static ColorProcessor cp;

	public PlotWindowColumn(String title1,String xLabel1,String yLabel1,float[] xValues1,int color){
		super(createImage(title1));
		p3=new PlotColumn(xLabel1,yLabel1,xValues1,color);
		savehist=false;
		imp.getCanvas().addMouseMotionListener(this);
	}
	
	public PlotWindowColumn(String title1,String xLabel1,String yLabel1,float[][] xValues1,int[] npts,int color){
		super(createImage(title1));
		p3=new PlotColumn(xLabel1,yLabel1,xValues1,npts,color);
		savehist=false;
		imp.getCanvas().addMouseMotionListener(this);
	}


	public PlotWindowColumn(String title1,PlotColumn plot1){
		super(createImage(title1));
		p3=plot1;
		savehist=false;
		imp.getCanvas().addMouseMotionListener(this);
	}

	public PlotWindowColumn(ImagePlus imp,PlotColumn p3){
		// turns the passed ImagePlus into a plot window
		// used for bioformats which wants to make the window for me
		super(imp);
		int width=PlotColumn.WIDTH+PlotColumn.LEFT_MARGIN+PlotColumn.RIGHT_MARGIN;
		int height=PlotColumn.HEIGHT+PlotColumn.TOP_MARGIN+PlotColumn.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		imp.setProcessor(cp);
		// cp=(ColorProcessor)this.imp.getProcessor();
		this.p3=p3;
		imp.getCanvas().addMouseMotionListener(this);
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
		int width=PlotColumn.WIDTH+PlotColumn.LEFT_MARGIN+PlotColumn.RIGHT_MARGIN;
		int height=PlotColumn.HEIGHT+PlotColumn.TOP_MARGIN+PlotColumn.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		return new ImagePlus(title1,cp);
	}

	public void setBinSize(float binsize){
		p3.setBinSize(binsize);
		updatePlot();
	}

	public float getBinSize(){
		return p3.getBinSize();
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

	/*public int[] getroiindices(){
		Roi roi=imp.getRoi();
		if(roi instanceof PolygonRoi){
			// Polygon poly=roi.getPolygon();
			// return p3.getpolygonindices(poly);
			return null;
		}else{
			Rectangle rect=imp.getProcessor().getRoi();
			return p3.getrectindices(rect);
		}
	}

	public void scaleroi(){
		Rectangle rect=imp.getProcessor().getRoi();
		if(rect!=null){
			p3.scalerect(rect);
			imp.killRoi();
			//updatePlot();
			//IJ.log("plot updated");
		}
	}*/

	public void updateSeries(float[] xValues1,int series,boolean rescale){
		p3.updateSeries(xValues1,series,rescale);
		updatePlot();
	}


	public void updatePlot(){
		cp=p3.getProcessor();
		imp.setProcessor(null,cp);
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
		coordinates = new Label("Density = intialized size large");
		coordinates.setFont(new Font("Monospaced", Font.PLAIN, 12));
		buttons.add(coordinates);
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
		//float temp=p3.getPlotCoords(x);
		//IJ.showStatus("X="+temp);
	}

	public void mouseMoved(MouseEvent e){
		update_density();
	}

	public void mouseDragged(MouseEvent e){
		update_density();
	}

	public void mouseClicked(MouseEvent e){
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
	}

	public void mouseReleased(MouseEvent e){
		update_density();
	}
	
	public void update_density(){
		Roi roi=imp.getRoi();
		if(roi!=null){
			//coordinates.setText("Density = "+jutils.formatted_string(p3.get_score(roi)));
		}else{
			coordinates.setText("Density = "+jutils.formatted_string(0.0f));
		}
	}

	void showList(){
		StringBuffer sb=new StringBuffer();
		StringBuffer headings=new StringBuffer();
		if(!savehist){
			float[][] tempxvals=p3.getXValues();
			int[] npts=p3.getNpts();
			int maxpts=p3.getmaxpts();
			int nser=p3.getNSeries();
			for(int i=0;i<nser;i++) headings.append("X"+(i+1)+"\t");
			headings.deleteCharAt(headings.length()-1);
			for(int j=0;j<maxpts;j++){
				for(int i=0;i<nser;i++){
					if(j<npts[i]) sb.append(""+tempxvals[i][j]+"\t");
					else sb.append("\t");
				}
				sb.deleteCharAt(sb.length()-1);
				sb.append("\n");
			}
			sb.deleteCharAt(sb.length()-1);
		}else{
			float[][] hist=p3.getHistogram();
			int nser=p3.getNSeries();
			for(int i=0;i<nser;i++) headings.append("X"+(i+1)+"\tFrequency"+(i+1)+"\t");
			headings.deleteCharAt(headings.length()-1);
			for(int i=0;i<hist[0].length;i++){
				for(int j=0;j<hist.length;j++){
					sb.append(""+hist[j][i]+"\t");
				}
				sb.deleteCharAt(sb.length()-1);
				sb.append("\n");
			}
			sb.deleteCharAt(sb.length()-1);
		}
		new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);
	}

	void saveAsText(){
		/*FileDialog fd=new FileDialog(this,"Save as Text...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();*/
		SaveDialog sd=new SaveDialog("Save as Text...",getTitle(),"txt");
		String name=sd.getFileName();
		String directory=sd.getDirectory();
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
		StringBuffer headings=new StringBuffer();
		if(!savehist){
			float[][] tempxvals=p3.getXValues();
			int[] npts=p3.getNpts();
			int maxpts=p3.getmaxpts();
			int nser=p3.getNSeries();
			for(int i=0;i<nser;i++) headings.append("X"+(i+1));
			pw.println(headings.toString());
			for(int j=0;j<maxpts;j++){
				StringBuffer sb=new StringBuffer();
				for(int i=0;i<nser;i++){
					if(j<npts[i]) sb.append(""+tempxvals[i]+"\t");
					else sb.append("\t");
				}
				sb.deleteCharAt(sb.length()-1);
				pw.println(sb.toString());
			}
		}else{
			float[][] hist=p3.getHistogram();
			int nser=p3.getNSeries();
			for(int i=0;i<nser;i++) headings.append("X"+(i+1)+"\tFrequency"+(i+1));
			pw.println(headings.toString());
			for(int i=0;i<hist[0].length;i++){
				StringBuffer sb=new StringBuffer();
				for(int j=0;j<hist.length;j++){
					sb.append(""+hist[j][i]+"\t");
				}
				sb.deleteCharAt(sb.length()-1);
				pw.println(sb.toString());
			}
		}
		pw.close();
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

	public void copyToClipboard(){
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
		if(!savehist){
			
			float[][] tempxvals=p3.getXValues();
			int[] npts=p3.getNpts();
			int maxpts=p3.getmaxpts();
			int nser=p3.getNSeries();
			for(int j=0;j<maxpts;j++){
				for(int i=0;i<nser;i++){
					if(j<npts[i]) sb.append(""+tempxvals[i][j]+"\t");
					else sb.append("\t");
				}
				sb.deleteCharAt(sb.length()-1);
				sb.append("\n");
			}
			sb.deleteCharAt(sb.length()-1);
		}else{
			float[][] tempxvals=p3.getHistogram();
			for(int i=0;i<tempxvals[0].length;i++){
				sb.append(""+tempxvals[0][i]+"\t"+tempxvals[1][i]+"\n");
			}
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
		gd.addCheckbox("Log x?",p3.getLogAxes()[0]);
		gd.addCheckbox("Log y?",p3.getLogAxes()[1]);
		gd.addStringField("x label",p3.getxLabel());
		gd.addStringField("y label",p3.getyLabel());
		boolean delsel=false;
		gd.addCheckbox("Delete Selected",delsel);
		boolean ascalex=false;
		gd.addCheckbox("AutoScale x",ascalex);
		boolean ascaley=false;
		gd.addCheckbox("AutoScale y",ascalex);
		boolean doscaleroi=false;
		gd.addCheckbox("Scale Roi",doscaleroi);
		gd.addNumericField("Magnification",p3.getmagnification(),5,10,null);
		gd.addNumericField("Bin Size",p3.getBinSizeUnits(),5,10,null);
		gd.addCheckbox("Save Histogram?",savehist);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		limits[0]=(float)gd.getNextNumber();
		limits[1]=(float)gd.getNextNumber();
		limits[2]=(float)gd.getNextNumber();
		limits[3]=(float)gd.getNextNumber();
		p3.setLimits(limits);
		//IJ.log(""+p3.getLimits()[2]);
		p3.setLogAxes(gd.getNextBoolean(),gd.getNextBoolean());
		//IJ.log(""+p3.getLimits()[2]);
		p3.setxLabel(gd.getNextString());
		p3.setyLabel(gd.getNextString());
		delsel=gd.getNextBoolean();
		ascalex=gd.getNextBoolean();
		ascaley=gd.getNextBoolean();
		doscaleroi=gd.getNextBoolean();
		p3.setmagnification((float)gd.getNextNumber());
		//IJ.log(""+p3.getLimits()[2]);
		float binsize=(float)gd.getNextNumber();
		//if(binsize!=p3.getBinSizeUnits()){
		//	p3.setBinSizeUnits(binsize);
		//}
		savehist=gd.getNextBoolean();
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
			//scaleroi();
		}
		//IJ.log(""+p3.getLimits()[2]);
		updatePlot();
	}

	public void lostOwnership(Clipboard clipboard,Transferable contents){
	}

	public void actionPerformed(ActionEvent e){
		Object b=e.getSource();
		if(b==list){
			showList();
		}else{
			if(b==save){
				//IJ.log("test");
				GenericDialog gd=new GenericDialog("Save Options");
				String[] savechoice={"Text","Plot Object"};
				gd.addChoice("File Type?",savechoice,savechoice[1]);
				//String[] binarytypechoice={"Float","Integer","Short"};
				//gd.addChoice("Binary File Type?",binarytypechoice,binarytypechoice[0]);
				gd.showDialog(); if(gd.wasCanceled()) return;
				int choiceindex=gd.getNextChoiceIndex();
				if(choiceindex==0){
					saveAsText();
				}else{
					saveAsObject();
				}
			}else{
				if(b==editbutton){
					editPlot();
				}else{
					if(b==copy){
						copyToClipboard();
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
	
	public float[] getYValues(int series){
		return p3.getYValues(series);
	}

	public float[] getErrors(int series,boolean upper){
		return p3.getErrors(series,upper);
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

	public float[] getLimits(){
		return p3.getLimits();
	}

	public boolean[] getLogAxes(){
		return p3.getLogAxes();
	}
	
	public String[] getAnnotations(){
		return p3.getAnnotations();
	}

	public ColorProcessor getProcessor(){
		return p3.getProcessor();
	}

	public PlotColumn getPlot(){
		return p3;
	}

}
