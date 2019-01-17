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
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.plugin.frame.Recorder;
import ij.process.ColorProcessor;
import ij.text.TextWindow;

import java.awt.Button;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Panel;
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
import java.io.PrintWriter;

public class PlotWindow3D extends ImageWindow implements ActionListener,ClipboardOwner{
	private Button list,save,copy,editbutton,selbutton,seldbutton,toimagebutton;
	private Button rotrightbutton,rotleftbutton,rotupbutton,rotdownbutton,rotcounterbutton,rotclockbutton,rotresetbutton;
	// private Label coordinates;
	private static String defaultDirectory=null;
	private static String title;
	public Plot3D p3;
	private static ColorProcessor cp;

	public PlotWindow3D(String title1,String xLabel1,String yLabel1,String zLabel1,float[][] xValues1,float[][] yValues1,float[][][] zValues1,Object npts1){
		super(createImage(title1));
		p3=new Plot3D(xLabel1,yLabel1,zLabel1,xValues1,yValues1,zValues1,npts1);
	}

	public PlotWindow3D(String title1,String xLabel1,String yLabel1,String zLabel1,float[] xValues1,float[] yValues1,float[][] zValues1){
		super(createImage(title1));
		p3=new Plot3D(xLabel1,yLabel1,zLabel1,xValues1,yValues1,zValues1);
	}

	public PlotWindow3D(String title1,String xLabel1,String yLabel1,String zLabel1,float[][] zValues1,int startxy){
		super(createImage(title1));
		p3=new Plot3D(xLabel1,yLabel1,zLabel1,zValues1,startxy);
	}

	public PlotWindow3D(String title1,String xLabel1,String yLabel1,String zLabel1,float[][][] zValues1,int startxy,Object npts1){
		super(createImage(title1));
		p3=new Plot3D(xLabel1,yLabel1,zLabel1,zValues1,startxy,npts1);
	}

	public PlotWindow3D(String title1,Plot3D p3){
		super(createImage(title1));
		this.p3=p3;
	}

	public PlotWindow3D(ImagePlus imp,Plot3D p3){
		// turns the passed ImagePlus into a plot window
		// used for bioformats which wants to make the window for me
		super(imp);
		int width=Plot3D.WIDTH+Plot3D.LEFT_MARGIN+Plot3D.RIGHT_MARGIN;
		int height=Plot3D.HEIGHT+Plot3D.TOP_MARGIN+Plot3D.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		imp.setProcessor(cp);
		// cp=(ColorProcessor)this.imp.getProcessor();
		this.p3=p3;
	}

	static ImagePlus createImage(String title1){
		int width=Plot3D.WIDTH+Plot3D.LEFT_MARGIN+Plot3D.RIGHT_MARGIN;
		int height=Plot3D.HEIGHT+Plot3D.TOP_MARGIN+Plot3D.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		title=title1;
		return new ImagePlus(title,cp);
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1,double zMin1,double zMax1){
		p3.setLimits(xMin1,xMax1,yMin1,yMax1,zMin1,zMax1);
		updatePlot();
	}

	public void setLimits(float[] limits){
		p3.setLimits(limits);
		updatePlot();
	}

	/** Sets the x-axis and y-axis range. */
	public void setLogAxes(boolean logx1,boolean logy1,boolean logz1){
		p3.setLogAxes(logx1,logy1,logz1);
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

	public void zautoscale(){
		p3.zautoscale();
		updatePlot();
	}

	public void updateSeries(float[] xValues1,float[] yValues1,float[][] zValues1,int series,boolean rescale){
		p3.updateSeries(xValues1,yValues1,zValues1,series,rescale);
		updatePlot();
	}

	public void updateSeries(float[][] zValues1,int series,boolean rescale){
		p3.updateSeries(zValues1,series,rescale);
		updatePlot();
	}

	public void deleteSeries(int series,boolean rescale){
		p3.deleteSeries(series,rescale);
		updatePlot();
	}

	public void addPoints(float[] xValues1,float[] yValues1,float[][] zValues1,boolean rescale){
		p3.addPoints(xValues1,yValues1,zValues1,rescale);
		updatePlot();
	}

	public void addPoints(float[][] zValues1,boolean rescale,int startxy){
		p3.addPoints(zValues1,rescale,startxy);
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
		selbutton=new Button("Select+");
		selbutton.addActionListener(this);
		buttons.add(selbutton);
		seldbutton=new Button("Select-");
		seldbutton.addActionListener(this);
		buttons.add(seldbutton);
		toimagebutton=new Button("To Image");
		toimagebutton.addActionListener(this);
		buttons.add(toimagebutton);
		Panel rotbuttons=new Panel();
		rotbuttons.setLayout(new FlowLayout(FlowLayout.RIGHT));
		rotrightbutton=new Button(">");
		rotrightbutton.addActionListener(this);
		rotbuttons.add(rotrightbutton);
		rotleftbutton=new Button("<");
		rotleftbutton.addActionListener(this);
		rotbuttons.add(rotleftbutton);
		rotupbutton=new Button("^");
		rotupbutton.addActionListener(this);
		rotbuttons.add(rotupbutton);
		rotdownbutton=new Button("v");
		rotdownbutton.addActionListener(this);
		rotbuttons.add(rotdownbutton);
		rotclockbutton=new Button("Clock+");
		rotclockbutton.addActionListener(this);
		rotbuttons.add(rotclockbutton);
		rotcounterbutton=new Button("Clock-");
		rotcounterbutton.addActionListener(this);
		rotbuttons.add(rotcounterbutton);
		rotresetbutton=new Button("Reset Rot.");
		rotresetbutton.addActionListener(this);
		rotbuttons.add(rotresetbutton);
		// coordinates = new Label("X=12345678, Y=12345678");
		// coordinates.setFont(new Font("Monospaced", Font.PLAIN, 12));
		// buttons.add(coordinates);
		add(buttons);
		add(rotbuttons);
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
		super.mouseMoved(x,y);
		/*
		 * float[] temp=getPlotCoords(x,y);
		 * coordinates.setText("X="+temp[0]+", Y="+temp[1]);
		 */
	}

	void showList(){
		if(p3 instanceof Traj3D){
			StringBuffer sb=new StringBuffer();
			StringBuffer headings=new StringBuffer();
			int maxpts=p3.getmaxxpts();
			float[][] tempxvals=p3.getXValues();
			float[][] tempyvals=p3.getYValues();
			float[][] tempzvals=p3.getZValues()[0];
			int[] npts=p3.getNpts()[0];
			int nser=p3.getNSeries();
			for(int i=0;i<nser;i++){
				headings.append("x"+i+"\ty"+i+"\tz"+i);
				if(i<nser-1)
					headings.append("\t");
			}
			for(int i=0;i<maxpts;i++){
				for(int j=0;j<nser;j++){
					if(i>=npts[j]) sb.append("NaN\tNaN\tNaN");
					else sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]+"\t"+tempzvals[j][i]);
					if(j<(nser-1))
						sb.append("\t");
				}
				sb.append("\n");
			}
			new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);
		}else{
			StringBuffer sb=new StringBuffer();
			StringBuffer headings=new StringBuffer();
			int tempmaxypts=p3.getmaxypts();
			float[][][] tempzvals=p3.getZValues();
			for(int j=0;j<tempmaxypts;j++){
				headings.append("y"+j);
				if(j<(tempmaxypts-1)){
					headings.append("\t");
				}
			}
			for(int i=0;i<p3.getNSeries();i++){
				for(int j=0;j<p3.getmaxxpts();j++){
					for(int k=0;k<tempmaxypts;k++){
						sb.append(""+tempzvals[i][j][k]);
						if(k<(tempmaxypts-1)){
							sb.append("\t");
						}
					}
					sb.append("\n");
				}
				sb.append("\n");
			}
			new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);
		}
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
		SaveDialog sd=new SaveDialog("Save as Text...",getTitle(),".txt");
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
		int maxypts=p3.getmaxypts();
		float[][][] zValues=p3.getZValues();
		for(int i=0;i<p3.getNSeries();i++){
			for(int j=0;j<p3.getmaxxpts();j++){
				StringBuffer sb=new StringBuffer();
				for(int k=0;k<maxypts;k++){
					sb.append(""+zValues[i][j][k]);
					if(k<(maxypts-1)){
						sb.append("\t");
					}
				}
				pw.println(sb.toString());
			}
			pw.println("\n");
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

	void saveToImage(){
		if(p3 instanceof Traj3D)
			return;
		float[] xvals=p3.getXValues(0);
		float pixsize=xvals[1]-xvals[0];
		int nseries=p3.getNSeries();
		int width=p3.getmaxxpts();
		int height=p3.getmaxypts();
		ImageStack stack=new ImageStack(width,height);
		for(int i=0;i<nseries;i++){
			float[] temp=new float[width*height];
			float[][] zvals=p3.getZValues(i);
			for(int j=0;j<height;j++){
				for(int k=0;k<width;k++){
					temp[k+j*width]=zvals[k][j];
				}
			}
			stack.addSlice("",temp);
		}
		ImagePlus imp10=new ImagePlus(title,stack);
		Calibration cal=imp10.getCalibration().copy();
		cal.pixelWidth=pixsize;
		cal.pixelHeight=pixsize;
		imp10.setCalibration(cal);
		imp10.show();
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
		StringBuffer headings=new StringBuffer();
		if(p3 instanceof Traj3D){
			int maxpts=p3.getmaxxpts();
			float[][] tempxvals=p3.getXValues();
			float[][] tempyvals=p3.getYValues();
			float[][] tempzvals=p3.getZValues()[0];
			int[] npts=p3.getNpts()[0];
			int nser=p3.getNSeries();
			for(int i=0;i<nser;i++){
				headings.append("x"+i+"\ty"+i+"\tz"+i);
				if(i<nser-1)
					headings.append("\t");
			}
			for(int i=0;i<maxpts;i++){
				for(int j=0;j<nser;j++){
					if(i>=npts[j]) sb.append("NaN\tNaN\tNaN");
					else sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]+"\t"+tempzvals[j][i]);
					if(j<(nser-1))
						sb.append("\t");
				}
				sb.append("\n");
			}
		}else{
			int tempmaxypts=p3.getmaxypts();
			float[][][] tempzvals=p3.getZValues();
			for(int j=0;j<tempmaxypts;j++){
				headings.append("y"+j);
				if(j<(tempmaxypts-1)){
					headings.append("\t");
				}
			}
			for(int i=0;i<p3.getNSeries();i++){
				for(int j=0;j<p3.getmaxxpts();j++){
					for(int k=0;k<tempmaxypts;k++){
						sb.append(""+tempzvals[i][j][k]);
						if(k<(tempmaxypts-1)){
							sb.append("\t");
						}
					}
					sb.append("\n");
				}
				sb.append("\n");
			}
		}
		String text=headings.toString()+"\n"+sb.toString();
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
		gd.addNumericField("z min",limits[4],5,10,null);
		gd.addNumericField("z max",limits[5],5,10,null);
		boolean[] logs=p3.getLogAxes();
		gd.addCheckbox("Log x?",logs[0]);
		gd.addCheckbox("Log y?",logs[1]);
		gd.addCheckbox("Log z?",logs[2]);
		gd.addStringField("x label",p3.getxLabel());
		gd.addStringField("y label",p3.getyLabel());
		gd.addStringField("z label",p3.getzLabel());
		boolean delsel=false;
		gd.addCheckbox("Delete Selected",delsel);
		boolean ascalex=false;
		gd.addCheckbox("AutoScale x",ascalex);
		boolean ascaley=false;
		gd.addCheckbox("AutoScale y",ascalex);
		boolean ascalez=false;
		gd.addCheckbox("AutoScale z",ascalez);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		limits[0]=(float)gd.getNextNumber();
		limits[1]=(float)gd.getNextNumber();
		limits[2]=(float)gd.getNextNumber();
		limits[3]=(float)gd.getNextNumber();
		limits[4]=(float)gd.getNextNumber();
		limits[5]=(float)gd.getNextNumber();
		p3.setLimits(limits);
		logs[0]=gd.getNextBoolean();
		logs[1]=gd.getNextBoolean();
		logs[2]=gd.getNextBoolean();
		p3.setLogAxes(logs[0],logs[1],logs[2]);
		p3.setxLabel(gd.getNextString());
		p3.setyLabel(gd.getNextString());
		p3.setzLabel(gd.getNextString());
		delsel=gd.getNextBoolean();
		ascalex=gd.getNextBoolean();
		ascaley=gd.getNextBoolean();
		ascalez=gd.getNextBoolean();
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
		if(ascalez){
			ascalez=false;
			p3.zautoscale();
		}
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
					saveAsObject();
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
							if(b==seldbutton){
								p3.selectSeries(p3.getSelected()-1);
								IJ.showStatus("Series "+(p3.getSelected()+1)+" of "+p3.getNSeries()+" Selected");
								updatePlot();
							}else{
								if(b==toimagebutton){
									saveToImage();
								}else{
									if(b==rotrightbutton){
										p3.setrotation(p3.getrotation()[0],p3.getrotation()[1],p3.getrotation()[2]+10.0);
										updatePlot();
									}else{
										if(b==rotleftbutton){
											p3.setrotation(p3.getrotation()[0],p3.getrotation()[1],p3.getrotation()[2]-10.0);
											updatePlot();
										}else{
											if(b==rotupbutton){
												p3.setrotation(p3.getrotation()[0]-10.0,p3.getrotation()[1],p3.getrotation()[2]);
												updatePlot();
											}else{
												if(b==rotdownbutton){
													p3.setrotation(p3.getrotation()[0]+10.0,p3.getrotation()[1],p3.getrotation()[2]);
													updatePlot();
												}else{
													if(b==rotclockbutton){
														p3.setrotation(p3.getrotation()[0],p3.getrotation()[1]-10.0,p3.getrotation()[2]);
														updatePlot();
													}else{
														if(b==rotcounterbutton){
															p3.setrotation(p3.getrotation()[0],p3.getrotation()[1]+10.0,p3.getrotation()[2]);
															updatePlot();
														}else{
															if(b==rotresetbutton){
																p3.setrotation(-60.0,0.0,-45.0);
																updatePlot();
															}
														}
													}
												}
											}
										}
									}
								}
							}
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

	public float[][][] getZValues(){
		return p3.getZValues();
	}

	public float[][] getZValues(int series){
		return p3.getZValues(series);
	}
	
	public String[] getAllLabels(){
		String[] temp={imp.getTitle(),p3.getxLabel(),p3.getyLabel(),p3.getzLabel()};
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

	public String getzLabel(){
		return p3.getzLabel();
	}

	public String[] getAnnotations(){
		return p3.getAnnotations();
	}
	
	public int[][] getNpts(){
		return p3.getNpts();
	}

	public int getNSeries(){
		return p3.getNSeries();
	}

	public float[] getLimits(){
		return p3.getLimits();
	}
	
	public boolean[] getLogAxes() {
		return p3.getLogAxes();
	}

	public int[] getShapes(){
		return p3.getShapes();
	}

	public int[] getColors(){
		return p3.getColors();
	}

	public Plot3D getPlot(){
		return p3;
	}

}
