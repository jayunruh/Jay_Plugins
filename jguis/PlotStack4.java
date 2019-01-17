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
import ij.process.ColorProcessor;
import ij.text.TextWindow;

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.Scrollbar;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;

/**
 * This class is an extended ImageWindow that displays line graphs. This class
 * was adapted from the PlotWindow class original to ImageJ Adaptations were
 * done by Jay Unruh (jru@stowers.org)
 */
public class PlotStack4 extends ImageWindow implements ActionListener,AdjustmentListener,ClipboardOwner{

	private Button list,save,copy,editbutton,extractbutton,selbutton,seldbutton;
	private Scrollbar sb;
	private Label slice;
	// private Label coordinates;
	private static String defaultDirectory=null;
	private static String title;
	public Plot4[] p3;
	public int currslice,nslices;
	private static ColorProcessor cp;

	public PlotStack4(String title1,Plot4 plot){
		super(createImage(title1));
		p3=new Plot4[1];
		p3[0]=plot;
		nslices=1;
		currslice=0;
	}

	public PlotStack4(String title1,Plot4[] plots){
		super(createImage(title1));
		p3=plots;
		nslices=plots.length;
		currslice=0;
	}

	public PlotStack4(String title1,String xLabel1,String yLabel1,float[][] xdata,float[][] ydata){
		super(createImage(title1));
		nslices=ydata.length;
		p3=new Plot4[nslices];
		for(int i=0;i<nslices;i++){
			p3[i]=new Plot4(xLabel1,yLabel1,xdata[i],ydata[i]);
		}
		currslice=0;
	}

	public PlotStack4(String title1,String xLabel1,String yLabel1,float[][][] xdata,float[][][] ydata){
		super(createImage(title1));
		nslices=ydata.length;
		p3=new Plot4[nslices];
		for(int i=0;i<nslices;i++){
			p3[i]=new Plot4(xLabel1,yLabel1,xdata[i],ydata[i],null);
		}
		currslice=0;
	}

	public PlotStack4(String title1,String xLabel1,String yLabel1,float[][][] ydata){
		super(createImage(title1));
		nslices=ydata.length;
		p3=new Plot4[nslices];
		for(int i=0;i<nslices;i++){
			p3[i]=new Plot4(xLabel1,yLabel1,ydata[i],null);
		}
		currslice=0;
	}

	public PlotStack4(String title1,String xLabel1,String yLabel1,float[] xdata,float[][][] ydata){
		super(createImage(title1));
		nslices=ydata.length;
		p3=new Plot4[nslices];
		float[][] tempxdata=new float[ydata[0].length][xdata.length];
		for(int i=0;i<ydata[0].length;i++)
			tempxdata[i]=xdata;
		for(int i=0;i<nslices;i++){
			p3[i]=new Plot4(xLabel1,yLabel1,tempxdata,ydata[i],null);
		}
		currslice=0;
	}

	public void setmagnification(float newmag){
		for(int i=0;i<nslices;i++){
			p3[i].setmagnification(newmag);
		}
		updatePlot();
	}

	public float getmagnification(){
		return p3[0].getmagnification();
	}

	static ImagePlus createImage(String title1){
		int width=Plot4.WIDTH+Plot4.LEFT_MARGIN+Plot4.RIGHT_MARGIN;
		int height=Plot4.HEIGHT+Plot4.TOP_MARGIN+Plot4.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		title=title1;
		return new ImagePlus(title,cp);
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1){
		p3[currslice].setLimits(xMin1,xMax1,yMin1,yMax1);
		updatePlot();
	}

	public void setLimitsAll(double xMin1,double xMax1,double yMin1,double yMax1){
		for(int i=0;i<nslices;i++){
			p3[i].setLimits(xMin1,xMax1,yMin1,yMax1);
		}
		updatePlot();
	}

	public void setLimits(float[] limits){
		p3[currslice].setLimits(limits);
		updatePlot();
	}

	/** Sets the x-axis and y-axis range. */
	public void setLogAxes(boolean logx1,boolean logy1){
		p3[currslice].setLogAxes(logx1,logy1);
		updatePlot();
	}

	public void setAllLogAxes(boolean logx1,boolean logy1){
		for(int i=0;i<nslices;i++){
			p3[i].setLogAxes(logx1,logy1);
		}
		updatePlot();
	}

	public void setSlice(int slice1){
		if(slice1>=nslices){
			slice1=nslices-1;
		}
		if(slice1<0){
			slice1=0;
		}
		currslice=slice1;
		slice.setText("Slice = "+(currslice+1));
		sb.setValue(currslice);
		updatePlot();
	}

	public void autoscale(){
		p3[currslice].autoscale();
		updatePlot();
	}

	public void xautoscale(){
		p3[currslice].xautoscale();
		updatePlot();
	}

	public void yautoscale(){
		p3[currslice].yautoscale();
		updatePlot();
	}

	public void autoscaleAll(){
		for(int i=0;i<nslices;i++){
			p3[i].autoscale();
		}
		updatePlot();
	}

	public void xautoscaleAll(){
		for(int i=0;i<nslices;i++){
			p3[i].xautoscale();
		}
		updatePlot();
	}

	public void yautoscaleAll(){
		for(int i=0;i<nslices;i++){
			p3[i].yautoscale();
		}
		updatePlot();
	}

	public void scaleroi(){
		Rectangle rect=imp.getProcessor().getRoi();
		if(rect!=null){
			p3[currslice].scalerect(rect);
			imp.killRoi();
			updatePlot();
		}
	}

	public void updatePlot(Plot4 plot,int plotnumber){
		p3[plotnumber]=plot;
		updatePlot();
	}

	public void deletePlot(int plotnumber){
		Plot4[] temp=new Plot4[nslices-1];
		for(int i=0;i<plotnumber;i++){
			temp[i]=p3[i];
		}
		for(int i=plotnumber;i<(nslices-1);i++){
			temp[i]=p3[i+1];
		}
		p3=temp;
		nslices--;
		sb.setMaximum(nslices);
		if(currslice==nslices){
			currslice--;
			sb.setValue(currslice);
		}
		updatePlot();
	}

	public void addPlot(Plot4 plot){
		Plot4[] temp=new Plot4[nslices+1];
		for(int i=0;i<nslices;i++){
			temp[i]=p3[i];
		}
		temp[nslices]=plot;
		nslices++;
		sb.setMaximum(nslices);
		updatePlot();
	}

	public void updatePlot(){
		cp=p3[currslice].getProcessor();
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
		extractbutton=new Button("Extract...");
		extractbutton.addActionListener(this);
		buttons.add(extractbutton);
		editbutton=new Button("Edit...");
		editbutton.addActionListener(this);
		buttons.add(editbutton);
		selbutton=new Button("Select+");
		selbutton.addActionListener(this);
		buttons.add(selbutton);
		seldbutton=new Button("Select-");
		seldbutton.addActionListener(this);
		buttons.add(seldbutton);
		Panel sbpanel=new Panel();
		sbpanel.setPreferredSize(new Dimension(Plot4.WIDTH+Plot4.LEFT_MARGIN+Plot4.RIGHT_MARGIN,20));
		sbpanel.setLayout(new BorderLayout());
		sb=new Scrollbar(Scrollbar.HORIZONTAL,0,1,0,nslices);
		sb.addAdjustmentListener(this);
		sbpanel.add("South",sb);
		slice=new Label("Slice = "+(currslice+1));
		buttons.add(slice);
		// coordinates = new Label("X=12345678, Y=12345678");
		// coordinates.setFont(new Font("Monospaced", Font.PLAIN, 12));
		// buttons.add(coordinates);
		add(buttons);
		add(sbpanel);
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
		StringBuffer sb=new StringBuffer();
		StringBuffer headings=new StringBuffer();
		int tempnseries=p3[currslice].getNSeries();
		float[][] tempxvals=p3[currslice].getXValues();
		float[][] tempyvals=p3[currslice].getYValues();
		for(int j=0;j<tempnseries;j++){
			headings.append("X"+(j+1)+"\t"+"Y"+(j+1));
			if(j<(tempnseries-1)){
				headings.append("\t");
			}
		}
		for(int i=0;i<p3[currslice].getmaxpts();i++){
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
		int tempnseries=p3[currslice].getNSeries();
		float[][] tempxvals=p3[currslice].getXValues();
		float[][] tempyvals=p3[currslice].getYValues();
		for(int i=0;i<p3[currslice].getmaxpts();i++){
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
		int[] tempnpts=p3[currslice].getNpts();
		float[][] tempyvals=p3[currslice].getYValues();
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
		setTitle(name);
		saveAsObject(directory+File.separator+name);
	}

	public void saveAsObject(String filename){
		p3[currslice].saveplot2file(filename);
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
		int tempnseries=p3[currslice].getNSeries();
		float[][] tempxvals=p3[currslice].getXValues();
		float[][] tempyvals=p3[currslice].getYValues();
		for(int i=0;i<p3[currslice].getmaxpts();i++){
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
		float[] limits=p3[currslice].getLimits();
		gd.addNumericField("x min",limits[0],5,10,null);
		gd.addNumericField("x max",limits[1],5,10,null);
		gd.addNumericField("y min",limits[2],5,10,null);
		gd.addNumericField("y max",limits[3],5,10,null);
		boolean[] logs=p3[currslice].getLogAxes();
		gd.addCheckbox("Log x?",logs[0]);
		gd.addCheckbox("Log y?",logs[1]);
		gd.addStringField("x label",p3[currslice].getxLabel());
		gd.addStringField("y label",p3[currslice].getyLabel());
		boolean delsel=false;
		gd.addCheckbox("Delete Selected",delsel);
		boolean ascalex=false;
		gd.addCheckbox("AutoScale x",ascalex);
		boolean ascaley=false;
		gd.addCheckbox("AutoScale y",ascalex);
		boolean doscaleroi=false;
		gd.addCheckbox("Scale Roi",doscaleroi);
		gd.addNumericField("Magnification",p3[currslice].getmagnification(),5,10,null);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		limits[0]=(float)gd.getNextNumber();
		limits[1]=(float)gd.getNextNumber();
		limits[2]=(float)gd.getNextNumber();
		limits[3]=(float)gd.getNextNumber();
		p3[currslice].setLimits(limits);
		logs[0]=gd.getNextBoolean();
		logs[1]=gd.getNextBoolean();
		p3[currslice].setLogAxes(logs[0],logs[1]);
		p3[currslice].setxLabel(gd.getNextString());
		p3[currslice].setyLabel(gd.getNextString());
		delsel=gd.getNextBoolean();
		ascalex=gd.getNextBoolean();
		ascaley=gd.getNextBoolean();
		doscaleroi=gd.getNextBoolean();
		setmagnification((float)gd.getNextNumber());
		if(delsel){
			delsel=false;
			p3[currslice].deleteSeries(p3[currslice].getSelected(),false);
		}
		if(ascalex){
			ascalex=false;
			p3[currslice].xautoscale();
		}
		if(ascaley){
			ascaley=false;
			p3[currslice].yautoscale();
		}
		if(doscaleroi){
			doscaleroi=false;
			scaleroi();
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
						if(b==extractbutton){
							GenericDialog gd2=new GenericDialog("Options");
							String[] choices={"Extract Single","Extract to Movie","Combine All","Extract All"};
							gd2.addChoice("Extraction_Options",choices,choices[0]);
							gd2.showDialog();
							if(gd2.wasCanceled()){
								return;
							}
							int index=gd2.getNextChoiceIndex();
							if(index==1){
								ImageStack is=new ImageStack(imp.getWidth(),imp.getHeight());
								for(int i=0;i<nslices;i++){
									is.addSlice("",p3[i].get8Processor());
								}
								new ImagePlus(this.getTitle()+" movie",is).show();
							}else{
								if(index==2){
									int totseries=0;
									int maxpts=0;
									for(int i=0;i<nslices;i++){
										totseries+=p3[i].getNSeries();
										int temp=p3[i].getmaxpts();
										if(temp>maxpts)
											maxpts=temp;
									}
									float[][] newxvals=new float[totseries][maxpts];
									float[][] newyvals=new float[totseries][maxpts];
									int[] npts=new int[totseries];
									int counter=0;
									for(int i=0;i<nslices;i++){
										int tempnser=p3[i].getNSeries();
										float[][] tempxvals=p3[i].getXValues();
										float[][] tempyvals=p3[i].getYValues();
										int[] tempnpts=p3[i].getNpts();
										for(int j=0;j<tempnser;j++){
											npts[counter]=tempnpts[j];
											System.arraycopy(tempxvals[j],0,newxvals[counter],0,npts[counter]);
											System.arraycopy(tempyvals[j],0,newyvals[counter],0,npts[counter]);
											counter++;
										}
									}
									new PlotWindow4(this.getPlotTitle()+"_combined",p3[0].getxLabel(),p3[0].getyLabel(),newxvals,newyvals,npts).draw();
								}else{
									if(index==3){
										for(int i=0;i<p3.length;i++) new PlotWindow4(this.getPlotTitle()+(i+1),p3[i]).draw();
									} else {
										new PlotWindow4(this.getPlotTitle(),p3[currslice]).draw();
									}

								}
							}
						}else{
							if(b==selbutton){
								p3[currslice].selectSeries(p3[currslice].getSelected()+1);
								updatePlot();
							}else{
								p3[currslice].selectSeries(p3[currslice].getSelected()-1);
								updatePlot();
							}
						}
					}
				}
			}
		}
	}

	public void adjustmentValueChanged(AdjustmentEvent e){
		if(e.getSource()==sb){
			setSlice(sb.getValue());
		}
	}

	public void selectSeries(int series){
		p3[currslice].selectSeries(series);
		updatePlot();
	}

	public float[][] getXValues(){
		return p3[currslice].getXValues();
	}

	public float[] getXValues(int series){
		return p3[currslice].getXValues(series);
	}

	public float[][] getYValues(){
		return p3[currslice].getYValues();
	}

	public float[] getYValues(int series){
		return p3[currslice].getYValues(series);
	}

	public String getPlotTitle(){
		return imp.getTitle();
	}

	public String getxLabel(){
		return p3[currslice].getxLabel();
	}

	public String getyLabel(){
		return p3[currslice].getyLabel();
	}

	public int[] getNpts(){
		return p3[currslice].getNpts();
	}

	public int getNSeries(){
		return p3[currslice].getNSeries();
	}

	public float[] getLimits(){
		return p3[currslice].getLimits();
	}

	public int[] getShapes(){
		return p3[currslice].getShapes();
	}

	public int[] getColors(){
		return p3[currslice].getColors();
	}

	public ColorProcessor getProcessor(){
		return p3[currslice].getProcessor();
	}

}
