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
import jalgs.jdist;
import jalgs.jstatistics;
import jalgs.jsim.rngs;

import java.awt.Button;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.Label;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.File;

/**
 * This class is an extended ImageWindow that displays 2D histograms. This class
 * was adapted from the PlotWindow class original to ImageJ Adaptations were
 * done by Jay Unruh (jru@stowers.org)
 */
public class PlotWindow2DHist extends ImageWindow implements ActionListener,ClipboardOwner,MouseMotionListener,MouseListener{

	private Button editbutton,list,copy;
	private Label coordinates;
	private static String title;
	public Plot2DHist p3;
	public boolean savehist,pearson;
	private static ColorProcessor cp;

	public PlotWindow2DHist(String title1,String xLabel1,String yLabel1,float[] xValues1,float[] yValues1,int[] lut1){
		super(createImage(title1));
		savehist=false;
		p3=new Plot2DHist(xLabel1,yLabel1,xValues1,yValues1,lut1);
		imp.getCanvas().addMouseMotionListener(this);
	}

	public PlotWindow2DHist(String title1,String xLabel1,String yLabel1,float[] hist2d,int width,int height,float startx,float endx,float starty,float endy,int[] lut1){
		super(createImage(title1));
		savehist=false;
		p3=new Plot2DHist(xLabel1,yLabel1,hist2d,width,height,startx,endx,starty,endy,lut1);
		imp.getCanvas().addMouseMotionListener(this);
	}

	public PlotWindow2DHist(String title1,Plot2DHist plot){
		super(createImage(title1));
		savehist=false;
		p3=plot;
		imp.getCanvas().addMouseMotionListener(this);
	}

	public PlotWindow2DHist(ImagePlus imp,Plot2DHist p3){
		// turns the passed ImagePlus into a plot window
		// used for bioformats which wants to make the window for me
		super(imp);
		int width=Plot2DHist.WIDTH+Plot2DHist.LEFT_MARGIN+Plot2DHist.RIGHT_MARGIN;
		int height=Plot2DHist.HEIGHT+Plot2DHist.TOP_MARGIN+Plot2DHist.BOTTOM_MARGIN;
		int[] temp=new int[width*height];
		for(int i=0;i<width*height;i++){
			temp[i]=0xffffffff;
		}
		cp=new ColorProcessor(width,height,temp);
		imp.setProcessor(cp);
		// cp=(ColorProcessor)this.imp.getProcessor();
		this.p3=p3;
	}

	public void setLut(int lut){
		p3.setLut(lut);
		updatePlot();
	}

	public void setBinSize(int binsize){
		p3.setBinSize(binsize);
		updatePlot();
	}

	public int getBinSize(){
		return p3.getBinSize();
	}

	public void setmagnification(float newmag){
		p3.setmagnification(newmag);
		updatePlot();
	}

	public float getmagnification(){
		return p3.getmagnification();
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
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1,double intMin1,double intMax1){
		p3.setLimits(xMin1,xMax1,yMin1,yMax1,intMin1,intMax1);
		updatePlot();
	}

	public void setLimits(float[] limits){
		p3.setLimits(limits);
		updatePlot();
	}
	
	public void setLogAxes(boolean logx,boolean logy){
		p3.setLogAxes(logx,logy);
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

	public void intautoscale(){
		p3.intautoscale();
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

	public int[] getroiindices(){
		Roi roi=imp.getRoi();
		if(roi!=null){
			return p3.getindices(roi);
		}
		return null;
	}

	public void updateData(float[] xValues1,float[] yValues1,boolean rescale){
		p3.updateData(xValues1,yValues1,rescale);
		updatePlot();
	}

	public void updatePlot(){
		update_density();
		cp=p3.getProcessor();
		imp.setProcessor(null,cp);
		imp.updateAndDraw();
	}

	/** Displays the plot. */
	public void draw(){
		Panel buttons=new Panel();
		// buttons.setLayout(new FlowLayout(FlowLayout.RIGHT));
		buttons.setLayout(null);
		buttons.setBounds(0,0,360,40);
		list=new Button(" Save ");
		list.setBounds(5,5,50,30);
		list.addActionListener(this);
		buttons.add(list);
		copy=new Button("Copy");
		copy.setBounds(5+50+5,5,50,30);
		copy.addActionListener(this);
		buttons.add(copy);
		editbutton=new Button("Edit");
		editbutton.setBounds(5+50+5+50+5,5,50,30);
		editbutton.addActionListener(this);
		buttons.add(editbutton);
		coordinates=new Label("Density = initialized size large");
		coordinates.setBounds(5+50+5+50+5+50+5,10,150,20);
		coordinates.setFont(new Font("Monospaced",Font.PLAIN,12));
		buttons.add(coordinates);
		add(buttons);
		// pack();
		// this.setSize(364,326);
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
			// IJ.log("test");
			if(pearson)
				coordinates.setText("Pearson = "+jutils.formatted_string(p3.get_pearson(roi)));
			else
				coordinates.setText("Density = "+jutils.formatted_string(p3.get_score(roi)));
		}else{
			if(pearson)
				coordinates.setText("Pearson = "+jutils.formatted_string(p3.get_pearson(null)));
			else
				coordinates.setText("Density = "+jutils.formatted_string(0.0f));
		}
	}

	public void showProgress(int pos,int end){
		IJ.showProgress(pos,end);
	}

	public void showLog(String text){
		IJ.log(text);
	}

	void showList(){
		StringBuffer sb=new StringBuffer();
		StringBuffer headings=new StringBuffer();
		if(!savehist){
			float[] tempxvals=p3.getXValues();
			float[] tempyvals=p3.getYValues();
			headings.append("X\tY");
			for(int i=0;i<tempxvals.length;i++){
				sb.append(""+tempxvals[i]+"\t"+tempyvals[i]+"\n");
			}
		}else{
			float[][] hist=p3.getHistogram();
			headings.append("X0");
			for(int i=1;i<hist.length;i++)
				headings.append("\tX"+i);
			for(int i=0;i<hist[0].length;i++){
				sb.append(""+hist[0][i]);
				for(int j=1;j<hist.length;j++){
					sb.append("\t"+hist[j][i]);
				}
				sb.append("\n");
			}
		}
		new TextWindow("Plot Values",headings.toString(),sb.toString(),200,400);
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
		if(!savehist){
			float[] tempxvals=p3.getXValues();
			float[] tempyvals=p3.getYValues();
			for(int i=0;i<tempxvals.length;i++){
				sb.append(""+tempxvals[i]+"\t"+tempyvals[i]+"\n");
			}
		}else{
			/*
			 * float[] limits=p3.getLimits(); boolean[] logs=p3.getLogAxes();
			 * float[] tempyvals=p3.getYValues(); float
			 * tbinsizey=((float)p3.getBinSize
			 * ()/(float)p3.histSize)*(limits[3]-limits[2]); if(logs[1]){ float
			 * logymin=(float)Math.log((double)limits[2]); float
			 * logymax=(float)Math.log((double)limits[3]);
			 * tbinsizey=(p3.getBinSize()/p3.histSize)*(logymax-logymin); }
			 */
			float[][] hist=p3.getHistogram();
			for(int i=0;i<hist.length;i++){
				sb.append(""+hist[i][0]);
				for(int j=1;j<hist[i].length;j++){
					sb.append("\t"+hist[i][j]);
				}
				sb.append("\n");
			}
		}
		String text=sb.toString();
		StringSelection contents=new StringSelection(text);
		systemClipboard.setContents(contents,this);
		IJ.showStatus(text.length()+" characters copied to Clipboard");
	}

	void saveAsObject(){
		/*FileDialog fd=new FileDialog(this,"Save as Plot Object...",FileDialog.SAVE);
		String temptitle=getTitle();
		if(!temptitle.endsWith(".pw2")){
			temptitle+=".pw2";
		}
		fd.setFile(temptitle);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
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

	void editPlot(){
		GenericDialog gd=new GenericDialog("Plot Options");
		float[] limits=p3.getLimits();
		gd.addNumericField("x min",limits[0],5,10,null);
		gd.addNumericField("x max",limits[1],5,10,null);
		gd.addNumericField("y min",limits[2],5,10,null);
		gd.addNumericField("y max",limits[3],5,10,null);
		gd.addNumericField("z min",limits[4],5,10,null);
		gd.addNumericField("z max",limits[5],5,10,null);
		gd.addStringField("x label",p3.getxLabel());
		gd.addStringField("y label",p3.getyLabel());
		gd.addCheckbox("Log x",p3.getLogAxes()[0]);
		gd.addCheckbox("Log y",p3.getLogAxes()[1]);
		boolean ascalex=false;
		gd.addCheckbox("AutoScale x",ascalex);
		boolean ascaley=false;
		gd.addCheckbox("AutoScale y",ascaley);
		boolean ascalez=false;
		gd.addCheckbox("AutoScale z",ascalez);
		boolean doscaleroi=false;
		gd.addCheckbox("Scale Roi",doscaleroi);
		gd.addNumericField("Magnification",p3.getmagnification(),5,10,null);
		gd.addNumericField("Bin Size",p3.getBinSize(),0,10,null);
		gd.addChoice("LUT",Plot2DHist.lutnames,Plot2DHist.lutnames[p3.getLut()]);
		gd.addCheckbox("Save Histogram?",savehist);
		gd.addCheckbox("Pearson?",pearson);
		gd.addCheckbox("Pearson_Analysis?",false);
		gd.addCheckbox("Quadrant_Analysis?",false);
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
		p3.setxLabel(gd.getNextString());
		p3.setyLabel(gd.getNextString());
		boolean tlogx=gd.getNextBoolean();
		boolean tlogy=gd.getNextBoolean();
		p3.setLogAxes(tlogx,tlogy);
		ascalex=gd.getNextBoolean();
		ascaley=gd.getNextBoolean();
		ascalez=gd.getNextBoolean();
		doscaleroi=gd.getNextBoolean();
		p3.setmagnification((float)gd.getNextNumber());
		p3.setBinSize((int)gd.getNextNumber());
		p3.setLut(gd.getNextChoiceIndex());
		savehist=gd.getNextBoolean();
		pearson=gd.getNextBoolean();
		if(gd.getNextBoolean()){
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addNumericField("#_of_iterations",1000,0);
			gd2.showDialog();
			if(!gd2.wasCanceled()){
				int iter=(int)gd2.getNextNumber();
				float[] vals=p3.pearson_analysis(iter,imp.getRoi(),this);
				new PlotWindowHist("Pearson_Analysis","pearson coef","frequency",vals,3).draw();
				float stdev=jstatistics.getstatistic("StDev",vals,null);
				float pearson=p3.get_pearson(imp.getRoi());
				IJ.log("Pearson = "+pearson);
				IJ.log("Randomized Pearson Stdev = "+stdev);
				double P=(new jdist()).tCumProb(100.0,Math.abs(pearson/stdev),false);
				IJ.log("Probability that distributions are random = "+(float)P);
			}
		}
		if(gd.getNextBoolean()){
			Rectangle roi=imp.getRoi().getBounds();
			float[] mins=p3.getPlotCoords(roi.x,roi.y+roi.height);
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addNumericField("X_Min",mins[0],5,15,null);
			gd2.addNumericField("Y_Min",mins[1],5,15,null);
			gd2.addNumericField("N_Trials_For_Errors",10000,0);
			gd2.addNumericField("Percentile_For_Errors",95.0,5,15,null);
			gd2.addCheckbox("Show Histograms",true);
			gd2.showDialog();
			if(!gd2.wasCanceled()){
				mins[0]=(float)gd2.getNextNumber();
				mins[1]=(float)gd2.getNextNumber();
				float[] scores=p3.quadrant_analysis(mins);
				IJ.log("Both Over = "+scores[3]);
				IJ.log("Only X Over = "+scores[1]);
				IJ.log("Only Y Over = "+scores[2]);
				IJ.log("Both under = "+scores[0]);
				float tot=scores[0]+scores[1]+scores[2]+scores[3];
				IJ.log("F Both Over = "+scores[3]/tot);
				IJ.log("F Only X Over = "+scores[1]/tot);
				IJ.log("F Only Y Over = "+scores[2]/tot);
				IJ.log("F Both under = "+scores[0]/tot);
				int ntrials=(int)gd2.getNextNumber();
				float percentile=(float)gd2.getNextNumber();
				float[] fscores={scores[0]/tot,scores[1]/tot,scores[2]/tot,scores[3]/tot};
				if(ntrials>1){
					Object[] errs=getQuadsErrs(fscores,(int)tot,ntrials,percentile);
					float[] errs2=(float[])errs[0];
					IJ.log("F Both Over Err = "+(errs2[3]-fscores[3]));
					IJ.log("F Only X Over Err = "+(errs2[1]-fscores[1]));
					IJ.log("F Only Y Over Err = "+(errs2[2]-fscores[2]));
					IJ.log("F Both under Err = "+(errs2[0]-fscores[0]));
					if(gd2.getNextBoolean()){
						float[][] simdata=(float[][])errs[1];
						String[] labels={"Both Under Histogram","X Over Histogram","Y Over Histogram","Both Over Histogram"};
						for(int i=0;i<4;i++){
							new PlotWindowHist(labels[i],"Counts","Sim Frequency",simdata[i],3).draw();
						}
					}
				}
			}
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
			p3.intautoscale();
		}
		if(doscaleroi){
			doscaleroi=false;
			scaleroi();
		}
		updatePlot();
	}
	
	private Object[] getQuadsErrs(float[] probs,int ntrials,int nsims,float percentile){
		float[][] simcounts=new float[4][nsims];
		rngs random=new rngs();
		for(int i=0;i<nsims;i++){
			float[] temp=simQuads(probs,ntrials,random);
			simcounts[0][i]=temp[0];
			simcounts[1][i]=temp[1];
			simcounts[2][i]=temp[2];
			simcounts[3][i]=temp[3];
		}
		float[] stats=new float[4];
		for(int i=0;i<4;i++){
			float temp1=percentile;
			stats[i]=jstatistics.getstatistic("Percentile",simcounts[i],new float[]{temp1});
			stats[i]/=ntrials;
		}
		return new Object[]{stats,simcounts};
	}
	
	private float[] simQuads(float[] probs,int ntrials,rngs random){
		float[] cumprobs={probs[0],probs[0]+probs[1],probs[0]+probs[1]+probs[2]};
		int[] counts=new int[4];
		for(int i=0;i<ntrials;i++){
			float rand=(float)random.unidev(1.0,0.0);
			if(rand<cumprobs[0]) counts[0]++;
			else if(rand<cumprobs[1]) counts[1]++;
			else if(rand<cumprobs[2]) counts[2]++;
			else counts[3]++;
		}
		return new float[]{counts[0],counts[1],counts[2],counts[3]};
	}

	public void lostOwnership(Clipboard clipboard,Transferable contents){
	}

	public void actionPerformed(ActionEvent e){
		Object b=e.getSource();
		if(b==editbutton){
			editPlot();
		}else{
			if(b==list){
				saveAsObject();
			}else{
				if(b==copy){
					copyToClipboard();
				}
			}
		}
	}

	public float[] getXValues(){
		return p3.getXValues();
	}

	public float[] getYValues(){
		return p3.getYValues();
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

	public int getNpts(){
		return p3.getNpts();
	}

	public float[] getLimits(){
		return p3.getLimits();
	}

	public ColorProcessor getProcessor(){
		return p3.getProcessor();
	}

	public Plot2DHist getPlot(){
		return p3;
	}

}
