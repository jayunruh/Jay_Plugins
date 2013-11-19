/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.gui.GenericDialog;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.datatransfer.*;
import javax.swing.*;

public class PlotPanel extends JPanel implements MouseListener,MouseMotionListener,ActionListener,ClipboardOwner{
	// this class implements a plot panel that can be added to any window
	// depends on the Plot4 class to draw the plot

	private static String defaultDirectory=null;
	private static String title;
	public Image plotimage;
	public Plot4 p3;
	private boolean zooming;
	private int startx,starty,currx,curry,panelwidth,panelheight;
	
	public static void launch_frame(PlotPanel pp,String title){
		final Frame f=new Frame(title);
		f.setLocation(10,10);
		f.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				Component[] comps=f.getComponents();
				for(int i=0;i<comps.length;i++){
					comps[i].setVisible(false);
				}
				f.dispose();
			}
		});

		f.setLayout(null);
		//Insets ins=f.getInsets();
		//pp.totalSize.height=AutoCorrFitWindow.H+ins.bottom+ins.top+65;
		//pp.totalSize.width=AutoCorrFitWindow.WR+ins.left+ins.right;
		//pp.setBounds(ins.top+5,ins.left+5,pp.totalSize.width,pp.totalSize.height);
		pp.setBounds(10,20,550,300);
		f.add(pp);
		f.pack();
		f.setResizable(false);
		//f.setSize(panel.totalSize);
		f.setSize(new Dimension(600,350));
		f.setVisible(true);
		pp.requestFocus();
	}

	public void init(String xLabel1,String yLabel1,Object xValues1,Object yValues1,Object npts1){
		zooming=false;
		if(yValues1 instanceof float[]){
			if(xValues1==null){
				p3=new Plot4(xLabel1,yLabel1,(float[])yValues1);
			}else{
				p3=new Plot4(xLabel1,yLabel1,(float[])xValues1,(float[])yValues1);
			}
		}
		if(yValues1 instanceof float[][]){
			if(xValues1==null){
				p3=new Plot4(xLabel1,yLabel1,(float[][])yValues1,npts1);
			}else{
				p3=new Plot4(xLabel1,yLabel1,(float[][])xValues1,(float[][])yValues1,npts1);
			}
		}
		panelwidth=Plot4.LEFT_MARGIN+Plot4.WIDTH+Plot4.RIGHT_MARGIN;
		panelheight=Plot4.TOP_MARGIN+Plot4.HEIGHT+Plot4.BOTTOM_MARGIN;
		addMouseMotionListener(this);
		addMouseListener(this);
		updatePlot();
	}
	
	public void init(String xLabel1,String yLabel1,float[] xVals1,float[] yVals1){
		init(xLabel1,yLabel1,xVals1,yVals1,null);
	}
	
	public void init(String xLabel1,String yLabel1,float[] yVals1){
		init(xLabel1,yLabel1,null,yVals1,null);
	}
	
	public void init(String xLabel1,String yLabel1,float[][] xVals1,float[][] yVals1,int[] npts1){
		init(xLabel1,yLabel1,xVals1,yVals1,npts1);
	}

	public void setBounds(int x,int y,int width,int height){
		panelwidth=width;
		panelheight=height;
		super.setBounds(x,y,width,height);
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1){
		p3.setLimits(xMin1,xMax1,yMin1,yMax1);
		updatePlot();
	}

	public void setLogAxes(boolean logx1,boolean logy1){
		p3.setLogAxes(logx1,logy1);
		updatePlot();
	}

	public void autoscale(){
		p3.autoscale();
		updatePlot();
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

	public void addPoints(float[] xValues1,float[] yValues1,boolean rescale){
		p3.addPoints(xValues1,yValues1,rescale);
		updatePlot();
	}

	public void addPoints(float[] yValues1,boolean rescale){
		p3.addPoints(yValues1,rescale);
		updatePlot();
	}

	public void changePlot(Plot4 plot1){
		p3=plot1;
		updatePlot();
	}

	void updatePlot(){
		plotimage=p3.getImage();
		repaint();
	}

	public Image getplotimage(){
		return plotimage;
	}

	public Plot4 getPlot(){
		return p3;
	}

	public void paint(Graphics g){
		// super.paint(g);
		g.setColor(Color.white);
		g.fillRect(0,0,panelwidth,panelheight);
		g.setColor(Color.black);
		g.clipRect(0,0,panelwidth,panelheight);
		g.drawImage(plotimage,0,0,this);
		if(zooming){
			g.setColor(Color.RED);
			int x,y,width,height;
			if(currx<startx){
				x=currx;
				width=startx-currx;
			}else{
				x=startx;
				width=currx-startx;
			}
			if(curry<starty){
				y=curry;
				height=starty-curry;
			}else{
				y=starty;
				height=curry-starty;
			}
			g.drawRect(x,y,width,height);
			g.setColor(Color.BLACK);
		}
	}

	public void update(Graphics g){
		paint(g);
	}

	public void mouseClicked(MouseEvent e){
	}

	public void mouseEntered(MouseEvent e){
	}

	public void mouseExited(MouseEvent e){
	}

	public void mousePressed(MouseEvent e){
		if(e.getButton()==MouseEvent.BUTTON3){
			handlepopup(e);
		}
		if(e.getButton()==MouseEvent.BUTTON1){
			startx=e.getX();
			starty=e.getY();
			float mag=p3.getmagnification();
			Rectangle frame=new Rectangle((int)(p3.LEFT_MARGIN*mag),(int)(p3.TOP_MARGIN*mag),(int)(p3.WIDTH*mag),(int)(p3.HEIGHT));
			if(frame.contains(startx,starty)){
				zooming=true;
			}
		}
	}

	public void mouseReleased(MouseEvent e){
		if(zooming){
			int endx=e.getX();
			int endy=e.getY();
			if(endx==startx||endy==starty){
				zooming=false;
				return;
			}
			int width=Math.abs(endx-startx);
			int height=Math.abs(endy-starty);
			int rx=startx;
			if(endx<startx){
				rx=endx;
			}
			int ry=starty;
			if(endy<starty){
				ry=endy;
			}
			Rectangle rect=new Rectangle(rx,ry,width,height);
			p3.scalerect(rect);

			updatePlot();
			zooming=false;
		}
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
		if(zooming){
			currx=e.getX();
			curry=e.getY();
			repaint();
		}
	}

	public void actionPerformed(ActionEvent e){
		String cmd=e.getActionCommand();
		if(cmd=="Autoscale"){
			autoscale();
		}
		if(cmd=="Edit"){
			editPlot();
		}
		if(cmd=="Copy"){
			copyToClipboard();
		}
		if(cmd=="List"){
			listData();
		}
		if(cmd=="Save As Text"){
			saveAsText();
		}
		if(cmd=="Save As Object"){
			saveAsObject();
		}
		// System.out.println("popup item selected = "+e.getActionCommand());
	}

	public void lostOwnership(Clipboard clipboard,Transferable contents){
	}

	void handlepopup(MouseEvent e){
		// System.out.println("popup selected");
		int x=e.getX();
		int y=e.getY();
		JPopupMenu popup=new JPopupMenu("");
		JMenuItem mi;
		mi=popup.add("Autoscale");
		mi.addActionListener(this);
		mi=popup.add("Edit");
		mi.addActionListener(this);
		mi=popup.add("Copy");
		mi.addActionListener(this);
		mi=popup.add("List");
		mi.addActionListener(this);
		mi=popup.add("Save As Text");
		mi.addActionListener(this);
		mi=popup.add("Save As Object");
		mi.addActionListener(this);
		add(popup);
		popup.show(this,x,y);
	}

	void saveAsText(){
		FileDialog fd=new FileDialog((Frame)this.getParent(),"Save as Text...",FileDialog.SAVE);
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
			System.out.println(""+e);
			return;
		}
		float[][] tempyvals=p3.getYValues();
		float[][] tempxvals=p3.getXValues();
		for(int i=0;i<p3.getmaxpts();i++){
			StringBuffer sb=new StringBuffer();
			for(int j=0;j<p3.getNSeries();j++){
				sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]);
				if(j<(p3.getNSeries()-1)){
					sb.append("\t");
				}
			}
			pw.println(sb.toString());
		}
		pw.close();
	}
	
	void saveAsObject(){
		FileDialog fd=new FileDialog((Frame)this.getParent(),"Save as Plot Object...",FileDialog.SAVE);
		if(defaultDirectory!=null)
			fd.setDirectory(defaultDirectory);
		String temptitle="plot.pw2";
		fd.setFile(temptitle);
		fd.show();
		String name=fd.getFile();
		String directory=fd.getDirectory();
		defaultDirectory=directory;
		fd.dispose();
		if(name==null||name==""||directory==null||directory=="")
			return;
		p3.saveplot2file(directory+File.separator+name);
	}

	void copyToClipboard(){
		Clipboard systemClipboard=null;
		try{
			systemClipboard=getToolkit().getSystemClipboard();
		}catch(Exception e){
			systemClipboard=null;
		}
		if(systemClipboard==null){
			System.out.println("error opening clipboard");
			return;
		}
		StringBuffer sb=new StringBuffer();
		float[][] tempyvals=p3.getYValues();
		float[][] tempxvals=p3.getXValues();
		for(int i=0;i<p3.getmaxpts();i++){
			for(int j=0;j<p3.getNSeries();j++){
				sb.append(""+tempxvals[j][i]+"\t"+tempyvals[j][i]);
				if(j<(p3.getNSeries()-1)){
					sb.append("\t");
				}
			}
			sb.append("\n");
		}
		String text=sb.toString();
		StringSelection contents=new StringSelection(text);
		systemClipboard.setContents(contents,this);
	}

	void listData(){
		float[][] tempyvals=p3.getYValues();
		float[][] tempxvals=p3.getXValues();
		int nseries=p3.getNSeries();
		int maxpts=p3.getmaxpts();
		Object[][] tabledata=new Object[maxpts][2*nseries];
		String[] columnlabels=new String[2*nseries];
		for(int i=0;i<nseries;i++){
			columnlabels[2*i]="x"+(i+1);
			columnlabels[2*i+1]="y"+(i+1);
			for(int j=0;j<maxpts;j++){
				tabledata[j][2*i]=new Float(tempxvals[i][j]);
				tabledata[j][2*i+1]=new Float(tempyvals[i][j]);
			}
		}
		TableDialog2.showDialog(null,null,"Plot Data",columnlabels,tabledata,null);
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
		updatePlot();
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
}
