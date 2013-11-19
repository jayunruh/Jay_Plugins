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

public class PlotPanel3D extends JPanel implements MouseListener,MouseMotionListener,ActionListener,ClipboardOwner{
	// this class implements a plot panel that can be added to any window
	// depends on the Plot3D class to draw the plot

	private static String defaultDirectory=null;
	public Image plotimage;
	public Plot3D p3;
	public int startxy,panelwidth,panelheight;

	public void init(String xLabel1,String yLabel1,String zLabel1,Object xValues1,Object yValues1,Object zValues1,Object npts1){
		startxy=0;
		if(zValues1 instanceof float[][]){
			if(xValues1==null||yValues1==null){
				p3=new Plot3D(xLabel1,yLabel1,zLabel1,(float[][])zValues1,startxy);
			}else{
				p3=new Plot3D(xLabel1,yLabel1,zLabel1,(float[])xValues1,(float[])yValues1,(float[][])zValues1);
			}
		}
		if(zValues1 instanceof float[][][]){
			if(xValues1==null||yValues1==null){
				p3=new Plot3D(xLabel1,yLabel1,zLabel1,(float[][][])zValues1,startxy,npts1);
			}else{
				p3=new Plot3D(xLabel1,yLabel1,zLabel1,(float[][])xValues1,(float[][])yValues1,(float[][][])zValues1,npts1);
			}
		}
		panelwidth=Plot3D.LEFT_MARGIN+Plot3D.WIDTH+Plot3D.RIGHT_MARGIN;
		panelheight=Plot3D.TOP_MARGIN+Plot3D.HEIGHT+Plot3D.BOTTOM_MARGIN;
		addMouseMotionListener(this);
		addMouseListener(this);
		updatePlot();
	}

	public void setBounds(int x,int y,int width,int height){
		panelwidth=width;
		panelheight=height;
		super.setBounds(x,y,width,height);
	}

	/** Sets the x-axis and y-axis range. */
	public void setLimits(double xMin1,double xMax1,double yMin1,double yMax1,double zMin1,double zMax1){
		p3.setLimits(xMin1,xMax1,yMin1,yMax1,zMin1,zMax1);
		updatePlot();
	}

	public void setLogAxes(boolean logx1,boolean logy1,boolean logz1){
		p3.setLogAxes(logx1,logy1,logz1);
		updatePlot();
	}

	public void autoscale(){
		p3.autoscale();
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

	public void addPoints(float[][] zValues1,boolean rescale){
		p3.addPoints(zValues1,rescale,startxy);
		updatePlot();
	}

	public void changePlot(Plot3D plot1){
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

	public Plot3D getPlot(){
		return p3;
	}

	public void paint(Graphics g){
		// super.paint(g);
		g.setColor(Color.white);
		g.fillRect(0,0,panelwidth,panelheight);
		g.setColor(Color.black);
		g.clipRect(0,0,panelwidth,panelheight);
		g.drawImage(plotimage,0,0,this);
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
	}

	public void mouseReleased(MouseEvent e){
	}

	public void mouseMoved(MouseEvent e){
	}

	public void mouseDragged(MouseEvent e){
	}

	public void actionPerformed(ActionEvent e){
		String cmd=e.getActionCommand();
		if(cmd=="Autoscale"){
			autoscale();
			return;
		}
		if(cmd=="Edit"){
			editPlot();
			return;
		}
		if(cmd=="Copy"){
			copyToClipboard();
			return;
		}
		if(cmd=="List"){
			listData();
			return;
		}
		if(cmd=="Save"){
			saveAsText();
			return;
		}

		double[] rotation=new double[3];
		if(cmd=="Rot. Right"){
			rotation[2]=10.0;
		}
		if(cmd=="Rot. Left"){
			rotation[2]=-10.0;
		}
		if(cmd=="Rot. Up"){
			rotation[0]=-10.0;
		}
		if(cmd=="Rot. Down"){
			rotation[0]=10.0;
		}
		if(cmd=="Rot. Clock"){
			rotation[1]=-10.0;
		}
		if(cmd=="Rot. Counter"){
			rotation[1]=10.0;
		}
		double[] currotation=p3.getrotation();
		currotation[0]+=rotation[0];
		currotation[1]+=rotation[1];
		currotation[2]+=rotation[2];
		p3.setrotation(currotation);
		updatePlot();
		// System.out.println("popup item selected = "+e.getActionCommand());
	}

	public void lostOwnership(Clipboard clipboard,Transferable contents){
	}

	void handlepopup(MouseEvent e){
		// System.out.println("popup selected");
		int x=e.getX();
		int y=e.getY();
		/*
		 * PopupMenu popup=new PopupMenu(""); MenuItem mi; mi=new
		 * MenuItem("Autoscale"); mi.addActionListener(this); popup.add(mi);
		 * mi=new MenuItem("Edit"); mi.addActionListener(this); popup.add(mi);
		 * mi=new MenuItem("Copy"); mi.addActionListener(this); popup.add(mi);
		 * mi=new MenuItem("Save"); mi.addActionListener(this); popup.add(mi);
		 */
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
		mi=popup.add("Save");
		mi.addActionListener(this);

		mi=popup.add("Rot. Right");
		mi.addActionListener(this);
		mi=popup.add("Rot. Left");
		mi.addActionListener(this);
		mi=popup.add("Rot. Up");
		mi.addActionListener(this);
		mi=popup.add("Rot. Down");
		mi.addActionListener(this);
		mi=popup.add("Rot. Clock");
		mi.addActionListener(this);
		mi=popup.add("Rot. Counter");
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
		int nseries=p3.getNSeries();
		int maxxpts=p3.getmaxxpts();
		int maxypts=p3.getmaxypts();
		float[][][] zValues=p3.getZValues();
		for(int i=0;i<nseries;i++){
			for(int j=0;j<maxxpts;j++){
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
		int nseries=p3.getNSeries();
		int maxxpts=p3.getmaxxpts();
		int maxypts=p3.getmaxypts();
		float[][][] zValues=p3.getZValues();
		for(int i=0;i<nseries;i++){
			for(int j=0;j<maxxpts;j++){
				for(int k=0;k<maxypts;k++){
					sb.append(""+zValues[i][j][k]);
					if(k<(maxypts-1)){
						sb.append("\t");
					}
				}
				sb.append("\n");
			}
			sb.append("\n");
		}
		String text=sb.toString();
		StringSelection contents=new StringSelection(text);
		systemClipboard.setContents(contents,this);
	}

	void listData(){
		float[][][] zValues=p3.getZValues();
		int maxxpts=p3.getmaxxpts();
		int maxypts=p3.getmaxypts();
		int nseries=p3.getNSeries();
		Object[][] tabledata=new Object[maxypts*nseries][maxxpts];
		String[] columnlabels=new String[maxxpts];
		for(int i=0;i<maxxpts;i++){
			columnlabels[i]="x"+i;
		}
		for(int i=0;i<nseries;i++){
			for(int j=0;j<maxxpts;j++){
				for(int k=0;k<maxypts;k++){
					tabledata[j][k+i*maxypts]=new Float(zValues[i][j][k]);
				}
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
		gd.addNumericField("z min",limits[4],5,10,null);
		gd.addNumericField("z max",limits[5],5,10,null);
		boolean[] logs=p3.getLogAxes();
		gd.addCheckbox("Log x?",logs[0]);
		gd.addCheckbox("Log y?",logs[1]);
		gd.addCheckbox("Log z?",logs[2]);
		gd.addStringField("x label",p3.getxLabel());
		gd.addStringField("y label",p3.getyLabel());
		gd.addStringField("z label",p3.getzLabel());
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
		ascalex=gd.getNextBoolean();
		ascaley=gd.getNextBoolean();
		ascalez=gd.getNextBoolean();
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

	public String getxLabel(){
		return p3.getxLabel();
	}

	public String getyLabel(){
		return p3.getyLabel();
	}

	public String getzLabel(){
		return p3.getzLabel();
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
}
