/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.text.TextWindow;
import jalgs.jdataio;
import jalgs.jdist;
import jalgs.jstatistics;
import jalgs.jfit.NLLSfit;
import jalgs.jfit.NLLSfit_v2;
import jalgs.jfit.NLLSfitinterface;
import jalgs.jfit.NLLSfitinterface_v2;
import jalgs.jfit.NLLSglobalfit;
import jalgs.jfit.NLLSglobalfit_v2;
import jalgs.jfit.monte_carlo_errors_v2;
import jalgs.jfit.support_plane_errors;
import jalgs.jfit.support_plane_errors_v2;

import java.awt.Button;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;

public class CrossCorrFitWindow extends Panel implements ActionListener,NLLSfitinterface_v2,ListSelectionListener,TableModelListener{

	public final static int H=750;
	public final static int WR=800;
	public Dimension totalSize=new Dimension();

	private Button fitavgbutton,fitglobalbutton,clearparamsbutton,undobutton,geterrorsbutton,savebutton,montecarlobutton;
	private Button editconsbutton;
	private PlotWindow4 pwavg,pwfit,pwavgtraj,pwfittraj;
	private int ncurves,nparams,npts,dispcurve,psfflag;
	private int[] nmeas,indices;
	private float[][][] corr,fit,undofit,xvals,trajectories;
	private float[][] avg,avgfit,avgs,vars,avgtraj;
	private Label globalc2label,dispcurvelabel,copylabel,betalabel;
	private boolean[] include;
	private TextField betaval;
	private TablePanel datapanel;
	private String[] names;
	private double[] intensity1,g01,c2,undoc2,intensity2,g02,g0cc,g0mincc;
	private double[][] globalparams,avgconstraints,undoparams;
	private double[][][] globalconstraints;
	private double[] avgparams;
	private double globalc2,beta,undoglobalc2,khz;
	private boolean checkc2,brightcorr,outerrors;
	private int[][] globalvflmatrix,undovflmatrix;
	private int[] avgfixes;
	private int[] fitplotcolors;
	private int[] avgplotcolors;
	private String[][] globalformulas,undoformulas;
	private String[] paramsnames;
	public double w0g,w0r;
	private TextWindow tw;
	NLLSglobalfit_v2 globalfitclass;
	NLLSfit_v2 fitclass;

	public static void launch_frame(CrossCorrFitWindow panel){
		final Frame f=new Frame("Cross Corr Analysis");
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
		Insets ins=f.getInsets();
		panel.totalSize.height=CrossCorrFitWindow.H+ins.bottom+ins.top+65;
		panel.totalSize.width=CrossCorrFitWindow.WR+ins.left+ins.right;
		panel.setBounds(ins.top+5,ins.left+5,panel.totalSize.width,panel.totalSize.height);
		f.add(panel);
		f.pack();
		f.setResizable(false);
		f.setSize(panel.totalSize);
		f.setVisible(true);
		panel.requestFocus();
	}
	
	public static void launch_from_file(String filename,boolean brightcorr1) {
		try{
			InputStream is=new BufferedInputStream(new FileInputStream(filename));
			final CrossCorrFitWindow cw=new CrossCorrFitWindow();
			cw.init_from_is(is,brightcorr1);
			is.close();
			CrossCorrFitWindow.launch_frame(cw);
		}catch(IOException e){
			return;
		}
	}
	
	public static boolean is_this(String filename) {
		try{
			InputStream is=new BufferedInputStream(new FileInputStream(filename));
			jdataio jdio=new jdataio();
			String temp=jdio.readstring(is);
			int type=jdio.readintelint(is);
			is.close();
			if(temp.equals("jfcs_file_type") && type==1) return true;
			else return false;
		}catch(IOException e){
			return false;
		}
	}
	
	public void init_from_is(InputStream is,boolean brightcorr1) {
		jdataio jdio=new jdataio();
		String temp=jdio.readstring(is);
		int type=jdio.readintelint(is);
		if(type!=1) {
			IJ.error("wrong file type");
		}
		int nseries=jdio.readintelint(is);
		int npts=jdio.readintelint(is);
		int trajpts=jdio.readintelint(is);
		float khz=jdio.readintelfloat(is);
		int psfflag=jdio.readintelint(is);
		boolean brightcorr=(jdio.readintelint(is)==1);
		boolean changebrightcorr=false;
		if(brightcorr1!=brightcorr) changebrightcorr=true;
		String[] names=new String[nseries];
		float[] xvals=new float[npts];
		float[][][] corr=new float[nseries][3][npts];
		float[][][] trajectories=new float[nseries][2][trajpts];
		float[][] avg=new float[3][nseries];
		float[][] var=new float[3][nseries];
		int[] selections=new int[nseries];
		for(int i=0;i<nseries;i++) {names[i]=jdio.readstring(is);}
		jdio.readintelfloatfile(is,npts,xvals);
		for(int i=0;i<nseries;i++) {
			for(int j=0;j<3;j++) {
				jdio.readintelfloatfile(is,npts,corr[i][j]);
			}
		}
		for(int i=0;i<nseries;i++) {
			for(int j=0;j<2;j++) {
				jdio.readintelfloatfile(is,trajpts,trajectories[i][j]);
			}
		}
		jdio.readintelfloatfile(is,nseries,avg[0]);
		jdio.readintelfloatfile(is,nseries,avg[1]);
		jdio.readintelfloatfile(is,nseries,avg[2]);
		jdio.readintelfloatfile(is,nseries,var[0]);
		jdio.readintelfloatfile(is,nseries,var[1]);
		jdio.readintelfloatfile(is,nseries,var[2]);
		jdio.readintelintfile(is,nseries,selections);
		//if we are changing the brightcorr status, do that here
		if(changebrightcorr) {
			if(brightcorr) { //here we undo brightcorr by dividing by avg
				for(int i=0;i<nseries;i++) {
					for(int j=0;j<npts;j++) {
						corr[i][0][j]/=khz*avg[0][i];
						corr[i][1][j]/=khz*avg[1][i];
						corr[i][2][j]/=khz*avg[2][i]; //is this right?
					}
				}
			} else { //here we apply brightcorr by multiplying by avg
				for(int i=0;i<nseries;i++) {
					for(int j=0;j<npts;j++) {
						corr[i][0][j]*=khz*avg[0][i];
						corr[i][1][j]*=khz*avg[1][i];
						corr[i][2][j]*=khz*avg[2][i];
					}
				}
			}
			brightcorr=!brightcorr;
		}
		init(names,corr,xvals,trajectories,avg,var,khz,psfflag,brightcorr);
		//now set up the selections
		for(int i=0;i<nseries;i++) {
			include[i]=(selections[i]==1);
			datapanel.table.setValueAt(new Boolean(include[i]),i,0);
		}
		updateavg();
		datapanel.table.setValueAt(""+(float)intensity1[ncurves],ncurves,2);
		datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
		pwavg.updateSeries(this.avg[0],0,true);
		pwavg.updateSeries(this.avg[1],1,true);
		pwavg.updateSeries(this.avg[2],2,true);
	}
	
	public void savefile() {
		SaveDialog sd=new SaveDialog("Save as JFCS Object...","My_Autocorr_File",".jfcs");
		String name=sd.getFileName();
		String directory=sd.getDirectory();
		if(name==null||name==""||directory==null||directory=="")
			return;
		if(!name.endsWith(".jfcs")){
			name+=".jfcs";
		}
		String dir2=directory.replace("\\","\\\\");
		try{
			OutputStream os=new BufferedOutputStream(new FileOutputStream(directory+File.separator+name));
			save2os(os);
			os.close();
		}catch(IOException e){
			return;
		}
	}
	
	public void save2os(OutputStream os) {
		jdataio jdio=new jdataio();
		jdio.writestring(os,"jfcs_file_type");
		jdio.writeintelint(os,1);
		jdio.writeintelint(os,ncurves);
		jdio.writeintelint(os,npts);
		float[][][] temptraj=trajectories;
		if(trajectories==null) {
			temptraj=new float[ncurves][2][10];
			for(int i=0;i<ncurves;i++) {
				for(int j=0;j<10;j++) {
					temptraj[i][0][j]=avgs[0][i];
					temptraj[i][1][j]=avgs[1][i];
				}
			}
		}
		jdio.writeintelint(os,temptraj[0][0].length);
		jdio.writeintelfloat(os,(float)khz);
		jdio.writeintelint(os,psfflag);
		jdio.writeintelint(os,brightcorr?1:0);
		for(int i=0;i<ncurves;i++) {jdio.writestring(os,names[i]);}
		jdio.writeintelfloatarray(os,xvals[0][0]);
		for(int i=0;i<ncurves;i++) {
			for(int j=0;j<3;j++) {
				jdio.writeintelfloatarray(os,corr[i][j]);
			}
		}
		for(int i=0;i<ncurves;i++) {
			jdio.writeintelfloatarray(os,temptraj[i][0]);
			jdio.writeintelfloatarray(os,temptraj[i][1]);
		}
		jdio.writeintelfloatarray(os,avgs[0]);
		jdio.writeintelfloatarray(os,avgs[1]);
		jdio.writeintelfloatarray(os,avgs[2]);
		jdio.writeintelfloatarray(os,vars[0]);
		jdio.writeintelfloatarray(os,vars[1]);
		jdio.writeintelfloatarray(os,vars[2]);
		for(int i=0;i<ncurves;i++) {
			if(include[i]) jdio.writeintelint(os,1);
			else jdio.writeintelint(os,0);
		}
	}

	public void init(String[] names1,float[][][] corr1,float[] xvals1,float[][][] trajectories1,float[][] avg1,float[][] var1,float khz1,int psfflag1,boolean brightcorr1){
		setLayout(null);
		names=names1;
		corr=corr1;
		psfflag=psfflag1;
		avgs=avg1;
		vars=var1;
		brightcorr=brightcorr1;
		khz=khz1;
		trajectories=trajectories1;
		ncurves=corr.length;
		nparams=11;
		npts=corr[0][0].length;

		include=new boolean[ncurves];
		intensity1=new double[ncurves+1];
		g01=new double[ncurves+1];
		intensity2=new double[ncurves+1];
		g02=new double[ncurves+1];
		g0cc=new double[ncurves+1];
		g0mincc=new double[ncurves+1];
		c2=new double[ncurves+1];
		nmeas=new int[ncurves+1];
		avg=new float[3][npts];
		indices=new int[ncurves];
		beta=0.0;

		getintbright();
		for(int i=0;i<ncurves;i++){
			include[i]=true;
			indices[i]=i;
		}
		updateavg();

		datapanel=new TablePanel();
		Object[][] tabledata=new Object[ncurves+1][9];
		String[] labels={"Include","Name","Ig (kHz)","G(0)g","Ir (kHz)","G(0)r","G(0)cc","min","chi^2"};
		if(brightcorr){
			labels[3]="Bg (kHz)";
			labels[5]="Br (kHz)";
			labels[6]="Bcc (kHz)";
		}
		for(int i=0;i<=ncurves;i++){
			if(i!=ncurves){
				tabledata[i][0]=new Boolean(true);
				tabledata[i][1]=""+names[i];
			} else {
				tabledata[i][1]="Avg";
			}
			tabledata[i][2]=""+intensity1[i];
			tabledata[i][3]=""+g01[i];
			tabledata[i][4]=""+intensity2[i];
			tabledata[i][5]=""+g02[i];
			tabledata[i][6]=""+g0cc[i];
			tabledata[i][7]=""+g0mincc[i];
			tabledata[i][8]=""+c2[i];
		}
		datapanel.init(labels,tabledata,null);
		datapanel.setBounds(10,30,600,750);
		// datapanel.setBounds(10,30,300,300);
		datapanel.table.getSelectionModel().addListSelectionListener(this);
		datapanel.table.getModel().addTableModelListener(this);

		int starty=60;
		int startx=10;
		int buttonsx=startx+30+210+50+50+50+50+50+50+90;

		fitavgbutton=new Button("Fit Avg");
		fitavgbutton.setBounds(buttonsx,starty-25,100,40);
		fitavgbutton.addActionListener(this);
		add(fitavgbutton);

		fitglobalbutton=new Button("Fit Global");
		fitglobalbutton.setBounds(buttonsx,starty-25+50,100,40);
		fitglobalbutton.addActionListener(this);
		add(fitglobalbutton);

		clearparamsbutton=new Button("Reset Fit Params");
		clearparamsbutton.setBounds(buttonsx,starty-25+50+50,100,40);
		clearparamsbutton.addActionListener(this);
		add(clearparamsbutton);

		checkc2=false;

		fitclass=new NLLSfit_v2(this,0.0001,50,0.1);
		globalfitclass=new NLLSglobalfit_v2(this,0.0001,50,0.1);
		avgfit=new float[3][npts];
		fit=new float[ncurves][3][npts];

		xvals=new float[ncurves][3][npts];
		for(int i=0;i<ncurves;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<npts;k++){
					xvals[i][j][k]=xvals1[k];
					fit[i][j][k]=0.0f;
				}
			}
		}

		globalc2label=new Label("Global chi^2 = "+(float)0.0);
		globalc2label.setBounds(buttonsx,starty-25+50+50+50,160,20);
		add(globalc2label);

		betalabel=new Label("Bleedthrough f");
		betalabel.setBounds(buttonsx,starty-25+50+50+50+30+30,100,20);
		add(betalabel);
		betaval=new TextField(""+(float)beta);
		betaval.setBounds(buttonsx+100,starty-25+50+50+50+30+30,40,20);
		betaval.addActionListener(this);
		add(betaval);

		beta=Double.parseDouble(betaval.getText());
		updatebeta();

		undobutton=new Button("Undo Global Fit");
		undobutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50,100,40);
		undobutton.addActionListener(this);
		add(undobutton);
		
		editconsbutton=new Button("Edit Constraints");
		editconsbutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50+50+50,100,40);
		editconsbutton.addActionListener(this);
		add(editconsbutton);

		geterrorsbutton=new Button("Get Errors");
		geterrorsbutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50+50+50+50,100,40);
		geterrorsbutton.addActionListener(this);
		add(geterrorsbutton);
		
		montecarlobutton=new Button("Monte Carlo");
		montecarlobutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50+50+50+50+50,100,40);
		montecarlobutton.addActionListener(this);
		add(montecarlobutton);
		
		savebutton=new Button("Save Analysis");
		savebutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50+50+50+50+50+50,100,40);
		savebutton.addActionListener(this);
		add(savebutton);

		copylabel=new Label("copyright 2009 Jay Unruh (jru@stowers.org) non-profit use only");
		copylabel.setBounds(10,790,400,20);
		add(copylabel);

		float[][] dumavgcorr=new float[3][avg[0].length];
		System.arraycopy(avg[0],0,dumavgcorr[0],0,avg[0].length);
		System.arraycopy(avg[1],0,dumavgcorr[1],0,avg[0].length);
		System.arraycopy(avg[2],0,dumavgcorr[2],0,avg[0].length);
		pwavg=new PlotWindow4("Avg","tau","G(tau)",xvals[0],dumavgcorr,null);
		pwavg.setLogAxes(true,false);
		pwavg.draw();
		for(int i=0;i<3;i++){
			pwavg.addPoints(xvals[0][i],new float[npts],true);
		}
		pwavg.addPoints(xvals[0][0],new float[npts],true);
		avgplotcolors=pwavg.getColors();
		avgplotcolors[0]=2;
		avgplotcolors[1]=3;
		avgplotcolors[2]=0;
		avgplotcolors[3]=1;
		avgplotcolors[4]=1;
		avgplotcolors[5]=1;
		avgplotcolors[6]=5;
		pwavg.setLogAxes(true,false);

		float[][] dumcorr=new float[3][corr[0][0].length];
		System.arraycopy(corr[dispcurve][0],0,dumcorr[0],0,corr[0][0].length);
		System.arraycopy(corr[dispcurve][1],0,dumcorr[1],0,corr[0][0].length);
		System.arraycopy(corr[dispcurve][2],0,dumcorr[2],0,corr[0][0].length);
		pwfit=new PlotWindow4("Selected Curve","tau","G(tau)",xvals[0],dumcorr,null);
		pwfit.setLogAxes(true,false);
		pwfit.draw();
		for(int i=0;i<3;i++){
			pwfit.addPoints(xvals[0][i],fit[dispcurve][i],true);
		}
		pwfit.addPoints(xvals[0][0],new float[npts],true);
		fitplotcolors=pwfit.getColors();
		fitplotcolors[0]=2;
		fitplotcolors[1]=3;
		fitplotcolors[2]=0;
		fitplotcolors[3]=1;
		fitplotcolors[4]=1;
		fitplotcolors[5]=1;
		fitplotcolors[6]=5;
		pwfit.setLogAxes(true,false);

		if(trajectories!=null){
			pwavgtraj=new PlotWindow4("Avg Traj","sample","Intensity",avgtraj,null);
			pwavgtraj.draw();
			int[] tempcolors=pwavgtraj.getColors();
			tempcolors[0]=2;
			tempcolors[1]=3;
			pwavgtraj.setLogAxes(false,false);

			float[][] dumtraj=new float[2][trajectories[0][0].length];
			System.arraycopy(trajectories[dispcurve][0],0,dumtraj[0],0,trajectories[0][0].length);
			System.arraycopy(trajectories[dispcurve][1],0,dumtraj[1],0,trajectories[0][0].length);
			pwfittraj=new PlotWindow4("Selected Traj","sample","Intensity",dumtraj,null);
			pwfittraj.draw();
			tempcolors=pwfittraj.getColors();
			tempcolors[0]=2;
			tempcolors[1]=3;
			pwfittraj.setLogAxes(false,false);
		}
		add(datapanel);
		resetparams();
	}

	public void setVisible(boolean b){
		if(!b){
			GenericDialog gd=new GenericDialog("Options");
			boolean closeplots=true;
			gd.addCheckbox("Close Plots?",closeplots);
			gd.showDialog();
			closeplots=gd.getNextBoolean();
			if(closeplots){
				if(pwavgtraj!=null&&pwavgtraj.isVisible()){
					pwavgtraj.close();
				}
				if(pwavg!=null&&pwavg.isVisible()){
					pwavg.close();
				}
				if(pwfit!=null&&pwfit.isVisible()){
					pwfit.close();
				}
				if(pwfittraj!=null&&pwfittraj.isVisible()){
					pwfittraj.close();
				}
			}
		}
		super.setVisible(b);
	}

	public void actionPerformed(ActionEvent e){
		beta=Double.parseDouble(betaval.getText());
		updatebeta();
		if(e.getSource()==fitavgbutton){
			fitavg();
		}
		if(e.getSource()==fitglobalbutton){
			fitglobal();
		}
		if(e.getSource()==clearparamsbutton){
			resetparams();
		}
		if(e.getSource()==undobutton){
			undoglobalfit();
		}
		if(e.getSource()==geterrorsbutton){
			geterrors();
		}
		if(e.getSource()==editconsbutton){
			showconstraintsdialog();
		}
		if(e.getSource()==savebutton){
			savefile();
		}
		if(e.getSource()==montecarlobutton){
			montecarlo();
		}
	}

	public void valueChanged(ListSelectionEvent e){
		if(e.getSource()==datapanel.table.getSelectionModel()){
			int selrow=datapanel.table.getSelectedRow();
			dispcurve=selrow;
			if(dispcurve==ncurves){
				dispcurve--;
			}
			for(int i=0;i<3;i++){
				pwfit.updateSeries(corr[dispcurve][i],i,true);
				pwfit.updateSeries(fit[dispcurve][i],3+i,true);
			}
			if(trajectories!=null){
				pwfittraj.updateSeries(trajectories[dispcurve][0],0,true);
				pwfittraj.updateSeries(trajectories[dispcurve][1],1,true);
			}
		}
	}

	public void tableChanged(TableModelEvent e){
		int row=e.getFirstRow();
		int column=e.getColumn();
		if(column==0){
			boolean valid=((Boolean)datapanel.table.getValueAt(row,0)).booleanValue();
			if(valid!=include[row]){
				beta=Double.parseDouble(betaval.getText());
				include[row]=valid;
				updateavg();
				updatebeta();
				datapanel.table.setValueAt(""+(float)intensity1[ncurves],ncurves,2);
				datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
				datapanel.table.setValueAt(""+(float)intensity2[ncurves],ncurves,4);
				datapanel.table.setValueAt(""+(float)g02[ncurves],ncurves,5);
				datapanel.table.setValueAt(""+(float)g0cc[ncurves],ncurves,6);
				datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
				for(int i=0;i<3;i++){
					pwavg.updateSeries(avg[i],i,true);
				}
				if(trajectories!=null){
					pwavgtraj.updateSeries(avgtraj[0],0,true);
					pwavgtraj.updateSeries(avgtraj[1],1,true);
				}
			}
		}
	}

	private void resetparams(){
		avgparams=new double[22];
		avgparams[0]=5.0;
		avgparams[1]=0.0;
		avgparams[2]=g01[ncurves];
		avgparams[3]=10.0;
		avgparams[4]=0.0;
		avgparams[5]=100.0;
		avgparams[6]=0.0;
		avgparams[7]=0.01;
		avgparams[8]=0.0;
		avgparams[9]=g02[ncurves];
		avgparams[10]=10.0;
		avgparams[11]=0.0;
		avgparams[12]=100.0;
		avgparams[13]=0.0;
		avgparams[14]=0.01;
		avgparams[15]=0.0;
		avgparams[16]=g0cc[ncurves];
		avgparams[17]=10.0;
		avgparams[18]=0.0;
		avgparams[19]=100.0;
		avgparams[20]=0.0;
		avgparams[21]=0.01;
		avgfixes=new int[22];
		avgfixes[0]=1;
		avgfixes[1]=0;
		avgfixes[2]=0;
		avgfixes[3]=0;
		avgfixes[4]=1;
		avgfixes[5]=1;
		avgfixes[6]=1;
		avgfixes[7]=1;
		avgfixes[8]=0;
		avgfixes[9]=0;
		avgfixes[10]=0;
		avgfixes[11]=1;
		avgfixes[12]=1;
		avgfixes[13]=1;
		avgfixes[14]=1;
		avgfixes[15]=0;
		avgfixes[16]=0;
		avgfixes[17]=0;
		avgfixes[18]=1;
		avgfixes[19]=1;
		avgfixes[20]=1;
		avgfixes[21]=1;

		double g01temp=Math.abs(g01[ncurves]);
		if(g01temp<0.01){
			g01temp=0.01;
		}
		double g02temp=Math.abs(g02[ncurves]);
		if(g02temp<0.01){
			g02temp=0.01;
		}
		double g0cctemp=Math.abs(g0cc[ncurves]);
		if(g0cctemp<0.01){
			g0cctemp=0.01;
		}
		double mult=10.0;

		avgconstraints=new double[2][22];
		avgconstraints[0][0]=1.0;
		avgconstraints[1][0]=20.0;
		avgconstraints[0][1]=-mult*g01temp;
		avgconstraints[0][2]=-mult*g01temp;
		avgconstraints[0][3]=0.01;
		avgconstraints[0][4]=-mult*g01temp;
		avgconstraints[0][5]=0.01;
		avgconstraints[0][6]=-2.0;
		avgconstraints[0][7]=0.1;
		avgconstraints[1][1]=mult*g01temp;
		avgconstraints[1][2]=mult*g01temp;
		avgconstraints[1][3]=100000.0;
		avgconstraints[1][4]=mult*g01temp;
		avgconstraints[1][5]=100000.0;
		avgconstraints[1][6]=2.0;
		avgconstraints[1][7]=1000.0;
		avgconstraints[0][8]=-mult*g02temp;
		avgconstraints[0][9]=-mult*g02temp;
		avgconstraints[0][10]=0.01;
		avgconstraints[0][11]=-mult*g02temp;
		avgconstraints[0][12]=0.01;
		avgconstraints[0][13]=-2.0;
		avgconstraints[0][14]=0.1;
		avgconstraints[1][8]=mult*g02temp;
		avgconstraints[1][9]=mult*g02temp;
		avgconstraints[1][10]=100000.0;
		avgconstraints[1][11]=mult*g02temp;
		avgconstraints[1][12]=100000.0;
		avgconstraints[1][13]=2.0;
		avgconstraints[1][14]=1000.0;
		avgconstraints[0][15]=-mult*g0cctemp;
		avgconstraints[0][16]=-mult*g0cctemp;
		avgconstraints[0][17]=0.01;
		avgconstraints[0][18]=-mult*g0cctemp;
		avgconstraints[0][19]=0.01;
		avgconstraints[0][20]=-2.0;
		avgconstraints[0][21]=0.1;
		avgconstraints[1][15]=mult*g0cctemp;
		avgconstraints[1][16]=mult*g0cctemp;
		avgconstraints[1][17]=100000.0;
		avgconstraints[1][18]=mult*g0cctemp;
		avgconstraints[1][19]=100000.0;
		avgconstraints[1][20]=2.0;
		avgconstraints[1][21]=1000.0;

		globalparams=new double[ncurves][22];
		globalvflmatrix=new int[ncurves][22];
		for(int i=0;i<ncurves;i++){
			globalparams[i][0]=5.0;
			globalparams[i][1]=0.0;
			globalparams[i][2]=g01[i];
			globalparams[i][3]=10.0;
			globalparams[i][4]=0.0;
			globalparams[i][5]=100.0;
			globalparams[i][6]=0.0;
			globalparams[i][7]=0.01;
			globalparams[i][8]=0.0;
			globalparams[i][9]=g02[i];
			globalparams[i][10]=10.0;
			globalparams[i][11]=0.0;
			globalparams[i][12]=100.0;
			globalparams[i][13]=0.0;
			globalparams[i][14]=0.01;
			globalparams[i][15]=0.0;
			globalparams[i][16]=g0cc[i];
			globalparams[i][17]=10.0;
			globalparams[i][18]=0.0;
			globalparams[i][19]=100.0;
			globalparams[i][20]=0.0;
			globalparams[i][21]=0.01;

			globalvflmatrix[i][0]=1;
			globalvflmatrix[i][1]=0;
			globalvflmatrix[i][2]=0;
			globalvflmatrix[i][3]=0;
			globalvflmatrix[i][4]=1;
			globalvflmatrix[i][5]=1;
			globalvflmatrix[i][6]=1;
			globalvflmatrix[i][7]=1;
			globalvflmatrix[i][8]=0;
			globalvflmatrix[i][9]=0;
			globalvflmatrix[i][10]=0;
			globalvflmatrix[i][11]=1;
			globalvflmatrix[i][12]=1;
			globalvflmatrix[i][13]=1;
			globalvflmatrix[i][14]=1;
			globalvflmatrix[i][15]=0;
			globalvflmatrix[i][16]=0;
			globalvflmatrix[i][17]=0;
			globalvflmatrix[i][18]=1;
			globalvflmatrix[i][19]=1;
			globalvflmatrix[i][20]=1;
			globalvflmatrix[i][21]=1;
		}

		globalconstraints=new double[2][ncurves][22];
		for(int i=0;i<ncurves;i++){
			globalconstraints[0][i][0]=1.0;
			globalconstraints[1][i][0]=20.0;
			globalconstraints[0][i][1]=-5.0*g01[ncurves];
			globalconstraints[0][i][2]=-5.0*g01[ncurves];
			globalconstraints[0][i][3]=0.01;
			globalconstraints[0][i][4]=-5.0*g01[ncurves];
			globalconstraints[0][i][5]=0.01;
			globalconstraints[0][i][6]=-2.0;
			globalconstraints[0][i][7]=0.1;
			globalconstraints[1][i][1]=5.0*g01[ncurves];
			globalconstraints[1][i][2]=5.0*g01[ncurves];
			globalconstraints[1][i][3]=100000.0;
			globalconstraints[1][i][4]=5.0*g01[ncurves];
			globalconstraints[1][i][5]=100000.0;
			globalconstraints[1][i][6]=2.0;
			globalconstraints[1][i][7]=1000.0;
			globalconstraints[0][i][8]=-5.0*g02[ncurves];
			globalconstraints[0][i][9]=-5.0*g02[ncurves];
			globalconstraints[0][i][10]=0.01;
			globalconstraints[0][i][11]=-5.0*g02[ncurves];
			globalconstraints[0][i][12]=0.01;
			globalconstraints[0][i][13]=-2.0;
			globalconstraints[0][i][14]=0.1;
			globalconstraints[1][i][8]=5.0*g02[ncurves];
			globalconstraints[1][i][9]=5.0*g02[ncurves];
			globalconstraints[1][i][10]=100000.0;
			globalconstraints[1][i][11]=5.0*g02[ncurves];
			globalconstraints[1][i][12]=100000.0;
			globalconstraints[1][i][13]=2.0;
			globalconstraints[1][i][14]=1000.0;
			globalconstraints[0][i][15]=-5.0*g0cc[ncurves];
			globalconstraints[0][i][16]=-5.0*g0cc[ncurves];
			globalconstraints[0][i][17]=0.01;
			globalconstraints[0][i][18]=-5.0*g0cc[ncurves];
			globalconstraints[0][i][19]=0.01;
			globalconstraints[0][i][20]=-2.0;
			globalconstraints[0][i][21]=0.1;
			globalconstraints[1][i][15]=5.0*g0cc[ncurves];
			globalconstraints[1][i][16]=5.0*g0cc[ncurves];
			globalconstraints[1][i][17]=100000.0;
			globalconstraints[1][i][18]=5.0*g0cc[ncurves];
			globalconstraints[1][i][19]=100000.0;
			globalconstraints[1][i][20]=2.0;
			globalconstraints[1][i][21]=1000.0;
		}

		for(int i=0;i<ncurves;i++){
			c2[i]=0.0;
			datapanel.table.setValueAt(""+(float)c2[i],i,8);
		}

		globalformulas=new String[ncurves][22];
		String[] tempnames={"z0/w0","offsetg","G01_g","td1_g(ms)","G02_g","td2_g(ms)","ftrip_g","ttrip_g(us)","offsetr","G01_r","td1_r(ms)","G02_r","td2_r(ms)","ftrip_r","ttrip_r(us)","offsetcc",
				"G01_cc","td1_cc(ms)","G02_cc","td2_cc(ms)","ftrip_cc","ttrip_cc(us)"};
		paramsnames=tempnames;

		undoparams=new double[ncurves][22];
		undovflmatrix=new int[ncurves][22];
		undoformulas=new String[ncurves][22];
		undoc2=new double[ncurves];
		undoglobalc2=0.0;
		undofit=new float[ncurves][3][npts];
	}

	private void fitavg(){
		if(showfitdialog()){
			double[] stats=new double[2];
			float[] tempdata=new float[3*npts];
			// float[] tempweights=new float[3*npts];
			for(int i=0;i<3;i++){
				for(int j=0;j<npts;j++){
					tempdata[i*npts+j]=avg[i][j];
				}
			}
			int tempmaxiter=fitclass.maxiter;
			if(checkc2){
				fitclass.maxiter=0;
			}
			float[] fit=fitclass.fitdata(avgparams,avgfixes,avgconstraints,tempdata,null,stats,false);
			fitclass.maxiter=tempmaxiter;
			for(int i=0;i<3;i++){
				for(int j=0;j<npts;j++){
					avgfit[i][j]=fit[i*npts+j];
				}
				pwavg.updateSeries(avgfit[i],i+3,true);
			}
			c2[ncurves]=stats[1];
			datapanel.table.setValueAt(""+(float)c2[ncurves],ncurves,8);
			g01[ncurves]=avgparams[2]+avgparams[4];
			g02[ncurves]=avgparams[9]+avgparams[11];
			g0cc[ncurves]=avgparams[16]+avgparams[18];
			datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
			datapanel.table.setValueAt(""+(float)g02[ncurves],ncurves,5);
			datapanel.table.setValueAt(""+(float)g0cc[ncurves],ncurves,6);
			updatebeta();
		}
	}

	private void fitglobal(){
		int nparams=22;
		int nsel=0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				nsel++;
			}
		}
		double[][] params=new double[nsel][nparams];
		String[][] tempformulas=new String[nsel][nparams];
		double[][][] constraints=new double[2][nsel][nparams];
		int[][] vflmatrix=new int[nsel][nparams];

		int counter=0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				for(int j=0;j<nparams;j++){
					params[counter][j]=globalparams[i][j];
					tempformulas[counter][j]=globalformulas[i][j];
					constraints[0][counter][j]=globalconstraints[0][i][j];
					constraints[1][counter][j]=globalconstraints[1][i][j];
					vflmatrix[counter][j]=globalvflmatrix[i][j];
				}
				counter++;
			}
			for(int j=0;j<nparams;j++){
				undoparams[i][j]=globalparams[i][j];
				undoformulas[i][j]=globalformulas[i][j];
				undovflmatrix[i][j]=globalvflmatrix[i][j];
			}
			for(int j=0;j<3;j++){
				for(int k=0;k<npts;k++){
					undofit[i][j][k]=fit[i][j][k];
				}
			}
			undoc2[i]=c2[i];
		}
		undoglobalc2=globalc2;
		if(showglobalfitdialog(params,tempformulas,vflmatrix)){
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<nparams;j++){
						globalparams[i][j]=params[counter][j];
						globalformulas[i][j]=tempformulas[counter][j];
						globalvflmatrix[i][j]=vflmatrix[counter][j];
					}
					counter++;
				}
			}
			double[] stats=new double[2];
			float[][] tempdata=new float[nsel][3*npts];
			// float[][] tempweights=new float[nsel][3*npts];
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<3;j++){
						for(int k=0;k<npts;k++){
							tempdata[counter][j*npts+k]=corr[i][j][k];
							// tempweights[counter][j+k*3]=weights[i][j][k];
						}
					}
					counter++;
				}
			}
			int tempmaxiter=globalfitclass.maxiter;
			if(checkc2){
				globalfitclass.changemaxiter(0);
			}
			double[] tempc2vals=new double[nsel];
			IJ.showStatus("Fitting Globally");
			float[][] tempfit=globalfitclass.fitdata(params,vflmatrix,tempformulas,paramsnames,constraints,tempdata,null,stats,tempc2vals,false);
			IJ.showStatus("Fit Complete");
			globalfitclass.changemaxiter(tempmaxiter);
			globalc2=stats[1];
			globalc2label.setText("Global chi^2 = "+(float)globalc2);
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<3;j++){
						for(int k=0;k<npts;k++){
							fit[i][j][k]=tempfit[counter][j*npts+k];
						}
						if(i==dispcurve){
							pwfit.updateSeries(fit[dispcurve][j],3+j,true);
						}
					}
					for(int j=0;j<nparams;j++){
						globalparams[i][j]=params[counter][j];
					}
					c2[i]=tempc2vals[counter];
					datapanel.table.setValueAt(""+(float)c2[i],i,8);
					g01[i]=globalparams[i][2]+globalparams[i][4];
					g02[i]=globalparams[i][9]+globalparams[i][11];
					g0cc[i]=globalparams[i][16]+globalparams[i][18];
					datapanel.table.setValueAt(""+(float)g01[i],i,3);
					datapanel.table.setValueAt(""+(float)g02[i],i,5);
					datapanel.table.setValueAt(""+(float)g0cc[i],i,6);
					counter++;
				}
			}
			updatebeta();
		}
	}
	
	private boolean showconstraintsdialog(){
		int nparams=avgparams.length;
		Object[][] tabledata=new Object[nparams][3];
		String[] columnlabels={"Parameters","Lower Limit","Upper Limit"};
		for(int i=0;i<nparams;i++){
			tabledata[i][0]=paramsnames[i];
			tabledata[i][1]=avgconstraints[0][i];
			tabledata[i][2]=avgconstraints[1][i];
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Avg Constraints",columnlabels,tabledata,null);
		if(retvals==null){
			return false;
		}
		for(int i=0;i<nparams;i++){
			avgconstraints[0][i]=((Double)retvals[i][1]).doubleValue();
			avgconstraints[1][i]=((Double)retvals[i][2]).doubleValue();
			for(int j=0;j<ncurves;j++){
				globalconstraints[0][j][i]=avgconstraints[0][i];
				globalconstraints[1][j][i]=avgconstraints[1][i];
			}
		}
		return true;
	}

	private boolean showfitdialog(){
		int nparams=avgparams.length;
		Object[][] tabledata=new Object[nparams+1][3];
		String[] columnlabels={"Parameters","Values","Fix?"};
		tabledata[0][0]="Check chi^2?";
		tabledata[0][1]=new Boolean(checkc2);
		tabledata[1][0]="z0/w0";
		tabledata[2][0]="offsetg";
		tabledata[3][0]="G01_g";
		tabledata[4][0]="td1_g(ms)";
		tabledata[5][0]="G02_g";
		tabledata[6][0]="td2_g(ms)";
		tabledata[7][0]="ftrip_g";
		tabledata[8][0]="ttrip_g(us)";
		tabledata[9][0]="offsetr";
		tabledata[10][0]="G01_r";
		tabledata[11][0]="td1_r(ms)";
		tabledata[12][0]="G02_r";
		tabledata[13][0]="td2_r(ms)";
		tabledata[14][0]="ftrip_r";
		tabledata[15][0]="ttrip_r(us)";
		tabledata[16][0]="offsetcc";
		tabledata[17][0]="G01_cc";
		tabledata[18][0]="td1_cc(ms)";
		tabledata[19][0]="G02_cc";
		tabledata[20][0]="td2_cc(ms)";
		tabledata[21][0]="ftrip_cc";
		tabledata[22][0]="ttrip_cc(us)";
		for(int i=0;i<nparams;i++){
			tabledata[i+1][1]=new Double(avgparams[i]);
			tabledata[i+1][2]=new Boolean(avgfixes[i]==1);
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Fit Parameters",columnlabels,tabledata,null);
		if(retvals==null){
			return false;
		}
		checkc2=((Boolean)retvals[0][1]).booleanValue();
		for(int i=0;i<nparams;i++){
			avgparams[i]=((Double)retvals[i+1][1]).doubleValue();
			avgfixes[i]=((Boolean)retvals[i+1][2]).booleanValue()?1:0;
		}
		return true;
	}

	private boolean showglobalfitdialog(double[][] params,String[][] formulas,int[][] vflmatrix){
		int nsets=params.length;
		int nparams=params[0].length;
		Object[][] tabledata=new Object[nparams+1][2+2*nsets];
		int[][] options=new int[nparams+1][2+2*nsets];
		String[] columnlabels=new String[2+2*nsets];
		columnlabels[0]="Parameters";
		columnlabels[1]="Linking";
		for(int i=0;i<nsets;i++){
			columnlabels[2*i+2]="values "+i;
			columnlabels[2*i+3]="link "+i;
		}
		String[] vfl={"Vary","Fix","Link","FLink"};
		String[] linking={"Custom","LinkAll","VaryAll"};
		tabledata[0][0]="Check chi^2?";
		tabledata[0][1]=new Boolean(checkc2);
		tabledata[1][0]="z0/w0";
		tabledata[2][0]="offsetg";
		tabledata[3][0]="G01_g";
		tabledata[4][0]="td1_g(ms)";
		tabledata[5][0]="G02_g";
		tabledata[6][0]="td2_g(ms)";
		tabledata[7][0]="ftrip_g";
		tabledata[8][0]="ttrip_g(us)";
		tabledata[9][0]="offsetr";
		tabledata[10][0]="G01_r";
		tabledata[11][0]="td1_r(ms)";
		tabledata[12][0]="G02_r";
		tabledata[13][0]="td2_r(ms)";
		tabledata[14][0]="ftrip_r";
		tabledata[15][0]="ttrip_r(us)";
		tabledata[16][0]="offsetcc";
		tabledata[17][0]="G01_cc";
		tabledata[18][0]="td1_cc(ms)";
		tabledata[19][0]="G02_cc";
		tabledata[20][0]="td2_cc(ms)";
		tabledata[21][0]="ftrip_cc";
		tabledata[22][0]="ttrip_cc(us)";
		for(int i=0;i<nparams;i++){
			tabledata[i+1][1]=linking;
			for(int j=0;j<nsets;j++){
				tabledata[i+1][2*j+2]=new String(""+(float)params[j][i]);
				tabledata[i+1][2*j+3]=vfl;
				options[i+1][2*j+3]=vflmatrix[j][i];
				if(vflmatrix[j][i]>2){
					tabledata[i+1][2*j+2]=formulas[j][i];
				}else{
					tabledata[i+1][2*j+2]=new String(""+(float)params[j][i]);
				}
			}
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Fit Parameters",columnlabels,tabledata,options);
		if(retvals==null){
			// IJ.showMessage("returned false");
			return false;
		}
		// IJ.showMessage("returned true");
		checkc2=((Boolean)retvals[0][1]).booleanValue();
		for(int i=0;i<nparams;i++){
			int linkindex=0;
			String retlink=(String)retvals[i+1][1];
			if(retlink.equals("LinkAll")){
				linkindex=1;
			}
			if(retlink.equals("VaryAll")){
				linkindex=2;
			}
			for(int j=0;j<nsets;j++){
				if(linkindex==0){
					String retvfl=(String)retvals[i+1][2*j+3];
					if(retvfl.equals("Vary")){
						vflmatrix[j][i]=0;
					}else{
						if(retvfl.equals("Fix")){
							vflmatrix[j][i]=1;
						}else{
							if(retvfl.equals("Link")){
								vflmatrix[j][i]=2;
							}else{
								vflmatrix[j][i]=3;
								// IJ.showMessage(""+j+" , "+i);
							}
						}
					}
				}else{
					String retvfl=(String)retvals[i+1][2*j+3];
					if(linkindex==1){
						if(retvfl.equals("Fix")&&j==0){
							vflmatrix[j][i]=1;
						}else{
							if(retvfl.equals("FLink")){
								vflmatrix[j][i]=3;
							}else{
								vflmatrix[j][i]=2;
							}
						}
					}else{
						if(retvfl.equals("FLink")){
							vflmatrix[j][i]=3;
						}else{
							vflmatrix[j][i]=0;
						}
					}
				}
				if(vflmatrix[j][i]<3){
					params[j][i]=Double.parseDouble((String)retvals[i+1][2*j+2]);
				}else{
					params[j][i]=0.0;
					formulas[j][i]=(String)retvals[i+1][2*j+2];
					// IJ.showMessage(formulas[j][i]);
				}
			}
		}
		return true;
	}

	private void undoglobalfit(){
		for(int i=0;i<ncurves;i++){
			for(int j=0;j<nparams;j++){
				globalparams[i][j]=undoparams[i][j];
				globalformulas[i][j]=undoformulas[i][j];
				globalvflmatrix[i][j]=undovflmatrix[i][j];
			}
			for(int j=0;j<3;j++){
				for(int k=0;k<npts;k++){
					fit[i][j][k]=undofit[i][j][k];
				}
				if(i==dispcurve){
					pwfit.updateSeries(fit[dispcurve][j],3+j,true);
				}
			}
			c2[i]=undoc2[i];
			datapanel.table.setValueAt(""+(float)c2[i],i,8);
		}
		globalc2=undoglobalc2;
		globalc2label.setText("Global chi^2 = "+(float)globalc2);
		updatebeta();
	}

	private void geterrors(){
		GenericDialog gd=new GenericDialog("Options");
		float conf=0.67f;
		gd.addNumericField("Confidence Limit",(int)(conf*100.0f),5,10,null);
		gd.addChoice("Error Parameter",paramsnames,paramsnames[0]);
		double spacing=0.01;
		gd.addNumericField("Chi^2 plot spacing (% of value)?",spacing*100.0,2,10,null);
		boolean globalerror=false;
		gd.addCheckbox("Global Fit Error?",globalerror);
		int dataset=0;
		gd.addNumericField("Data Set (for Global Error)",dataset,0);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		conf=0.01f*(float)gd.getNextNumber();
		int paramindex=gd.getNextChoiceIndex();
		spacing=0.01*gd.getNextNumber();
		globalerror=gd.getNextBoolean();
		dataset=(int)gd.getNextNumber();

		if(globalerror){
			support_plane_errors_v2 erclass=new support_plane_errors_v2(this,0.0001,50,true,0.1);
			int[] erindeces={paramindex,dataset};
			// need to set up all the matrices
			int nsel=0;
			int nparams=22;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					nsel++;
				}
			}
			double[][] params=new double[nsel][nparams];
			String[][] tempformulas=new String[nsel][nparams];
			double[][][] constraints=new double[2][nsel][nparams];
			int[][] vflmatrix=new int[nsel][nparams];

			float[][] tempdata=new float[nsel][3*npts];

			int nfit=0;
			int counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<nparams;j++){
						params[counter][j]=globalparams[i][j];
						tempformulas[counter][j]=globalformulas[i][j];
						constraints[0][counter][j]=globalconstraints[0][i][j];
						constraints[1][counter][j]=globalconstraints[1][i][j];
						vflmatrix[counter][j]=globalvflmatrix[i][j];
						if(vflmatrix[counter][j]==0||(j==0&&vflmatrix[counter][j]==2)){
							nfit++;
						}
					}
					for(int j=0;j<3;j++){
						for(int k=0;k<npts;k++){
							tempdata[counter][k+j*npts]=corr[i][j][k];
						}
					}
					counter++;
				}
			}
			int dofnum=3*npts*nsel-(nfit-1)-1;
			int dofden=3*npts*nsel-nfit-1;
			// double flim=FLimit(dofnum,dofden,(double)conf);
			double flim=(new jdist()).FLimit(dofnum,dofden,conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN&&flim<1.0){
				IJ.showMessage("Invalid Limiting F Value");
				return;
			}
			double truespacing=Math.abs(params[erindeces[1]][erindeces[0]]*spacing);
			double[][] c2plot=erclass.geterrorsglobal(params,vflmatrix,tempformulas,paramsnames,constraints,tempdata,null,flim,truespacing,erindeces);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			int templength=c2plot[0].length;
			float[][] c2plotf=new float[2][templength-1];
			for(int i=0;i<(templength-1);i++){
				c2plotf[0][i]=(float)c2plot[0][i+1];
				c2plotf[1][i]=(float)c2plot[1][i+1];
			}
			new PlotWindow4("c2 plot",paramsnames[paramindex]+"["+dataset+"]","Chi^2",c2plotf[0],c2plotf[1]).draw();
		}else{
			support_plane_errors_v2 erclass=new support_plane_errors_v2(this,0.0001,50,false,0.1);
			int errindex=paramindex;

			float[] tempdata=new float[npts*3];
			for(int i=0;i<3;i++){
				for(int j=0;j<npts;j++){
					tempdata[j+i*npts]=avg[i][j];
				}
			}

			int nfit=0;
			for(int i=0;i<7;i++){
				if(avgfixes[i]==0){
					nfit++;
				}
			}
			int dofnum=3*npts-(nfit-1)-1;
			int dofden=3*npts-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN&&flim<1.0){
				IJ.showMessage("Invalid Limiting F Value");
				return;
			}
			double truespacing=Math.abs(avgparams[errindex]*spacing);
			double[][] c2plot=erclass.geterrors(avgparams,avgfixes,avgconstraints,tempdata,null,flim,truespacing,errindex);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			int templength=c2plot[0].length;
			float[][] c2plotf=new float[2][templength-1];
			for(int i=0;i<(templength-1);i++){
				c2plotf[0][i]=(float)c2plot[0][i+1];
				c2plotf[1][i]=(float)c2plot[1][i+1];
			}
			new PlotWindow4("c2 plot",paramsnames[errindex],"Chi^2",c2plotf[0],c2plotf[1]).draw();
		}
	}
	
	private void montecarlo(){
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("#_of_trials",500,0);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		int ntrials=(int)gd.getNextNumber();
		float[] tempdata=new float[3*npts];
		// float[] tempweights=new float[3*npts];
		for(int i=0;i<3;i++){
			for(int j=0;j<npts;j++){
				tempdata[i*npts+j]=avg[i][j];
			}
		}
		monte_carlo_errors_v2 mc=new monte_carlo_errors_v2(this,0.0001,50,false,0.1);
		outerrors=true;
		float[] weights=null;
		boolean useweights=false;
		/*if(useweights&&ninclude>1){
			weights=new float[npts];
			for(int j=0;j<npts;j++){
				weights[j]=1.0f/(errs[j]*errs[j]);
			}
		}*/
		StringBuffer sb=new StringBuffer();
		sb.append("Trial\t");
		for(int i=0;i<paramsnames.length;i++){
			if(avgfixes[i]==0)
				sb.append(paramsnames[i]+"\t");
		}
		sb.append("chi^2");
		tw=new TextWindow("Monte Carlo Results",sb.toString(),"",400,400);
		outerrors=true;
		double[][] errors=mc.geterrors(avgparams,avgfixes,avgconstraints,tempdata,weights,ntrials);
		sb=new StringBuffer();
		sb.append("StDev\t");
		for(int i=0;i<errors.length;i++){
			float[] ferr=new float[errors[0].length];
			for(int j=0;j<ferr.length;j++)
				ferr[j]=(float)errors[i][j];
			float stdev=jstatistics.getstatistic("StDev",ferr,null);
			sb.append(""+stdev);
			if(i<(errors.length-1))
				sb.append("\t");
		}
		tw.append(sb.toString());
		outerrors=false;
	}

	private void getintbright(){
		for(int i=0;i<ncurves;i++){
			intensity1[i]=khz*avgs[0][i];
			intensity2[i]=khz*avgs[1][i];
			g01[i]=(vars[0][i]-avgs[0][i])/(avgs[0][i]*avgs[0][i]);
			g02[i]=(vars[1][i]-avgs[1][i])/(avgs[1][i]*avgs[1][i]);
			g0cc[i]=(vars[2][i])/(avgs[0][i]*avgs[1][i]);
			g0mincc[i]=(g01[i]*beta*intensity1[i])/intensity2[i];
			if(brightcorr){
				g01[i]*=khz*avgs[0][i];
				g02[i]*=khz*avgs[1][i];
				g0cc[i]*=khz*Math.sqrt(avgs[0][i]*avgs[1][i]);
				g0mincc[i]=g01[i]*beta*Math.sqrt(intensity1[i]/intensity2[i]);
			}
		}
	}

	private void updateavg(){
		avg=new float[3][npts];
		if(trajectories!=null){
			avgtraj=new float[2][trajectories[0][0].length];
		}
		int ninclude=0;
		double tempavg=0.0;
		double tempavg2=0.0;
		double tempvar=0.0;
		double tempvar2=0.0;
		double tempvar3=0.0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				ninclude++;
				for(int j=0;j<3;j++){
					for(int k=0;k<npts;k++){
						avg[j][k]+=corr[i][j][k];
					}
				}
				if(trajectories!=null){
					for(int j=0;j<2;j++){
						for(int k=0;k<trajectories[0][0].length;k++){
							avgtraj[j][k]+=trajectories[i][j][k];
						}
					}
				}
				tempavg+=avgs[0][i];
				tempavg2+=avgs[1][i];
				tempvar+=vars[0][i];
				tempvar2+=vars[1][i];
				tempvar3+=vars[2][i];
			}
		}
		for(int j=0;j<3;j++){
			for(int k=0;k<npts;k++){
				avg[j][k]/=ninclude;
			}
		}
		if(trajectories!=null){
			for(int j=0;j<2;j++){
				for(int k=0;k<trajectories[0][0].length;k++){
					avgtraj[j][k]/=ninclude;
				}
			}
		}
		tempavg/=ninclude;
		tempavg2/=ninclude;
		tempvar/=ninclude;
		tempvar2/=ninclude;
		tempvar3/=ninclude;
		intensity1[ncurves]=khz*tempavg;
		intensity2[ncurves]=khz*tempavg2;
		g01[ncurves]=(tempvar-tempavg)/(tempavg*tempavg);
		g02[ncurves]=(tempvar2-tempavg2)/(tempavg2*tempavg2);
		g0cc[ncurves]=tempvar3/(tempavg*tempavg2);
		g0mincc[ncurves]=(g01[ncurves]*beta*intensity1[ncurves])/intensity2[ncurves];
		if(brightcorr){
			g01[ncurves]*=khz*tempavg;
			g02[ncurves]*=khz*tempavg2;
			g0cc[ncurves]*=khz*Math.sqrt(tempavg*tempavg2);
			g0mincc[ncurves]=g01[ncurves]*beta*Math.sqrt(intensity1[ncurves]/intensity2[ncurves]);
		}
	}

	public void updatebeta(){
		for(int i=0;i<=ncurves;i++){
			g0mincc[i]=(g01[i]*beta*intensity1[i])/intensity2[i];
			if(brightcorr){
				g0mincc[i]=(g01[i]*beta)*Math.sqrt(intensity1[i]/intensity2[i]);
			}
			datapanel.table.setValueAt(""+(float)g0mincc[i],i,7);
		}
		update_bleed_curves();
	}

	public void update_bleed_curves(){
		if(pwfit!=null&&pwavg!=null){
			float[] temp=new float[npts];
			for(int j=0;j<npts;j++){
				if(!brightcorr){
					temp[j]=(float)(fit[dispcurve][0][j]*beta*intensity1[dispcurve]/intensity2[dispcurve]);
				}else{
					temp[j]=(float)(fit[dispcurve][0][j]*beta*Math.sqrt(intensity1[dispcurve]/intensity2[dispcurve]));
				}
			}
			pwfit.updateSeries(temp,6,true);
			for(int j=0;j<npts;j++){
				if(!brightcorr){
					temp[j]=(float)(avgfit[0][j]*beta*intensity1[ncurves]/intensity2[ncurves]);
				}else{
					temp[j]=(float)(avgfit[0][j]*beta*(float)Math.sqrt(intensity1[ncurves]/intensity2[ncurves]));
				}
			}
			pwavg.updateSeries(temp,6,true);
		}
	}

	public double[] fitfunc(double[] params){
		// params are
		// r,baselineg,G01g,td1g,G02g,td2g,ftripg,ttripg,baseliner,G01r,td1r,G02r,td2r,ftripr,ttripr,baselinecc,G01cc,td1cc,G02cc,td2cc,ftripcc,ttripcc;
		double[] fit=new double[3*npts];
		for(int i=0;i<3*npts;i++) {
    		int curveindex=(int)((double)i/(double)(3*npts));
    		int remvar=i-curveindex*3*npts;
    		int acindex=(int)((double)remvar/(double)(npts));
    		remvar-=npts*acindex;
    		double t_td1=xvals[0][0][remvar]/(params[acindex*7+3]/1000.0);
    		double t_td2=xvals[0][0][remvar]/(params[acindex*7+5]/1000.0);
    		double factor=1.0/(params[0]*params[0]);
    		double temp1=1.0/(1.0+t_td1);
    		boolean scanning=(psfflag==3);
    		if(scanning) psfflag=0;
    		if(psfflag==2)
    			temp1=Math.sqrt(temp1);
    		double temp2=0.0;
    		if(params[acindex*7+4]!=0.0){
    			temp2=1.0/(1.0+t_td2);
    			if(psfflag==2)
    				temp2=Math.sqrt(temp2);
    		}
    		if(psfflag!=1){
    			temp1*=1.0/Math.sqrt(1.0+factor*t_td1);
    			if(params[acindex*7+4]!=0.0){
    				temp2*=1.0/Math.sqrt(1.0+factor*t_td2);
    			}
    		}
    		if(scanning){
    			psfflag=3;
    			float tpix=xvals[0][0][1]-xvals[0][0][0]; //assume the x axis is linear along the time coordinate
    			float x=xvals[0][0][remvar]/tpix;
    			if(w0g==0.0) w0g=5.0; if(w0r==0.0) w0r=w0g*1.16; //this assumes 40 nm pixels and diffraction limited
    			float w0=(float)w0g;
    			if(curveindex==1) w0=(float)w0r;
    			if(curveindex==2) w0=(float)Math.sqrt(w0g*w0g+w0r*w0r);
    			temp1*=sfcs(x,w0,(float)t_td1);
    			temp2*=sfcs(x,w0,(float)t_td2);
    		}
    		double temp3=0.0;
    		if(params[acindex*7+6]!=0.0){
    			temp3=params[acindex*7+6]*Math.exp(-(double)xvals[0][0][remvar]/(params[acindex*7+7]/1000000.0));
    		}
    		fit[i]= params[acindex*7+1]+(params[acindex*7+2]*temp1+params[acindex*7+4]*temp2)*(1.0+temp3);
		}
		return fit;
	}
	
	public double sfcs(float x,float w0,float t_td){
		float retval=(float)Math.exp(x*x/(w0*w0+w0*w0*t_td));
		return retval;
	}

	public void showresults(String results){
		if(!outerrors)
			IJ.log(results);
		else
			tw.append(results);
	}

}
