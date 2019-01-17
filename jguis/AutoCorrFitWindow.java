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
import ij.plugin.frame.Recorder;
import ij.text.TextWindow;
import jalgs.jdataio;
import jalgs.jdist;
import jalgs.jstatistics;
import jalgs.jfit.NLLSfit_v2;
import jalgs.jfit.NLLSfitinterface_v2;
import jalgs.jfit.NLLSglobalfit_v2;
import jalgs.jfit.monte_carlo_errors_v2;
import jalgs.jfit.support_plane_errors_v2;

import java.awt.Button;
import java.awt.Checkbox;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Insets;
import java.awt.Label;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
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

public class AutoCorrFitWindow extends Panel implements ActionListener,NLLSfitinterface_v2,ListSelectionListener,TableModelListener,ItemListener{
	public final static int H=775;
	public final static int WR=600;
	public Dimension totalSize=new Dimension();

	private Button fitavgbutton,fitglobalbutton,clearparamsbutton,undobutton,geterrorsbutton,outbutton,montecarlobutton,savebutton;
	private Button editconsbutton;
	private Checkbox useweightscheck;
	private PlotWindow4 pwavg,pwfit,pwavgtraj,pwfittraj;
	private int ncurves,nparams,npts,dispcurve,psfflag,ninclude;
	private int[] nmeas,indices;
	private float[][] corr,fit,undofit,xvals,trajectories;
	private float[] avg,errs,avgfit,avgs,vars,avgtraj;
	private Label globalc2label,copylabel;
	private TablePanel datapanel;
	private boolean[] include;
	private String[] names;
	private double[] intensity1,g01,c2,undoc2;
	private double[][] globalparams,avgconstraints,undoparams;
	private double[][][] globalconstraints;
	private double[] avgparams;
	private double globalc2,undoglobalc2,khz;
	private boolean checkc2,brightcorr,outerrors,useweights;
	private int[][] globalvflmatrix,undovflmatrix;
	private int[] avgfixes;
	private int[] fitplotcolors;
	private int[] avgplotcolors;
	private String[][] globalformulas,undoformulas;
	private String[] paramsnames;
	private TextWindow tw;
	NLLSglobalfit_v2 globalfitclass;
	NLLSfit_v2 fitclass;

	public static void launch_frame(AutoCorrFitWindow panel){
		final Frame f=new Frame("Auto Corr Analysis");
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
		panel.totalSize.height=AutoCorrFitWindow.H+ins.bottom+ins.top+65;
		panel.totalSize.width=AutoCorrFitWindow.WR+ins.left+ins.right;
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
			final AutoCorrFitWindow cw=new AutoCorrFitWindow();
			cw.init_from_is(is,brightcorr1);
			is.close();
			AutoCorrFitWindow.launch_frame(cw);
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
			if(temp.equals("jfcs_file_type") && type==0) return true;
			else return false;
		}catch(IOException e){
			return false;
		}
	}
	
	public void init_from_is(InputStream is,boolean brightcorr1) {
		jdataio jdio=new jdataio();
		String temp=jdio.readstring(is);
		int type=jdio.readintelint(is);
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
		float[][] corr=new float[nseries][npts];
		float[][] trajectories=new float[nseries][trajpts];
		float[] avg=new float[nseries];
		float[] var=new float[nseries];
		int[] selections=new int[nseries];
		for(int i=0;i<nseries;i++) {names[i]=jdio.readstring(is);}
		jdio.readintelfloatfile(is,npts,xvals);
		for(int i=0;i<nseries;i++) jdio.readintelfloatfile(is,npts,corr[i]);
		for(int i=0;i<nseries;i++) jdio.readintelfloatfile(is,trajpts,trajectories[i]);
		jdio.readintelfloatfile(is,nseries,avg);
		jdio.readintelfloatfile(is,nseries,var);
		jdio.readintelintfile(is,nseries,selections);
		//if we are changing the brightcorr status, do that here
		if(changebrightcorr) {
			if(brightcorr) { //here we undo brightcorr by dividing by avg
				for(int i=0;i<nseries;i++) {
					for(int j=0;j<npts;j++) corr[i][j]/=khz*avg[i];
				}
			} else { //here we apply brightcorr by multiplying by avg
				for(int i=0;i<nseries;i++) {
					for(int j=0;j<npts;j++) corr[i][j]*=khz*avg[i];
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
		pwavg.updateSeries(this.avg,0,true);
		float[][] temperrs=new float[2][];
		temperrs[0]=errs;
		temperrs[1]=new float[errs.length];
		pwavg.addErrors(temperrs);
		if(trajectories!=null){
			pwavgtraj.updateSeries(avgtraj,0,true);
		}
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
		jdio.writeintelint(os,0);
		jdio.writeintelint(os,ncurves);
		jdio.writeintelint(os,npts);
		float[][] temptraj=trajectories;
		if(trajectories==null) {
			temptraj=new float[ncurves][10];
			for(int i=0;i<ncurves;i++) {
				for(int j=0;j<10;j++) {
					temptraj[i][j]=avgs[i];
				}
			}
		}
		jdio.writeintelint(os,temptraj[0].length);
		jdio.writeintelfloat(os,(float)khz);
		jdio.writeintelint(os,psfflag);
		jdio.writeintelint(os,brightcorr?1:0);
		for(int i=0;i<ncurves;i++) {jdio.writestring(os,names[i]);}
		jdio.writeintelfloatarray(os,xvals[0]);
		for(int i=0;i<ncurves;i++) {jdio.writeintelfloatarray(os,corr[i]);}
		for(int i=0;i<ncurves;i++) {jdio.writeintelfloatarray(os,temptraj[i]);}
		jdio.writeintelfloatarray(os,avgs);
		jdio.writeintelfloatarray(os,vars);
		for(int i=0;i<ncurves;i++) {
			if(include[i]) jdio.writeintelint(os,1);
			else jdio.writeintelint(os,0);
		}
	}

	public void init(String[] names1,float[][] corr1,float[] xvals1,float[][] trajectories1,float[] avg1,float[] var1,float khz1,int psfflag1,boolean brightcorr1){
		setLayout(null);
		names=names1;
		corr=corr1;
		psfflag=psfflag1;
		avgs=avg1;
		vars=var1;
		trajectories=trajectories1;
		brightcorr=brightcorr1;
		khz=khz1;
		ncurves=corr.length;
		nparams=11;
		npts=corr[0].length;
		useweights=false;
		if(ncurves>1)
			useweights=true;

		include=new boolean[ncurves];
		intensity1=new double[ncurves+1];
		g01=new double[ncurves+1];
		c2=new double[ncurves+1];
		nmeas=new int[ncurves+1];
		avg=new float[npts];
		indices=new int[ncurves];

		getintbright();
		for(int i=0;i<ncurves;i++){
			include[i]=true;
			indices[i]=i;
		}
		updateavg();

		datapanel=new TablePanel();
		Object[][] tabledata=new Object[ncurves+1][5];
		String[] labels={"Include","Name","I (kHz)","G(0)","chi^2"};
		if(brightcorr){
			labels[3]="B (kHz)";
		}
		for(int i=0;i<=ncurves;i++){
			if(i!=ncurves){
				tabledata[i][0]=new Boolean(true);
				tabledata[i][1]=""+names[i];
			}else{
				tabledata[i][1]="Avg";
			}
			tabledata[i][2]=""+intensity1[i];
			tabledata[i][3]=""+g01[i];
			tabledata[i][4]=""+c2[i];
		}
		datapanel.init(labels,tabledata,null);
		datapanel.setBounds(10,30,420,750);
		datapanel.table.getSelectionModel().addListSelectionListener(this);
		datapanel.table.getModel().addTableModelListener(this);

		int starty=60;
		int startx=10;

		int buttonsx=startx+30+210+50+50+90;

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
		// fitclass=new NLLSfit(this,50);
		globalfitclass=new NLLSglobalfit_v2(this,0.0001,50,0.1);
		avgfit=new float[npts];
		fit=new float[ncurves][npts];

		xvals=new float[ncurves][npts];
		for(int i=0;i<ncurves;i++){
			for(int k=0;k<npts;k++){
				xvals[i][k]=xvals1[k];
				fit[i][k]=0.0f;
			}
		}

		globalc2label=new Label("Global chi^2 = "+(float)0.0);
		globalc2label.setBounds(buttonsx,starty-25+50+50+50,160,20);
		add(globalc2label);

		dispcurve=0;

		undobutton=new Button("Undo Global Fit");
		undobutton.setBounds(buttonsx,starty-25+50+50+50+30,100,40);
		undobutton.addActionListener(this);
		add(undobutton);

		geterrorsbutton=new Button("Get Errors");
		geterrorsbutton.setBounds(buttonsx,starty-25+50+50+50+30+50,100,40);
		geterrorsbutton.addActionListener(this);
		add(geterrorsbutton);

		outbutton=new Button("Plot All");
		outbutton.setBounds(buttonsx,starty-25+50+50+50+30+50+50,100,40);
		outbutton.addActionListener(this);
		add(outbutton);

		montecarlobutton=new Button("Monte Carlo");
		montecarlobutton.setBounds(buttonsx,starty-25+50+50+50+30+50+50+50,100,40);
		montecarlobutton.addActionListener(this);
		add(montecarlobutton);

		useweightscheck=new Checkbox("Use Weights");
		useweightscheck.setBounds(buttonsx,starty-25+50+50+50+30+50+50+50+50,100,20);
		useweightscheck.setState(useweights);
		useweightscheck.addItemListener(this);
		add(useweightscheck);

		editconsbutton=new Button("Edit Constraints");
		editconsbutton.setBounds(buttonsx,starty-25+50+50+50+30+50+50+50+50+30,100,40);
		editconsbutton.addActionListener(this);
		add(editconsbutton);
		
		savebutton=new Button("Save Analysis");
		savebutton.setBounds(buttonsx,starty-25+50+50+50+30+50+50+50+50+30+50,100,40);
		savebutton.addActionListener(this);
		add(savebutton);

		copylabel=new Label("copyright 2009 Jay Unruh (jru@stowers.org)");
		copylabel.setBounds(10,810,250,20);
		add(copylabel);

		if(brightcorr){
			pwavg=new PlotWindow4("Avg","tau","B (kHz)",xvals[0],avg);
		}else{
			pwavg=new PlotWindow4("Avg","tau","G(tau)",xvals[0],avg);
		}
		pwavg.setLogAxes(true,false);
		pwavg.draw();
		pwavg.addPoints(xvals[0],new float[npts],true);
		float[][] temperrs=new float[2][];
		temperrs[0]=errs;
		temperrs[1]=new float[errs.length];
		pwavg.addErrors(temperrs);

		float[] dumcorr=new float[corr[0].length];
		System.arraycopy(corr[dispcurve],0,dumcorr,0,corr[0].length);
		if(brightcorr){
			pwfit=new PlotWindow4("Selected Curve","tau","B (kHz)",xvals[0],dumcorr);
		}else{
			pwfit=new PlotWindow4("Selected Curve","tau","G(tau)",xvals[0],dumcorr);
		}
		pwfit.setLogAxes(true,false);
		pwfit.draw();
		pwfit.addPoints(xvals[0],new float[npts],true);

		if(trajectories!=null){
			pwavgtraj=new PlotWindow4("Avg Traj","sample","Intensity (kHz)",avgtraj);
			pwavgtraj.draw();

			float[] dumtraj=new float[trajectories[0].length];
			System.arraycopy(trajectories[dispcurve],0,dumtraj,0,trajectories[0].length);
			pwfittraj=new PlotWindow4("Selected Traj","sample","Intensity (kHz)",dumtraj);
			pwfittraj.draw();
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
		if(e.getSource()==outbutton){
			outall();
		}
		if(e.getSource()==montecarlobutton){
			montecarlo();
		}
		if(e.getSource()==editconsbutton){
			showconstraintsdialog();
		}
		if(e.getSource()==savebutton){
			savefile();
		}
	}

	public void valueChanged(ListSelectionEvent e){
		if(e.getSource()==datapanel.table.getSelectionModel()){
			int selrow=datapanel.table.getSelectedRow();
			dispcurve=selrow;
			if(dispcurve==ncurves){
				dispcurve--;
			}
			pwfit.updateSeries(corr[dispcurve],0,true);
			pwfit.updateSeries(fit[dispcurve],1,true);
			if(trajectories!=null){
				pwfittraj.updateSeries(trajectories[dispcurve],0,true);
			}
		}
	}

	public void tableChanged(TableModelEvent e){
		int row=e.getFirstRow();
		int column=e.getColumn();
		if(column==0){
			boolean valid=((Boolean)datapanel.table.getValueAt(row,0)).booleanValue();
			if(valid!=include[row]){
				include[row]=valid;
				updateavg();
				datapanel.table.setValueAt(""+(float)intensity1[ncurves],ncurves,2);
				datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
				pwavg.updateSeries(avg,0,true);
				float[][] temperrs=new float[2][];
				temperrs[0]=errs;
				temperrs[1]=new float[errs.length];
				pwavg.addErrors(temperrs);
				if(trajectories!=null){
					pwavgtraj.updateSeries(avgtraj,0,true);
				}
			}
		}
	}

	public void itemStateChanged(ItemEvent e){
		if(e.getSource()==useweightscheck){
			useweights=useweightscheck.getState();
			updateavg();
			float[][] temperrs=new float[2][];
			temperrs[0]=errs;
			temperrs[1]=new float[errs.length];
			pwavg.addErrors(temperrs);
		}
	}

	private void resetparams(){
		// params are r,baseline,G01,td1,G02,td2,ftrip,ttrip;
		avgparams=new double[8];
		avgparams[0]=5.0;
		avgparams[1]=0.0;
		avgparams[2]=g01[ncurves];
		avgparams[3]=10.0;
		avgparams[4]=0.0;
		avgparams[5]=100.0;
		avgparams[6]=0.0;
		avgparams[7]=250.0;
		avgfixes=new int[8];
		avgfixes[0]=1;
		avgfixes[1]=0;
		avgfixes[2]=0;
		avgfixes[3]=0;
		avgfixes[4]=1;
		avgfixes[5]=1;
		avgfixes[6]=1;
		avgfixes[7]=1;

		double g01temp=Math.abs(g01[ncurves]);
		if(g01temp<0.01){
			g01temp=0.01;
		}
		double mult=10.0;

		avgconstraints=new double[2][8];
		avgconstraints[0][0]=1.0;
		avgconstraints[1][0]=20.0;
		avgconstraints[0][1]=-mult*g01temp;
		avgconstraints[0][2]=0.0;
		avgconstraints[0][3]=0.01;
		avgconstraints[0][4]=0.0;
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

		globalparams=new double[ncurves][8];
		globalvflmatrix=new int[ncurves][8];
		for(int i=0;i<ncurves;i++){
			globalparams[i][0]=5.0;
			globalparams[i][1]=0.0;
			globalparams[i][2]=g01[i];
			globalparams[i][3]=10.0;
			globalparams[i][4]=0.0;
			globalparams[i][5]=100.0;
			globalparams[i][6]=0.0;
			globalparams[i][7]=0.01;
			globalvflmatrix[i][0]=1;
			globalvflmatrix[i][1]=0;
			globalvflmatrix[i][2]=0;
			globalvflmatrix[i][3]=0;
			globalvflmatrix[i][4]=1;
			globalvflmatrix[i][5]=1;
			globalvflmatrix[i][6]=1;
			globalvflmatrix[i][7]=1;
		}

		globalconstraints=new double[2][ncurves][8];
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
		}

		for(int i=0;i<ncurves;i++){
			c2[i]=0.0;
			datapanel.table.setValueAt(""+(float)c2[i],i,4);
		}

		globalformulas=new String[ncurves][8];
		String[] tempnames={"z0/w0","offset","G01","td1(ms)","G02","td2(ms)","ftrip","ttrip(us)"};
		paramsnames=tempnames;

		undoparams=new double[ncurves][8];
		undovflmatrix=new int[ncurves][8];
		undoformulas=new String[ncurves][8];
		undoc2=new double[ncurves];
		undoglobalc2=0.0;
		undofit=new float[ncurves][npts];
	}

	private void fitavg(){
		if(showfitdialog()){
			double[] stats=new double[2];
			float[] tempdata=new float[npts];
			float[] weights=null;
			if(useweights)
				weights=new float[npts];
			for(int j=0;j<npts;j++){
				tempdata[j]=avg[j];
				if(useweights)
					weights[j]=1.0f/(errs[j]*errs[j]);
			}
			int tempmaxiter=fitclass.maxiter;
			if(checkc2){
				fitclass.maxiter=0;
			}
			float[] fit=fitclass.fitdata(avgparams,avgfixes,avgconstraints,tempdata,weights,stats,false);
			fitclass.maxiter=tempmaxiter;
			for(int j=0;j<npts;j++){
				avgfit[j]=fit[j];
			}
			pwavg.updateSeries(avgfit,1,true);
			c2[ncurves]=stats[1];
			datapanel.table.setValueAt(""+(float)c2[ncurves],ncurves,4);
			g01[ncurves]=avgparams[2]+avgparams[4];
			datapanel.table.setValueAt(""+(float)g01[ncurves],ncurves,3);
		}
	}

	private void fitglobal(){
		int nparams=8;
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
			for(int k=0;k<npts;k++){
				undofit[i][k]=fit[i][k];
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
			float[][] tempdata=new float[nsel][npts];
			// float[][] tempweights=new float[nsel][npts];
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int k=0;k<npts;k++){
						tempdata[counter][k]=corr[i][k];
						// tempweights[counter][k]=weights[i][k];
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
					for(int k=0;k<npts;k++){
						fit[i][k]=tempfit[counter][k];
					}
					if(i==dispcurve){
						pwfit.updateSeries(fit[dispcurve],1,true);
					}
					for(int j=0;j<nparams;j++){
						globalparams[i][j]=params[counter][j];
					}
					c2[i]=tempc2vals[counter];
					datapanel.table.setValueAt(""+(float)c2[i],i,4);
					g01[i]=globalparams[i][2]+globalparams[i][4];
					datapanel.table.setValueAt(""+(float)g01[i],i,3);
					counter++;
				}
			}
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
		tabledata[2][0]="offset";
		tabledata[3][0]="G01";
		tabledata[4][0]="td1(ms)";
		tabledata[5][0]="G02";
		tabledata[6][0]="td2(ms)";
		tabledata[7][0]="ftrip";
		tabledata[8][0]="ttrip(us)";
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
		tabledata[2][0]="offset";
		tabledata[3][0]="G01";
		tabledata[4][0]="td1(ms)";
		tabledata[5][0]="G02";
		tabledata[6][0]="td2(ms)";
		tabledata[7][0]="ftrip";
		tabledata[8][0]="ttrip(us)";
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
			for(int k=0;k<npts;k++){
				fit[i][k]=undofit[i][k];
			}
			if(i==dispcurve){
				pwfit.updateSeries(fit[dispcurve],1,true);
			}
			c2[i]=undoc2[i];
			datapanel.table.setValueAt(""+(float)c2[i],i,4);
		}
		globalc2=undoglobalc2;
		globalc2label.setText("Global chi^2 = "+(float)globalc2);
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
			int nparams=8;
			support_plane_errors_v2 erclass=new support_plane_errors_v2(this,0.0001,50,true,0.1);
			int[] erindeces={paramindex,dataset};
			// need to set up all the matrices
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
			float[][] tempdata=new float[nsel][npts];

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
					for(int j=0;j<npts;j++){
						tempdata[counter][j]=corr[i][j];
					}
					counter++;
				}
			}
			int dofnum=npts*nsel-(nfit-1)-1;
			int dofden=npts*nsel-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN&&flim<1.0){
				IJ.showMessage("Invalid Limiting F Value");
				return;
			}
			double truespacing=Math.abs(params[erindeces[1]][erindeces[0]]*spacing);
			double[][] c2plot=erclass.geterrorsglobal(params,vflmatrix,tempformulas,paramsnames,constraints,tempdata,null,flim,truespacing,erindeces);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			IJ.log("upper error = "+(c2plot[1][0]-params[erindeces[1]][erindeces[0]])+" lower error = "+(params[erindeces[1]][erindeces[0]]-c2plot[0][0]));
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
			int nfit=0;
			for(int i=0;i<7;i++){
				if(avgfixes[i]==0){
					nfit++;
				}
			}
			int dofnum=npts-(nfit-1)-1;
			int dofden=npts-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN&&flim<1.0){
				IJ.showMessage("Invalid Limiting F Value");
				return;
			}
			double truespacing=Math.abs(avgparams[errindex]*spacing);
			double[][] c2plot=erclass.geterrors(avgparams,avgfixes,avgconstraints,avg,null,flim,truespacing,errindex);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			IJ.log("upper error = "+(c2plot[1][0]-avgparams[errindex])+" lower error = "+(avgparams[errindex]-c2plot[0][0]));
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
		monte_carlo_errors_v2 mc=new monte_carlo_errors_v2(this,0.0001,50,false,0.1);
		outerrors=true;
		float[] weights=null;
		if(useweights&&ninclude>1){
			weights=new float[npts];
			for(int j=0;j<npts;j++){
				weights[j]=1.0f/(errs[j]*errs[j]);
			}
		}
		StringBuffer sb=new StringBuffer();
		sb.append("Trial\t");
		for(int i=0;i<paramsnames.length;i++){
			if(avgfixes[i]==0)
				sb.append(paramsnames[i]+"\t");
		}
		sb.append("chi^2");
		tw=new TextWindow("Monte Carlo Results",sb.toString(),"",400,400);
		outerrors=true;
		double[][] errors=mc.geterrors(avgparams,avgfixes,avgconstraints,avg,weights,ntrials);
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
			intensity1[i]=khz*avgs[i];
			g01[i]=(vars[i]-avgs[i])/(avgs[i]*avgs[i]);
			if(brightcorr){
				g01[i]*=khz*avgs[i];
			}
		}
	}

	private void updateavg(){
		avg=new float[npts];
		errs=new float[npts];
		if(trajectories!=null){
			avgtraj=new float[trajectories[0].length];
		}
		ninclude=0;
		double tempavg=0.0;
		double tempvar=0.0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				ninclude++;
				for(int k=0;k<npts;k++){
					avg[k]+=corr[i][k];
					errs[k]+=corr[i][k]*corr[i][k];
				}
				if(trajectories!=null){
					for(int k=0;k<trajectories[0].length;k++){
						avgtraj[k]+=trajectories[i][k];
					}
				}
				tempavg+=avgs[i];
				tempvar+=vars[i];
			}
		}
		for(int k=0;k<npts;k++){
			avg[k]/=ninclude;
			if(useweights&&ninclude>1){
				errs[k]/=ninclude;
				errs[k]-=avg[k]*avg[k];
				errs[k]*=((float)ninclude)/((float)(ninclude-1));
				errs[k]=(float)Math.sqrt(errs[k]);
				errs[k]/=(float)Math.sqrt(ninclude);
			}else{
				errs[k]=0.0f;
			}
		}
		if(trajectories!=null){
			for(int k=0;k<trajectories[0].length;k++){
				avgtraj[k]/=ninclude;
			}
		}
		tempavg/=ninclude;
		tempvar/=ninclude;
		intensity1[ncurves]=khz*tempavg;
		g01[ncurves]=(tempvar-tempavg)/(tempavg*tempavg);
		if(brightcorr){
			g01[ncurves]*=khz*tempavg;
		}
	}

	private void outall(){
		ninclude=0;
		for(int i=0;i<ncurves;i++){
			if(include[i])
				ninclude++;
		}
		float[][] curves=new float[ninclude][];
		float[][] outx=new float[ninclude][];
		int counter=0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				curves[counter]=corr[i];
				outx[counter]=xvals[0];
				counter++;
			}
		}
		new PlotWindow4("All Selected Curves","tau (s)","G(tau)",outx,curves,null).draw();
	}

	public double[] fitfunc(double[] params){
		double[] fit=new double[npts];
		// params are r,baseline,G01,td1,G02,td2,ftrip,ttrip;
		for(int i=0;i<npts;i++){
			int indvar=i;
			double t_td1=xvals[0][indvar]/(params[3]/1000.0);
			double t_td2=xvals[0][indvar]/(params[5]/1000.0);
			double factor=1.0/(params[0]*params[0]);
			double temp1=1.0/(1.0+t_td1);
			if(psfflag==2)
				temp1=Math.sqrt(temp1);
			double temp2=0.0;
			if(params[4]!=0.0){
				temp2=1.0/(1.0+t_td2);
				if(psfflag==2)
					temp2=Math.sqrt(temp2);
			}
			if(psfflag!=1){
				temp1*=1.0/Math.sqrt(1.0+factor*t_td1);
				if(params[4]!=0.0){
					temp2*=1.0/Math.sqrt(1.0+factor*t_td2);
				}
			}
			double temp3=0.0;
			if(params[6]!=0.0){
				temp3=(params[6]/(1-params[6]))*Math.exp(-(double)xvals[0][indvar]/(params[7]/1000000.0));
			}
			fit[i]=params[1]+(params[2]*temp1+params[4]*temp2)*(1.0+temp3);
		}
		return fit;
	}

	public void showresults(String results){
		if(!outerrors)
			IJ.log(results);
		else
			tw.append(results);
	}

}
