/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import javax.swing.*;
import java.io.*;
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;
import java.awt.event.*;
import ij.io.*;

public class analysis_pch implements PlugIn {
	//this plugin is a gui for fitting photon counting histograms singly or globally for 
	//3D Gaussian or Gaussian Lorentzian squared point spread functions
	//(solution confocal, solution two photon)
	//copyright 2009 Jay Unruh, Stowers Institute for Medical Research

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		double sfreq=20000.0;
		gd.addNumericField("Sampling Frequency?",sfreq,1,10,null);
		String[] psfchoice={"3D Gaussian","Gaus-Lorentz^2","2D Gaussian"};
		gd.addChoice("PSF Type?",psfchoice,psfchoice[0]);
		String[] filetypechoice={"Confocor 3 raw","Short binary trajectory","PlotWindow trajectory","Ascii Text File"};
		gd.addChoice("File Type?",filetypechoice,filetypechoice[0]);
		boolean showtraj=false;
		gd.addCheckbox("Show Trajectories?",showtraj);
		int trajbin=20;
		gd.addNumericField("Bin Trajectory By?",trajbin,0,10,null);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=gd.getNextNumber();
		int psfflag=gd.getNextChoiceIndex();
		int fileflag=gd.getNextChoiceIndex();
		showtraj=gd.getNextBoolean();
		trajbin=(int)gd.getNextNumber();
		int nfiles=0;
		Object[] histograms=null;
		int max=0;
		String[] names=null;
		Object[] trajectories=null;
		int[] trajlengths=null;
		if(fileflag<2){
			jdataio ioclass=new jdataio();
			File[] filearray=ioclass.openfiles(OpenDialog.getDefaultDirectory(),IJ.getInstance());
			if(filearray.length==0){return;}
			String dir=filearray[0].getAbsolutePath();
			int sepindex=dir.lastIndexOf(File.separator);
			String newdir=dir.substring(0,sepindex+1);
			OpenDialog.setDefaultDirectory(newdir);
			nfiles=filearray.length;
			if(nfiles>25){nfiles=25;}
			histograms=new Object[nfiles];
			trajectories=new Object[nfiles];
			trajlengths=new int[nfiles];
			names=new String[nfiles+1];
			names[nfiles]="avg";
			for(int i=0;i<nfiles;i++){
				try{
					names[i]=filearray[i].getName();
					int length=(int)(((double)filearray[i].length()-128.0)/4.0);
					int length2=(int)(((double)filearray[i].length())/2.0);
					InputStream instream=new BufferedInputStream(new FileInputStream(filearray[i]));
					if(fileflag==0){
						int[] pmdata=new int[length];
						if(!ioclass.skipstreambytes(instream,128)){showioerror(); instream.close(); return;}
						if(!ioclass.readintelintfile(instream,length,pmdata)){showioerror(); instream.close(); return;}
						histograms[i]=(new pmodeconvert()).pm2pch(pmdata,sfreq,20000000);
						if(showtraj){
							trajectories[i]=(new pmodeconvert()).pm2tm(pmdata,sfreq/(double)trajbin,20000000);
							trajlengths[i]=((float[])trajectories[i]).length;
						}
					} else {
						float[] tmdata=new float[length2];
						if(!ioclass.readintelshortfile(instream,length2,tmdata)){showioerror(); instream.close(); return;}
						histograms[i]=(new pmodeconvert()).create_histogram(tmdata);
						if(showtraj){
							int tempnewlength=(int)(length/trajbin);
							float[] binned=new float[tempnewlength];
							for(int j=0;j<tempnewlength;j++){
								for(int k=0;k<trajbin;k++){
									binned[j]+=tmdata[j*trajbin+k];
								}
							}
							trajectories[i]=binned;
							trajlengths[i]=binned.length;
						}
					}
					if(((float[])histograms[i]).length>max){max=((float[])histograms[i]).length;}
					instream.close();
				} catch(IOException e){
					showioerror();
					return;
				}
			}
		} else {
			if(fileflag==2){
				ImageWindow iw=WindowManager.getCurrentWindow();
				float[][] trajectories2=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
				float[][] tempxvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
				sfreq=1.0/((double)tempxvals[0][1]);
				nfiles=trajectories2.length;
				if(nfiles>25){nfiles=25;}
				names=new String[nfiles+1];
				names[nfiles]="avg";
				histograms=new Object[nfiles];
				showtraj=false;
				for(int i=0;i<nfiles;i++){
					names[i]="trajectory "+(i+1);
					histograms[i]=(new pmodeconvert()).create_histogram(trajectories2[i]);
					if(((float[])histograms[i]).length>max){max=((float[])histograms[i]).length;}
				}
			} else {
				GenericDialog gd2=new GenericDialog("Ascii Import Options");
				boolean includex=false;
				gd2.addCheckbox("Include tab delim x column?",includex);
				int skiphead=0;
				gd2.addNumericField("Header Rows",skiphead,0);
				gd2.addCheckbox("Is Trajectory?",false);
				gd2.showDialog(); if(gd2.wasCanceled()){return;}
				includex=gd2.getNextBoolean();
				skiphead=(int)gd2.getNextNumber();
				boolean istraj=gd2.getNextBoolean();
				jdataio ioclass=new jdataio();
				File[] filearray=ioclass.openfiles(OpenDialog.getDefaultDirectory(),IJ.getInstance());
				if(filearray.length==0){return;}
				String dir=filearray[0].getAbsolutePath();
				int sepindex=dir.lastIndexOf(File.separator);
				String newdir=dir.substring(0,sepindex+1);
				OpenDialog.setDefaultDirectory(newdir);
				nfiles=filearray.length;
				if(nfiles>25){nfiles=25;}
				histograms=new Object[nfiles];
				names=new String[nfiles+1];
				names[nfiles]="avg";
				showtraj=false;
				for(int i=0;i<nfiles;i++){
					try{
						names[i]=filearray[i].getName();
						BufferedReader d=new BufferedReader(new FileReader(filearray[i]));
						int flines=0;
						while(d.readLine()!=null){flines++;}
						d.close();
						if(!istraj && flines>4096) flines=4096;
						d=new BufferedReader(new FileReader(filearray[i]));
						String[] temphist=new String[flines];
						int counter=0;
						do{
							temphist[counter]=d.readLine();
							//IJ.log(temphist[counter]);
							counter++;
						}while((temphist[counter-1]!=null && temphist[counter-1]!="") && counter<flines);
						float[] temphist2=new float[counter-1-skiphead];
						if(!istraj && (counter-1)>max){max=(counter-1);}
						if(!includex){
							for(int j=skiphead;j<(counter-1);j++){
								temphist2[j-skiphead]=(float)Integer.parseInt(temphist[j]);
							}
						} else {
							for(int j=skiphead;j<(counter-1);j++){
								int tabindex=temphist[j].indexOf('\t');
								temphist2[j-skiphead]=(float)Integer.parseInt(temphist[j].substring(tabindex+1));
							}
						}
						d.close();
						if(istraj){
							histograms[i]=(new pmodeconvert()).create_histogram(temphist2);
							int tempmax=((float[])histograms[i]).length+1;
							if(tempmax>max) max=tempmax;
						} else {
							histograms[i]=temphist2;
						}
					}
					catch(IOException e){
						showioerror();
						return;
					}
				}
			}	
		}
		float[][] pch=new float[nfiles][max];
		for(int i=0;i<nfiles;i++){
			System.arraycopy(((float[])histograms[i]),0,pch[i],0,((float[])histograms[i]).length);
		}

		float[][] disptraj=null;
		float[][] trajx=null;
		if(showtraj){
			int trajmax=0;
			for(int i=0;i<nfiles;i++){if(trajlengths[i]>trajmax){trajmax=trajlengths[i];}}
			disptraj=new float[nfiles][trajmax];
			trajx=new float[nfiles][trajmax];
			double deltat=(double)trajbin/sfreq;
			for(int j=0;j<nfiles;j++){
				System.arraycopy((float[])trajectories[j],0,disptraj[j],0,trajlengths[j]);
				for(int i=0;i<trajmax;i++){
					trajx[j][i]=(float)((double)i*deltat);
				}
			}
		}

		final PCHFitWindow cw = new PCHFitWindow();
		cw.init(names,pch,psfflag,disptraj,trajx,showtraj);

		final  Frame f = new Frame("PCH Analysis");
		f.setLocation(300,50);
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				f.dispose();
			}
		});

		f.add(cw);
		f.pack();
		f.setResizable(false);
		Insets ins = f.getInsets();
		cw.totalSize.height = PCHFitWindow.H + ins.bottom + ins.top + 65;
		cw.totalSize.width  = PCHFitWindow.WR + ins.left + ins.right;
		f.setSize(cw.totalSize);
		f.setVisible(true);
		cw.requestFocus();
		
	}

	private void showioerror(){
		IJ.showMessage("Error in file io");
	}

}

class PCHFitWindow extends Panel implements ActionListener, ItemListener, NLLSfitinterface {

	public final static int H = 600;
	public final static int WR = 650;
	public Dimension totalSize = new Dimension();

	private Choice dispcurvechoice;
	private Button fitavgbutton,fitglobalbutton,clearparamsbutton,undobutton,geterrorsbutton;
	private PlotWindow4 pwcurves,pwavg,pwfit,pwtraj;
	private int psfflag,ncurves,nparams,npts,dispcurve;
	private int[] nphotons,indices;
	private float[][] pch,fit,undofit,xvals,weights,traj,trajx;
	private float[] avg,avgfit,avgweights;
	private Label namelabel,intlabel,brightlabel,c2label,globalc2label,dispcurvelabel,copylabel,n_b_label,nlabel;
	private Checkbox[] checkarray;
	private boolean[] include;
	private TextField[] namearray,intarray,earray,narray,c2array;
	private String[] names;
	private double[] intensity,bright,c2,undoc2,number;
	private double[][] globalparams,avgconstraints,undoparams;
	private double[][][] globalconstraints;
	private double[] avgparams;
	private double globalc2,undoglobalc2;
	private boolean checkc2,showtraj;
	private int[][] globalvflmatrix,undovflmatrix;
	private String[][] globalformulas,undoformulas;
	private String[] paramsnames;
	private int[] avgfixes;
	NLLSglobalfit globalfitclass;
	NLLSfit fitclass;
	pch pchfunc;

	void init(String[] names1,float[][] pch1,int psfflag1,float[][] traj1,float[][] trajx1,boolean showtraj1){
		setLayout(null);
		names=names1;
		pch=pch1;
		psfflag=psfflag1;
		ncurves=pch.length;
		nparams=6;
		npts=pch[0].length;
		traj=traj1;
		trajx=trajx1;
		showtraj=showtraj1;

		checkarray=new Checkbox[ncurves];
		include=new boolean[ncurves];
		namearray=new TextField[ncurves+1];
		intarray=new TextField[ncurves+1];
		intensity=new double[ncurves+1];
		earray=new TextField[ncurves+1];
		narray=new TextField[ncurves+1];
		bright=new double[ncurves+1];
		number=new double[ncurves+1];
		c2array=new TextField[ncurves+1];
		c2=new double[ncurves+1];
		nphotons=new int[ncurves+1];
		avg=new float[npts];
		indices=new int[ncurves];

		getintbright();
		for(int i=0;i<ncurves;i++){include[i]=true; indices[i]=i;}
		updateavg();

		int starty=60;
		int startx=10;
		int yinc=25;
		for(int i=0;i<=ncurves;i++){
			if(i!=ncurves){
				checkarray[i]=new Checkbox("",include[i]);
				checkarray[i].setBounds(startx,starty+i*yinc,20,20);
				checkarray[i].addItemListener(this);
				add(checkarray[i]);
			}

			namearray[i]=new TextField(names[i]);
			namearray[i].setBounds(startx+30,starty+i*yinc,200,20);
			add(namearray[i]);

			intarray[i]=new TextField(""+(float)intensity[i]);
			intarray[i].setBounds(startx+30+210,starty+i*yinc,40,20);
			add(intarray[i]);

			earray[i]=new TextField(""+(float)bright[i]);
			earray[i].setBounds(startx+30+210+50,starty+i*yinc,40,20);
			add(earray[i]);

			narray[i]=new TextField(""+(float)number[i]);
			narray[i].setBounds(startx+30+210+50+50,starty+i*yinc,40,20);
			add(narray[i]);
			
			c2[i]=0.0;
			c2array[i]=new TextField(""+(float)c2[i]);
			c2array[i].setBounds(startx+30+210+50+50+50,starty+i*yinc,80,20);
			add(c2array[i]);
		}

		namelabel=new Label("Filename");
		namelabel.setBounds(startx+30,starty-25,100,20);
		add(namelabel);

		intlabel=new Label("<I>");
		intlabel.setBounds(startx+30+210,starty-25,40,20);
		add(intlabel);

		brightlabel=new Label("<e>");
		brightlabel.setBounds(startx+30+210+50,starty-25,40,20);
		add(brightlabel);

		nlabel=new Label("<N>");
		nlabel.setBounds(startx+30+210+50+50,starty-25,40,20);
		add(nlabel);

		c2label=new Label("chi^2");
		c2label.setBounds(startx+30+210+50+50+50,starty-25,80,20);
		add(c2label);

		int buttonsx=startx+30+210+50+50+50+90;

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

		fitclass=new NLLSfit(this,0.0001,50,0.1);
		globalfitclass=new NLLSglobalfit(this,0.0001,50,0.1);
		pchfunc=new pch((int)((double)npts*1.5),psfflag);
		avgfit=new float[npts];
		fit=new float[ncurves][npts];

		xvals=new float[ncurves][npts];
		for(int i=0;i<ncurves;i++){
			for(int j=0;j<npts;j++){
				xvals[i][j]=(float)j;
				fit[i][j]=1.0f;
			}
		}

		globalc2label=new Label("Global chi^2 = "+(float)0.0);
		globalc2label.setBounds(buttonsx,starty-25+50+50+50,140,20);
		add(globalc2label);

		dispcurvelabel=new Label("Display Fit #");
		dispcurvelabel.setBounds(buttonsx,starty-25+50+50+50+30,70,20);
		add(dispcurvelabel);

		dispcurvechoice=new Choice();
		for(int i=0;i<ncurves;i++){dispcurvechoice.add(""+(i+1));}
		dispcurve=0;
		dispcurvechoice.select(0);
		dispcurvechoice.setBounds(buttonsx+80,starty-25+50+50+50+30,40,20);
		dispcurvechoice.addItemListener(this);
		add(dispcurvechoice);

		undobutton=new Button("Undo Global Fit");
		undobutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50,100,40);
		undobutton.addActionListener(this);
		add(undobutton);

		geterrorsbutton=new Button("Get Errors");
		geterrorsbutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50+50,100,40);
		geterrorsbutton.addActionListener(this);
		add(geterrorsbutton);

		copylabel=new Label("copyright 2009 Jay Unruh (jru@stowers.org) non-profit use only");
		copylabel.setBounds(10,630,400,20);
		add(copylabel);

		n_b_label=new Label("N and B Analysis");
		n_b_label.setBounds(250,10,100,20);
		add(n_b_label);

		pwavg=new PlotWindow4("Avg","k","Frequency",xvals[0],avg);
		pwavg.setLogAxes(false,true);
		pwavg.draw();
		pwavg.addPoints(xvals[0],new float[npts],true);

		pwfit=new PlotWindow4("Selected Curve","k","Frequency",xvals[0],pch[dispcurve]);
		pwfit.setLogAxes(false,true);
		pwfit.draw();
		pwfit.addPoints(xvals[0],fit[dispcurve],true);
		pwfit.selectSeries(0);

		if(showtraj){
			pwtraj=new PlotWindow4("Selected Trajectory","time (s)","Intensity",trajx[0],traj[0]);
			pwtraj.draw();
		} else {
			pwtraj=null;
		}

		resetparams();
	}

	public void paint(Graphics g){
	}

	public void update(Graphics g)
	{
		paint(g);
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
	}

	public void itemStateChanged(ItemEvent e){
		boolean avgchanged=false;
		for(int i=0;i<ncurves;i++){
			if(checkarray[i].getState()!=include[i]){
				avgchanged=true;
				if(include[i]){
					include[i]=false;
					//pwcurves.deleteSeries(indices[i],true);
					for(int j=0;j<ncurves;j++){
						if(indices[j]>indices[i]){
							indices[j]--;
						}
					}
					indices[i]=-1;
				} else {
					include[i]=true;
					//pwcurves.addPoints(xvals[0],pch[i],true);
					//reorder the curves skipping those that aren't included
					int counter=0;
					for(int j=0;j<ncurves;j++){
						if(include[j]){
							//pwcurves.updateSeries(pch[j],counter,true);
							indices[j]=counter;
							counter++;
						}
					}
				}
			}
		}
		if(avgchanged){
			updateavg();
			intarray[ncurves].setText(""+(float)intensity[ncurves]);
			earray[ncurves].setText(""+(float)bright[ncurves]);
			narray[ncurves].setText(""+(float)number[ncurves]);
			pwavg.updateSeries(avg,0,true);
		}
		if(e.getSource()==dispcurvechoice){
			dispcurve=dispcurvechoice.getSelectedIndex();
			pwfit.updateSeries(pch[dispcurve],0,true);
			pwfit.updateSeries(fit[dispcurve],1,true);
			if(showtraj){
				pwtraj.updateSeries(trajx[dispcurve],traj[dispcurve],0,true);
			}
		}
	}

	private void resetparams(){
		avgparams=new double[7];
		avgparams[0]=0.0; avgparams[1]=bright[ncurves]; avgparams[2]=number[ncurves]; avgparams[3]=0.1; avgparams[4]=0.0; avgparams[5]=0.1; avgparams[6]=0.0;
		avgfixes=new int[7];
		avgfixes[0]=1; avgfixes[1]=0; avgfixes[2]=0; avgfixes[3]=1; avgfixes[4]=1; avgfixes[5]=1; avgfixes[6]=1;

		avgconstraints=new double[2][7];
		avgconstraints[0][0]=0.0001; avgconstraints[0][1]=0.0001; avgconstraints[0][2]=0.0; avgconstraints[0][3]=0.0001; avgconstraints[0][4]=0.0; avgconstraints[0][5]=0.0001; avgconstraints[0][6]=0.0; 
		avgconstraints[1][0]=intensity[ncurves]; avgconstraints[1][1]=100.0*bright[ncurves]; avgconstraints[1][2]=100.0*number[ncurves]; avgconstraints[1][3]=100.0*bright[ncurves];
		avgconstraints[1][4]=100.0*number[ncurves]; avgconstraints[1][5]=100.0*bright[ncurves]; avgconstraints[1][6]=100.0*number[ncurves];

		globalparams=new double[ncurves][7];
		globalvflmatrix=new int[ncurves][7];
		for(int i=0;i<ncurves;i++){
			globalparams[i][0]=0.0; globalparams[i][1]=bright[i]; globalparams[i][2]=number[i]; globalparams[i][3]=0.1; globalparams[i][4]=0.0; globalparams[i][5]=0.1; globalparams[i][6]=0.0;
			globalvflmatrix[i][0]=1; globalvflmatrix[i][1]=0; globalvflmatrix[i][2]=0; globalvflmatrix[i][3]=1; globalvflmatrix[i][4]=1; globalvflmatrix[i][5]=1; globalvflmatrix[i][6]=1;
		}

		globalconstraints=new double[2][ncurves][7];
		for(int i=0;i<ncurves;i++){
			globalconstraints[0][i][0]=0.0001; globalconstraints[0][i][1]=0.0001; globalconstraints[0][i][2]=0.0; globalconstraints[0][i][3]=0.0001; globalconstraints[0][i][4]=0.0; globalconstraints[0][i][5]=0.0001; globalconstraints[0][i][6]=0.0; 
			globalconstraints[1][i][0]=intensity[i]; globalconstraints[1][i][1]=100.0*bright[i]; globalconstraints[1][i][2]=100.0*number[i]; globalconstraints[1][i][3]=100.0*bright[i]; globalconstraints[1][i][4]=100.0*number[i];
			globalconstraints[1][i][5]=100.0*bright[i]; globalconstraints[1][i][6]=100.0*number[i];
		}

		for(int i=0;i<ncurves;i++){
			c2[i]=0.0;
			c2array[i].setText(""+(float)c2[i]);
		}

		globalformulas=new String[ncurves][7];
		String[] tempnames={"background","e1","N1","e2","N2","e3","N3"};
		paramsnames=tempnames;

		undoparams=new double[ncurves][7];
		undovflmatrix=new int[ncurves][7];
		undoformulas=new String[ncurves][7];
		undoc2=new double[ncurves];
		undoglobalc2=0.0;
		undofit=new float[ncurves][npts];
	}

	private void fitavg(){
		if(showfitdialog()){
			double[] stats=new double[2];
			float[] tempdata=new float[npts];
			for(int i=0;i<avg.length;i++){tempdata[i]=(float)((double)avg[i]/(double)nphotons[ncurves]);}
			int tempmaxiter=fitclass.maxiter;
			if(checkc2){fitclass.maxiter=0;}
			float[] fit=fitclass.fitdata(avgparams,avgfixes,avgconstraints,tempdata,avgweights,stats,false);
			fitclass.maxiter=tempmaxiter;
			for(int i=0;i<avg.length;i++){avgfit[i]=fit[i]*(float)nphotons[ncurves];}
			pwavg.updateSeries(avgfit,1,true);
			c2[ncurves]=stats[1];
			c2array[ncurves].setText(""+(float)c2[ncurves]);
		}
	}

	private void fitglobal(){
		int nsel=0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){nsel++;}
		}
		double[][] params=new double[nsel][7];
		String[][] tempformulas=new String[nsel][7];
		double[][][] constraints=new double[2][nsel][7];
		int[][] vflmatrix=new int[nsel][7];

		int counter=0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				for(int j=0;j<7;j++){
					params[counter][j]=globalparams[i][j];
					tempformulas[counter][j]=globalformulas[i][j];
					constraints[0][counter][j]=globalconstraints[0][i][j];
					constraints[1][counter][j]=globalconstraints[1][i][j];
					vflmatrix[counter][j]=globalvflmatrix[i][j];
				}
				counter++;
			}
			for(int j=0;j<7;j++){
				undoparams[i][j]=globalparams[i][j];
				undoformulas[i][j]=globalformulas[i][j];
				undovflmatrix[i][j]=globalvflmatrix[i][j];
			}
			for(int j=0;j<npts;j++){
				undofit[i][j]=fit[i][j];
			}
			undoc2[i]=c2[i];
		}
		undoglobalc2=globalc2;
		if(showglobalfitdialog(params,tempformulas,vflmatrix)){
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<7;j++){
						globalparams[i][j]=params[counter][j];
						globalformulas[i][j]=tempformulas[counter][j];
						globalvflmatrix[i][j]=vflmatrix[counter][j];
					}
					counter++;
				}
			}
			double[] stats=new double[2];
			float[][] tempdata=new float[nsel][npts];
			float[][] tempweights=new float[nsel][npts];
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<npts;j++){
						tempdata[counter][j]=(float)((double)pch[i][j]/(double)nphotons[i]);
						tempweights[counter][j]=weights[i][j];
					}
					counter++;
				}
			}
			int tempmaxiter=globalfitclass.maxiter;
			if(checkc2){globalfitclass.changemaxiter(0);}
			double[] tempc2vals=new double[nsel];
			IJ.showStatus("Fitting Globally");
			float[][] tempfit=globalfitclass.fitdata(params,vflmatrix,tempformulas,paramsnames,constraints,tempdata,tempweights,stats,tempc2vals,false);
			IJ.showStatus("Fit Complete");
			globalfitclass.changemaxiter(tempmaxiter);
			globalc2=stats[1];
			globalc2label.setText("Global chi^2 = "+(float)globalc2);
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<npts;j++){
						fit[i][j]=tempfit[counter][j]*(float)nphotons[i];
					}
					for(int j=0;j<7;j++){
						globalparams[i][j]=params[counter][j];
					}
					c2[i]=tempc2vals[counter];
					c2array[i].setText(""+(float)c2[i]);
					counter++;
				}
			}
			pwfit.updateSeries(fit[dispcurve],1,true);
		}
	}

	private boolean showfitdialog(){
		int nparams=avgparams.length;
		Object[][] tabledata=new Object[nparams+1][3];
		String[] columnlabels={"Parameters","Values","Fix?"};
		tabledata[0][0]="Check chi^2?"; tabledata[0][1]=new Boolean(checkc2);
		tabledata[1][0]="background"; tabledata[2][0]="e1"; tabledata[3][0]="N1"; tabledata[4][0]="e2"; tabledata[5][0]="N2"; tabledata[6][0]="e3"; tabledata[7][0]="N3";
		for(int i=0;i<nparams;i++){
			tabledata[i+1][1]=new Double(avgparams[i]);
			tabledata[i+1][2]=new Boolean(avgfixes[i]==1);
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Fit Parameters",columnlabels,tabledata,null);
		if(retvals==null){return false;}
		checkc2=((Boolean)retvals[0][1]).booleanValue();
		for(int i=0;i<nparams;i++){
			avgparams[i]=((Double)retvals[i+1][1]).doubleValue();
			avgfixes[i]=((Boolean)retvals[i+1][2]).booleanValue() ? 1:0;
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
		tabledata[1][0]="background"; tabledata[2][0]="e1"; tabledata[3][0]="N1"; tabledata[4][0]="e2"; tabledata[5][0]="N2"; tabledata[6][0]="e3"; tabledata[7][0]="N3";
		for(int i=0;i<nparams;i++){
			tabledata[i+1][1]=linking;
			for(int j=0;j<nsets;j++){
				tabledata[i+1][2*j+2]=new Double(params[j][i]);
				tabledata[i+1][2*j+3]=vfl;
				options[i+1][2*j+3]=vflmatrix[j][i];
				if(vflmatrix[j][i]>2){
					tabledata[i+1][2*j+2]=formulas[j][i];
				} else {
					tabledata[i+1][2*j+2]=new String(""+(float)params[j][i]);
				}
			}
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Fit Parameters",columnlabels,tabledata,options);
		if(retvals==null){return false;}
		checkc2=((Boolean)retvals[0][1]).booleanValue();
		for(int i=0;i<nparams;i++){
			int linkindex=0;
			String retlink=(String)retvals[i+1][1];
			if(retlink.equals("LinkAll")){linkindex=1;}
			if(retlink.equals("VaryAll")){linkindex=2;}
			for(int j=0;j<nsets;j++){
				if(linkindex==0){
					String retvfl=(String)retvals[i+1][2*j+3];
					if(retvfl.equals("Vary")){
						vflmatrix[j][i]=0;
					} else{
						if(retvfl.equals("Fix")){
							vflmatrix[j][i]=1;
						} else {
							if(retvfl.equals("Link")){
								vflmatrix[j][i]=2;
							} else {
								vflmatrix[j][i]=3;
								//IJ.showMessage(""+j+" , "+i);
							}
						}
					}
				} else {
					String retvfl=(String)retvals[i+1][2*j+3];
					if(linkindex==1){
						if(retvfl.equals("Fix") && j==0){vflmatrix[j][i]=1;}
						else{
							if(retvfl.equals("FLink")){vflmatrix[j][i]=3;}
							else{vflmatrix[j][i]=2;}
						}
					} else {
						if(retvfl.equals("FLink")){vflmatrix[j][i]=3;}
						else{vflmatrix[j][i]=0;}
					}
				}
				if(vflmatrix[j][i]<3){
					params[j][i]=Double.parseDouble((String)retvals[i+1][2*j+2]);
				} else{
					params[j][i]=0.0; formulas[j][i]=(String)retvals[i+1][2*j+2];
					//IJ.showMessage(formulas[j][i]);
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
			for(int j=0;j<npts;j++){
				fit[i][j]=undofit[i][j];
			}
			if(i==dispcurve){
				pwfit.updateSeries(fit[dispcurve],1,true);
			}
			c2[i]=undoc2[i];
			c2array[i].setText(""+(float)c2[i]);
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
		gd.showDialog(); if(gd.wasCanceled()){return;}
		conf=0.01f*(float)gd.getNextNumber();
		int paramindex=(int)gd.getNextChoiceIndex();
		spacing=0.01*gd.getNextNumber();
		globalerror=gd.getNextBoolean();
		dataset=(int)gd.getNextNumber();
		

		if(globalerror){
			support_plane_errors erclass=new support_plane_errors(this,0.0001,50,true,0.1);
			int[] erindeces={paramindex,dataset};
			//need to set up all the matrices
			int nsel=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){nsel++;}
			}
			double[][] params=new double[nsel][7];
			String[][] tempformulas=new String[nsel][7];
			double[][][] constraints=new double[2][nsel][7];
			int[][] vflmatrix=new int[nsel][7];
			float[][] tempdata=new float[nsel][npts];
			float[][] tempweights=new float[nsel][npts];

			int nfit=0;
			int counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<7;j++){
						params[counter][j]=globalparams[i][j];
						tempformulas[counter][j]=globalformulas[i][j];
						constraints[0][counter][j]=globalconstraints[0][i][j];
						constraints[1][counter][j]=globalconstraints[1][i][j];
						vflmatrix[counter][j]=globalvflmatrix[i][j];
						if(vflmatrix[counter][j]==0 || (j==0 && vflmatrix[counter][j]==2)){nfit++;}
					}
					for(int j=0;j<npts;j++){
						tempdata[counter][j]=(float)((double)pch[i][j]/(double)nphotons[i]);
						tempweights[counter][j]=weights[i][j];
					}
					counter++;
				}
			}
			int dofnum=npts*nsel-(nfit-1)-1;
			int dofden=npts*nsel-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,(double)conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN && flim < 1.0){IJ.showMessage("Invalid Limiting F Value"); return;}
			double truespacing=Math.abs(params[erindeces[1]][erindeces[0]]*spacing);
			double[][] c2plot=erclass.geterrorsglobal(params,vflmatrix,tempformulas,paramsnames,constraints,tempdata,tempweights,flim,truespacing,erindeces);
			IJ.log("upper limit = "+c2plot[1][0]+" lower limit = "+c2plot[0][0]);
			int templength=c2plot[0].length;
			float[][] c2plotf=new float[2][templength-1];
			for(int i=0;i<(templength-1);i++){
				c2plotf[0][i]=(float)c2plot[0][i+1];
				c2plotf[1][i]=(float)c2plot[1][i+1];
			}
			new PlotWindow4("c2 plot",paramsnames[paramindex]+"["+dataset+"]","Chi^2",c2plotf[0],c2plotf[1]).draw();
		} else {
			support_plane_errors erclass=new support_plane_errors(this,0.0001,50,false,0.1);
			int errindex=paramindex;
			float[] tempdata=new float[npts];
			for(int i=0;i<avg.length;i++){tempdata[i]=(float)((double)avg[i]/(double)nphotons[ncurves]);}
			int nfit=0;
			for(int i=0;i<7;i++){
				if(avgfixes[i]==0){
					nfit++;
				}
			}
			int dofnum=npts-(nfit-1)-1;
			int dofden=npts-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,(double)conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN && flim < 1.0){IJ.showMessage("Invalid Limiting F Value"); return;}
			double truespacing=Math.abs(avgparams[errindex]*spacing);
			double[][] c2plot=erclass.geterrors(avgparams,avgfixes,avgconstraints,tempdata,avgweights,flim,truespacing,errindex);
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

	private void getintbright(){
		weights=new float[ncurves][npts];
		for(int i=0;i<ncurves;i++){
			for(int j=0;j<npts;j++){nphotons[i]+=pch[i][j];}
			double tempavg=0.0; double tempavg2=0.0;
			for(int j=0;j<npts;j++){
				double normed=(double)pch[i][j]/(double)nphotons[i];
				if(pch[i][j]>0.0f){
					weights[i][j]=(float)((double)nphotons[i]/(normed*(1.0f-normed)));
				} else {
					weights[i][j]=1.0f;
				}
				tempavg+=(double)j*normed;
				tempavg2+=(double)j*(double)j*normed;
			}
			tempavg2-=tempavg*tempavg;
			tempavg2/=tempavg;
			bright[i]=(tempavg2-1.0);
		if(psfflag==0){
			bright[i]/=0.3536;
		} else {
			if(psfflag==1){
				bright[i]/=0.078;
			} else {
				bright[i]/=0.5;
			}
		}
			intensity[i]=(float)tempavg;
			number[i]=intensity[i]/bright[i];
		}
	}

	private void updateavg(){
		nphotons[ncurves]=0;
		avg=new float[npts];
		avgweights=new float[npts];
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				for(int j=0;j<npts;j++){
					avg[j]+=pch[i][j];
					nphotons[ncurves]+=pch[i][j];
				}
			}
		}
		double tempavg=0.0; double tempavg2=0.0;
		for(int i=0;i<npts;i++){
			double normed=(double)avg[i]/(double)nphotons[ncurves];
			if(avg[i]>0.0f){
				avgweights[i]=(float)((double)nphotons[ncurves]/(normed*(1.0f-normed)));
			} else {
				avgweights[i]=1.0f;
			}
			tempavg+=(double)i*normed;
			tempavg2+=(double)i*(double)i*normed;
		}
		tempavg2-=tempavg*tempavg;
		tempavg2/=tempavg;
		bright[ncurves]=(tempavg2-1.0);
		if(psfflag==0){
			bright[ncurves]/=0.3536;
		} else {
			if(psfflag==1){
				bright[ncurves]/=0.078;
			} else {
				bright[ncurves]/=0.5;
			}
		}
		intensity[ncurves]=(float)tempavg;
		number[ncurves]=intensity[ncurves]/bright[ncurves];
	}

	public double fitfunc(double[] params,int indvar){
		double[] brightnesses=new double[3];
		double[] numbers=new double[3];
		double background=params[0];
		double nmax=0.0;
		for(int i=0;i<3;i++){
			brightnesses[i]=params[2*i+1];
			numbers[i]=params[2*i+2];
			if(numbers[i]>nmax){nmax=numbers[i];}
		}
		int nlength=20;
		if(nmax>15.0f){nlength=(int)(nmax*1.3);}
		return pchfunc.getpch(indvar,brightnesses,numbers,background,nlength);
	}

	public void showresults(String results){
		IJ.log(results);
	}
			
}
