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
import java.text.*;
import java.awt.datatransfer.*;

public class analysis_pch_2D implements PlugIn {
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
		boolean ch2green=true;
		gd.addCheckbox("Ch2 is green?",ch2green);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		sfreq=gd.getNextNumber();
		int psfflag=gd.getNextChoiceIndex();
		int fileflag=gd.getNextChoiceIndex();
		ch2green=gd.getNextBoolean();
		int nfiles=0;
		Object[] histograms=null;
		int xmax=0;
		int ymax=0;
		String[] names=null;
		if(fileflag<2){
			jdataio ioclass=new jdataio();
			File[] filearray=ioclass.openfiles(OpenDialog.getDefaultDirectory(),IJ.getInstance());
			if(filearray.length==0){return;}
			String dir=filearray[0].getAbsolutePath();
			int sepindex=dir.lastIndexOf(File.separator);
			String newdir=dir.substring(0,sepindex+1);
			OpenDialog.setDefaultDirectory(newdir);
			nfiles=filearray.length/2;
			if(nfiles>25){nfiles=25;}
			histograms=new Object[nfiles];
			names=organize_c3_files(filearray);
			for(int i=0;i<nfiles;i++){
				try{
					int length1=(int)(((double)filearray[2*i].length()-128.0)/4.0);
					int length2=(int)(((double)filearray[2*i+1].length()-128.0)/4.0);
					int length3=(int)(((double)filearray[2*i].length())/2.0);
					int length4=(int)(((double)filearray[2*i+1].length())/2.0);
					InputStream instream=new BufferedInputStream(new FileInputStream(filearray[2*i]));
					InputStream instream2=new BufferedInputStream(new FileInputStream(filearray[2*i+1]));
					if(fileflag==0){
						int[] pmdata=new int[length1];
						int[] pmdata2=new int[length2];
						if(!ioclass.skipstreambytes(instream,128)){showioerror(); instream.close(); return;}
						if(!ioclass.skipstreambytes(instream2,128)){showioerror(); instream2.close(); return;}
						if(!ioclass.readintelintfile(instream,length1,pmdata)){showioerror(); instream.close(); return;}
						if(!ioclass.readintelintfile(instream2,length2,pmdata2)){showioerror(); instream2.close(); return;}
						if(ch2green){
							histograms[i]=(new pmodeconvert()).pm2pch(pmdata2,pmdata,sfreq,20000000);
						} else {
							histograms[i]=(new pmodeconvert()).pm2pch(pmdata,pmdata2,sfreq,20000000);
						}
					} else {
						float[] tmdata=new float[length3];
						float[] tmdata2=new float[length4];
						if(!ioclass.readintelshortfile(instream,length3,tmdata)){showioerror(); instream.close(); return;}
						if(!ioclass.readintelshortfile(instream2,length4,tmdata2)){showioerror(); instream2.close(); return;}
						if(ch2green){
							histograms[i]=(new pmodeconvert()).create_2Dhistogram(tmdata2,tmdata);
						} else {
							histograms[i]=(new pmodeconvert()).create_2Dhistogram(tmdata,tmdata2);
						}
					}
					if(((float[][])histograms[i]).length>xmax){xmax=((float[][])histograms[i]).length;}
					if(((float[][])histograms[i])[0].length>ymax){ymax=((float[][])histograms[i])[0].length;}
					instream.close();
					instream2.close();
				} catch(IOException e){
					showioerror();
					return;
				}
			}
		} else {
			if(fileflag==2){
				ImageWindow iw=WindowManager.getCurrentWindow();
				float[][] trajectories=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
				float[][] tempxvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
				sfreq=1.0/((double)tempxvals[0][1]);
				nfiles=trajectories.length/2;
				if(nfiles>25){nfiles=25;}
				names=new String[nfiles+1];
				names[nfiles]="avg";
				histograms=new Object[nfiles];
				for(int i=0;i<nfiles;i++){
					names[i]="trajectory "+(i+1);
					if(ch2green){
						histograms[i]=(new pmodeconvert()).create_2Dhistogram(trajectories[2*i+1],trajectories[2*i]);
					} else {
						histograms[i]=(new pmodeconvert()).create_2Dhistogram(trajectories[2*i],trajectories[2*i+1]);
					}
					if(((float[][])histograms[i]).length>xmax){xmax=((float[][])histograms[i]).length;}
					if(((float[][])histograms[i])[0].length>ymax){ymax=((float[][])histograms[i])[0].length;}
				}
			} else {
				//here we read tab delimited lines from files
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
				for(int i=0;i<nfiles;i++){
					try{
						names[i]=filearray[i].getName();
						BufferedReader d=new BufferedReader(new FileReader(filearray[i]));
						String[] lines=new String[256];
						int counter=0;
						do{
							lines[counter]=d.readLine();
							counter++;
						}while((lines[counter-1]!=null && lines[counter-1]!="") && counter<256);
						int numcolumns=0;
						for(int j=0;j<counter-1;j++){
							int temp=getncolumns(lines[j]);
							if(temp>numcolumns){numcolumns=temp;}
						}
						float[][] temphist2=null;
						if(ch2green){
							temphist2=new float[numcolumns][counter-1];
						} else {
							temphist2=new float[counter-1][numcolumns];
						}
						for(int k=0;k<counter-1;k++){
							float[] temp=tab_delim2float(lines[k]);
							for(int j=0;j<numcolumns;j++){
								if(ch2green){
									temphist2[j][k]=temp[j];
								} else {
									temphist2[k][j]=temp[j];
								}
							}
						}							
						histograms[i]=temphist2;
						d.close();
					}
					catch(IOException e){
						showioerror();
						return;
					}
				}
				for(int i=0;i<nfiles;i++){
					if(((float[][])histograms[i]).length>xmax){xmax=((float[][])histograms[i]).length;}
					if(((float[][])histograms[i])[0].length>ymax){ymax=((float[][])histograms[i])[0].length;}
				}
			}
		}
		//note that here x is green and y is red
		float[][][] pch=new float[nfiles][xmax][ymax];
		for(int i=0;i<nfiles;i++){
			for(int j=0;j<((float[][])histograms[i]).length;j++){
				for(int k=0;k<((float[][])histograms[i])[j].length;k++){
					pch[i][j][k]=((float[][])histograms[i])[j][k];
				}
			}
		}

		final PCH2DFitWindow cw = new PCH2DFitWindow();
		cw.init(names,pch,psfflag);

		final  Frame f = new Frame("PCH 2D Analysis");
		f.setLocation(10,10);
		f.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				f.dispose();
			}
		});

		f.add(cw);
		f.pack();
		f.setResizable(false);
		Insets ins = f.getInsets();
		cw.totalSize.height = PCH2DFitWindow.H + ins.bottom + ins.top + 65;
		cw.totalSize.width  = PCH2DFitWindow.WR + ins.left + ins.right;
		f.setSize(cw.totalSize);
		f.setVisible(true);
		cw.requestFocus();
		
	}

	private int getncolumns(String input){
		int counter=0;
		int index=-1;
		do{
			index=input.indexOf('\t',index+1);
			counter++;
		}while(index>=0);
		return counter;
	}

	private float[] tab_delim2float(String line){
		int numcolumns=getncolumns(line);
		String templine=line.toString();
		float[] tempdata=new float[numcolumns];
		int tabindex=-1;
		for(int i=0;i<(numcolumns-1);i++){
			tabindex=templine.indexOf('\t');
			if(tabindex<0){break;}
			tempdata[i]=Float.parseFloat(templine.substring(0,tabindex));
			templine=templine.substring(tabindex+1);
		}
		if(tabindex>=0){tempdata[numcolumns-1]=Float.parseFloat(templine.substring(0,templine.length()));}
		if(numcolumns==1){
			tempdata[0]=Float.parseFloat(templine.substring(0,templine.length()));
		}
		return tempdata;
	}
		

	private String[] organize_c3_files(File[] filearray){
		int length=filearray.length;
		int[] assign=new int[length];
		for(int i=0;i<length;i++){
			assign[i]=-1;
		}
		int counter=0;
		File[] temp=new File[length];
		String[] outnames=new String[length/2+1];
		outnames[length/2]="avg";
		for(int i=0;i<length;i++){
			if(assign[i]<0){
				String bait=filearray[i].getName();
				String baittrunc=bait.substring(0,bait.length()-5);
				for(int j=i+1;j<length;j++){
					if(assign[j]<0){
						String target=filearray[j].getName();
						if(target.substring(0,target.length()-5).equals(baittrunc)){
							assign[j]=counter;
							assign[i]=counter;
							if(bait.charAt(bait.length()-5)>target.charAt(target.length()-5)){
								temp[2*counter]=filearray[i];
								temp[2*counter+1]=filearray[j];
								outnames[counter]=target.substring(0,target.length()-8);
							} else {
								temp[2*counter]=filearray[j];
								temp[2*counter+1]=filearray[i];
								outnames[counter]=bait.substring(0,bait.length()-8);
							}
							counter++;
							break;
						}
					}
				}
			}
		}
		filearray=temp;
		return outnames;
	}

	private void showioerror(){
		IJ.showMessage("Error in file io");
	}

}

class PCH2DFitWindow extends Panel implements ActionListener, ItemListener, NLLSfitinterface,ClipboardOwner {

	public final static int H = 750;
	public final static int WR = 950;
	public Dimension totalSize = new Dimension();

	private Choice dispcurvechoice;
	private Button fitavgbutton,fitglobalbutton,clearparamsbutton,undobutton,geterrorsbutton;
	private PlotWindow3D pwavg,pwfit;
	private int psfflag,ncurves,nparams,xpts,ypts,dispcurve;
	private int[] nmeas,indices;
	private float[][][] pch,fit,undofit,xvals,weights;
	private float[][] avg,avgfit,avgweights;
	private Label namelabel,intlabel,brightlabel,c2label,globalc2label,dispcurvelabel,copylabel,n_b_label,nlabel,int2label,bright2label,n2label,brightcclabel,brightccminlabel,betalabel;
	private Checkbox[] checkarray;
	private boolean[] include;
	private TextField[] namearray,int1array,e1array,n1array,int2array,e2array,n2array,eccarray,eminccarray,c2array;
	private TextField betaval;
	private String[] names;
	private double[] intensity1,bright1,c2,undoc2,number1,intensity2,bright2,number2,brightcc,brightmincc;
	private double[][] globalparams,avgconstraints,undoparams;
	private double[][][] globalconstraints;
	private double[] avgparams;
	private double globalc2,beta,undoglobalc2;
	private boolean checkc2;
	private int[][] globalvflmatrix,undovflmatrix;
	private int[] avgfixes;
	private int[] fitplotshapes;
	private int[] avgplotshapes;
	private String[][] globalformulas,undoformulas;
	private String[] paramsnames;
	NLLSglobalfit globalfitclass;
	NLLSfit fitclass;
	pch2D pchfunc;

	void init(String[] names1,float[][][] pch1,int psfflag1){
		setLayout(null);
		names=names1;
		pch=pch1;
		psfflag=psfflag1;
		ncurves=pch.length;
		nparams=11;
		xpts=pch[0].length;
		ypts=pch[0][0].length;

		checkarray=new Checkbox[ncurves];
		include=new boolean[ncurves];
		namearray=new TextField[ncurves+1];
		int1array=new TextField[ncurves+1];
		intensity1=new double[ncurves+1];
		e1array=new TextField[ncurves+1];
		n1array=new TextField[ncurves+1];
		bright1=new double[ncurves+1];
		number1=new double[ncurves+1];
		int2array=new TextField[ncurves+1];
		intensity2=new double[ncurves+1];
		e2array=new TextField[ncurves+1];
		n2array=new TextField[ncurves+1];
		bright2=new double[ncurves+1];
		number2=new double[ncurves+1];
		eccarray=new TextField[ncurves+1];
		brightcc=new double[ncurves+1];
		eminccarray=new TextField[ncurves+1];
		brightmincc=new double[ncurves+1];
		c2array=new TextField[ncurves+1];
		c2=new double[ncurves+1];
		nmeas=new int[ncurves+1];
		avg=new float[xpts][ypts];
		indices=new int[ncurves];
		beta=0.05;

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

			int1array[i]=new TextField(""+(float)intensity1[i]);
			int1array[i].setBounds(startx+30+210,starty+i*yinc,40,20);
			add(int1array[i]);

			e1array[i]=new TextField(""+(float)bright1[i]);
			e1array[i].setBounds(startx+30+210+50,starty+i*yinc,40,20);
			add(e1array[i]);

			n1array[i]=new TextField(""+(float)number1[i]);
			n1array[i].setBounds(startx+30+210+50+50,starty+i*yinc,40,20);
			add(n1array[i]);

			int2array[i]=new TextField(""+(float)intensity2[i]);
			int2array[i].setBounds(startx+30+210+50+50+50,starty+i*yinc,40,20);
			add(int2array[i]);

			e2array[i]=new TextField(""+(float)bright2[i]);
			e2array[i].setBounds(startx+30+210+50+50+50+50,starty+i*yinc,40,20);
			add(e2array[i]);

			n2array[i]=new TextField(""+(float)number2[i]);
			n2array[i].setBounds(startx+30+210+50+50+50+50+50,starty+i*yinc,40,20);
			add(n2array[i]);

			eccarray[i]=new TextField(""+(float)brightcc[i]);
			eccarray[i].setBounds(startx+30+210+50+50+50+50+50+50,starty+i*yinc,40,20);
			add(eccarray[i]);

			eminccarray[i]=new TextField(""+(float)brightmincc[i]);
			eminccarray[i].setBounds(startx+30+210+50+50+50+50+50+50+50,starty+i*yinc,40,20);
			add(eminccarray[i]);
			
			c2[i]=0.0;
			c2array[i]=new TextField(""+(float)c2[i]);
			c2array[i].setBounds(startx+30+210+50+50+50+50+50+50+50+50,starty+i*yinc,80,20);
			add(c2array[i]);
		}

		namelabel=new Label("Filename");
		namelabel.setBounds(startx+30,starty-25,100,20);
		add(namelabel);

		intlabel=new Label("<Ig>");
		intlabel.setBounds(startx+30+210,starty-25,40,20);
		add(intlabel);

		brightlabel=new Label("<eg>");
		brightlabel.setBounds(startx+30+210+50,starty-25,40,20);
		add(brightlabel);

		nlabel=new Label("<Ng>");
		nlabel.setBounds(startx+30+210+50+50,starty-25,40,20);
		add(nlabel);

		int2label=new Label("<Ir>");
		int2label.setBounds(startx+30+210+50+50+50,starty-25,40,20);
		add(int2label);

		bright2label=new Label("<er>");
		bright2label.setBounds(startx+30+210+50+50+50+50,starty-25,40,20);
		add(bright2label);

		n2label=new Label("<Nr>");
		n2label.setBounds(startx+30+210+50+50+50+50+50,starty-25,40,20);
		add(n2label);

		brightcclabel=new Label("<ecc>");
		brightcclabel.setBounds(startx+30+210+50+50+50+50+50+50,starty-25,40,20);
		add(brightcclabel);

		brightccminlabel=new Label("min");
		brightccminlabel.setBounds(startx+30+210+50+50+50+50+50+50+50,starty-25,40,20);
		add(brightccminlabel);

		c2label=new Label("chi^2");
		c2label.setBounds(startx+30+210+50+50+50+50+50+50+50+50,starty-25,80,20);
		add(c2label);

		int buttonsx=startx+30+210+50+50+50+50+50+50+50+50+90;

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
		pchfunc=new pch2D((int)((double)xpts*1.5),(int)((double)ypts*1.5),psfflag);
		avgfit=new float[xpts][ypts];
		fit=new float[ncurves][xpts][ypts];

		xvals=new float[ncurves][xpts][ypts];
		for(int i=0;i<ncurves;i++){
			for(int j=0;j<xpts;j++){
				for(int k=0;k<ypts;k++){
					xvals[i][j][k]=(float)k;
					fit[i][j][k]=1.0f;
				}
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

		betalabel=new Label("Bleedthrough f");
		betalabel.setBounds(buttonsx,starty-25+50+50+50+30+30,100,20);
		add(betalabel);
		betaval=new TextField(""+(float)beta);
		betaval.setBounds(buttonsx+110,starty-25+50+50+50+30+30,40,20);
		betaval.addActionListener(this);
		add(betaval);

		beta=Double.parseDouble(betaval.getText());
		updatebeta();

		undobutton=new Button("Undo Global Fit");
		undobutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50,100,40);
		undobutton.addActionListener(this);
		add(undobutton);

		geterrorsbutton=new Button("Get Errors");
		geterrorsbutton.setBounds(buttonsx,starty-25+50+50+50+30+30+50+50,100,40);
		geterrorsbutton.addActionListener(this);
		add(geterrorsbutton);

		copylabel=new Label("copyright 2009 Jay Unruh (jru@stowers.org) non-profit use only");
		copylabel.setBounds(10,790,400,20);
		add(copylabel);

		n_b_label=new Label("N and B Analysis");
		n_b_label.setBounds(250,10,100,20);
		add(n_b_label);

		pwavg=new PlotWindow3D("Avg","kg","kr","Frequency",avg,0);
		pwavg.setLogAxes(false,false,true);
		pwavg.draw();
		pwavg.addPoints(avgfit,true,0);
		float[] temp=pwavg.getLimits(); temp[4]=1.0f; pwavg.setLimits(temp);

		float[][] temppch=new float[xpts][ypts];
		for(int i=0;i<xpts;i++){
			System.arraycopy(pch[dispcurve][i],0,temppch[i],0,ypts);
		}
		pwfit=new PlotWindow3D("Selected Curve","kg","kr","Frequency",temppch,0);
		pwfit.setLogAxes(false,false,true);
		pwfit.draw();
		pwfit.addPoints(fit[dispcurve],true,0);
		float[] temp2=pwfit.getLimits(); temp2[4]=1.0f; pwfit.setLimits(temp2);

		resetparams();
		repaint();
	}

	public void paint(Graphics g){
	}

	public void update(Graphics g)
	{
		paint(g);
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
	}

	public void lostOwnership(java.awt.datatransfer.Clipboard clipboard, Transferable contents) {}

	public void itemStateChanged(ItemEvent e){
		beta=Double.parseDouble(betaval.getText());
		updatebeta();
		boolean avgchanged=false;
		for(int i=0;i<ncurves;i++){
			if(checkarray[i].getState()!=include[i]){
				avgchanged=true;
				if(include[i]){
					include[i]=false;
				} else {
					include[i]=true;
				}
			}
		}
		if(avgchanged){
			updateavg();
			updatebeta();
			int1array[ncurves].setText(""+(float)intensity1[ncurves]);
			e1array[ncurves].setText(""+(float)bright1[ncurves]);
			n1array[ncurves].setText(""+(float)number1[ncurves]);
			int2array[ncurves].setText(""+(float)intensity2[ncurves]);
			e2array[ncurves].setText(""+(float)bright2[ncurves]);
			n2array[ncurves].setText(""+(float)number2[ncurves]);
			eccarray[ncurves].setText(""+(float)brightcc[ncurves]);
			pwavg.updateSeries(avg,0,true);
			float[] temp=pwavg.getLimits(); temp[4]=1.0f; pwavg.setLimits(temp);
		}
		if(e.getSource()==dispcurvechoice){
			dispcurve=dispcurvechoice.getSelectedIndex();
			pwfit.updateSeries(pch[dispcurve],0,true);
			pwfit.updateSeries(fit[dispcurve],1,true);
			float[] temp=pwfit.getLimits(); temp[4]=1.0f; pwfit.setLimits(temp);
		}
	}

	private void resetparams(){
		avgparams=new double[11];
		avgparams[0]=0.0; avgparams[1]=0.0; avgparams[2]=bright1[ncurves]; avgparams[3]=bright1[ncurves]*beta; avgparams[4]=number1[ncurves];
		avgfixes=new int[11];
		avgfixes[0]=1; avgfixes[1]=1; avgfixes[2]=1; avgfixes[3]=1; avgfixes[4]=0;
		for(int i=0;i<2;i++){
			avgparams[3*i+5]=0.1; avgfixes[3*i+5]=1;
			avgparams[3*i+6]=0.1; avgfixes[3*i+6]=1;
			avgparams[3*i+7]=0.0; avgfixes[3*i+7]=1;
		}
		avgparams[5]=0.0; avgparams[6]=bright2[ncurves]; avgparams[7]=number2[ncurves];
		avgfixes[7]=0;

		avgconstraints=new double[2][11];
		avgconstraints[0][0]=0.0001; avgconstraints[1][0]=intensity1[ncurves]; avgconstraints[0][1]=0.0001; avgconstraints[1][1]=intensity2[ncurves];
		for(int i=0;i<3;i++){
			avgconstraints[0][3*i+2]=0.0001; avgconstraints[0][3*i+3]=0.0001; avgconstraints[0][3*i+4]=0.0;
			avgconstraints[1][3*i+2]=100.0*bright1[ncurves]; avgconstraints[1][3*i+3]=100.0*bright2[ncurves]; avgconstraints[1][3*i+4]=50.0*(number1[ncurves]+number2[ncurves]);
		}

		globalparams=new double[ncurves][11];
		globalvflmatrix=new int[ncurves][11];
		for(int i=0;i<ncurves;i++){
			globalparams[i][0]=0.0; globalparams[i][1]=0.0; globalparams[i][2]=bright1[i]; globalparams[i][3]=beta*bright1[i]; globalparams[i][4]=number1[i];
			globalvflmatrix[i][0]=1; globalvflmatrix[i][1]=1; globalvflmatrix[i][2]=2; globalvflmatrix[i][3]=2; globalvflmatrix[i][4]=0;
			for(int j=0;j<2;j++){
				globalparams[i][3*j+5]=0.1; globalvflmatrix[i][3*j+5]=1;
				globalparams[i][3*j+6]=0.1; globalvflmatrix[i][3*j+6]=1;
				globalparams[i][3*j+7]=0.0; globalvflmatrix[i][3*j+7]=1;
			}
			globalparams[i][5]=0.0; globalparams[i][6]=bright2[i]; globalparams[i][7]=number2[i];
			globalvflmatrix[i][6]=2; globalvflmatrix[i][7]=0;
		}
		globalvflmatrix[0][2]=1;
		globalvflmatrix[0][3]=1;
		globalvflmatrix[0][6]=1;

		globalconstraints=new double[2][ncurves][11];
		for(int j=0;j<ncurves;j++){
			globalconstraints[0][j][0]=0.0001; globalconstraints[1][j][0]=intensity1[j]; globalconstraints[0][j][1]=0.0001; globalconstraints[1][j][1]=intensity2[j];
			for(int i=0;i<3;i++){
				globalconstraints[0][j][3*i+2]=0.0001; globalconstraints[0][j][3*i+3]=0.0001; globalconstraints[0][j][3*i+4]=0.0;
				globalconstraints[1][j][3*i+2]=100.0*bright1[j]; globalconstraints[1][j][3*i+3]=100.0*bright2[j]; globalconstraints[1][j][3*i+4]=50.0*(number1[j]+number2[j]);
			}
		}

		for(int i=0;i<ncurves;i++){
			c2[i]=0.0;
			c2array[i].setText(""+(float)c2[i]);
		}

		globalformulas=new String[ncurves][11];
		String[] tempnames={"backgroundg","backgroundr","e1g","e1r","N1","e2g","e2r","N2","e3g","e3r","N3"};
		paramsnames=tempnames;

		undoparams=new double[ncurves][11];
		undovflmatrix=new int[ncurves][11];
		undoformulas=new String[ncurves][11];
		undoc2=new double[ncurves];
		undoglobalc2=0.0;
		undofit=new float[ncurves][xpts][ypts];
	}


	private void fitavg(){
		if(showfitdialog()){
			double[] stats=new double[2];
			float[] tempdata=new float[xpts*ypts];
			float[] tempweights=new float[xpts*ypts];
			for(int i=0;i<xpts;i++){
				for(int j=0;j<ypts;j++){
					tempdata[i+j*xpts]=(float)((double)avg[i][j]/(double)nmeas[ncurves]);
					tempweights[i+j*xpts]=avgweights[i][j];
				}
			}
			int tempmaxiter=fitclass.maxiter;
			if(checkc2){fitclass.maxiter=0;}
			float[] fit=fitclass.fitdata(avgparams,avgfixes,avgconstraints,tempdata,tempweights,stats,false);
			fitclass.maxiter=tempmaxiter;
			for(int i=0;i<xpts;i++){
				for(int j=0;j<ypts;j++){
					avgfit[i][j]=fit[i+j*xpts]*(float)nmeas[ncurves];
				}
			}
			pwavg.updateSeries(avgfit,1,false);
			float[] temp=pwavg.getLimits(); temp[4]=1.0f; pwavg.setLimits(temp);
			c2[ncurves]=stats[1];
			c2array[ncurves].setText(""+(float)c2[ncurves]);
		}
	}

	private void fitglobal(){
		int nparams=11;
		int nsel=0;
		for(int i=0;i<ncurves;i++){
			if(include[i]){nsel++;}
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
			for(int j=0;j<xpts;j++){
				for(int k=0;k<ypts;k++){
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
			float[][] tempdata=new float[nsel][xpts*ypts];
			float[][] tempweights=new float[nsel][xpts*ypts];
			counter=0;
			for(int i=0;i<ncurves;i++){
				if(include[i]){
					for(int j=0;j<xpts;j++){
						for(int k=0;k<ypts;k++){
							tempdata[counter][j+k*xpts]=(float)((double)pch[i][j][k]/(double)nmeas[i]);
							tempweights[counter][j+k*xpts]=weights[i][j][k];
						}
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
					for(int j=0;j<xpts;j++){
						for(int k=0;k<ypts;k++){
							fit[i][j][k]=tempfit[counter][j+xpts*k]*(float)nmeas[i];
						}
					}
					if(i==dispcurve){
						pwfit.updateSeries(fit[dispcurve],1,true);
					}
					for(int j=0;j<nparams;j++){
						globalparams[i][j]=params[counter][j];
					}
					c2[i]=tempc2vals[counter];
					c2array[i].setText(""+(float)c2[i]);
					counter++;
				}
			}
			float[] temp=pwfit.getLimits(); temp[4]=1.0f; pwfit.setLimits(temp);
		}
	}

	private boolean showfitdialog(){
		int nparams=avgparams.length;
		Object[][] tabledata=new Object[nparams+1][3];
		String[] columnlabels={"Parameters","Values","Fix?"};
		tabledata[0][0]="Check chi^2?"; tabledata[0][1]=new Boolean(checkc2);
		tabledata[1][0]="backgroundg"; tabledata[2][0]="backgroundr";
		tabledata[3][0]="e1g"; tabledata[4][0]="e1r"; tabledata[5][0]="N1"; tabledata[6][0]="e2g"; tabledata[7][0]="e2r"; tabledata[8][0]="N2";
		tabledata[9][0]="e3g"; tabledata[10][0]="e3r"; tabledata[11][0]="N3";
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
		tabledata[1][0]="backgroundg"; tabledata[2][0]="backgroundr";
		tabledata[3][0]="e1g"; tabledata[4][0]="e1r"; tabledata[5][0]="N1"; tabledata[6][0]="e2g"; tabledata[7][0]="e2r"; tabledata[8][0]="N2";
		tabledata[9][0]="e3g"; tabledata[10][0]="e3r"; tabledata[11][0]="N3";
		for(int i=0;i<nparams;i++){
			tabledata[i+1][1]=linking;
			for(int j=0;j<nsets;j++){
				tabledata[i+1][2*j+2]=new String(""+(float)params[j][i]);
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
		if(retvals==null){
			//IJ.showMessage("returned false");
			return false;
		}
		//IJ.showMessage("returned true");
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
			for(int j=0;j<xpts;j++){
				for(int k=0;k<ypts;k++){
					fit[i][j][k]=undofit[i][j][k];
				}

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
			int nparams=11;
			for(int i=0;i<ncurves;i++){
				if(include[i]){nsel++;}
			}
			double[][] params=new double[nsel][nparams];
			String[][] tempformulas=new String[nsel][nparams];
			double[][][] constraints=new double[2][nsel][nparams];
			int[][] vflmatrix=new int[nsel][nparams];

			float[][] tempdata=new float[nsel][xpts*ypts];
			float[][] tempweights=new float[nsel][xpts*ypts];

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
						if(vflmatrix[counter][j]==0 || (j==0 && vflmatrix[counter][j]==2)){nfit++;}
					}
					for(int j=0;j<xpts;j++){
						for(int k=0;k<ypts;k++){
							tempdata[counter][j+k*xpts]=(float)((double)pch[i][j][k]/(double)nmeas[i]);
							tempweights[counter][j+k*xpts]=weights[i][j][k];
						}
					}
					counter++;
				}
			}
			int dofnum=xpts*ypts*nsel-(nfit-1)-1;
			int dofden=xpts*ypts*nsel-nfit-1;
			//double flim=FLimit(dofnum,dofden,(double)conf);
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

			float[] tempdata=new float[xpts*ypts];
			float[] tempweights=new float[xpts*ypts];
			for(int i=0;i<xpts;i++){
				for(int j=0;j<ypts;j++){
					tempdata[i+j*xpts]=(float)((double)avg[i][j]/(double)nmeas[ncurves]);
					tempweights[i+j*xpts]=avgweights[i][j];
				}
			}

			int nfit=0;
			for(int i=0;i<7;i++){
				if(avgfixes[i]==0){
					nfit++;
				}
			}
			int dofnum=xpts*ypts-(nfit-1)-1;
			int dofden=xpts*ypts-nfit-1;
			double flim=(new jdist()).FLimit(dofnum,dofden,(double)conf);
			IJ.log("FLimit = "+(float)flim);
			if(flim==Double.NaN && flim < 1.0){IJ.showMessage("Invalid Limiting F Value"); return;}
			double truespacing=Math.abs(avgparams[errindex]*spacing);
			double[][] c2plot=erclass.geterrors(avgparams,avgfixes,avgconstraints,tempdata,tempweights,flim,truespacing,errindex);
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
		weights=new float[ncurves][xpts][ypts];
		for(int i=0;i<ncurves;i++){
			nmeas[i]=0;
			for(int j=0;j<xpts;j++){
				for(int k=0;k<ypts;k++){nmeas[i]+=(int)pch[i][j][k];}
			}
			double tempavg=0.0; double tempavg2=0.0; double temp2avg=0.0; double temp2avg2=0.0; double tempccavg=0.0;
			for(int j=0;j<xpts;j++){
				for(int k=0;k<ypts;k++){
					double normed=(double)pch[i][j][k]/(double)nmeas[i];
					if(pch[i][j][k]>0.0f){
						weights[i][j][k]=(float)((double)nmeas[i]/(normed*(1.0f-normed)));
					} else {
						weights[i][j][k]=1.0f;
					}
					tempavg+=normed*(double)j;
					tempavg2+=normed*(double)j*(double)j;
					temp2avg+=normed*(double)k;
					temp2avg2+=normed*(double)k*(double)k;
					tempccavg+=normed*(double)k*(double)j;
				}
			}
			tempccavg-=tempavg*temp2avg;
			brightcc[i]=tempccavg/Math.sqrt(tempavg*temp2avg);
			tempavg2-=tempavg*tempavg;
			tempavg2/=tempavg;
			bright1[i]=(tempavg2-1.0);
			temp2avg2-=temp2avg*temp2avg;
			temp2avg2/=temp2avg;
			bright2[i]=(temp2avg2-1.0);
			intensity1[i]=tempavg;
			intensity2[i]=temp2avg;
			if(psfflag==0){
				bright1[i]/=0.3536;
				bright2[i]/=0.3536;
				brightcc[i]/=0.3536;
			} else {
				if(psfflag==1){
					bright1[i]/=0.078;
					bright2[i]/=0.078;
					brightcc[i]/=0.078;
				} else {
					bright1[i]/=0.5;
					bright2[i]/=0.5;
					brightcc[i]/=0.5;
				}
			}
			number1[i]=intensity1[i]/bright1[i];
			number2[i]=intensity2[i]/bright2[i];
			brightmincc[i]=(bright1[i]*beta)*Math.sqrt(intensity1[i]/intensity2[i]);
		}
	}

	private void updateavg(){
		nmeas[ncurves]=0;
		avg=new float[xpts][ypts];
		avgweights=new float[xpts][ypts];
		for(int i=0;i<ncurves;i++){
			if(include[i]){
				for(int j=0;j<xpts;j++){
					for(int k=0;k<ypts;k++){
						avg[j][k]+=pch[i][j][k];
						nmeas[ncurves]+=(int)pch[i][j][k];
					}
				}
			}
		}
		double tempavg=0.0; double tempavg2=0.0; double temp2avg=0.0; double temp2avg2=0.0; double tempccavg=0.0;
		for(int i=0;i<xpts;i++){
			for(int j=0;j<ypts;j++){
				double normed=(double)avg[i][j]/(double)nmeas[ncurves];
				avgweights[i][j]=(float)((double)nmeas[ncurves]/(normed*(1.0f-normed)));
				if(avg[i][j]>0.0f){
					avgweights[i][j]=(float)((double)nmeas[ncurves]/(normed*(1.0f-normed)));
				} else {
					avgweights[i][j]=1.0f;
				}
				tempavg+=(double)i*normed;
				tempavg2+=(double)i*(double)i*normed;
				temp2avg+=(double)j*normed;
				temp2avg2+=(double)j*(double)j*normed;
				tempccavg+=(double)i*(double)j*normed;
			}
		}
		tempccavg-=tempavg*temp2avg;
		brightcc[ncurves]=tempccavg/Math.sqrt(tempavg*temp2avg);
		tempavg2-=tempavg*tempavg;
		tempavg2/=tempavg;
		bright1[ncurves]=(tempavg2-1.0);
		temp2avg2-=temp2avg*temp2avg;
		temp2avg2/=temp2avg;
		bright2[ncurves]=(temp2avg2-1.0);
		intensity1[ncurves]=tempavg;
		intensity2[ncurves]=temp2avg;
		if(psfflag==0){
			bright1[ncurves]/=0.3536;
			bright2[ncurves]/=0.3536;
			brightcc[ncurves]/=0.3536;
		} else {
			if(psfflag==1){
				bright1[ncurves]/=0.078;
				bright2[ncurves]/=0.078;
				brightcc[ncurves]/=0.078;
			} else {
				bright1[ncurves]/=0.5;
				bright2[ncurves]/=0.5;
				brightcc[ncurves]/=0.5;
			}
		}
		number1[ncurves]=intensity1[ncurves]/bright1[ncurves];
		number2[ncurves]=intensity2[ncurves]/bright2[ncurves];
		brightmincc[ncurves]=(bright1[ncurves]*beta)*Math.sqrt(intensity1[ncurves]/intensity2[ncurves]);
	}

	public void updatebeta(){
		for(int i=0;i<=ncurves;i++){
			brightmincc[i]=(bright1[i]*beta)/Math.sqrt(intensity1[i]/intensity2[i]);
			eminccarray[i].setText(""+(float)brightmincc[i]);
		}
	}

	public double fitfunc(double[] params,int indvar){
		double[] brightnesses1=new double[3];
		double[] brightnesses2=new double[3];
		double[] numbers=new double[3];
		double background1=params[0];
		double background2=params[1];
		double nmax=0.0;
		for(int i=0;i<3;i++){
			brightnesses1[i]=params[3*i+2];
			brightnesses2[i]=params[3*i+3];
			numbers[i]=params[3*i+4];
			if(numbers[i]>nmax){nmax=numbers[i];}
		}
		int nlength=20;
		if(nmax>15.0f){nlength=(int)(nmax*1.3);}
		int yvar=(int)(indvar/xpts);
		int xvar=indvar-yvar*xpts;
		return pchfunc.getpch(xvar,yvar,brightnesses1,brightnesses2,numbers,background1,background2,nlength);
	}

	public void showresults(String results){
		IJ.showStatus(results);
		IJ.log(results);
	}
			
}
