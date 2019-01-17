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
import java.io.*;
import jguis.*;
import jalgs.*;
import jalgs.jsim.*;
import java.lang.reflect.*;
import ij.measure.*;
import ij.io.*;

public class simfluc_scanning_jru_v2 implements PlugIn,simflucinterface {
	float w0,z0,boxsize,yflowrate,adcoff,readstdev,gain,exptime;
	float gridsize,jumpprob,dfraction;
	float[] D,bkd,kbleach;
	float[][] bright,trate;
	int[] num,bsteps;
	int frames,boxpixels,noiseindex,boundaryindex,confineindex,dimensionindex;
	boolean startcenter,restrict,bleached_reappear,transitions,tpe,batch,bleachbyphoton;
	int nspecies,nchannels,batchruns;

	public void run(String arg) {
		init_options();
		tpe=false;
		bleachbyphoton=false;
		if(!settings_dialog()){return;}
		if(!particle_dialog()){return;}
		if(transitions){
			if(!transitions_dialog()){return;}
		} else {
			trate=new float[nspecies][nspecies];
		}
		simfluc_scanning ss=new simfluc_scanning(this);
		if(tpe){
			ss.tpe=true;
		}
		simfluc_camera sc=new simfluc_camera(this);
		ss.bleachbyphoton=bleachbyphoton;
		sc.bleachbyphoton=bleachbyphoton;
		String batchdir=null;
		if(batch){
			DirectoryChooser dc=new DirectoryChooser("Choose Ouput Dir");
			batchdir=dc.getDirectory();
			if(batchdir==null) return;
		}
		for(int batchindex=0;batchindex<batchruns;batchindex++){
		if(dimensionindex==0){
			float[][] data=ss.do_simfluc_point(boxsize,boxpixels,w0,z0,frames,exptime,bkd,yflowrate,confineindex,startcenter,restrict,bleached_reappear,
				noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,adcoff,readstdev,num,bright,D,trate,kbleach,bsteps);

			float[] xvals=new float[frames];
			float xinc=exptime*(float)(1.0e-6);
			for(int i=0;i<frames;i++){
				xvals[i]=xinc*(float)i;
			}
			PlotWindow4 pw=new PlotWindow4("Point Simulation","time","Intensity",xvals,data[0]);
			pw.draw();
			for(int i=1;i<data.length;i++){
				pw.addPoints(xvals,data[i],true);
			}
			if(batch){
				pw.saveAsObject(batchdir+File.separator+"simtraj"+batchindex+".pw");
				pw.close();
			}
		}

		if(dimensionindex==1){
			float[][] data=ss.do_simfluc_line(boxsize,boxpixels,w0,z0,frames,exptime,bkd,yflowrate,confineindex,startcenter,restrict,bleached_reappear,
				noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,adcoff,readstdev,num,bright,D,trate,kbleach,bsteps);

			ImageStack stack=new ImageStack(boxpixels,frames);
			for(int i=0;i<data.length;i++){
				stack.addSlice("",data[i]);
			}
			ImagePlus imp=new ImagePlus("Raster Line Simulation",stack);
			Calibration cal=imp.getCalibration().copy();
			cal.setUnit("um");
			cal.pixelWidth=(double)boxsize/(double)boxpixels;
			cal.pixelHeight=cal.pixelWidth;
			imp.setCalibration(cal);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(data.length,1,1);
			if(data.length>1){
				new CompositeImage(imp,CompositeImage.COMPOSITE).show();
			} else {
				imp.show();
			}
			if(batch){
				new FileSaver(imp).saveAsTiff(batchdir+File.separator+"simcarpet"+batchindex+".tif");
				imp.close();
			}
		}

		if(dimensionindex==2){
			float[][][] data=ss.do_simfluc_image(boxsize,boxpixels,w0,z0,frames,exptime,bkd,yflowrate,confineindex,startcenter,restrict,bleached_reappear,
				noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,adcoff,readstdev,num,bright,D,trate,kbleach,bsteps);

			ImageStack stack=new ImageStack(boxpixels,boxpixels);
			for(int i=0;i<frames;i++){
				for(int j=0;j<data.length;j++){
					stack.addSlice("",data[j][i]);
				}
			}
			ImagePlus imp=new ImagePlus("Raster Image Simulation",stack);
			Calibration cal=imp.getCalibration().copy();
			cal.setUnit("um");
			cal.pixelWidth=(double)boxsize/(double)boxpixels;
			cal.pixelHeight=cal.pixelWidth;
			imp.setCalibration(cal);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(data.length,1,frames);
			if(data.length>1){
				new CompositeImage(imp,CompositeImage.COMPOSITE).show();
			} else {
				imp.show();
			}
			if(batch){
				new FileSaver(imp).saveAsTiff(batchdir+File.separator+"simimage"+batchindex+".tif");
				imp.close();
			}
		}

		if(dimensionindex==3){
			float[][] data=sc.do_simfluc_line(boxsize,boxpixels,w0,z0,frames,exptime,bkd,yflowrate,confineindex,startcenter,restrict,bleached_reappear,
				noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,adcoff,readstdev,num,bright,D,trate,kbleach,bsteps);

			ImageStack stack=new ImageStack(boxpixels,frames);
			for(int i=0;i<data.length;i++){
				stack.addSlice("",data[i]);
			}
			ImagePlus imp=new ImagePlus("Line Simulation",stack);
			Calibration cal=imp.getCalibration().copy();
			cal.setUnit("um");
			cal.pixelWidth=(double)boxsize/(double)boxpixels;
			cal.pixelHeight=cal.pixelWidth;
			imp.setCalibration(cal);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(data.length,1,1);
			if(data.length>1){
				new CompositeImage(imp,CompositeImage.COMPOSITE).show();
			} else {
				imp.show();
			}
			if(batch){
				new FileSaver(imp).saveAsTiff(batchdir+File.separator+"simcarpet"+batchindex+".tif");
				imp.close();
			}
		}

		if(dimensionindex==4){
			float[][][] data=sc.do_simfluc_lineimage(boxsize,boxpixels,w0,z0,frames,exptime,bkd,yflowrate,confineindex,startcenter,restrict,bleached_reappear,
				noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,adcoff,readstdev,num,bright,D,trate,kbleach,bsteps);

			ImageStack stack=new ImageStack(boxpixels,boxpixels);
			for(int i=0;i<frames;i++){
				for(int j=0;j<data.length;j++){
					stack.addSlice("",data[j][i]);
				}
			}
			ImagePlus imp=new ImagePlus("Line Scan Simulation",stack);
			Calibration cal=imp.getCalibration().copy();
			cal.setUnit("um");
			cal.pixelWidth=(double)boxsize/(double)boxpixels;
			cal.pixelHeight=cal.pixelWidth;
			imp.setCalibration(cal);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(data.length,1,frames);
			if(data.length>1){
				new CompositeImage(imp,CompositeImage.COMPOSITE).show();
			} else {
				imp.show();
			}
			if(batch){
				new FileSaver(imp).saveAsTiff(batchdir+File.separator+"simimage"+batchindex+".tif");
				imp.close();
			}
		}

		if(dimensionindex==5){
			float[][][] data=sc.do_simfluc_image(boxsize,boxpixels,w0,z0,frames,exptime,bkd,yflowrate,confineindex,startcenter,restrict,bleached_reappear,
				noiseindex,boundaryindex,gridsize,dfraction,jumpprob,gain,adcoff,readstdev,num,bright,D,trate,kbleach,bsteps);

			ImageStack stack=new ImageStack(boxpixels,boxpixels);
			for(int i=0;i<frames;i++){
				for(int j=0;j<data.length;j++){
					stack.addSlice("",data[j][i]);
				}
			}
			ImagePlus imp=new ImagePlus("Image Simulation",stack);
			Calibration cal=imp.getCalibration().copy();
			cal.setUnit("um");
			cal.pixelWidth=(double)boxsize/(double)boxpixels;
			cal.pixelHeight=cal.pixelWidth;
			imp.setCalibration(cal);
			imp.setOpenAsHyperStack(true);
			imp.setDimensions(data.length,1,frames);
			if(data.length>1){
				new CompositeImage(imp,CompositeImage.COMPOSITE).show();
			} else {
				imp.show();
			}
			if(batch){
				new FileSaver(imp).saveAsTiff(batchdir+File.separator+"simimage"+batchindex+".tif");
				imp.close();
			}
		}
		if(!batch) break;
		}
		set_options();
	}

	public void showprogress(int i,int end){
		IJ.showProgress(i,end);
	}

	public void showmessage(String message){
		IJ.log(message);
	}

	private boolean particle_dialog(){
		String[] columnlabels=new String[nchannels+4];
		columnlabels[0]="#"; columnlabels[nchannels+1]="D (um^2/s)"; columnlabels[nchannels+2]="kbleach (1/s or 1/p)"; columnlabels[nchannels+3]="Bleach Steps";
		for(int i=0;i<nchannels;i++){columnlabels[i+1]="Bright"+(i+1)+" (cpsm)";}
		Object[][] tabledata=new Object[nspecies][nchannels+4];
		for(int i=0;i<nspecies;i++){
			tabledata[i][0]=new Integer(num[i]);
			for(int j=0;j<nchannels;j++){
				tabledata[i][j+1]=new Float(bright[i][j]);
			}
			tabledata[i][nchannels+1]=new Float(D[i]); tabledata[i][nchannels+2]=new Float(kbleach[i]); tabledata[i][nchannels+3]=new Integer(bsteps[i]);
		}
		Object[][] retvals=TableDialog2.showDialog(null,null,"Particle Options",columnlabels,tabledata,null);
		if(retvals==null){return false;}
		for(int i=0;i<nspecies;i++){
			num[i]=((Double)retvals[i][0]).intValue();
			for(int j=0;j<nchannels;j++){
				bright[i][j]=((Double)retvals[i][j+1]).floatValue();
			}
			D[i]=((Double)retvals[i][nchannels+1]).floatValue(); kbleach[i]=((Double)retvals[i][nchannels+2]).floatValue(); bsteps[i]=((Double)retvals[i][nchannels+3]).intValue();
		}
		return true;
	}

	private boolean transitions_dialog(){
		String[] columnlabels=new String[nspecies+1];
		Object[][] tabledata=new Object[nspecies][nspecies+1];
		for(int i=0;i<nspecies;i++){
			columnlabels[i+1]="k"+(i+1)+" (1/s)";
			tabledata[i][0]="k"+(i+1)+" (1/s)";
			for(int j=0;j<nspecies;j++){
				if(i==j){
					tabledata[i][j+1]="NA";
				} else {
					tabledata[i][j+1]=new Float(trate[i][j]);
				}
			}
		}
		Object[][] retvals=jguis.TableDialog2.showDialog(null,null,"Transition Options",columnlabels,tabledata,null);
		if(retvals==null){return false;}
		for(int i=0;i<nspecies;i++){
			for(int j=0;j<nspecies;j++){
				if(i!=j){
					trate[i][j]=((Double)retvals[i][j+1]).floatValue();
				} else {
					trate[i][j]=0.0f;
				}
			}
		}
		return true;
	}

	private boolean settings_dialog(){
		GenericDialog gd=new GenericDialog("Settings");
		gd.addNumericField("Box size (um)",boxsize,5,10,null);
		gd.addNumericField("Box pixels",boxpixels,0);
		gd.addNumericField("w0 (um)",w0,5,10,null);
		gd.addNumericField("z0 (um)",z0,5,10,null);
		gd.addNumericField("# of frames",frames,0,10,null);
		gd.addNumericField("Exposure time (us)",exptime,5,10,null);
		for(int i=0;i<nchannels;i++){
			gd.addNumericField("Ch"+(i+1)+" Background (cps)",bkd[i],5,10,null);
		}
		gd.addNumericField("X flow rate (um/s)",yflowrate,5,10,null);
		String[] confine_options={"None","xy plane","xz plane","1D xz"};
		gd.addChoice("Plane Confinement",confine_options,confine_options[confineindex]);
		String[] measure_options={"Point","Raster Line","Raster Image","Line","Line Scan Image","Image"};
		gd.addChoice("Scanning Mode",measure_options,measure_options[dimensionindex]);
		gd.addCheckbox("TPE Excitation",false);
		gd.addCheckbox("Allow Transitions b/tw Particles?",transitions);
		gd.addCheckbox("Start with molecules in center?",startcenter);
		gd.addCheckbox("Grid of restricted diffusion?",restrict);
		gd.addCheckbox("Bleached Molecules Reappear?",bleached_reappear);
		gd.addCheckbox("Photon_Budget_Bleaching",bleachbyphoton);
		String[] noiseoptions={"Poisson","Analog","None"};
		gd.addChoice("Noise Options",noiseoptions,noiseoptions[noiseindex]);
		String[] boundaryoptions={"Periodic","Reflective","New Particle"};
		gd.addChoice("Boundary Conditions",boundaryoptions,boundaryoptions[boundaryindex]);
		gd.addCheckbox("Batch Mode",false);
		gd.addNumericField("#_of_Batches",4,0);
		gd.showDialog(); if(gd.wasCanceled()){return false;}
		boxsize=(float)gd.getNextNumber();
		boxpixels=(int)gd.getNextNumber();
		w0=(float)gd.getNextNumber();
		z0=(float)gd.getNextNumber();
		frames=(int)gd.getNextNumber();
		exptime=(float)gd.getNextNumber();
		for(int i=0;i<nchannels;i++){
			bkd[i]=(float)gd.getNextNumber();
		}
		yflowrate=(float)gd.getNextNumber();
		confineindex=gd.getNextChoiceIndex();
		dimensionindex=gd.getNextChoiceIndex();
		tpe=gd.getNextBoolean();
		transitions=gd.getNextBoolean();
		startcenter=gd.getNextBoolean();
		restrict=gd.getNextBoolean();
		bleached_reappear=gd.getNextBoolean();
		bleachbyphoton=gd.getNextBoolean();
		noiseindex=gd.getNextChoiceIndex();
		boundaryindex=gd.getNextChoiceIndex();
		batch=gd.getNextBoolean();
		batchruns=(int)gd.getNextNumber();
		if(noiseindex==1){
			GenericDialog gd2=new GenericDialog("Noise Settings");
			gd2.addNumericField("Offset",adcoff,5,10,null);
			gd2.addNumericField("Read St Dev",readstdev,5,10,null);
			gd2.addNumericField("Gain",gain,5,10,null);
			gd2.showDialog(); if(gd2.wasCanceled()){return false;}
			adcoff=(float)gd2.getNextNumber();
			readstdev=(float)gd2.getNextNumber();
			gain=(float)gd2.getNextNumber();
		}
		if(restrict){
			GenericDialog gd2=new GenericDialog("Restriction Settings");
			gd2.addNumericField("Grid Size (um)",gridsize,5,10,null);
			gd2.addNumericField("Jump Probability",jumpprob,5,10,null);
			gd2.addNumericField("D Multiplier",dfraction,5,10,null);
			gd2.showDialog(); if(gd2.wasCanceled()){return false;}
			gridsize=(float)gd2.getNextNumber();
			jumpprob=(float)gd2.getNextNumber();
			dfraction=(float)gd2.getNextNumber();
		}
		return true;
	}

	void init_options(){
		String dir=System.getProperty("user.home");
		try{
			File f=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"simfluc_scanning_jru_v2.jrn");
			BufferedReader b=new BufferedReader(new FileReader(f));
			nspecies=(int)getstringnum(b.readLine());
			nchannels=(int)getstringnum(b.readLine());

			D=new float[nspecies];
			bright=new float[nspecies][nchannels];
			trate=new float[nspecies][nspecies];
			bkd=new float[nchannels];
			num=new int[nspecies];
			kbleach=new float[nspecies];
			bsteps=new int[nspecies];

			for(int i=0;i<nspecies;i++){
				num[i]=(int)getstringnum(b.readLine());
				D[i]=getstringnum(b.readLine());
				for(int j=0;j<nchannels;j++){
					bright[i][j]=getstringnum(b.readLine());
				}
				for(int j=0;j<nspecies;j++){
					trate[i][j]=getstringnum(b.readLine());
				}
				kbleach[i]=getstringnum(b.readLine());
				bsteps[i]=(int)getstringnum(b.readLine());
			}

			boxsize=getstringnum(b.readLine());
			boxpixels=(int)getstringnum(b.readLine());
			w0=getstringnum(b.readLine());
			z0=getstringnum(b.readLine());
			frames=(int)getstringnum(b.readLine());
			exptime=getstringnum(b.readLine());
			for(int i=0;i<nchannels;i++){
				bkd[i]=getstringnum(b.readLine());
			}
			yflowrate=getstringnum(b.readLine());
			confineindex=(int)getstringnum(b.readLine());
			dimensionindex=(int)getstringnum(b.readLine());
			transitions= ((int)getstringnum(b.readLine())==1) ? true : false;
			startcenter= ((int)getstringnum(b.readLine())==1) ? true : false;
			noiseindex=(int)getstringnum(b.readLine());
			boundaryindex=(int)getstringnum(b.readLine());
			adcoff=getstringnum(b.readLine());
			readstdev=getstringnum(b.readLine());
			gain=getstringnum(b.readLine());
			restrict= ((int)getstringnum(b.readLine())==1) ? true : false;
			bleached_reappear= ((int)getstringnum(b.readLine())==1) ? true : false;
			gridsize=getstringnum(b.readLine());
			jumpprob=getstringnum(b.readLine());
			dfraction=getstringnum(b.readLine());
			b.close();
		}
		catch(IOException e){
			nspecies=6;
			nchannels=2;

			D=new float[nspecies];
			bright=new float[nspecies][nchannels];
			trate=new float[nspecies][nspecies];
			bkd=new float[nchannels];
			num=new int[nspecies];
			kbleach=new float[nspecies];
			bsteps=new int[nspecies];

			num[0]=100;
			D[0]=40.0f;
			bright[0][0]=10000.0f;
			bsteps[0]=1;
			boxsize=3.2f;
			boxpixels=64;
			w0=0.17f;
			z0=0.7f;
			frames=135000;
			exptime=10.0f;
			yflowrate=0.0f;
			confineindex=0;
			dimensionindex=0;
			transitions=false;
			startcenter=false;
			noiseindex=0;
			boundaryindex=0;
			adcoff=1000.0f;
			readstdev=25.0f;
			gain=100.0f;
			restrict=false;
			bleached_reappear=false;
			gridsize=0.5f;
			jumpprob=1.0f;
			dfraction=1.0f;
		}
		return;
	}
	
	void set_options(){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){a.mkdir();}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"simfluc_scanning_jru_v2.jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			d.write(""+nspecies+" number_of_species\n");
			d.write(""+nchannels+" number_of_channels\n");
			for(int i=0;i<nspecies;i++){
				d.write(""+num[i]+" number_of_particle"+(i+1)+"\n");
				d.write(""+D[i]+" D"+(i+1)+"(um^2/s)\n");
				for(int j=0;j<nchannels;j++){
					d.write(""+bright[i][j]+" ch"+(j+1)+"_brightness"+(i+1)+"(cpsm)\n");
				}
				for(int j=0;j<nspecies;j++){
					d.write(""+trate[i][j]+" "+(i+1)+"_to_"+(j+1)+"_transition_rate(1/s)\n");
				}
				d.write(""+kbleach[i]+" bleach_rate"+(i+1)+"(1/s)\n");
				d.write(""+bsteps[i]+" bleach_steps"+(i+1)+"\n");
			}
			d.write(""+boxsize+" Box_size(um)\n");
			d.write(""+boxpixels+" Image_size(pixels)\n");
			d.write(""+w0+" Beam_radial_waist(um)\n");
			d.write(""+z0+" Beam_axial_waist(um)\n");
			d.write(""+frames+" #_of_frames\n");
			d.write(""+exptime+" exposure_time(us)\n");
			for(int i=0;i<nchannels;i++){
				d.write(""+bkd[i]+" Ch"+(i+1)+"_background(cps)\n");
			}
			d.write(""+yflowrate+" y_flow_rate(um/s)\n");
			d.write(""+confineindex+" confinement_options\n");
			d.write(""+dimensionindex+" scan_options\n");
			d.write(""+(transitions ? 1 : 0)+" allow_transitions?\n");
			d.write(""+(startcenter ? 1 : 0)+" initialize_molecules_at_center?\n");
			d.write(""+noiseindex+" noise_type(0_poisson,1_analog,2_none)\n");
			d.write(""+boundaryindex+" boundary_condition(1_periodic,2_reflective,3_newparticle)\n");
			d.write(""+adcoff+" ADC_offset\n");
			d.write(""+readstdev+" Read_noise_stdev\n");
			d.write(""+gain+" Analog_gain\n");
			d.write(""+(restrict ? 1 : 0)+" rectricted_grid_diffusion?\n");
			d.write(""+(bleached_reappear ? 1 : 0)+" bleached_molecules_reappear?\n");
			d.write(""+gridsize+" grid_size(um)\n");
			d.write(""+jumpprob+" restricted_jump_prob\n");
			d.write(""+dfraction+" restricted_D_fraction\n");
			d.close();
		}
		catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}

	private float getstringnum(String s){
		int whiteindex=s.indexOf(" ");
		if(whiteindex==-1){whiteindex=s.indexOf("\t");}
		return Float.parseFloat(s.substring(0,whiteindex));
	}

}
