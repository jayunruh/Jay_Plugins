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
import jalgs.jfit.*;

public class linear_unmixing_jru_v1 implements PlugIn {
	int species,mainoptionsindex,totspecies,maxlength;
	boolean subback;
	String[] speciesnames;
	float[][] speciesspectra;
	float[] backgroundspectrum;
	int[] speciesindices;

	public void run(String arg) {
		//this plugin performs linear unmixing and in so doing finds the intensity of multiple components in frequency resolved image
		//the image could be frequency resolved in color or in temporal fourier frequency (lifetime or photoactivation)
		//a set of basis spectra as well as a background basis spectrum can be "read in" from plots
		//the default values are stored in a journal file in "ImageJ_defaults" in the users home directory
		init_options();
		GenericDialog gd=new GenericDialog("Options");
		String[] mainoptions={"Read_Reference_Spectrum","Unmix","Update_Reference_Spectrum","Update_Background_Spectrum","Delete_Reference_Spectrum","Delete_all_Reference_Spectra","Show_Reference_Spectrum","Show_Background_Spectrum","Carpet_Unmix"};
		gd.addChoice("What_would_you_like_to_do?",mainoptions,mainoptions[mainoptionsindex]);
		gd.addCheckbox("Use_Plot_Spectra",false);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		mainoptionsindex=gd.getNextChoiceIndex();
		boolean plot=gd.getNextBoolean();
		if(mainoptionsindex==1){
			float[][] tempspectra2=speciesspectra;
			int maxspecies=totspecies;
			String[] tempnames=speciesnames;
			int[] tempindices=speciesindices;
			if(plot){
				String[] plotlabels={"Spectra Plot"};
				ImageWindow[] iw=jutils.selectPlots(false,1,plotlabels);
				tempspectra2=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
				maxspecies=tempspectra2.length;
				species=maxspecies;
				tempnames=new String[maxspecies];
				tempindices=new int[maxspecies];
				for(int i=0;i<maxspecies;i++){tempindices[i]=i; tempnames[i]="series"+i;}
				speciesindices=tempindices;
				if(species>maxspecies) species=maxspecies;
				maxlength=tempspectra2[0].length;
				//pad the background if necessary to achieve the correct length
				if(backgroundspectrum==null) backgroundspectrum=new float[maxlength];
				if(backgroundspectrum.length<maxlength) padlength();
			}
			GenericDialog gd2=new GenericDialog("More Options");
			gd2.addNumericField("How_Many_Species_are_Present",species,0);
			gd2.addCheckbox("Subtract_Background?",subback);
			gd2.addCheckbox("Output_Residuals?",false);
			gd2.addCheckbox("Output_chi^2?",false);
			gd2.addCheckbox("Truncate_Negative_Values",true);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			species=(int)gd2.getNextNumber();
			subback=gd2.getNextBoolean();
			boolean outres=gd2.getNextBoolean();
			boolean outc2=gd2.getNextBoolean();
			boolean truncate=gd2.getNextBoolean();
			if(species>maxspecies || species==0){IJ.showMessage("Not enough reference spectra"); return;}
			GenericDialog gd3=new GenericDialog("Pick Reference Spectra");
			for(int i=0;i<species;i++){gd3.addChoice("Species_"+(i+1),tempnames,tempnames[tempindices[i]]);}
			gd3.showDialog(); if(gd3.wasCanceled()){return;}
			for(int i=0;i<species;i++){tempindices[i]=gd3.getNextChoiceIndex();}
			if(!plot) for(int i=0;i<species;i++){speciesindices[i]=tempindices[i];}
			//for(int i=0;i<species;i++){speciesindices[i]=tempindices[i];}
			ImagePlus imp=WindowManager.getCurrentImage();
			ImageStack stack=imp.getStack();
			int nch=imp.getNChannels();
			int slices=imp.getNSlices();
			if(nch>maxlength){nch=maxlength;}
			int width=imp.getWidth();
			int height=imp.getHeight();
			ImageStack resstack=new ImageStack(width,height);
			ImageStack c2stack=new ImageStack(width,height);
			ImageStack residstack=new ImageStack(width,height);
			for(int i=0;i<slices;i++){
				for(int j=0;j<species;j++){
					resstack.addSlice(tempnames[tempindices[j]]+i,new float[width*height]);
				}
				if(outc2){
					c2stack.addSlice("",new float[width*height]);
				}
				if(outres){
					for(int j=0;j<nch;j++){
						residstack.addSlice("",new float[width*height]);
					}
				}
			}
			float[][] tempspectra=new float[species][];
			for(int i=0;i<species;i++){
				tempspectra[i]=tempspectra2[speciesindices[i]];
			}
			linleastsquares lls=new linleastsquares(tempspectra,false,0,nch-1);
			boolean isfloat=(stack.getPixels(1) instanceof float[]);
			for(int j=0;j<slices;j++){
				for(int i=0;i<width*height;i++){
					float[] spectrum=new float[nch];
					if(isfloat){
						if(subback){
							for(int k=0;k<nch;k++){spectrum[k]=((float[])stack.getPixels(j*nch+k+1))[i]-backgroundspectrum[k];}
						} else {
							for(int k=0;k<nch;k++){spectrum[k]=((float[])stack.getPixels(j*nch+k+1))[i];}
						}
					} else {
						if(subback){
							for(int k=0;k<nch;k++){spectrum[k]=(float)(((short[])stack.getPixels(j*nch+k+1))[i]&0xffff)-backgroundspectrum[k];}
						} else {
							for(int k=0;k<nch;k++){spectrum[k]=(float)(((short[])stack.getPixels(j*nch+k+1))[i]&0xffff);}
						}
					}
					double[] contributions=lls.fitdata(spectrum,null);
					if(truncate){for(int k=0;k<species;k++){if(contributions[k]<0.0) contributions[k]=0.0;}}
					for(int k=0;k<species;k++){((float[])resstack.getPixels(j*species+k+1))[i]=(float)contributions[k];}
					if(outres){
						float[] res=lls.get_fresid(contributions,spectrum,null);
						for(int k=0;k<nch;k++){((float[])residstack.getPixels(j*nch+k+1))[i]=(float)res[k];}
					}
					if(outc2){
						((float[])c2stack.getPixels(j+1))[i]=(float)lls.get_c2(contributions,spectrum,null);
					}
				}
				IJ.showProgress(j,slices);
			}
			ImagePlus imp5=new ImagePlus("Unmixed Stack",resstack);
			imp5.setDimensions(species,slices,1);
			imp5.setOpenAsHyperStack(true);
			imp5.copyScale(imp);
			if(species>1){
				new CompositeImage(imp5,CompositeImage.COLOR).show();
			} else {
				imp5.show();
			}
			if(outc2){
				new ImagePlus("chi^2",c2stack).show();
			}
			if(outres){
				ImagePlus imp6=new ImagePlus("Residuals Stack",residstack);
				imp6.setDimensions(nch,slices,1);
				imp6.setOpenAsHyperStack(true);
				imp6.show();
			}
		}
		if(mainoptionsindex==0){
			GenericDialog gd2=new GenericDialog("Options");
			String tempname="";
			gd2.addStringField("Species_Name?",tempname);
			gd2.addCheckbox("Normalize_Spectrum?",true);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			tempname=gd2.getNextString();
			boolean norm=gd2.getNextBoolean();
			ImagePlus imp=WindowManager.getCurrentImage();
			ImageWindow iw=imp.getWindow();
			int selected=((Integer)jutils.runPW4VoidMethod(iw,"getSelected")).intValue();
			if(selected<0) selected=0;
			float[] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[selected];
			float[][] tempobject=new float[totspecies+1][];
			if(speciesspectra!=null){System.arraycopy(speciesspectra,0,tempobject,0,totspecies);}
			tempobject[totspecies]=new float[yvals.length];
			if(norm){
				double sum=0.0;
				for(int i=0;i<yvals.length;i++){
					sum+=(float)yvals[i];
				}
				for(int i=0;i<yvals.length;i++){
					tempobject[totspecies][i]=yvals[i]/(float)sum;
				}
			} else {
				for(int i=0;i<yvals.length;i++){
					tempobject[totspecies][i]=yvals[i];
				}
			}
			speciesspectra=tempobject;
			tempobject=null;
			String[] tempstring=new String[totspecies+1];
			if(speciesnames!=null && speciesnames.length>0){System.arraycopy(speciesnames,0,tempstring,0,totspecies);}
			tempstring[totspecies]=tempname;
			speciesnames=tempstring;
			tempstring=null;
			/*int[] tempint=new int[totspecies+1];
			if(speciesindices!=null && speciesindices.length>0){System.arraycopy(speciesindices,0,tempint,0,totspecies);}
			speciesindices=tempint;
			tempint=null;*/
			totspecies++;
			if(yvals.length>maxlength){maxlength=yvals.length; padlength();}
			if(yvals.length<maxlength){padlength();}
		}
		if(mainoptionsindex==2){
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Species_to_edit?",speciesnames,speciesnames[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int tempindex=gd2.getNextChoiceIndex();
			ImagePlus imp=WindowManager.getCurrentImage();
			ImageWindow iw=imp.getWindow();
			int selected=((Integer)jutils.runPW4VoidMethod(iw,"getSelected")).intValue();
			if(selected<0) selected=0;
			float[] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[selected];
			speciesspectra[tempindex]=new float[yvals.length];
			for(int i=0;i<yvals.length;i++){
				speciesspectra[tempindex][i]=yvals[i];
			}
			if(yvals.length>maxlength){maxlength=yvals.length; padlength();}
			if(yvals.length<maxlength){padlength();}
		}
		if(mainoptionsindex==3){
			ImagePlus imp=WindowManager.getCurrentImage();
			ImageWindow iw=imp.getWindow();
			float[] yvals=((float[][])jutils.runPW4VoidMethod(iw,"getYValues"))[0];
			backgroundspectrum=new float[yvals.length];
			for(int i=0;i<yvals.length;i++){
				backgroundspectrum[i]=yvals[i];
			}
		}
		if(mainoptionsindex==4){
			if(speciesspectra==null){return;}
			GenericDialog gd2=new GenericDialog("Options");
			gd2.addChoice("Species_to_delete?",speciesnames,speciesnames[0]);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int tempindex=gd2.getNextChoiceIndex();
			float[][] tempobject=new float[totspecies-1][];
			String[] tempstring=new String[totspecies-1];
			int[] tempint=new int[totspecies-1];
			for(int i=0;i<tempindex;i++){
				tempobject[i]=speciesspectra[i];
				tempstring[i]=speciesnames[i];
				tempint[i]=speciesindices[i];
			}
			for(int i=(tempindex+1);i<totspecies;i++){
				tempobject[i-1]=speciesspectra[i];
				tempstring[i-1]=speciesnames[i];
				tempint[i-1]=speciesindices[i];
			}
			speciesspectra=tempobject;
			tempobject=null;
			speciesnames=tempstring;
			tempstring=null;
			speciesindices=tempint;
			tempint=null;
			totspecies--;
		}
		if(mainoptionsindex==5){
			totspecies=0;
			mainoptionsindex=0;
			species=0;
			maxlength=0;
			subback=false;
		}
		if(mainoptionsindex==6){
			GenericDialog gd3=new GenericDialog("Options");
			gd3.addNumericField("Number of Spectra to Show?",1,0);
			gd3.showDialog(); if(gd3.wasCanceled()){return;}
			int nshow=(int)gd3.getNextNumber();
			GenericDialog gd2=new GenericDialog("Options");
			for(int i=0;i<nshow;i++){
				gd2.addChoice("Species_"+(i+1)+"_to_show?",speciesnames,speciesnames[0]);
			}
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			int[] indices=new int[nshow];
			for(int i=0;i<nshow;i++){
				indices[i]=gd2.getNextChoiceIndex();
			}
			String title=speciesnames[indices[0]]+" Spectrum";
			if(nshow>1){
				title="Unmixing Spectra";
			}
			PlotWindow4 pw=new PlotWindow4(title,"Spectral Unit","Intensity",speciesspectra[indices[0]]);
			pw.draw();
			for(int i=1;i<nshow;i++){
				pw.addPoints(speciesspectra[indices[i]],true);
			}
		}
		if(mainoptionsindex==7){
			new PlotWindow4("Background Spectrum","Spectral Unit","Intensity",backgroundspectrum).draw();
		}
		if(mainoptionsindex==8){
			float[][] tempspectra2=speciesspectra;
			int maxspecies=totspecies;
			String[] tempnames=speciesnames;
			int[] tempindices=speciesindices;
			if(plot){
				String[] plotlabels={"Spectra Plot"};
				ImageWindow[] iw=jutils.selectPlots(false,1,plotlabels);
				tempspectra2=(float[][])jutils.runPW4VoidMethod(iw[0],"getYValues");
				maxspecies=tempspectra2.length;
				tempnames=new String[maxspecies];
				for(int i=0;i<maxspecies;i++){tempindices[i]=i; tempnames[i]="series"+i;}
				if(species>maxspecies) species=maxspecies;
			}
			GenericDialog gd2=new GenericDialog("More Options");
			gd2.addNumericField("How_Many_Species_are_Present",species,0);
			gd2.addCheckbox("Subtract_Background?",subback);
			gd2.addCheckbox("Output_Residuals?",false);
			gd2.addCheckbox("Output_chi^2?",false);
			gd2.addCheckbox("Truncate_Negative_Values",true);
			gd2.showDialog(); if(gd2.wasCanceled()){return;}
			species=(int)gd2.getNextNumber();
			subback=gd2.getNextBoolean();
			boolean outres=gd2.getNextBoolean();
			boolean outc2=gd2.getNextBoolean();
			boolean truncate=gd2.getNextBoolean();
			if(species>maxspecies || species==0){IJ.showMessage("Not enough reference spectra"); return;}
			GenericDialog gd3=new GenericDialog("Pick Reference Spectra");
			for(int i=0;i<species;i++){gd3.addChoice("Species_"+(i+1),tempnames,tempnames[tempindices[i]]);}
			gd3.showDialog(); if(gd3.wasCanceled()){return;}
			for(int i=0;i<species;i++){tempindices[i]=gd3.getNextChoiceIndex();}
			if(!plot) for(int i=0;i<species;i++){speciesindices[i]=tempindices[i];}
			ImagePlus imp=WindowManager.getCurrentImage();
			int nch=imp.getWidth();
			int width=nch;
			int height=imp.getHeight();
			if(nch>maxlength){nch=maxlength;}
			float[] carpet=(float[])imp.getProcessor().convertToFloat().getPixels();
			float[][] resstack=new float[species][height];
			float[] c2stack=new float[height];
			float[] residstack=new float[nch*height];
			float[][] tempspectra=new float[species][];
			for(int i=0;i<species;i++){
				tempspectra[i]=tempspectra2[tempindices[i]];
			}
			linleastsquares lls=new linleastsquares(tempspectra,false,0,nch-1);
			for(int i=0;i<height;i++){
				float[] spectrum=new float[nch];
				if(subback){
					for(int k=0;k<nch;k++){spectrum[k]=carpet[i*width+k]-backgroundspectrum[k];}
				} else {
					for(int k=0;k<nch;k++){spectrum[k]=carpet[i*width+k];}
				}
				double[] contributions=lls.fitdata(spectrum,null);
				if(truncate){for(int k=0;k<species;k++){if(contributions[k]<0.0) contributions[k]=0.0;}}
				for(int k=0;k<species;k++){resstack[k][i]=(float)contributions[k];}
				if(outres){
					float[] res=lls.get_fresid(contributions,spectrum,null);
					for(int k=0;k<nch;k++){residstack[k+i*nch]=(float)res[k];}
				}
				if(outc2){
					c2stack[i]=(float)lls.get_c2(contributions,spectrum,null);
				}
			}
			new PlotWindow4("Unmixed Profiles","Channel","Intensity",resstack,null).draw();
			if(outc2){
				new PlotWindow4("c2 Profiles","Channel","c2",c2stack).draw();
			}
			if(outres){
				ImagePlus imp6=new ImagePlus("Residuals Carpet",new FloatProcessor(nch,height,residstack,null));
			}
		}
		set_options();
	}

	private void padlength(){
		for(int i=0;i<totspecies;i++){
			float[] tempfloat=speciesspectra[i];
			if(tempfloat.length<maxlength){
				float[] tempfloat2=new float[maxlength];
				for(int j=0;j<tempfloat.length;j++){
					tempfloat2[i]=tempfloat[i];
				}
				tempfloat=tempfloat2;
			}
		}
		float[] tempfloat2=new float[maxlength];
		if(backgroundspectrum!=null){
			System.arraycopy(backgroundspectrum,0,tempfloat2,0,backgroundspectrum.length);
			backgroundspectrum=tempfloat2;
		} else {
			backgroundspectrum=new float[maxlength];
		}
	}

	private void init_options(){
		String dir=System.getProperty("user.home");
		try{
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"linear_unmixing_jru_v1.jrn");
			BufferedReader d=new BufferedReader(new FileReader(b));
			String temp1=d.readLine();
			if(temp1==null){
				totspecies=0;
				mainoptionsindex=0;
				species=0;
				maxlength=0;
				subback=false;
				return;
			}
			totspecies=Integer.parseInt(temp1);
			mainoptionsindex=Integer.parseInt(d.readLine());
			species=Integer.parseInt(d.readLine());
			maxlength=Integer.parseInt(d.readLine());
			subback = (Integer.parseInt(d.readLine())==0) ? false : true;
			speciesindices=new int[species];
			for(int i=0;i<species;i++){
				speciesindices[i]=Integer.parseInt(d.readLine());
			}
			speciesnames=new String[totspecies];
			speciesspectra=new float[totspecies][];
			for(int i=0;i<totspecies;i++){
				speciesnames[i]=d.readLine();	
				speciesspectra[i]=new float[maxlength];
				for(int j=0;j<maxlength;j++){
					speciesspectra[i][j]=Float.parseFloat(d.readLine());
				}
			}
			backgroundspectrum=new float[maxlength];
			for(int i=0;i<maxlength;i++){
				backgroundspectrum[i]=Float.parseFloat(d.readLine());
			}
			d.close();
		}
		catch(IOException e){
			totspecies=0;
			mainoptionsindex=0;
			species=0;
			maxlength=0;
			subback=false;
		}
		return;
	}
	
	private void set_options(){
		String dir=System.getProperty("user.home");
		try{
			File a=new File(dir+File.separator+"ImageJ_defaults");
			if(!a.exists()){a.mkdir();}
			File b=new File(dir+File.separator+"ImageJ_defaults"+File.separator+"linear_unmixing_jru_v1.jrn");
			BufferedWriter d=new BufferedWriter(new FileWriter(b));
			d.write(""+totspecies+"\n");
			d.write(""+mainoptionsindex+"\n");
			d.write(""+species+"\n");
			d.write(""+maxlength+"\n");
			d.write(""+(subback ? 1:0)+"\n");
			for(int i=0;i<species;i++){
				d.write(""+speciesindices[i]+"\n");
			}
			for(int i=0;i<totspecies;i++){
				d.write(speciesnames[i]+"\n");
				float[] tempspec=(float[])speciesspectra[i];
				for(int j=0;j<tempspec.length;j++){
					d.write(""+tempspec[j]+"\n");
				}
				for(int j=tempspec.length;j<maxlength;j++){
					d.write("0.0\n");
				}
			}
			if(backgroundspectrum==null){backgroundspectrum=new float[maxlength];}
			for(int i=0;i<maxlength;i++){
				d.write(""+backgroundspectrum[i]+"\n");
			}
			d.close();
		}
		catch(IOException e){
			IJ.showMessage("error writing file");
			return;
		}
		return;
	}
	

}
