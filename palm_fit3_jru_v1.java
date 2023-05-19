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
import java.awt.event.*;
import java.awt.image.*;
import ij.plugin.*;
import jalgs.*;
import jguis.*;
import jalgs.jseg.*;
import jalgs.jfit.*;

public class palm_fit3_jru_v1 implements PlugIn{
	double stdev;
	int fitsize,fitsized2;
	boolean[] fitmask;
	int maxiterations=10;

	public void run(String arg) {
		//this plugin finds peaks corresponding to single molecules within an image
		//it then fits them to gaussians and outputs an image with pixel values corresponding
		//to molecule intensity
		//this version uses an expidited version of nonlinear least squares for analog intensity fitting

		ImagePlus imp = WindowManager.getCurrentImage();
		if(imp.getStack().getSize()<2){
			final PALMWindow cw = new PALMWindow();
			cw.init((float[])imp.getProcessor().getPixels());
			final  Frame f = new Frame("Interactive PALM Plot");
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
			cw.totalSize.height = PALMWindow.H + ins.bottom + ins.top + 65;
			cw.totalSize.width  = PALMWindow.WR + ins.left + ins.right;
			f.setSize(cw.totalSize);
			f.setVisible(true);
			cw.requestFocus();
			return;
		}
		int width = imp.getWidth();
		int height = imp.getHeight();

		boolean morefids=false;
		Rectangle[] fidrect=new Rectangle[10];
		int numfids=-1;
		int totfidsize=0;
		boolean usefids=true;
		do{
			(new WaitForUserDialog("Select Fiducial with Rectangle and Click OK")).show();
			numfids++;
			fidrect[numfids]=(imp.getProcessor()).getRoi();
			if(fidrect[numfids].width==width && fidrect[numfids].height==height){
				IJ.showStatus("No Fiducials Selected");
				morefids=false;
				usefids=false;
				numfids--;
				break;
			}
			totfidsize+=(fidrect[numfids].width*fidrect[numfids].height);
			GenericDialog gd100=new GenericDialog("Options");
			gd100.addCheckbox("Add Another Fiducial?",morefids);
			gd100.showDialog(); 
			if(gd100.wasCanceled()){morefids=false;}
			else{morefids=gd100.getNextBoolean();}
		}while(morefids);
		numfids++;
		
		float avg,max;
		double sum=0.0;
		ImageProcessor ip = imp.getProcessor();
		float[] pixels=(float[])(ip.convertToFloat()).getPixels();

		max=pixels[0];
		int[] fidx=new int[numfids];
		int[] fidy=new int[numfids];
		float[] fidmax=new float[numfids];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				boolean infid=false;
				for(int k=0;k<numfids;k++){
					if(fidrect[k].contains((float)j,(float)i)){
						infid=true;
						if(pixels[j+i*width]>fidmax[k]){
							fidmax[k]=pixels[j+i*width];
							fidx[k]=j;
							fidy[k]=i;
						}
					}
				}
				if(!infid){
					sum+=(double)pixels[j+i*width];
					if(pixels[j+i*width]>max){max=pixels[j+i*width];}
				}
			}
		}
		avg=(float)(sum/(double)(width*height-totfidsize));

		//IJ.showMessage(""+fidmax+" , "+fidx+" , "+fidy);
		
		GenericDialog gd = new GenericDialog("Options");
		float edge_buffer=10.0f;
		gd.addNumericField("Edge Buffer (pixels)",edge_buffer,0);
		float threshhold=avg+(max-avg)/4.0f;
		gd.addNumericField("Threshhold",threshhold,2);
		float baseline=avg;
		gd.addNumericField("Baseline",baseline,3);
		int minarea=4;
		gd.addNumericField("Min Pixels in Blob",minarea,0);
		int maxarea=1000;
		gd.addNumericField("Max Pixels in Blob",maxarea,0);
		float min_int=10.0f;
		gd.addNumericField("Minimum Intensity",min_int,3);
		int maxblobs=1000;
		gd.addNumericField("Maximum # of Blobs per Image",maxblobs,0);
		stdev=1.62;
		gd.addNumericField("PSF St Dev",stdev,10,15,"");
		float blob_separation=6.0f*(float)stdev;
		gd.addNumericField("Minimum Blob Separation",blob_separation,0);
		int maxon=50;
		gd.addNumericField("Maximum Frames On",maxon,0);
		float toler=4.0f;
		gd.addNumericField("Pixel tolerance for same molecule",toler,10,15,"");
		float G=4.5f;
		gd.addNumericField("Camera Gain",G,10,15,"");
		int fidavg=50;
		gd.addNumericField("# Frames to Avg Fiducials",fidavg,0,15,null);
		boolean showfits=false;
		gd.addCheckbox("Show Molecules and Fits?",showfits);
		gd.addCheckbox("Plot Window Drift Correction",false);
		gd.showDialog();
		if(gd.wasCanceled()){return;}
		edge_buffer=(float)gd.getNextNumber();
		threshhold=(float)gd.getNextNumber();
		baseline=(float)gd.getNextNumber();
		minarea=(int)gd.getNextNumber();
		maxarea=(int)gd.getNextNumber();
		min_int=(float)gd.getNextNumber();
		maxblobs=(int)gd.getNextNumber();
		stdev=(double)gd.getNextNumber();
		blob_separation=(float)gd.getNextNumber();
		maxon=(int)gd.getNextNumber();
		toler=(float)gd.getNextNumber();
		G=(float)gd.getNextNumber();
		fidavg=(int)gd.getNextNumber();
		showfits=gd.getNextBoolean();
		boolean pwdrift=gd.getNextBoolean();
		float S=G+G;
		float[] criteria={(float)minarea,(float)maxarea,blob_separation,(float)maxblobs,threshhold,blob_separation,edge_buffer};
		findblobs fb=new findblobs(width,height,criteria);
		ImageStack stack=imp.getStack();
		int slices=stack.getSize();
		fitsized2=(int)(blob_separation/Math.pow(2.0,1.5));
		fitsize=fitsized2+fitsized2;

		fitmask=new boolean[fitsize*fitsize];
		for(int i=0;i<fitsize;i++){
			for(int j=0;j<fitsize;j++){
				float tempr=(float)Math.sqrt((double)((i-fitsized2)*(i-fitsized2)+(j-fitsized2)*(j-fitsized2)));
				if(tempr>fitsized2){fitmask[j+i*fitsize]=true;}
				else{fitmask[j+i*fitsize]=false;}
			}
		}

		fit_gaussian fg=new fit_gaussian((float)stdev,fitsize,0.001f);

		//now fit the fiducial, averaging every fidavg frames
		int fidslices=(int)((float)slices/(float)fidavg);
		float[][] finalfidcoords=new float[2][fidslices];
		float[][] expfidcoords=new float[2][slices];
		if(usefids){
			float[][][] fidcoords=new float[numfids][3][fidslices];
			IJ.showStatus("Fitting Fiducials");
			//IJ.log("Frame , Fiducial x, Fiducial y");
			//ImageStack fidstack=new ImageStack(fitsize,fitsize);
			//ImageStack fidfitstack=new ImageStack(fitsize,fitsize);
			for(int i=0;i<fidslices;i++){
				float[][] tempfid=new float[numfids][fitsize*fitsize];
				for(int j=0;j<fidavg;j++){
					float[] data=(float[])stack.getProcessor(i*fidavg+j+1).convertToFloat().getPixels();
					for(int l=0;l<fitsize;l++){
						for(int m=0;m<fitsize;m++){
							if(!fitmask[m+l*fitsize]){
								for(int n=0;n<numfids;n++){
									tempfid[n][m+l*fitsize]+=data[m+fidx[n]-fitsized2+width*(l+fidy[n]-fitsized2)];
								}
							}
						}
					}
				}
				for(int j=0;j<numfids;j++){
					float[] params={(float)fitsized2,(float)fitsized2,baseline*(float)fidavg,(float)fidavg*(fidmax[j]-baseline)};
					float[] fitstats=new float[2];
					//fidstack.addSlice("",tempfid);
					float[] subfit=fg.do_fit_gaussian2(params,fitstats,tempfid[j],fitmask,S,maxiterations);
					//fidfitstack.addSlice("",subfit);
					fidcoords[j][0][i]=params[0];
					fidcoords[j][1][i]=params[1];
					fidcoords[j][2][i]=params[3];
				}
				IJ.showProgress(i,fidslices);
			}
			int selfid=0;
			float startx=fidcoords[0][0][0];
			float starty=fidcoords[0][1][0];
			float startamp=fidcoords[0][2][0];
			for(int i=1;i<numfids;i++){
				if(fidcoords[i][2][0]>startamp){
					startamp=fidcoords[i][2][0];
					startx=fidcoords[i][0][0];
					starty=fidcoords[i][1][0];
					selfid=i;
				}
			}
			for(int i=0;i<fidslices;i++){
				if(fidcoords[selfid][2][i]<0.25f*startamp){
					int oldselfid=selfid;
					for(int j=0;j<numfids;j++){
						if(j!=selfid){
							if(fidcoords[j][2][i]>fidcoords[selfid][2][i]){
								selfid=j;
							}
						}
					}
					if(selfid!=oldselfid){
						startx+=fidcoords[selfid][0][i-1]-fidcoords[oldselfid][0][i-1];
						starty+=fidcoords[selfid][1][i-1]-fidcoords[oldselfid][1][i-1];
					}
				}
				finalfidcoords[0][i]=fidcoords[selfid][0][i]-startx;
				finalfidcoords[1][i]=fidcoords[selfid][1][i]-starty;
			}
			new PlotWindow4("Fiducial Drift","Time","Drift (pixels)",finalfidcoords,null).draw();
			fidcoords=null;
							
			for(int i=0;i<slices;i++){
				float place=(float)i/(float)fidavg;
				int prev=(int)place;
				float rem=place-(float)prev;
				if((prev+1)<fidslices){
					expfidcoords[0][i]=rem*(finalfidcoords[0][prev+1]-finalfidcoords[0][prev])+finalfidcoords[0][prev];
					expfidcoords[1][i]=rem*(finalfidcoords[1][prev+1]-finalfidcoords[1][prev])+finalfidcoords[1][prev];
				} else {
					expfidcoords[0][i]=finalfidcoords[0][fidslices-1];
					expfidcoords[1][i]=finalfidcoords[1][fidslices-1];
				}
			}
		} else {
			if(pwdrift){
				ImageWindow[] iw=jutils.selectPlots(false,1);
				float[] xvals=((float[][])jutils.runPW4VoidMethod(iw[0],"getXValues"))[0];
				float[] yvals=((float[][])jutils.runPW4VoidMethod(iw[0],"getYValues"))[0];
				for(int i=0;i<slices;i++){
					expfidcoords[0][i]=xvals[i]-xvals[0];
					expfidcoords[1][i]=yvals[i]-yvals[0];
				}
			}
		}

		if(edge_buffer<(float)fitsized2){edge_buffer=(float)fitsized2;}
		float[][] molecules=new float[6][500000];
		//the molecules array contains xc,yc,maxint,c2,oxc,oyc
		int[][] intmolecules=new int[4][500000];
		//the int molecules array contains x0,y0,onframe,onframes
		float[] baselines=new float[500000];
		//x0 and y0 are the first found coordinates
		boolean[] dead=new boolean[500000]; //this is a flag for whether or not this molecule has died
		Object[] fitarrays=new Object[10000]; //this contains the summed images of molecules

		int totmolecules=0;
		ImageStack molecule_stack=new ImageStack(fitsize,fitsize);
		ImageStack fit_stack=new ImageStack(fitsize,fitsize);
		int counter=0;
		//IJ.log("# , iterations , chisquared , x , y , Int");
		IJ.showStatus("Fitting Molecules");
		for(int i=0;i<slices;i++){
			float[] mask=new float[width*height];
			float[] data=(float[])stack.getProcessor(i+1).convertToFloat().getPixels();
			float[][] stats=fb.dofindblobs2(data,mask);
			for(int j=0;j<stats.length;j++){
				//check if each molecule exists
				boolean exists=false;
				for(int k=0;k<totmolecules;k++){
					if(!dead[k]){
						if(Math.abs(molecules[4][k]-stats[j][0])<toler){
							if(Math.abs(molecules[5][k]-stats[j][1])<toler){
								//the molecule exists, add its image to the existing one
								molecules[2][k]+=stats[j][2]-baseline;
								molecules[0][k]+=(stats[j][2]-baseline)*stats[j][0];
								molecules[1][k]+=(stats[j][2]-baseline)*stats[j][1];
								intmolecules[3][k]++;
								exists=true;
								float[] temp=(float[])fitarrays[k%10000];
								for(int l=0;l<fitsize;l++){
									for(int m=0;m<fitsize;m++){
										if(!fitmask[m+l*fitsize]){temp[m+l*fitsize]+=data[m+intmolecules[0][k]+width*(l+intmolecules[1][k])];}
									}
								}
								exists=true;
							}
						}
						if((i-intmolecules[2][k])>=maxon){
							//the molecule just died, fit its image to a gaussian
							dead[k]=true;
							molecules[0][k]/=molecules[2][k];
							molecules[1][k]/=molecules[2][k];
							float[] params={molecules[0][k]-(float)intmolecules[0][k],molecules[1][k]-(float)intmolecules[1][k],baseline*(float)intmolecules[3][k],molecules[2][k]};
							float[] fitstats=new float[2];
							float[] subdata=(float[])fitarrays[k%10000];
							molecule_stack.addSlice("",(Object)subdata);
							float[] subfit=fg.do_fit_gaussian2(params,fitstats,subdata,fitmask,S,maxiterations);
							fit_stack.addSlice("",(Object)subfit);
							molecules[0][k]=(float)params[0]+(float)intmolecules[0][k]-expfidcoords[0][i];
							molecules[1][k]=(float)params[1]+(float)intmolecules[1][k]-expfidcoords[1][i];
							molecules[2][k]=(float)params[3];
							molecules[3][k]=(float)fitstats[1];
							baselines[k]=params[2]/(float)intmolecules[3][k];
							//IJ.log(""+counter+" ,"+(int)fitstats[0]+" , "+fitstats[1]+" , "+molecules[0][k]+" , "+molecules[1][k]+" , "+molecules[2][k]);
							counter++;
						}
					}
				}
				if(!exists){
					if(!isfiducial(stats[j][0],stats[j][1],fidrect,numfids)){
						//this is a new molecule, initialize it
						totmolecules++;
						intmolecules[0][totmolecules-1]=(int)stats[j][0]-fitsized2;
						intmolecules[1][totmolecules-1]=(int)stats[j][1]-fitsized2;
						intmolecules[2][totmolecules-1]=i;
						intmolecules[3][totmolecules-1]=1;
						molecules[0][totmolecules-1]=stats[j][0]*(stats[j][2]-baseline);
						molecules[1][totmolecules-1]=stats[j][1]*(stats[j][2]-baseline);
						molecules[2][totmolecules-1]=stats[j][2]-baseline;
						molecules[4][totmolecules-1]=stats[j][0];
						molecules[5][totmolecules-1]=stats[j][1];
						float[] temp=new float[fitsize*fitsize];
						for(int l=0;l<fitsize;l++){
							for(int m=0;m<fitsize;m++){
								if(!fitmask[m+l*fitsize]){temp[m+l*fitsize]=data[m+intmolecules[0][totmolecules-1]+width*(l+intmolecules[1][totmolecules-1])];}
							}
						}
						fitarrays[(totmolecules-1)%10000]=temp;
					}
				}
			}
			IJ.showProgress(i,slices);
			if(IJ.escapePressed()){return;}
		}

		//analyze all of the molecules that haven't died yet
		for(int k=0;k<totmolecules;k++){
			if(!dead[k]){
				molecules[0][k]/=molecules[2][k];
				molecules[1][k]/=molecules[2][k];
				float[] params={molecules[0][k]-(float)intmolecules[0][k],molecules[1][k]-(float)intmolecules[1][k],baseline*(float)intmolecules[3][k],molecules[2][k]};
				float[] fitstats=new float[2];
				float[] subdata=(float[])fitarrays[k%10000];
				molecule_stack.addSlice("",(Object)subdata);
				float[] subfit=fg.do_fit_gaussian2(params,fitstats,subdata,fitmask,S,maxiterations);
				fit_stack.addSlice("",(Object)subfit);
				molecules[0][k]=(float)params[0]+(float)intmolecules[0][k]-expfidcoords[0][slices-1];
				molecules[1][k]=(float)params[1]+(float)intmolecules[1][k]-expfidcoords[1][slices-1];
				molecules[2][k]=(float)params[3];
				molecules[3][k]=(float)fitstats[1];
				baselines[k]=params[2]/(float)intmolecules[3][k];
				//IJ.log(""+counter+" , "+(int)fitstats[0]+" , "+fitstats[1]+" , "+molecules[0][k]+" , "+molecules[1][k]+" , "+molecules[2][k]);
				counter++;
			}
		}
		if(showfits){
			(new ImagePlus("Molecules",molecule_stack)).show();
			(new ImagePlus("Fits",fit_stack)).show();
		}

		float[] output=new float[5*totmolecules+5];
		output[0]=(float)width;
		if(height>width){
			output[0]=(float)height;
		}
		output[2]=(float)stdev;
		output[3]=(float)G;
		double totbase=0.0;
		for(int i=0;i<totmolecules;i++){
			output[5*i+5]=molecules[0][i];
			output[5*i+5+1]=molecules[1][i];
			output[5*i+5+2]=molecules[2][i];
			output[5*i+5+3]=molecules[3][i];
			output[5*i+5+4]=(float)intmolecules[2][i];
			totbase+=baselines[i];
		}
		output[1]=(float)(totbase/(double)totmolecules);
		//IJ.log("Baseline = "+output[1]);
		new ImagePlus("Molecule Data",new FloatProcessor(5,totmolecules+1,output,null)).show();

		final PALMWindow cw = new PALMWindow();
		cw.init(output);
		PALMWindow.launch_frame(cw);
	}

	private boolean isfiducial(float x,float y,Rectangle[] fidrect,int numfids){
		for(int i=0;i<numfids;i++){
			if(fidrect[i].contains(x,y)){
				return true;
			}
		}
		return false;
	}

	private float[] short2float(short[] data){
		float[] temp=new float[data.length];
		for(int i=0;i<data.length;i++){
			temp[i]=(float)(data[i]&0xffff);
		}
		return temp;
	}
}

