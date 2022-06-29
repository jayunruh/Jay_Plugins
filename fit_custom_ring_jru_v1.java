/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
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
import jalgs.*;
import jalgs.jfit.*;
import jguis.*;
import ij.plugin.frame.*;
import ij.text.*;

public class fit_custom_ring_jru_v1 implements PlugIn, NLLSfitinterface_v2 {
	int xpts,ypts,zpts,profwidth,ngaus;
	float zratio;
	gausfunc gf;

	//this plugin fits and realigns a 3D ring in the selected channel
	//the input is four xy positions moving clockwise from 12 oclock
	//we then fit the orthogonal and diagonal xz profiles to refine those four points constrained by the 3D ring equation
	//finally we realign so that all four points are in the central plane
	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		xpts=imp.getWidth(); ypts=imp.getHeight();
		ImageStack stack=imp.getStack();
		//zpts=stack.getSize();
		zratio=(float)jutils.get_zratio(imp);
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Z_Ratio",zratio,5,15,null);
		gd.addNumericField("X_Stdev",0.95,5,15,null);
		gd.addNumericField("Z_Stdev",4.0,5,15,null);
		gd.addCheckbox("Calibrate?",true);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		zratio=(float)gd.getNextNumber();
		float startstdev=(float)gd.getNextNumber();
		float startzstdev=(float)gd.getNextNumber();
		boolean cal=gd.getNextBoolean();
		float psize=(float)jutils.get_psize(imp);
		RoiManager rman=RoiManager.getInstance();
		if(rman==null){
			IJ.error("need point rois");
		}
		Roi[] rois=rman.getRoisAsArray();
		int slices=imp.getNSlices();
		zpts=slices;
		int frames=imp.getNFrames();
		int channels=imp.getNChannels();
		int currframe=imp.getT();
		int currchan=imp.getC();
		int currslice=imp.getZ();
		Object[] stack2=jutils.get3DZSeries(stack,currchan-1,currframe-1,frames,slices,channels);
		//Object[] stack2=jutils.stack2array(stack);
		//duplicate the stack to get a float array
		float[][] fstack=algutils.get_region2(stack2,0,0,xpts,ypts,xpts,ypts);
		ngaus=rois.length; //this should always be 4

		//first get the horizontal and vertical positions
		int vxc1=rois[0].getBounds().x;
		int vyc1=rois[0].getBounds().y;
		int vxc2=rois[2].getBounds().x;
		int vyc2=rois[2].getBounds().y;
		float vxc=0.5f*((float)vxc1+(float)vxc2);
		float vyc=0.5f*((float)vyc1+(float)vyc2);
		float vd=(float)(vyc1-vyc2);
		int hxc1=rois[1].getBounds().x;
		int hyc1=rois[1].getBounds().y;
		int hxc2=rois[3].getBounds().x;
		int hyc2=rois[3].getBounds().y;
		float hxc=0.5f*((float)hxc1+(float)hxc2);
		float hyc=0.5f*((float)hyc1+(float)hyc2);
		float hd=(float)(hxc1-hxc2);
		float maxd=(float)Math.max(hd,vd);
		float xc=vxc;
		float yc=hyc;
		//the profiles move in a counter-clockwise direction
		//next get the vertical xz profile
		int xzwidth=2;
		float[] temp=profiler.get2DLineProfile(new float[]{xc,yc-maxd,xc,yc+maxd},fstack[0],1,0,xzwidth,xpts,ypts);
		profwidth=temp.length;
		float[] vxzprof=new float[profwidth*zpts];
		System.arraycopy(temp,0,vxzprof,0,profwidth);
		for(int i=1;i<zpts;i++){
			temp=profiler.get2DLineProfile(new float[]{xc,yc-maxd,xc,yc+maxd},fstack[i],1,0,xzwidth,xpts,ypts);
			System.arraycopy(temp,0,vxzprof,i*profwidth,profwidth);
		}
		//and now the horizontal xz profile
		float[] hxzprof=new float[profwidth*zpts];
		for(int i=0;i<zpts;i++){
			temp=profiler.get2DLineProfile(new float[]{xc-maxd,yc,xc+maxd,yc},fstack[i],1,0,xzwidth,xpts,ypts);
			System.arraycopy(temp,0,hxzprof,i*profwidth,profwidth);
		}
		//and now the downward diagonal profile
		//need different x and y distances
		float dist2=maxd*0.70711f; //(mult by cos(45 deg))
		float[] d1xzprof=new float[profwidth*zpts];
		for(int i=0;i<zpts;i++){
			temp=profiler.get2DLineProfile(new float[]{xc-dist2,yc-dist2,xc+dist2,yc+dist2},fstack[i],1,0,xzwidth,xpts,ypts);
			if(temp.length<profwidth) System.arraycopy(temp,0,d1xzprof,i*profwidth,temp.length);
			else System.arraycopy(temp,0,d1xzprof,i*profwidth,profwidth);
		}
		//finally the upward diagonal profile
		float[] d2xzprof=new float[profwidth*zpts];
		for(int i=0;i<zpts;i++){
			temp=profiler.get2DLineProfile(new float[]{xc-dist2,yc+dist2,xc+dist2,yc-dist2},fstack[i],1,0,xzwidth,xpts,ypts);
			if(temp.length<profwidth) System.arraycopy(temp,0,d2xzprof,i*profwidth,temp.length);
			else System.arraycopy(temp,0,d2xzprof,i*profwidth,profwidth);
		}
		float[][] profiles=new float[][]{vxzprof,d1xzprof,hxzprof,d2xzprof};
		new ImagePlus("Profiles",jutils.array2stack(profiles,profwidth,zpts)).show();
		//concatenate the profiles for fitting
		float[] fitdata=new float[4*profiles[0].length];
		for(int i=0;i<4;i++) System.arraycopy(profiles[i],0,fitdata,i*profiles[i].length,profiles[i].length);
		//now initialize the fit parameters
		//params are 0baseline, 1radius, 2xc, 3yc, 4zc, 5zenith angle, 6azimuth angle, 7stdev, 8zstdev, and 4 amplitudes and amp ratios
		String[] paramsnames={"baseline","radius","xc","yc","zc","zenith","azimuth","stdev","zstdev","ampv","amprv","amp45","ampr45","amph","amprh","amp135","ampr135"};
		double[][] zstats=initZParams(profiles,maxd,3);
		double[] params=new double[17];
		float[] mins=jstatistics.getspectrum("Min",profiles,null);
		float[] maxs=jstatistics.getspectrum("Max",profiles,null);
		float max=jstatistics.getstatistic("Max",maxs,null);
		params[0]=jstatistics.getstatistic("Min",mins,null);
		params[1]=0.5*(double)maxd;
		IJ.log("guess_rad = "+params[1]);
		//xc,yc, zenith, and azimuth values are initialized to 0
		params[4]=zstats[profiles.length*2][0];
		params[7]=(double)startstdev;
		params[8]=(double)startzstdev;
		for(int i=0;i<4;i++){
			params[9+2*i]=0.5*(zstats[i*2][1]+zstats[i*2+1][1]);
			params[10+2*i]=zstats[i*2][1]/zstats[i*2+1][1];
		}
		//now set the constraints
		double[][] constraints=new double[2][17];
		constraints[0][0]=-max;  constraints[1][0]=max;
		constraints[0][1]=0.25*(double)maxd; constraints[1][1]=1.5*(double)maxd;
		constraints[0][2]=-0.2*params[1]; constraints[1][2]=0.2*params[1];
		constraints[0][3]=-0.2*params[1]; constraints[1][3]=0.2*params[1];
		constraints[0][4]=params[4]-zratio; constraints[1][4]=params[4]+zratio;
		constraints[0][5]=0.0; constraints[1][5]=0.25*3.14159265; //tilt angle is between 0 and 45 deg
		constraints[0][6]=-Math.PI; constraints[1][6]=Math.PI; //orientation (azimuth) is from -180 to +180
		constraints[0][7]=0.5*params[7]; constraints[1][7]=1.5*params[7];
		constraints[0][8]=0.5*params[8]; constraints[1][8]=1.5*params[8];
		for(int i=0;i<4;i++){
			constraints[0][9+2*i]=0.5*params[9+2*i]; constraints[1][9+2*i]=1.5*params[9+2*i];
			constraints[0][10+2*i]=0.5; constraints[1][10+2*i]=2.0;
		}
		int[] fixes=new int[17];
		fixes[6]=1; //fix the azimuth (orientation) at different values and fit each time
		gf=new gausfunc();
		//set up the fit
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0.0001,10,0.1);
		double[] fitstats=new double[2];
		double[] origparams=params.clone();
		double[] bestparams=null;
		double[] bestfitstats=null;
		float[] bestfit=null;
		double bestc2=-1.0;
		int bestiter=0;
		//cycle through the azimuth options and pick the best one
		for(double az=-Math.PI;az<Math.PI; az+=Math.PI/10.0){
			params=origparams.clone();
			params[6]=az;
			float[] tempfit=fitclass.fitdata(params,fixes,constraints,fitdata,null,fitstats,true);
			if(bestc2<0.0 || fitstats[1]<bestc2){
				bestc2=fitstats[1];
				bestiter=(int)fitstats[0];
				bestfit=tempfit;
				bestparams=params;
			}
		}
		float[][] fit2=new float[4][];
		for(int i=0;i<4;i++){
			fit2[i]=algutils.convert_arr_float(algutils.get_subarray(bestfit,i*zpts*profwidth,zpts*profwidth));
		}
		new ImagePlus("Fit_Profiles",jutils.array2stack(fit2,profwidth,zpts)).show();
		//now report the fit values
		TextWindow tw=jutils.selectTable("Custom Ring Parameters");
		if(tw==null) tw=new TextWindow("Custom Ring Parameters","name\tc2\titerations\t"+table_tools.print_string_array(paramsnames),"",400,200);
		tw.append(imp.getTitle()+"\t"+(float)bestc2+"\t"+bestiter+"\t"+table_tools.print_double_array(bestparams));
		//float[][] fitsim=simRing(bestparams,xc,yc);
		//new ImagePlus("3D Ring Sim",jutils.array2stack(fitsim,xpts,ypts)).show();
		//now we need to rotate the ring so that it is flat and the best phase is pointing upwards
		//start by getting the amplitude profile linearly along the ring in ring angle coordinates
		float[] ampprofile=getThetaAmpProfile(params,40);
		double phase=findPhase(algutils.convert_arr_double(ampprofile));
		double orthphase=wrap(phase+Math.PI/2.0);
		double oppphase=wrap(phase+Math.PI); 
		double orthphase2=wrap(phase+1.5*Math.PI);
		//use 4 unique points to define the transformation: rotate clockwise around the ring from the best phase (this way the centroid is the center of the ring)
		float[][] aligncoords=new float[3][4];
		double[] dtemp=get3DRing(phase,bestparams[1],bestparams[2],bestparams[3],bestparams[4],bestparams[5],bestparams[6]);
		aligncoords[0][0]=(float)dtemp[0]+(float)xc; aligncoords[1][0]=(float)dtemp[1]+(float)yc; aligncoords[2][0]=(float)dtemp[2];
		dtemp=get3DRing(orthphase,bestparams[1],bestparams[2],bestparams[3],bestparams[4],bestparams[5],bestparams[6]);
		aligncoords[0][1]=(float)dtemp[0]+(float)xc; aligncoords[1][1]=(float)dtemp[1]+(float)yc; aligncoords[2][1]=(float)dtemp[2];
		dtemp=get3DRing(oppphase,bestparams[1],bestparams[2],bestparams[3],bestparams[4],bestparams[5],bestparams[6]);
		aligncoords[0][2]=(float)dtemp[0]+(float)xc; aligncoords[1][2]=(float)dtemp[1]+(float)yc; aligncoords[2][2]=(float)dtemp[2];
		dtemp=get3DRing(orthphase2,bestparams[1],bestparams[2],bestparams[3],bestparams[4],bestparams[5],bestparams[6]);
		aligncoords[0][3]=(float)dtemp[0]+(float)xc; aligncoords[1][3]=(float)dtemp[1]+(float)yc; aligncoords[2][3]=(float)dtemp[2];
		IJ.log(table_tools.print_float_array(aligncoords));
		
		//we are transforming the coordinates to versions centered at 0,0,0
		float[][] refcoords={{0.0f,(float)bestparams[1],0.0f,-(float)bestparams[1]},{(float)bestparams[1],0.0f,-(float)bestparams[1],0.0f},{0.0f,0.0f,0.0f,0.0f}};
		
		IJ.log(table_tools.print_float_array(refcoords));
		Object[] trans=(new jreg()).scaled_rotation_fiducials(refcoords,aligncoords);
		//float[] t=(float[])trans[1];
		float[][] r=(float[][])trans[0];
		IJ.log(table_tools.print_float_array(r));
		//now map realigned image coordinates back to the original image and interpolate to realign
		float[] newcenter={0.5f*(float)xpts,0.5f*(float)ypts,0.5f*(float)zpts*zratio};
		float[][] t3D=new float[zpts*channels][xpts*ypts];
		//float[][] tsim=new float[zpts][xpts*ypts];
		for(int i=0;i<channels;i++){
			Object[] srcpix=jutils.get3DZSeries(stack,i,0,frames,slices,channels);
			for(int j=0;j<zpts;j++){
				for(int k=0;k<ypts;k++){
					for(int l=0;l<xpts;l++){
						float[] coords={(float)l-newcenter[0],(float)k-newcenter[1],zratio*(float)j-newcenter[2]};
						temp=matmult(coords,r);
						temp[0]+=(float)(bestparams[2]+(double)xc); temp[1]+=(float)(bestparams[3]+(double)yc); temp[2]=(temp[2]+(float)bestparams[4])/zratio;
						t3D[j*channels+i][l+k*xpts]=interpolation.interp3D(srcpix,xpts,ypts,temp[0],temp[1],temp[2]);
						//if(i==0) tsim[j][l+k*xpts]=interpolation.interp3D(fitsim,xpts,ypts,temp[0],temp[1],temp[2]);
					}
				}
			}
		}
		//new ImagePlus("3D Ring Sim Realigned",jutils.array2stack(tsim,xpts,ypts)).show();
		ImagePlus reimp=new ImagePlus("Realigned",jutils.array2stack(t3D,xpts,ypts));
		reimp.setOpenAsHyperStack(true);
		reimp.setDimensions(channels,zpts,1);
		new CompositeImage(reimp,CompositeImage.COLOR).show();
		//firstly, rotate the ring about the 
		//for testing purposes rotate the ring at a specified tilt angle
		/*ImageStack teststack=new ImageStack(profwidth,zpts);
		params[5]=1.0;
		for(double az=-Math.PI;az<Math.PI; az+=Math.PI/20.0){
			params[6]=az;
			double[] fit=fitfunc(params);
			for(int i=0;i<4;i++){
				teststack.addSlice(""+az,algutils.convert_arr_float(algutils.get_subarray(fit,i*zpts*profwidth,zpts*profwidth)));
			}
		}
		ImagePlus testimp=new ImagePlus("Fit_Profiles",teststack);
		testimp.setDimensions(1,4,teststack.getSize()/4);
		testimp.setOpenAsHyperStack(true);
		testimp.show();*/
	}

	public float[] matmult(float[] vec,float[][] mat){
		float[] newvec=new float[vec.length];
		for(int i=0;i<vec.length;i++){
			for(int j=0;j<vec.length;j++){
				newvec[i]+=vec[j]*mat[i][j];
			}
		}
		return newvec;
	}

	public double[][] getAmpProfile(double[] params){
		//calculate the amplitude profile
		double avgamp=0.25*(params[9]+params[11]+params[13]+params[15]);
		double[] ampprofile=new double[8];
		double[] anglevals=new double[8];
		for(int i=0;i<4;i++){
			ampprofile[i+4]=2.0*params[9+2*i]/(1.0+params[10+2*i]);
			ampprofile[i]=2.0*params[9+2*i]-ampprofile[i+4];
			ampprofile[i]/=avgamp; ampprofile[i+4]/=avgamp;
			anglevals[i]=2.0*Math.PI*(double)i/8.0;
			anglevals[i+4]=2.0*Math.PI*(double)(i+4)/8.0;
		}
		return new double[][]{anglevals,ampprofile};
	}

	public double findPhase(double[] ampprofile1){
		float[] ampprofile=algutils.convert_arr_float(ampprofile1);
		int len=ampprofile.length;
		linleastsquares lls=new linleastsquares();
		double minc2=-1.0;
		double minphase=0.0;
		for(double phase=0.0;phase<2.0*Math.PI;phase+=Math.PI/20.0){
			double[] cosvals=new double[len];
			for(int i=0;i<len;i++) cosvals[i]=0.5+0.5*Math.cos(2.0*Math.PI*(double)i/(double)len);
			double[] ampoff=lls.get_amp_offset(cosvals,ampprofile,true);
			double c2=lls.get_amp_offset_c2(cosvals,ampprofile,ampoff);
			if(minc2<0.0 || c2<minc2){
				minc2=c2;
				minphase=phase;
			}
		}
		return minphase;
	}

	public double[][] initZParams(float[][] profiles,float dist,int width){
		//here we get the max z position and the baseline and amplitude for each point of each xz profile
		int start1=(int)(0.5f*(float)profwidth-0.5f*dist-0.5f*(float)width);
		int start2=(int)(0.5f*(float)profwidth+0.5f*dist-0.5f*(float)width);
		double[][] stats=new double[profiles.length*2+2][2];
		for(int i=0;i<profiles.length;i++){
			float[][] zprofiles=new float[4][zpts];
			for(int j=0;j<zpts;j++){
				for(int k=0;k<width;k++){
					float temp1=profiles[i][start1+k+j*profwidth];
					zprofiles[0][j]+=temp1;
					float temp2=profiles[i][start2+k+j*profwidth];
					zprofiles[1][j]+=temp2;
					if(temp1>zprofiles[2][j]) zprofiles[2][j]=temp1;
					if(temp2>zprofiles[3][j]) zprofiles[3][j]=temp2;
				}
			}
			float maxint1=0.0f; int maxpos1=0; float amp1=0.0f;
			float maxint2=0.0f; int maxpos2=0; float amp2=0.0f;
			for(int j=0;j<zpts;j++){
				if(zprofiles[0][j]>maxint1){maxint1=zprofiles[0][j]; maxpos1=j;}
				if(zprofiles[2][j]>amp1) amp1=zprofiles[2][j];
				if(zprofiles[1][j]>maxint2){maxint2=zprofiles[1][j]; maxpos2=j;}
				if(zprofiles[3][j]>amp2) amp2=zprofiles[3][j];
			}
			stats[i*2][0]=zratio*(float)maxpos1;
			stats[i*2][1]=amp1;
			stats[i*2+1][0]=zratio*(float)maxpos2;
			stats[i*2+1][1]=amp2;
		}
		//calculate the average best z positions
		double avgz=0.0;
		for(int i=0;i<profiles.length*2;i++){avgz+=stats[i][0];}
		avgz/=(double)(profiles.length*2);
		stats[profiles.length*2][0]=avgz;
		return stats;
	}

	public float[][] simRing(double[] params,float xshift,float yshift){
		float[][] simpix=new float[zpts][xpts*ypts];
		//work our way around the ring simulating gaussians as we go
		double stdevz=params[8]/zratio;
		//double stdev=params[7];
		double stdev=0.5;
		//IJ.log(""+stdev);
		//IJ.log(""+stdevz);
		int nangles=40;
		float[] tampprofile=getThetaAmpProfile(params,nangles);
		int counter=0;
		for(double t=0.0; t<2.0*Math.PI;t+=2.0*Math.PI/nangles){
			double[] coords=get3DRing(t,params[1],params[2],params[3],params[4],params[5],params[6]);
			//float amp=tampprofile[counter];
			float amp=1.0f;
			//IJ.log(""+coords[0]+"\t"+coords[1]+"\t"+coords[2]+"\t"+amp);
			gf.draw_3D_func(simpix,coords[0]+(double)xshift,coords[1]+(double)yshift,coords[2]/zratio,xpts,ypts,stdev,stdevz,amp);
			counter++;
		}
		return simpix;
	}

	public float[] getThetaAmpProfile(double[] params,int nangles){
		//our crossing points aren't perfectly spaced around the ring
		//we need to interpolate to map these back to our ring theta values
		//note that the crossings matrix jumps back and forth across the ring while the ampprofile goes in counterclockwise order
		double[][] crossings=getCrossings(params[1],params[2],params[3],params[4],params[5],params[6]);
		double[][] ampprofile=getAmpProfile(params);
		double[] thetavals=new double[crossings.length];
		for(int i=0;i<4;i++) {thetavals[i]=crossings[2*i][4]; thetavals[4+i]=crossings[2*i+1][4];}
		float[] tampprofile=new float[nangles];
		int counter=0;
		for(double theta=0.0;theta<2.0*Math.PI;theta+=(2.0*Math.PI)/nangles){
			//for each theta value, find the crossing point above and below it
			int above=0; int below=0;
			double abovedist; double belowdist; double tabovedist; double tbelowdist;
			if(thetavals[0]>theta){abovedist=thetavals[0]-theta; belowdist=2.0*Math.PI-abovedist;}
			else{belowdist=theta-thetavals[0]; abovedist=2.0*Math.PI-belowdist;}
			for(int i=1;i<8;i++){
				if(thetavals[i]>theta){tabovedist=thetavals[i]-theta; tbelowdist=2.0*Math.PI-tabovedist;}
				else{tbelowdist=theta-thetavals[i]; tabovedist=2.0*Math.PI-tbelowdist;}
				if(tabovedist<abovedist){abovedist=tabovedist; above=i;}
				if(tbelowdist<belowdist){belowdist=tbelowdist; below=i;}
			}
			//now interpolate between these values to find the amplitude
			double fval=belowdist/(abovedist+belowdist);
			tampprofile[counter]=(float)(fval*(ampprofile[1][above]-ampprofile[1][below])+ampprofile[1][below]);
			counter++;
		}
		return tampprofile;
	}

	public double[] fitfunc(double[] params){
		//here we are fitting 4 xc profiles (vert, 45 deg, hor, 135 deg) to a 3D ring
		//params are 0baseline, 1radius, 2xc, 3yc, 4zc, 5zenith angle, 6azimuth angle, 7stdev, 8zstdev, and 4 amplitudes and amp ratios
		double zc=params[4]/zratio;
		double stdevz=params[8]/zratio;
		double[] func2=new double[zpts*profwidth*4];
		double[][] crossings=getCrossings(params[1],params[2],params[3],params[4],params[5],params[6]);
		for(int i=0;i<4;i++){
			float[] func=new float[zpts*profwidth];
			double xc1=-crossings[2*i][3]+0.5*(double)profwidth;
			double xc2=crossings[2*i+1][3]+0.5*(double)profwidth;
			double zc1=crossings[2*i][2]/zratio;
			double zc2=crossings[2*i+1][2]/zratio;
			double a2=2.0*params[9+2*i]/(1.0+params[10+2*i]);
			double a1=2.0*params[9+2*i]-a2;
			gf.draw_2D_func(func,xc1,zc1,profwidth,zpts,params[7],stdevz,(float)a1);
			gf.draw_2D_func(func,xc2,zc2,profwidth,zpts,params[7],stdevz,(float)a2);
			for(int j=0;j<zpts*profwidth;j++){
				func2[j+i*zpts*profwidth]=(double)func[j]+params[0];
			}
		}
		return func2;
	}

	public double[] get3DRing(double t,double r,double xc,double yc,double zc,double zenith, double azimuth){
		double x=xc-r*Math.cos(t)*Math.sin(azimuth)+r*Math.sin(t)*Math.cos(zenith)*Math.cos(azimuth);
		double y=yc+r*Math.cos(t)*Math.cos(azimuth)+r*Math.sin(t)*Math.cos(zenith)*Math.sin(azimuth);
		double z=zc-r*Math.sin(t)*Math.sin(zenith);
		return new double[]{x,y,z};
	}

	public double[][] getCrossings(double r,double xc,double yc,double zc,double zenith,double azimuth){
		//here we get all the coordinates where a ring crosses the xz, yz and diagonal z planes
		//start by getting the angle along the ring at which these crossings happen
		//first the yz plane (x=0)
		double b=Math.cos(zenith)*Math.cos(azimuth);
		double a1=Math.sin(azimuth)/b;
		double c1=xc/(r*b);
		double temp=Math.sqrt(1+a1*a1-c1*c1);
		double yztheta1=Math.atan2(-c1-a1*temp,a1*c1-temp); 
		double yztheta2=Math.atan2(-c1+a1*temp,a1*c1+temp);
		yztheta1=wrap(yztheta1); yztheta2=wrap(yztheta2);
		//now the xz plane (y=0)
		b=Math.cos(zenith)*Math.sin(azimuth);
		double xztheta1=0.0; double xztheta2=0.0;
		if(b==0.0f){
			double s=yc/r;
			xztheta1=Math.atan2(Math.sqrt(1.0-s*s),s);
			xztheta2=-xztheta1;
			xztheta1=wrap(xztheta1); xztheta2=wrap(xztheta2);
		} else {
			a1=Math.cos(azimuth)/b;
			c1=yc/(r*b);
			temp=Math.sqrt(1+a1*a1-c1*c1);
			xztheta1=Math.atan2(c1+a1*temp,a1*c1-temp);
			xztheta2=Math.atan2(c1-a1*temp,a1*c1+temp);
			xztheta1=wrap(xztheta1); xztheta2=wrap(xztheta2);
		}
		//now the upward diagonal plane (x=y)
		b=Math.cos(zenith)*Math.sin(azimuth)-Math.cos(zenith)*Math.cos(azimuth);
		a1=(Math.sin(azimuth)+Math.cos(azimuth))/b;
		c1=(xc-yc)/(r*b);
		temp=Math.sqrt(1+a1*a1-c1*c1);
		double d1ztheta1=Math.atan2(c1+a1*temp,a1*c1-temp);
		double d1ztheta2=Math.atan2(c1-a1*temp,a1*c1+temp);
		d1ztheta1=wrap(d1ztheta1); d1ztheta2=wrap(d1ztheta2);
		//now the downward diagonal plane (x=-y)
		b=-Math.cos(zenith)*Math.sin(azimuth)-Math.cos(zenith)*Math.cos(azimuth);
		a1=(Math.sin(azimuth)-Math.cos(azimuth))/b;
		c1=(xc+yc)/(r*b);
		temp=Math.sqrt(1+a1*a1-c1*c1);
		double d2ztheta1=Math.atan2(c1+a1*temp,a1*c1-temp);
		double d2ztheta2=Math.atan2(c1-a1*temp,a1*c1+temp);
		d2ztheta1=wrap(d2ztheta1); d2ztheta2=wrap(d2ztheta2);
		double[] angles=new double[]{yztheta1,yztheta2,d1ztheta1,d1ztheta2,xztheta1,xztheta2,d2ztheta1,d2ztheta2};
		//IJ.log(table_tools.print_double_array(angles));
		double[][] crossings=new double[8][4];
		//for each crossing, list x, y, z, d (xy distance from origin), theta
		for(int i=0;i<8;i++){
			double[] coords=get3DRing(angles[i],r,xc,yc,zc,zenith,azimuth);
			crossings[i]=new double[]{coords[0],coords[1],coords[2],0.0,angles[i]};
			crossings[i][3]=Math.sqrt(coords[0]*coords[0]+coords[1]*coords[1]);
		}
		//need to make sure the crossings are in the right order
		//start with the yz crossings
		if(crossings[0][1]>crossings[1][1]){double[] temp1=crossings[1]; crossings[1]=crossings[0]; crossings[0]=temp1;}
		//now the xz crossings
		if(crossings[2][0]>crossings[3][0]){double[] temp1=crossings[3]; crossings[3]=crossings[2]; crossings[2]=temp1;}
		//now the upward diagonal
		if(crossings[4][0]>crossings[5][0]){double[] temp1=crossings[5]; crossings[5]=crossings[4]; crossings[4]=temp1;}
		//and finally the downward diagonal
		if(crossings[6][1]<crossings[7][1]){double[] temp1=crossings[7]; crossings[7]=crossings[6]; crossings[6]=temp1;}
		return crossings;
	}

	public double wrap(double angle){
		if(angle<0.0) return angle+2.0*Math.PI;
		if(angle>=2.0*Math.PI) return angle-2.0*Math.PI;
		return angle;
	}

	public void showresults(String results){
		IJ.log(results);
	}

}
