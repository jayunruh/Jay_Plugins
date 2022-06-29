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
import ij.plugin.frame.*;
import jalgs.*;
import jalgs.jsim.*;
import jalgs.jfft.*;
import jalgs.jfit.*;
import jguis.*;

public class sim_exo_endo_jru_v2 implements PlugIn {
	//this is the full simulation including membrane insertion
	//here we confine exocytosis in the cap to five patches
	//endocytosis is excluded from these patches
	//the diffusion coefficient is slower inside the patches
	int size,maxsize,center,exosize,endosize,intarea,exocount,endocount,totarea,npatches,dsteps,patchsize;
	float windrad,exorate,endorate,endorate2,D,Dpatch,T,tot;
	float windrad2,exorate2,D2,D2patch;
	float intconc,pwind,endoconcratio,exoconcratio,membconc,minconc,maxconc,kbind,kunbind,boundconc;
	float ksyn,kdeg,ksyn2,kdeg2;
	float[][] profile,mask,bound,D2profile;
	float[][] patchcoords;
	boolean usepatches;
	rngs random;
	float[] fftconvprofile,fftconvprofilepatch;
	convolution2D conv2D;
	gauss_convolution gc;

	public void run(String arg) {
		maxsize=150;
		size=114;
		center=75;
		exosize=2;
		endosize=1;
		float mintime=30.0f; //total time in minutes
		float tottime=mintime*60.0f; //total time in seconds
		T=0.1f; //frametime in seconds (changed to 0.1 for compatibility with FRAP analysis)
		Dpatch=0.0043f;
		D=0.011f;
		usepatches=true; //tells whether to use the patches for endo/exocytosis
		endoconcratio=1.0f;
		exoconcratio=10.0f;
		float intarearatio=1.0f;
		boolean convolve=true; //if true we convolve with a confocal psf
		patchsize=1; //this is the patch half size used for exocytosis
		float patchsize2=0.5f; //this is the patch size value used for the diffusion coefficient
		int patchdist=7; //distance of radial patches from center
		int patchdist2=13;
		int frames=(int)(tottime/T);
		float psize=0.088f; //pixel size in um
		float perimeter=15.71f;
		float diameter=perimeter/(float)Math.PI;
		windrad=0.13f*perimeter; //window radius in um
		float Awind=(float)Math.PI*windrad*windrad;
		float Amemb=(float)Math.PI*(float)size*(float)size/4.0f;
		windrad2=windrad/psize;
		float windstdev2=windrad2/1.177f;
		endorate=1.0f/0.6f; //total endocytosis rate in vessicles/sec from Lew et al
		endorate2=endorate*T; //vessicles/frame
		exorate=endorate/4.0f; //total exocytosis rate in vessicles/sec
		exorate2=exorate*T; //vessicles/frame
		dsteps=1;
		npatches=19;
		float patchangle1=(float)Math.toRadians(60.0f);
		float patchangle2=(float)Math.toRadians(30.0f);
		//patchcoords=new int[npatches][2]; for(int i=0;i<npatches;i++){patchcoords[i][0]=center; patchcoords[i][1]=center;}
		//patchcoords[0][0]-=patchdist; patchcoords[1][0]+=patchdist; patchcoords[2][1]-=patchdist; patchcoords[3][1]+=patchdist;
		patchcoords=new float[npatches][2];
		for(int i=0;i<6;i++){
			patchcoords[i]=rotate_vector(patchdist,(float)i*patchangle1,center,center);
		}
		for(int i=0;i<12;i++){
			patchcoords[i+6]=rotate_vector(patchdist2,(float)i*patchangle2,center,center);
		}
		patchcoords[18]=new float[]{center,center};
		intarea=(int)(Math.PI*(float)size*(float)size);
		intconc=1.0f;
		float ratio=3.0f;
		pwind=ratio*Awind/(Amemb+(ratio-1.0f)*Awind); //probability of endocytosis in the window
		//kbind=0.00001f; //rate of binding to non diffusible non endocytosible region in proteins/pixel/frame
		//kunbind=0.0001f; //rate of unbinding to non diffusible non endo or exocytosible region in proteins/pixel/frame
		kbind=0.0f; kunbind=0.0f;
		//ksyn=0.00016f; //rate of synthesis in fraction of total molecules per second
		//ksyn=0.016f;
		ksyn=0.0f;
		kdeg=0.0f; //rate of degradation in per seconds
		D2=D*T/(psize*psize);
		D2patch=Dpatch*T/(psize*psize);
		profile=new float[maxsize][maxsize];
		//bound=new float[maxsize][maxsize];
		mask=new float[maxsize][maxsize];
		totarea=0;
		float startconc=1.0f;
		for(int i=0;i<maxsize;i++){for(int j=0;j<maxsize;j++){profile[i][j]=startconc;}}
		for(int i=(center-size/2);i<(center-size/2+size);i++){
			for(int j=(center-size/2);j<(center-size/2+size);j++){
				double rad=Math.sqrt((i-center)*(i-center)+(j-center)*(j-center));
				if(rad<=size/2){
					totarea++;
					mask[i][j]=1.0f;
				}
			}
		}
		D2profile=getDprofile(maxsize,maxsize,patchcoords,patchsize2,D2,D2patch);
		intarea=(int)(intarearatio*totarea);
		totarea+=intarea;
		tot=intconc*(float)intarea+startconc*(float)(totarea-intarea);
		ksyn2=ksyn*T/(tot*(float)intarea); //here is the synthesis rate in units of molecules per frame per unit area for the internal membrane
		kdeg2=kdeg*T; //rate of degradation in per frame
		float[] outprofile=new float[maxsize*frames];
		float[] areaprofile=new float[frames];
		float[] intareaprofile=new float[frames];
		float[] intconcprofile=new float[frames];
		float[] membconcprofile=new float[frames];
		float[] totconcprofile=new float[frames];
		float[] ptratio=new float[frames];
		ImagePlus imp=new ImagePlus("Simulation",new FloatProcessor(maxsize,maxsize,new float[maxsize*maxsize],null));
		//ImagePlus imp2=new ImagePlus("Mask",new FloatProcessor(maxsize,maxsize,new float[maxsize*maxsize],null));
		imp.show();
		//imp2.show();
		ImageStack simstack=new ImageStack(maxsize,maxsize);
		random=new rngs();
		float[] pixels=(float[])imp.getProcessor().getPixels();
		//float[] maskpix=(float[])imp2.getProcessor().getPixels();
		boolean dispframe=false;
		endocount=0; exocount=0;
		for(int i=0;i<frames;i++){
			step_time();
			int temp=0;
			float temp2=0.0f;
			float temp3=0.0f;
			minconc=1000000.0f;
			maxconc=0.0f;
			if((i%100)==0){dispframe=true; for(int j=0;j<(maxsize*maxsize);j++){pixels[j]=0.0f;}}
			for(int j=0;j<maxsize;j++){
				if(!convolve && mask[j][center]>0.5f) outprofile[i*maxsize+j]=profile[j][center];
				for(int k=0;k<maxsize;k++){
					if(mask[j][k]>0.5f){
						temp++;
						temp2+=profile[j][k];
						//temp3+=bound[j][k];
						if(!convolve){
							if(profile[j][k]<minconc) minconc=profile[j][k];
							if(profile[j][k]>maxconc) maxconc=profile[j][k];
						}
						if(dispframe){pixels[j+k*maxsize]=profile[j][k];}
					}
				}
			}
			if(convolve){
				float[] conv=convolve_psf(0.20f/psize,0.8f/psize);
				for(int j=0;j<maxsize;j++){
					if(mask[size/2][j]>0.5f){
						outprofile[i*maxsize+j]=conv[j];
						if(conv[j]<minconc) minconc=conv[j];
						if(conv[j]>maxconc) maxconc=conv[j];
					}
				}
			}
			areaprofile[i]=(float)temp;
			membconc=temp2/(float)temp;
			//intarea=totarea-temp;
			//intconc=(tot-temp2)/(float)intarea;
			intareaprofile[i]=(float)intarea;
			intconcprofile[i]=intconc;
			membconcprofile[i]=membconc;
			totconcprofile[i]=(membconc*(float)temp+intconc*(float)intarea)/((float)temp+(float)intarea);
			ptratio[i]=maxconc/minconc;
			if(dispframe){
				imp.updateAndDraw();
				simstack.addSlice("",pixels.clone());
			}
			corr_conc(profile,tot,membconc*(float)temp+intconc*(float)intarea,temp);
			if(dispframe) IJ.showStatus(""+i*T+" sec out of "+tottime+" sec");
			dispframe=false;
			IJ.showProgress(i,frames);
			if(IJ.escapePressed()){break;}
		}
		//IJ.log("# endo: "+endocount+" ,# exo: "+exocount);
		PlotWindow4 pwarea=new PlotWindow4("Area Profiles","frame","area",areaprofile);
		pwarea.addPoints(intareaprofile,true);
		pwarea.draw();
		float[] limits=pwarea.getLimits();
		limits[2]=0.0f;
		pwarea.setLimits(limits);
		PlotWindow4 pwconc=new PlotWindow4("Conc Profiles","frame","conc",membconcprofile);
		pwconc.addPoints(intconcprofile,true);
		pwconc.addPoints(totconcprofile,true);
		pwconc.draw();
		new PlotWindow4("Peak to Trough Profile","frame","ptratio",ptratio).draw();
		new ImagePlus("Simulation Movie",simstack).show();
		new ImagePlus("Simulation Carpet",new FloatProcessor(maxsize,frames,outprofile,null)).show();
		IJ.run("bin image jru v1", "bin_type=[Spatial Bin] bin_by_x?=1 bin_by_y?=100");
	}

	public void corr_conc(float[][] profile,float oldtot,float tot,int membarea){
		float corrval=tot-oldtot;
		float corrval2=corrval/(float)membarea;
		for(int i=0;i<maxsize;i++){
			for(int j=0;j<maxsize;j++){
				profile[i][j]-=corrval2;
			}
		}
	}

	public float[][] getDprofile(int width,int height,float[][] patchcoords,float patchsize,float Dlow,float Dhigh){
		float[][] Dprofile=new float[width][height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				Dprofile[i][j]=Dlow;
				for(int k=0;k<patchcoords.length;k++){
					float xdist2=(float)Math.abs((float)i-patchcoords[k][0]);
					float ydist2=(float)Math.abs((float)j-patchcoords[k][1]);
					double dist2=(double)(xdist2*xdist2+ydist2*ydist2);
					float tempval=2.0f*(float)patchsize;
					tempval=(float)Math.pow(tempval*tempval/((float)dist2+tempval*tempval),0.8);
					Dprofile[i][j]+=(Dhigh-Dlow)*tempval;
					if(Dprofile[i][j]<Dhigh) Dprofile[i][j]=Dhigh;
					if(Dprofile[i][j]>Dlow) Dprofile[i][j]=Dlow;
				}
			}
		}
		//make the D matrix so it contains 32 unique values or less
		float nlevels=32.0f;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				int scaled=(int)(nlevels*(Dprofile[i][j]-Dhigh)/(Dlow-Dhigh));
				float scaled2=((float)scaled)/nlevels;
				Dprofile[i][j]=Dhigh+scaled2*(Dlow-Dhigh);
			}
		}
		return Dprofile;
	}

	public float[] rotate_vector(float r,float theta,float xcenter,float ycenter){
		return new float[]{xcenter+r*(float)Math.cos(theta),ycenter+r*(float)Math.sin(theta)};
	}

	public float[] convolve_psf(float w02,float z02){
		int xhalfsize=(int)(w02*1.5f);
		int zhalfsize=(int)(z02*1.5f);
		gausfunc gf=new gausfunc();
		float[] zprofile=gf.get_norm_func(-zhalfsize,2*zhalfsize,1.0,z02*0.5f);
		float[] xprofile=gf.get_norm_func(-xhalfsize,2*xhalfsize,1.0,w02*0.5f);
		float[] conv=new float[maxsize];
		for(int i=xhalfsize;i<(maxsize-xhalfsize);i++){
			for(int j=i-xhalfsize;j<(i+xhalfsize);j++){
				for(int k=center-zhalfsize;k<(center+zhalfsize);k++){
					conv[i]+=profile[j][k]*xprofile[j-i+xhalfsize]*zprofile[k-center+zhalfsize];
				}
			}
		}
		return conv;
	}

	public void step_time(){
		//exocytosis only occurs in the "window"
		//endocytosis occurs faster in the "window"
		if(random.expdev(1.0/exorate2)<1.0){
			exocount++;
			//we have an exocytic event
			int[] temp=get_exo_coords();
			float exoconc=exoconcratio*intconc;
			if(exoconc>0.0f){
				if(insert_membrane(temp[0],temp[1],exoconc)){
					float inttot=intconc*(float)intarea;
					intarea-=exosize*exosize;
					inttot-=(float)(exosize*exosize)*exoconc;
					intconc=inttot/intarea;
				} else {
					IJ.log("exo failed");
				}
			} else {IJ.log("exo failed");}
		}
		if(random.expdev(1.0/endorate2)<1.0){
			endocount++;
			//we have an endocytic event
			float endoconc=minconc;
			float left=0.0f;
			if(random.unidev(1.0,0.0)<pwind){
				//endocytose inside the window outside patches
				int[] temp=get_endo_coords2();
				endoconc=endoconcratio*profile[temp[0]][temp[1]];
				//if(endoconc>1.0f) endoconc=1.0f;
				left=remove_membrane(temp[0],temp[1],endoconc);
				if(left==endoconc) IJ.log("endo failed");

			} else {
				//endocytose outside the window
				int[] temp=get_endo_coords();
				endoconc=endoconcratio*profile[temp[0]][temp[1]];
				//if(endoconc>1.0f) endoconc=1.0f;
				left=remove_membrane(temp[0],temp[1],endoconc);
				if(left==endoconc) IJ.log("endo failed");
			}
			float inttot=(float)intarea*intconc;
			intarea+=endosize*endosize;
			inttot+=(float)(endosize*endosize)*endoconc-left;
			intconc=inttot/(float)intarea;
		}
		fill_outside();
		handle_multi_diffusion(dsteps);
		//handle_binding();
		//handle_synthesis();
	}

	public int[] get_endo_coords(){
		//this is outside the window
		int[] temp={(int)random.unidev((double)maxsize-1.0-Double.MIN_VALUE,1.0),(int)random.unidev((double)maxsize-1.0-Double.MIN_VALUE,1.0)};
		double length=Math.sqrt((temp[0]-0.5*(double)maxsize)*(temp[0]-0.5*(double)maxsize)+(temp[1]-0.5*(double)maxsize)*(temp[1]-0.5*(double)maxsize));
		while(length<windrad2 || !inside_mask(temp[0],temp[1])){
			int[] temp2={(int)random.unidev((double)maxsize-1.0-Double.MIN_VALUE,1.0),(int)random.unidev((double)maxsize-1.0-Double.MIN_VALUE,1.0)};
			temp=temp2;
			length=Math.sqrt((temp[0]-0.5*(double)maxsize)*(temp[0]-0.5*(double)maxsize)+(temp[1]-0.5*(double)maxsize)*(temp[1]-0.5*(double)maxsize));
		}
		return temp;
	}

	public int[] get_endo_coords2(){
		//this is inside the window, outside the patches
		double[] temp=random.random_circle(windrad2);
		int[] temp2={(int)temp[0]+center,(int)temp[1]+center};
		if(usepatches){
			while(inpatch(temp2[0],temp2[1])){
				temp=random.random_circle(windrad2);
				temp2[0]=(int)temp[0]+center; temp2[1]=(int)temp[1]+center;
			}
		}
		return temp2;
	}

	public int[] get_exo_coords(){
		//this is inside the window, inside the patches
		double[] temp=random.random_circle(windrad2);
		int[] temp2={(int)temp[0]+center,(int)temp[1]+center};
		if(usepatches){
			while(!inpatch(temp2[0],temp2[1])){
				temp=random.random_circle(windrad2);
				temp2[0]=(int)temp[0]+center; temp2[1]=(int)temp[1]+center;
			}
		}
		//double[] temp=random.random_circle(6.3f); //a central circle of exocytosis
		//int[] temp2={(int)temp[0]+center,(int)temp[1]+center};
		return temp2;
	}

	public boolean inpatch(int x,int y){
		for(int i=0;i<npatches;i++){
			if(Math.abs(x-patchcoords[i][0])<=patchsize && Math.abs(y-patchcoords[i][1])<=patchsize) return true;
			//if(x==patchcoords[i][0] && y==patchcoords[i][1]) return true;
		}
		//double dist=Math.sqrt((x-center)*(x-center)+(y-center)*(y-center));
		//if(dist<6.3) return true;
		return false;
	}

	public boolean inwindownotpatch(int x,int y){
		float length=(float)Math.sqrt((x-center)*(x-center)+(y-center)*(y-center));
		if(length<windrad2){
			if(!inpatch(x,y)) return true;
		}
		return false;
	}

	public void fill_outside(){
		float avg_edge=0.0f;
		int edgecounter=0;
		for(int i=1;i<(maxsize-1);i++){
			for(int j=1;j<(maxsize-1);j++){
				if(mask[i][j]>0.5f){
					if(!inside_mask(i,j)){
						avg_edge+=profile[i][j];
						edgecounter++;
					}
				}
			}
		}
		avg_edge/=(float)edgecounter;
		//IJ.log(""+avg_edge);
		for(int i=0;i<maxsize;i++){
			for(int j=0;j<maxsize;j++){
				if(mask[i][j]<=0.5f){
					profile[i][j]=avg_edge;
				}
			}
		}
	}

	public boolean inside_mask(int x,int y){
		//decide if we are inside the mask including a one pixel border
		if(mask[x][y]<0.5f){ return false;}
		else{
			if(mask[x-1][y-1]<0.5f) return false;
			if(mask[x][y-1]<0.5f) return false;
			if(mask[x+1][y-1]<0.5f) return false;
			if(mask[x-1][y]<0.5f) return false;
			if(mask[x+1][y]<0.5f) return false;
			if(mask[x-1][y+1]<0.5f) return false;
			if(mask[x][y+1]<0.5f) return false;
			if(mask[x+1][y+1]<0.5f) return false;
			return true;
		}
	}

	public void handle_multi_diffusion(int substeps){
		//here we convolve the profile with a gaussian for in patch and out of patch components to obtain multicomponent diffusion results	
		if(gc==null){
			gc=new gauss_convolution(0,1.0f,1.0f,D2profile,maxsize,maxsize);
			//new ImagePlus("Diffusion Profile1",new FloatProcessor((float[][])gc.profiles[0])).show();
			//new ImagePlus("Diffusion Profile2",new FloatProcessor((float[][])gc.profiles[1])).show();
		}
		profile=(float[][])gc.steptime(profile);
	}

	public void handle_binding(){
		int start=(int)center-(int)windrad2;
		int end=start+2*(int)windrad2;
		for(int i=start;i<end;i++){
			for(int j=start;j<end;j++){
				double tempr=Math.sqrt((j-center)*(j-center)+(i-center)*(i-center));
				if((float)tempr<windrad2){
					float dpdt=kunbind*bound[i][j]-kbind*profile[i][j];
					float dbdt=-dpdt;
					profile[i][j]+=dpdt;
					bound[i][j]+=dbdt;
				}
			}
		}
	}

	public void handle_synthesis(){
		//the internal membrane goes up at a constant rate
		intconc+=ksyn2;
		//degradation occurs everywhere;
		if(kdeg2>0.0f){
			intconc-=kdeg2*intconc;
			for(int i=0;i<maxsize;i++){
				for(int j=0;j<maxsize;j++){
					if(mask[i][j]>0.5f){
						profile[i][j]-=profile[i][j]*kdeg2;
					}
				}
			}
		}
	}

	public boolean insert_membrane(int x,int y,float value){
		//int temparea=0;
		//for(int i=0;i<maxsize;i++){for(int j=0;j<maxsize;j++){if(mask[i][j]>0.5f){temparea++;}}};
		//exocytosis is always assumed to insert 4 pixels
		float angle=(float)random.unidev(2.0*Math.PI,0.0);
		float[] tempvec=rotatevector((float)x,(float)y,-angle);
		int newx=(int)tempvec[0]; int newy=(int)tempvec[1];
		if(newx<2 || newx>=(maxsize-2) || newy<2 || newy>=(maxsize-2)){return false;}
		profile=interpolation.rotate_image(profile,angle,center,center);
		mask=interpolation.rotate_image(mask,angle,center,center);
		//have to account for all exocytic directions to get an isotropic effect
		float[] newprofile=profile[newx-1].clone();
		float[] newprofile1=profile[newx].clone();
		float[] newmask=mask[newx-1].clone();
		float[] newmask1=mask[newx].clone();
		int direction=0;
		if(newy<center) direction=-1;
		if(newy>center) direction=1;
		if(direction==0){
			//here we shift equally in both directions:
			for(int j=0;j<(newy-1);j++){profile[newx-1][j]=newprofile[j+1]; mask[newx-1][j]=newmask[j+1];}
			profile[newx-1][newy-1]=value;mask[newx-1][newy-1]=1.0f; profile[newx-1][newy]=value; mask[newx-1][newy]=1.0f;
			for(int j=(newy+1);j<maxsize;j++){profile[newx-1][j]=newprofile[j-1]; mask[newx-1][j]=newmask[j-1];}
			for(int j=0;j<(newy-1);j++){profile[newx][j]=newprofile[j+1]; mask[newx][j]=newmask[j+1];}
			profile[newx][newy-1]=value; mask[newx][newy-1]=1.0f; profile[newx][newy]=value; mask[newx][newy]=1.0f;
			for(int j=(newy+1);j<maxsize;j++){profile[newx][j]=newprofile[j-1]; mask[newx][j]=newmask[j-1];}
		} else {
			if(direction==-1){
				//here we just shift up
				for(int j=0;j<(newy-2);j++){profile[newx-1][j]=newprofile[j+2]; mask[newx-1][j]=newmask[j+2];}
				profile[newx-1][newy-2]=value;mask[newx-1][newy-2]=1.0f; profile[newx-1][newy-1]=value; mask[newx-1][newy-1]=1.0f;
				for(int j=0;j<(newy-2);j++){profile[newx][j]=newprofile[j+2]; mask[newx][j]=newmask[j+2];}
				profile[newx][newy-2]=value; mask[newx][newy-2]=1.0f; profile[newx][newy-1]=value; mask[newx][newy-1]=1.0f;
			} else {
				//here we just shift down
				profile[newx-1][newy]=value;mask[newx-1][newy]=1.0f; profile[newx-1][newy+1]=value; mask[newx-1][newy+1]=1.0f;
				for(int j=(newy+2);j<maxsize;j++){profile[newx-1][j]=newprofile[j-2]; mask[newx-1][j]=newmask[j-2];}
				profile[newx][newy]=value; mask[newx][newy]=1.0f; profile[newx][newy+1]=value; mask[newx][newy+1]=1.0f;
				for(int j=(newy+2);j<maxsize;j++){profile[newx][j]=newprofile[j-2]; mask[newx][j]=newmask[j-2];}
			}
		}
		profile=interpolation.rotate_image(profile,-angle,center,center);
		mask=interpolation.rotate_image(mask,-angle,center,center);
		renorm_mask();
		//int temparea2=0;
		//for(int i=0;i<maxsize;i++){for(int j=0;j<maxsize;j++){if(mask[i][j]>0.5f) temparea2++;}};
		//IJ.log(""+(temparea2-temparea));
		return true;
	}

	public float remove_membrane(int x,int y,float value){
		//int temparea=0;
		//for(int i=0;i<maxsize;i++){for(int j=0;j<maxsize;j++){if(mask[i][j]>0.5f){temparea++;}}};
		//endocytosis is always assumed to consume 1 pixel
		float angle=(float)random.unidev(2.0*Math.PI,0.0);
		float[] tempvec=rotatevector((float)x,(float)y,-angle);
		int newx=(int)tempvec[0]; int newy=(int)tempvec[1];
		if(newx<1 || newx>=(maxsize-1) || newy<1 || newy>=(maxsize-1)){return value;}
		profile=interpolation.rotate_image(profile,angle,center,center);
		mask=interpolation.rotate_image(mask,angle,center,center);
		//int dither=random.jrnd.nextBoolean()?1:0;
		int dither=1;
		if(newy>center) dither=0;
		int perimeter=8;
		float excess=profile[newx][newy]-value;
		excess/=(float)perimeter;
		float min=0.0f;
		float retval=0.0f;
		//redistribute the excess protein around the endocytic site
		profile[newx-1][newy-1]+=excess; min=Math.min(profile[newx-1][newy-1],min);
		profile[newx-1][newy]+=excess; min=Math.min(profile[newx-1][newy],min);
		profile[newx-1][newy+1]+=excess; min=Math.min(profile[newx-1][newy+1],min);
		profile[newx][newy-1]+=excess; min=Math.min(profile[newx][newy-1],min);
		profile[newx][newy+1]+=excess; min=Math.min(profile[newx][newy+1],min);
		profile[newx+1][newy-1]+=excess; min=Math.min(profile[newx+1][newy-1],min);
		profile[newx+1][newy]+=excess; min=Math.min(profile[newx+1][newy],min);
		profile[newx+1][newy+1]+=excess; min=Math.min(profile[newx+1][newy+1],min);
		if(min<0.0){
			//correct for negative concentrations
			profile[newx-1][newy-1]-=min;
			profile[newx-1][newy]-=min;
			profile[newx-1][newy+1]-=min;
			profile[newx][newy-1]-=min;
			profile[newx][newy+1]-=min;
			profile[newx+1][newy-1]-=min;
			profile[newx+1][newy]-=min;
			profile[newx+1][newy+1]-=min;
			retval=(float)perimeter*(-min);
		}
		//here is a vertical deletion
		float[] newprofile=profile[newx].clone();
		float[] newmask=mask[newx].clone();
		//if dither is 1, we shift down
		if(dither==1){
			//we shift down
			for(int j=1;j<newy;j++){
				profile[newx][j+1]=newprofile[j]; mask[newx][j+1]=newmask[j];
			}
		} else {
			//we shift up
			for(int j=(newy+1);j<maxsize;j++){
				profile[newx][j-1]=newprofile[j]; mask[newx][j-1]=newmask[j];
			}
		}
		if(dither==1){mask[newx][0]=0.0f; profile[newx][0]=0.0f;}
		else{mask[newx][maxsize-1]=0.0f; profile[newx][maxsize-1]=0.0f;}
		profile=interpolation.rotate_image(profile,-angle,center,center);
		mask=interpolation.rotate_image(mask,-angle,center,center);
		renorm_mask();
		//int temparea2=0;
		//for(int i=0;i<maxsize;i++){for(int j=0;j<maxsize;j++){if(mask[i][j]>0.5f) temparea2++;}};
		//IJ.log(""+(temparea2-temparea));
		return retval;
	}

	public float[] rotatevector(float x,float y,float angle){
		float tempx=x-(float)center;
		float tempy=y-(float)center;
		float newx=(float)(Math.cos(-angle)*(double)tempx-Math.sin(-angle)*(double)tempy+(double)center);
		float newy=(float)(Math.sin(-angle)*(double)tempx+Math.cos(-angle)*(double)tempy+(double)center);
		float[] temp={newx,newy};
		return temp;
	}

	public float[][] arrcopy(float[][] input){
		float[][] output=new float[input.length][];
		for(int i=0;i<input.length;i++){
			output[i]=input[i].clone();
		}
		return output;
	}

	public float[][] bool2float(boolean[][] input){
		float[][] output=new float[input.length][input[0].length];
		for(int i=0;i<input.length;i++){
			for(int j=0;j<input[0].length;j++){if(input[i][j]) output[i][j]=1.0f;}
		}
		return output;
	}

	public boolean[][] float2bool(float[][] input){
		boolean[][] output=new boolean[input.length][input[0].length];
		for(int i=0;i<input.length;i++){
			for(int j=0;j<input[0].length;j++){output[i][j]=(input[i][j]>0.5f);}
		}
		return output;
	}

	public void renorm_mask(){
		for(int i=0;i<maxsize;i++){
			for(int j=0;j<maxsize;j++){
				if(mask[i][j]>0.5f) mask[i][j]=1.0f;
				else mask[i][j]=0.0f;
			}
		}
	}

}
