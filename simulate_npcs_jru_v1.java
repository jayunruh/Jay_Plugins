/*******************************************************************************
 * Copyright (c) 2018 Jay Unruh, Stowers Institute for Medical Research.
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
import jguis.*;
import jalgs.jfit.*;
import jalgs.jsim.*;
import jalgs.jseg.*;
import ij.text.*;

public class simulate_npcs_jru_v1 implements PlugIn {
	//this plugin simulates random points on the surface of a sphere convolved with a 3D psf
	//the image is 128 x 128 x 128
	public rngs random;

	public void run(String arg) {
		GenericDialog gd=new GenericDialog("Options");
		gd.addNumericField("Sphere_Radius (nm)",1000.0,5,15,null);
		gd.addNumericField("Min_npc_dist (nm)",0.0,5,15,null);
		gd.addNumericField("Number_of_points",100,0);
		gd.addNumericField("PSF_FWHM (nm)",100.0,5,15,null);
		gd.addNumericField("PSF_Z_FWHM (nm)",300.0,5,15,null);
		gd.addNumericField("Pixel_Size (nm)",40.0,5,15,null);
		gd.addNumericField("Max_Intensity (photons)",100.0,5,15,null);
		gd.addCheckbox("Add_Noise",false);
		gd.addNumericField("Read_Noise_Stdev",40.0f,5,15,null);
		gd.addNumericField("Gain (mutiplier)",50.0f,5,15,null);
		gd.addCheckbox("Add_Spb (2nd channel)",false);
		gd.addNumericField("Spb_Max_Intensity (photons)",100.0,5,15,null);
		gd.addNumericField("Spb_Separation (nm)",180.0f,5,15,null);
		gd.addCheckbox("Spb_on_top",false);
		gd.showDialog(); if(gd.wasCanceled()) return;
		float srad=(float)gd.getNextNumber();
		float minnpcdist=(float)gd.getNextNumber();
		int npts=(int)gd.getNextNumber();
		float xystdev=((float)gd.getNextNumber())/2.35f;
		float zstdev=((float)gd.getNextNumber())/2.35f;
		float psize=(float)gd.getNextNumber();
		float amp=(float)gd.getNextNumber();
		boolean addnoise=gd.getNextBoolean();
		float readstdev=(float)gd.getNextNumber();
		float gain=(float)gd.getNextNumber();
		boolean addspb=gd.getNextBoolean();
		float spbamp=(float)gd.getNextNumber();
		float spbsep=(float)gd.getNextNumber();
		boolean spbtop=gd.getNextBoolean();
		boolean table=true;
		random=new rngs();
		float[][] spherepts=makeSphere(npts,srad,64.0f*psize,64.0f*psize,64.0f*psize,minnpcdist);
		plotPoints("Npcs",spherepts,psize,table);
		Object[] stack=new Object[128];
		for(int i=0;i<128;i++) stack[i]=new float[128*128];
		drawPoints(stack,spherepts,psize,psize,xystdev,zstdev,128,128,amp);
		if(addnoise) addNoise(stack,readstdev,gain);
		if(addspb){
			Object[] spbstack=new Object[128];
			for(int i=0;i<128;i++) spbstack[i]=new float[128*128];
			float[][] spbpts=null;
			if(!spbtop) spbpts=makeSpb(spbsep,srad,64.0f*psize,64.0f*psize,64.0f*psize);
			else spbpts=makeSpbTop(spbsep,srad,64.0f*psize,64.0f*psize,64.0f*psize);
			plotPoints("SPB",spbpts,psize,table);
			drawPoints(spbstack,spbpts,psize,psize,xystdev,zstdev,128,128,spbamp);
			if(addnoise) addNoise(spbstack,readstdev,gain);
			//now interleave the stacks
			Object[] stack2=new Object[256];
			for(int i=0;i<128;i++){
				stack2[2*i]=stack[i];
				stack2[2*i+1]=spbstack[i];
			}
			ImagePlus imp=jutils.create_hyperstack("Sim NPCs",jutils.array2stack(stack2,128,128),1,128,2,true,null);
			jutils.set_psize(imp,psize/1000.0f);
			imp.show();
		} else {
			ImagePlus imp=new ImagePlus("Sim NPCs",jutils.array2stack(stack,128,128));
			jutils.set_psize(imp,psize/1000.0f);
			imp.show();
		}
	}

	public void addNoise(Object[] stack,float readstdev,float gain){
		for(int i=0;i<stack.length;i++){
			float[] pix=(float[])stack[i];
			for(int j=0;j<pix.length;j++){
				int photons=random.poidev((double)pix[j]);
				float measured=gain*(float)photons;
				measured+=random.gasdev(0.0,readstdev);
				pix[j]=measured;
			}
		}
	}

	public void plotPoints(String name,float[][] coords,float psize,boolean table){
		float[][] xpts=new float[coords.length][1];
		float[][] ypts=new float[coords.length][1];
		float[][] zpts=new float[coords.length][1];
		TextWindow tw=null;
		if(table) tw=new TextWindow("Simulated Coordinates: "+name,"point\tx\ty\tz","",400,200);
		for(int i=0;i<coords.length;i++){
			xpts[i][0]=coords[i][0];
			xpts[i][0]/=psize; //convert to pixels
			ypts[i][0]=coords[i][1];
			ypts[i][0]/=psize;
			zpts[i][0]=coords[i][2];
			zpts[i][0]/=psize;
			if(table) tw.append(""+(i+1)+"\t"+xpts[i][0]+"\t"+ypts[i][0]+"\t"+zpts[i][0]);
		}
		Traj3D traj=new Traj3D("x","y","z",xpts,ypts,zpts,null);
		int[] shapes=traj.getShapes();
		for(int i=0;i<coords.length;i++){
			shapes[i]=1;
		}
		new PlotWindow3D("Sphere Positions: "+name,traj).draw();
		//now create a pseudoimage with points masked
		float[][] stack=new float[128][128*128];
		for(int i=0;i<coords.length;i++){
			int xpos=(int)Math.round(xpts[i][0]);
			if(xpos<0) xpos=0; if(xpos>=128) xpos=128;
			int ypos=(int)Math.round(ypts[i][0]);
			if(ypos<0) ypos=0; if(ypos>=128) ypos=128;
			int zpos=(int)Math.round(zpts[i][0]);
			if(xpos<0) xpos=0; if(xpos>=128) xpos=128;
			stack[zpos][xpos+ypos*128]=1.0f;
		}
		new ImagePlus("Sphere Positions Image: "+name,jutils.array2stack(stack,128,128)).show();
	}

	public float[][] makeSphere(int nparticles,float radius,float xc, float yc, float zc,float minnpcdist){
		float[][] coords=new float[nparticles][3];
		for(int i=0;i<nparticles;i++){
			double[] dcoords=random.random_sphere(radius);
			coords[i][0]=(float)dcoords[0]+xc;
			coords[i][1]=(float)dcoords[1]+yc;
			coords[i][2]=(float)dcoords[2]+zc;
			if(minnpcdist>0.0f && i>0){
				while(mindist(coords,coords[i],i)<minnpcdist){
					dcoords=random.random_sphere(radius);
					coords[i][0]=(float)dcoords[0]+xc;
					coords[i][1]=(float)dcoords[1]+yc;
					coords[i][2]=(float)dcoords[2]+zc;
				}
			}
		}
		return coords;
	}

	public float mindist(float[][] coords,float[] query,int nparticles){
		float dist=(float)Math.sqrt((query[0]-coords[0][0])*(query[0]-coords[0][0])+(query[1]-coords[0][1])*(query[1]-coords[0][1])+(query[2]-coords[0][2])*(query[0]-coords[0][2]));
		for(int i=1;i<nparticles;i++){
			float tdist=(float)Math.sqrt((query[0]-coords[i][0])*(query[0]-coords[i][0])+(query[1]-coords[i][1])*(query[1]-coords[i][1])+(query[2]-coords[i][2])*(query[0]-coords[i][2]));
			if(tdist<dist) dist=tdist;
		}
		return dist;
	}

	public float[][] makeSpb(float bridgedist,float radius,float xc,float yc,float zc){
		//this makes a mother and daughter spb with bridgedist between them
		double[] dcoords=random.random_sphere(radius);
		float[][] coords=new float[2][3];
		coords[0][0]=(float)dcoords[0]+xc;
		coords[0][1]=(float)dcoords[1]+yc;
		coords[0][2]=(float)dcoords[2]+zc;
		//now that we have the mother coordinate, we need to put the daughter in a random cone
		//centered on the mother vector
		//need to convert the bridge distance to an angle distance in radians
		double angledist=Math.atan(bridgedist/radius);
		//IJ.log(""+angledist);
		double[] normvec={dcoords[0]/(double)radius,dcoords[1]/(double)radius,dcoords[2]/(double)radius};
		//calculate a random cone vector
		/*double[] conevec=random_cone_vector(normvec,angledist);
		coords[1][0]=-(float)conevec[0]*radius+xc;
		coords[1][1]=-(float)conevec[1]*radius+yc;
		coords[1][2]=-(float)conevec[2]*radius+zc;*/
		//alternatively, can always rotate about the vector z axis plane
		/*double[] rotaxisvector=new double[3];
		rotaxisvector[0]=normvec[1]; rotaxisvector[1]=-normvec[0];
		double rval=Math.sqrt(rotaxisvector[0]*rotaxisvector[0]+rotaxisvector[1]*rotaxisvector[1]); //normalize the dot product
		rotaxisvector[0]/=rval;
		rotaxisvector[1]/=rval;*/
		//or rotate about a random orthogonal vector
		double[] randomvec=random.random_sphere(1.0);
		double[] rotaxisvector=measure_object.crossProd(randomvec,normvec);
		double rval=Math.sqrt(rotaxisvector[0]*rotaxisvector[0]+rotaxisvector[1]*rotaxisvector[1]+rotaxisvector[2]*rotaxisvector[2]);
		rotaxisvector[0]/=rval;
		rotaxisvector[1]/=rval;
		rotaxisvector[2]/=rval;
		rotate_vector(normvec,rotaxisvector,angledist);
		coords[1][0]=(float)normvec[0]*radius+xc;
		coords[1][1]=(float)normvec[1]*radius+yc;
		coords[1][2]=(float)normvec[2]*radius+zc;
		float spbdist=(float)Math.sqrt((coords[1][0]-coords[0][0])*(coords[1][0]-coords[0][0])+(coords[1][1]-coords[0][1])*(coords[1][1]-coords[0][1])+(coords[1][2]-coords[0][2])*(coords[1][2]-coords[0][2]));
		//IJ.log("simulated spb dist "+spbdist);
		return coords;
	}

	public float[][] makeSpbTop(float bridgedist,float radius,float xc,float yc,float zc){
		//this makes a mother and daughter spb with bridgedist between them
		//in this version we put the spb on the top of the nucleus pointing "up" the y axis
		//double[] dcoords=random.random_sphere(radius);
		double[] dcoords={0.0,0.0,(double)radius};
		float[][] coords=new float[2][3];
		coords[0][0]=(float)dcoords[0]+xc;
		coords[0][1]=(float)dcoords[1]+yc;
		coords[0][2]=(float)dcoords[2]+zc;
		//now that we have the mother coordinate, we need to put the daughter in a random cone
		//centered on the mother vector
		//need to convert the bridge distance to an angle distance in radians
		double angledist=Math.atan(bridgedist/radius);
		//IJ.log(""+angledist);
		double[] normvec={dcoords[0]/(double)radius,dcoords[1]/(double)radius,dcoords[2]/(double)radius};
		//in this case we rotate about the x axis
		rotate_vector(normvec,new double[]{1.0,0.0,0.0},angledist);
		coords[1][0]=(float)normvec[0]*radius+xc;
		coords[1][1]=(float)normvec[1]*radius+yc;
		coords[1][2]=(float)normvec[2]*radius+zc;
		float spbdist=(float)Math.sqrt((coords[1][0]-coords[0][0])*(coords[1][0]-coords[0][0])+(coords[1][1]-coords[0][1])*(coords[1][1]-coords[0][1])+(coords[1][2]-coords[0][2])*(coords[1][2]-coords[0][2]));
		//IJ.log("simulated spb dist "+spbdist);
		return coords;
	}

	public double[] random_cone_vector(double[] center_vector,double theta){
		//here theta is in radians: the half arc size of the cone
		//first calculate the random cone vector about the z axis
		double[] randvector=new double[3];
		randvector[2]=Math.cos(theta);
		double axialr=Math.sin(theta);
		double angle=random.unidev(2.0*Math.PI,0.0); //this is the radial angle
		randvector[0]=axialr*Math.cos(angle);
		randvector[1]=axialr*Math.sin(angle);
		//now find the vector orthogonal to the center_vector-zaxis plane (dot product)--have to normalize it
		double[] rotaxisvector=new double[3];
		rotaxisvector[0]=center_vector[1];
		rotaxisvector[1]=-center_vector[0];
		double rval=Math.sqrt(rotaxisvector[0]*rotaxisvector[0]+rotaxisvector[1]*rotaxisvector[1]);
		rotaxisvector[0]/=rval;
		rotaxisvector[1]/=rval;
		//now find the angle between the z axis and the center vector (acos of dot product)
		double rotangle=Math.acos(center_vector[2]);
		//rotate the cone vector about the orthogonal vector
		rotate_vector(randvector,rotaxisvector,rotangle);
		return randvector;
	}

	public void rotate_vector(double[] invector,double[] rot_axis,double theta){
		//this subroutine rotates a vector about rot_axis according to Rodrigues' formula
		double xtemp,ytemp,ztemp,xrot,yrot,zrot,cos_theta,sin_theta;
		cos_theta=Math.cos(theta);
		sin_theta=Math.sin(theta);
		xtemp=invector[0]; ytemp=invector[1]; ztemp=invector[2]; xrot=rot_axis[0]; yrot=rot_axis[1]; zrot=rot_axis[2];
		invector[0]=((1.0-cos_theta)*xrot*xrot+cos_theta)*xtemp+((1-cos_theta)*yrot*xrot-zrot*sin_theta)*ytemp+((1.0-cos_theta)*zrot*xrot+yrot*sin_theta)*ztemp;
		invector[1]=((1.0-cos_theta)*xrot*yrot+zrot*sin_theta)*xtemp+((1-cos_theta)*yrot*yrot+cos_theta)*ytemp+((1.0-cos_theta)*zrot*yrot-xrot*sin_theta)*ztemp;
		invector[2]=((1.0-cos_theta)*xrot*zrot-yrot*sin_theta)*xtemp+((1-cos_theta)*yrot*zrot+xrot*sin_theta)*ytemp+((1.0-cos_theta)*zrot*zrot+cos_theta)*ztemp;
		return;
	}

	public void drawPoints(Object[] stack,float[][] coords,float psize,float zsize,float xystdev,float zstdev,int width,int height,float amp){
		gausfunc gf=new gausfunc();
		float xystdevpix=xystdev/psize;
		float zstdevpix=zstdev/zsize;
		for(int i=0;i<coords.length;i++){
			float xpos=coords[i][0]/psize;
			float ypos=coords[i][1]/psize;
			float zpos=coords[i][2]/zsize;
			gf.draw_3D_func(stack,xpos,ypos,zpos,width,height,xystdevpix,zstdevpix,amp);
		}
		return;
	}

}
