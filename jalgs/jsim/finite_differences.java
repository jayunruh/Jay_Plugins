/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

import jalgs.algutils;

public class finite_differences{
	public int boundarycond; // these are periodic, reflective, infinite, and
	// zero
	public int dsteps; // number of substeps for diffusion
	public int width,height;
	public float[][] D2;
	public float pixsize,T;
	public boolean oneD;

	// here we have methods to simulate diffusion with finite differences in one
	// or two dimensions

	public finite_differences(int boundarycond,int dsteps,float pixsize,float T,float D,int width,int height){
		this.boundarycond=boundarycond;
		this.dsteps=dsteps;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		D2=new float[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				D2[i][j]=D*T/(pixsize*pixsize);
			}
		}
		if(height==1)
			oneD=true;
	}

	public finite_differences(int boundarycond,int dsteps,float pixsize,float T,float[][] D,int width,int height){
		this.boundarycond=boundarycond;
		this.dsteps=dsteps;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		D2=new float[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				D2[i][j]=D[i][j]*T/(pixsize*pixsize);
			}
		}
		if(height==1)
			oneD=true;
	}

	public finite_differences(int boundarycond,int dsteps,float pixsize,float T,float[] D,int width,int height){
		this.boundarycond=boundarycond;
		this.dsteps=dsteps;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		D2=new float[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				D2[i][j]=D[j+i*width]*T/(pixsize*pixsize);
			}
		}
		if(height==1)
			oneD=true;
	}

	public Object steptime(Object profile){
		if(oneD)
			return handle_1D_diffusion((float[])profile);
		else{
			if(profile instanceof float[][]){
				return handle_2D_diffusion((float[][])profile);
			}else{
				return handle_2D_diffusion((float[])profile);
			}
		}
	}

	public float[][] handle_2D_diffusion(float[][] profile){
		float[][] temp=algutils.clone_multidim_array(profile);
		for(int i=0;i<dsteps;i++){
			diff2D(temp,1.0f/dsteps);
		}
		return temp;
	}

	public float[] handle_2D_diffusion(float[] profile){
		float[] temp=profile.clone();
		for(int i=0;i<dsteps;i++){
			diff2D(temp,1.0f/dsteps);
		}
		return temp;
	}

	private void diff2D(float[][] profile,float multiplier){
		float[][] temp=algutils.clone_multidim_array(profile);
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++)
				temp[i][j]*=D2[i][j]*multiplier;
		}
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				float derx=temp[i][j+1]-2.0f*temp[i][j]+temp[i][j-1];
				float dery=temp[i+1][j]-2.0f*temp[i][j]+temp[i-1][j];
				profile[i][j]+=(derx+dery);
			}
		}
		// these are periodic, reflective, fixed, and zero
		if(boundarycond==0){
			// upper left corder
			float derx=temp[0][1]-2.0f*temp[0][0]+temp[0][width-1];
			float dery=temp[1][0]-2.0f*temp[0][0]+temp[height-1][0];
			profile[0][0]+=(derx+dery);
			// lower left corner
			derx=temp[height-1][1]-2.0f*temp[height-1][0]+temp[height-1][width-1];
			dery=temp[0][0]-2.0f*temp[height-1][0]+temp[height-2][0];
			profile[height-1][0]+=(derx+dery);
			// upper right corner
			derx=temp[0][0]-2.0f*temp[0][width-1]+temp[0][width-2];
			dery=temp[1][width-1]-2.0f*temp[0][width-1]+temp[height-1][width-1];
			profile[0][width-1]+=(derx+dery);
			// lower right corner
			derx=temp[height-1][0]-2.0f*temp[height-1][width-1]+temp[height-1][width-2];
			dery=temp[0][width-1]-2.0f*temp[height-1][width-1]+temp[height-2][width-1];
			profile[height-1][width-1]+=(derx+dery);
			// right and left edges
			for(int i=1;i<(height-1);i++){
				derx=temp[i][1]-2.0f*temp[i][0]+temp[i][width-1];
				dery=temp[i+1][0]-2.0f*temp[i][0]+temp[i-1][0];
				profile[i][0]+=(derx+dery);
				derx=temp[i][0]-2.0f*temp[i][width-1]+temp[i][width-2];
				dery=temp[i+1][width-1]-2.0f*temp[i][width-1]+temp[i-1][width-1];
				profile[i][width-1]+=(derx+dery);
			}
			// top and bottom edges
			for(int i=1;i<(width-1);i++){
				derx=temp[0][i+1]-2.0f*temp[0][i]+temp[0][i-1];
				dery=temp[1][i]-2.0f*temp[0][i]+temp[height-1][i];
				profile[0][i]+=(derx+dery);
				derx=temp[height-1][i+1]-2.0f*temp[height-1][i]+temp[height-1][i-1];
				dery=temp[0][i]-2.0f*temp[height-1][i]+temp[height-2][i];
				profile[height-1][i]+=(derx+dery);
			}
			return;
		}
		if(boundarycond==1){
			// upper left corder
			float derx=2.0f*(temp[0][1]-temp[0][0]);
			float dery=2.0f*(temp[1][0]-temp[0][0]);
			profile[0][0]+=(derx+dery);
			// lower left corner
			derx=2.0f*(temp[height-1][1]-temp[height-1][0]);
			dery=2.0f*(temp[height-2][0]-temp[height-1][0]);
			profile[height-1][0]+=(derx+dery);
			// upper right corner
			derx=2.0f*(temp[0][width-2]-temp[0][width-1]);
			dery=2.0f*(temp[1][width-1]-temp[0][width-1]);
			profile[0][width-1]+=(derx+dery);
			// lower right corner
			derx=2.0f*(temp[height-1][width-2]-temp[height-1][width-1]);
			dery=2.0f*(temp[height-2][width-1]-temp[height-1][width-1]);
			profile[height-1][width-1]+=(derx+dery);
			// right and left edges
			for(int i=1;i<(height-1);i++){
				derx=2.0f*(temp[i][1]-temp[i][0]);
				dery=temp[i+1][0]-2.0f*temp[i][0]+temp[i-1][0];
				profile[i][0]+=(derx+dery);
				derx=2.0f*(temp[i][width-2]-temp[i][width-1]);
				dery=temp[i+1][width-1]-2.0f*temp[i][width-1]+temp[i-1][width-1];
				profile[i][width-1]+=(derx+dery);
			}
			// top and bottom edges
			for(int i=1;i<(width-1);i++){
				derx=temp[0][i+1]-2.0f*temp[0][i]+temp[0][i-1];
				dery=2.0f*(temp[1][i]-temp[0][i]);
				profile[0][i]+=(derx+dery);
				derx=temp[height-1][i+1]-2.0f*temp[height-1][i]+temp[height-1][i-1];
				dery=2.0f*(temp[height-2][i]-temp[height-1][i]);
				profile[height-1][i]+=(derx+dery);
			}
			return;
		}
		if(boundarycond==2){
			return;
		}
		if(boundarycond==3){
			// upper left corder
			float derx=temp[0][1]-2.0f*temp[0][0];
			float dery=temp[1][0]-2.0f*temp[0][0];
			profile[0][0]+=(derx+dery);
			// lower left corner
			derx=temp[height-1][1]-2.0f*temp[height-1][0];
			dery=temp[0][0]-2.0f*temp[height-1][0];
			profile[height-1][0]+=(derx+dery);
			// upper right corner
			derx=temp[0][0]-2.0f*temp[0][width-1];
			dery=temp[0][width-1]-2.0f*temp[0][width-1];
			profile[0][width-1]+=(derx+dery);
			// lower right corner
			derx=temp[height-1][0]-2.0f*temp[height-1][width-1];
			dery=temp[0][width-1]-2.0f*temp[height-1][width-1];
			profile[height-1][width-1]+=(derx+dery);
			// right and left edges
			for(int i=1;i<(height-1);i++){
				derx=temp[i][1]-2.0f*temp[i][0];
				dery=temp[i+1][0]-2.0f*temp[i][0];
				profile[i][0]+=(derx+dery);
				derx=temp[i][0]-2.0f*temp[i][width-1];
				dery=temp[i+1][width-1]-2.0f*temp[i][width-1];
				profile[i][width-1]+=(derx+dery);
			}
			// top and bottom edges
			for(int i=1;i<(width-1);i++){
				derx=temp[0][i+1]-2.0f*temp[0][i];
				dery=temp[1][i]-2.0f*temp[0][i];
				profile[0][i]+=(derx+dery);
				derx=temp[height-1][i+1]-2.0f*temp[height-1][i];
				dery=temp[0][i]-2.0f*temp[height-1][i];
				profile[height-1][i]+=(derx+dery);
			}
			return;
		}
	}

	private void diff2D(float[] profile,float multiplier){
		float[][] temp=algutils.array2multidim(profile,width,height);
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				float derx=calc_update(D2[i][j-1],D2[i][j],D2[i][j+1],temp[i][j-1],temp[i][j],temp[i][j+1])*multiplier;
				float dery=calc_update(D2[i-1][j],D2[i][j],D2[i+1][j],temp[i-1][j],temp[i][j],temp[i+1][j])*multiplier;
				profile[j+i*width]+=(derx+dery);
			}
		}
		// these are periodic, reflective, fixed, and zero
		if(boundarycond==0){
			// upper left corder
			float derx=calc_update(D2[0][width-1],D2[0][0],D2[0][1],temp[0][width-1],temp[0][0],temp[0][1])*multiplier;
			float dery=calc_update(D2[height-1][0],D2[0][0],D2[1][0],temp[height-1][0],temp[0][0],temp[1][0])*multiplier;
			profile[0]+=(derx+dery);
			// lower left corner
			derx=calc_update(D2[height-1][width-1],D2[height-1][0],D2[height-1][1],temp[height-1][width-1],temp[height-1][0],temp[height-1][1])*multiplier;
			dery=calc_update(D2[height-1][0],D2[0][0],D2[1][0],temp[height-1][0],temp[0][0],temp[1][0])*multiplier;
			profile[width*(height-1)]+=(derx+dery);
			// upper right corner
			derx=calc_update(D2[0][width-2],D2[0][width-1],D2[0][0],temp[0][width-2],temp[0][width-1],temp[0][0])*multiplier;
			dery=calc_update(D2[height-1][width-1],D2[0][width-1],D2[1][width-1],temp[height-1][width-1],temp[0][width-1],temp[1][width-1])*multiplier;
			profile[width-1]+=(derx+dery);
			// lower right corner
			derx=calc_update(D2[height-1][width-2],D2[height-1][width-1],D2[height-1][0],temp[height-1][width-2],temp[height-1][width-1],temp[height-1][0])*multiplier;
			dery=calc_update(D2[height-2][width-1],D2[height-1][width-1],D2[0][width-1],temp[height-2][width-1],temp[height-1][width-1],temp[0][width-1])*multiplier;
			profile[width*(height-1)+width-1]+=(derx+dery);
			// right and left edges
			for(int i=1;i<(height-1);i++){
				derx=calc_update(D2[i][width-1],D2[i][0],D2[i][1],temp[i][width-1],temp[i][0],temp[i][1])*multiplier;
				dery=calc_update(D2[i-1][0],D2[i][0],D2[i+1][0],temp[i-1][0],temp[i][0],temp[i+1][0])*multiplier;
				profile[i*width]+=(derx+dery);
				derx=calc_update(D2[i][width-2],D2[i][width-1],D2[i][0],temp[i][width-2],temp[i][width-1],temp[i][0])*multiplier;
				dery=calc_update(D2[i-1][width-1],D2[i][width-1],D2[i+1][width-1],temp[i-1][width-1],temp[i][width-1],temp[i+1][width-1])*multiplier;
				profile[i*width+width-1]+=(derx+dery);
			}
			// top and bottom edges
			for(int i=1;i<(width-1);i++){
				derx=calc_update(D2[0][i-1],D2[0][i],D2[0][i+1],temp[0][i-1],temp[0][i],temp[0][i+1])*multiplier;
				dery=calc_update(D2[height-1][i],D2[0][i],D2[1][i],temp[height-1][i],temp[0][i],temp[1][i])*multiplier;
				profile[i]+=(derx+dery);
				derx=calc_update(D2[height-1][i-1],D2[height-1][i],D2[height-1][i+1],temp[height-1][i-1],temp[height-1][i],temp[height-1][i+1])*multiplier;
				dery=calc_update(D2[height-2][i],D2[height-1][i],D2[0][0],temp[height-2][i],temp[height-1][i],temp[0][i])*multiplier;
				profile[(height-1)*width+i]+=(derx+dery);
			}
			return;
		}
		if(boundarycond==1){
			// upper left corder
			float derx=2.0f*(temp[0][1]-temp[0][0]);
			float dery=2.0f*(temp[1][0]-temp[0][0]);
			profile[0]+=(derx+dery);
			// lower left corner
			derx=2.0f*(temp[height-1][1]-temp[height-1][0]);
			dery=2.0f*(temp[height-2][0]-temp[height-1][0]);
			profile[(height-1)*width]+=(derx+dery);
			// upper right corner
			derx=2.0f*(temp[0][width-2]-temp[0][width-1]);
			dery=2.0f*(temp[1][width-1]-temp[0][width-1]);
			profile[width-1]+=(derx+dery);
			// lower right corner
			derx=2.0f*(temp[height-1][width-2]-temp[height-1][width-1]);
			dery=2.0f*(temp[height-2][width-1]-temp[height-1][width-1]);
			profile[height*width-1]+=(derx+dery);
			// right and left edges
			for(int i=1;i<(height-1);i++){
				derx=2.0f*(temp[i][1]-temp[i][0]);
				dery=temp[i+1][0]-2.0f*temp[i][0]+temp[i-1][0];
				profile[i*width]+=(derx+dery);
				derx=2.0f*(temp[i][width-2]-temp[i][width-1]);
				dery=temp[i+1][width-1]-2.0f*temp[i][width-1]+temp[i-1][width-1];
				profile[i*width+width-1]+=(derx+dery);
			}
			// top and bottom edges
			for(int i=1;i<(width-1);i++){
				derx=temp[0][i+1]-2.0f*temp[0][i]+temp[0][i-1];
				dery=2.0f*(temp[1][i]-temp[0][i]);
				profile[i]+=(derx+dery);
				derx=temp[height-1][i+1]-2.0f*temp[height-1][i]+temp[height-1][i-1];
				dery=2.0f*(temp[height-2][i]-temp[height-1][i]);
				profile[(height-1)*width+i]+=(derx+dery);
			}
			return;
		}
		if(boundarycond==2){
			return;
		}
		if(boundarycond==3){
			// upper left corder
			float derx=temp[0][1]-2.0f*temp[0][0];
			float dery=temp[1][0]-2.0f*temp[0][0];
			profile[0]+=(derx+dery);
			// lower left corner
			derx=temp[height-1][1]-2.0f*temp[height-1][0];
			dery=temp[0][0]-2.0f*temp[height-1][0];
			profile[(height-1)*width]+=(derx+dery);
			// upper right corner
			derx=temp[0][0]-2.0f*temp[0][width-1];
			dery=temp[1][width-1]-2.0f*temp[0][width-1];
			profile[width-1]+=(derx+dery);
			// lower right corner
			derx=temp[height-1][0]-2.0f*temp[height-1][width-1];
			dery=temp[0][width-1]-2.0f*temp[height-1][width-1];
			profile[(height-1)*width+width-1]+=(derx+dery);
			// right and left edges
			for(int i=1;i<(height-1);i++){
				derx=temp[i][1]-2.0f*temp[i][0];
				dery=temp[i+1][0]-2.0f*temp[i][0];
				profile[i*width]+=(derx+dery);
				derx=temp[i][0]-2.0f*temp[i][width-1];
				dery=temp[i+1][width-1]-2.0f*temp[i][width-1];
				profile[i*width+width-1]+=(derx+dery);
			}
			// top and bottom edges
			for(int i=1;i<(width-1);i++){
				derx=temp[0][i+1]-2.0f*temp[0][i];
				dery=temp[1][i]-2.0f*temp[0][i];
				profile[i]+=(derx+dery);
				derx=temp[height-1][i+1]-2.0f*temp[height-1][i];
				dery=temp[0][i]-2.0f*temp[height-1][i];
				profile[(height-1)*width+i]+=(derx+dery);
			}
			return;
		}
	}

	private void diff1D(float[] profile,float multiplier){
		float[] temp=profile.clone();
		for(int i=0;i<temp.length;i++)
			temp[i]*=D2[0][i]*multiplier;
		for(int i=1;i<(width-1);i++){
			float der=calc_update(D2[0][i-1],D2[0][i],D2[0][i+1],temp[i-1],temp[i],temp[i+1])*multiplier;
			profile[i]+=der;
		}
		// these are periodic, reflective, fixed, and zero
		if(boundarycond==0){
			float der=calc_update(D2[0][width-1],D2[0][0],D2[0][1],temp[width-1],temp[0],temp[1])*multiplier;
			profile[0]+=der;
			der=calc_update(D2[0][width-2],D2[0][width-1],D2[0][0],temp[width-2],temp[width-1],temp[0])*multiplier;
			profile[width-1]+=der;
			return;
		}
		if(boundarycond==1){
			float der=calc_update(D2[0][1],D2[0][0],D2[0][1],temp[1],temp[0],temp[1])*multiplier;
			profile[0]+=der;
			der=calc_update(D2[0][width-2],D2[0][width-1],D2[0][width-2],temp[width-2],temp[width-1],temp[width-2])*multiplier;
			profile[width-1]+=der;
			return;
		}
		if(boundarycond==2){
			return;
		}
		if(boundarycond==3){
			float der=calc_update(0.0f,D2[0][0],D2[0][1],0.0f,temp[0],temp[1])*multiplier;
			profile[0]+=der;
			der=calc_update(D2[0][width-2],D2[0][width-1],0.0f,temp[width-2],temp[width-1],0.0f)*multiplier;
			profile[width-1]+=der;
			return;
		}
	}

	private float calc_update(float d1,float d2,float d3,float c1,float c2,float c3){
		// return 0.5f*((d3+d2)*c3-(d3+d2+d2+d1)*c2+(d2+d1)*c1);
		return c3-2.0f*c2+c1+0.25f*(d3-d1)*(c3-c1);
	}

	public float[] handle_1D_diffusion(float[] profile){
		float[] temp=profile.clone();
		for(int i=0;i<dsteps;i++){
			diff1D(temp,1.0f/dsteps);
		}
		return temp;
	}
}
