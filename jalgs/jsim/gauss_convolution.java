/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

import jalgs.jfit.gausfunc;

public class gauss_convolution{
	public int boundarycond; // these are periodic, reflective, infinite, and
	// zero
	public int width,height;
	public int[][] Dmap;
	public float[] D2,widths;
	public int[] sizes;
	public Object[] profiles;
	public float pixsize,T,widthmult;
	public boolean oneD;

	// here we have methods to simulate diffusion with finite differences in one
	// or two dimensions

	public gauss_convolution(int boundarycond,float pixsize,float T,float D,int width,int height){
		this.boundarycond=boundarycond;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		D2=new float[1];
		D2[0]=D*T/(pixsize*pixsize);
		Dmap=new int[height][width];
		if(height==1) oneD=true;
		widthmult=3.0f;
		init_profiles();
	}

	public gauss_convolution(int boundarycond,float pixsize,float T,float[][] D,int width,int height){
		this.boundarycond=boundarycond;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		init_D_lookup(D);
		if(height==1)
			oneD=true;
		widthmult=3.0f;
		init_profiles();
	}

	public gauss_convolution(int boundarycond,float pixsize,float T,float[] D,int[][] Dmap,int width,int height){
		// here we initialize with the Dmap already generated
		this.boundarycond=boundarycond;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		this.Dmap=Dmap;
		D2=new float[D.length];
		for(int i=0;i<D.length;i++)
			D2[i]=D[i]*T/(pixsize*pixsize);
		if(height==1)
			oneD=true;
		widthmult=3.0f;
		init_profiles();
	}

	public gauss_convolution(int boundarycond,float pixsize,float T,float[] D,int width,int height){
		this.boundarycond=boundarycond;
		this.pixsize=pixsize;
		this.T=T;
		this.width=width;
		this.height=height;
		init_D_lookup(D);
		if(height==1)
			oneD=true;
		widthmult=3.0f;
		init_profiles();
	}

	public void init_D_lookup(float[][] D){
		int nDs=1;
		float[] temp=new float[width*height];
		temp[0]=D[0][0];
		Dmap=new int[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if((j+i)>0){
					// search through the lookup for matching D values
					boolean found=false;
					for(int k=0;k<nDs;k++){
						if(D[i][j]==temp[k]){
							Dmap[i][j]=k;
							found=true;
							break;
						}
					}
					if(!found){
						Dmap[i][j]=nDs;
						temp[nDs]=D[i][j];
						nDs++;
					}
				}
			}
		}
		D2=new float[nDs];
		for(int i=0;i<nDs;i++){
			D2[i]=temp[i]*T/(pixsize*pixsize);
		}
	}

	public void init_D_lookup(float[] D){
		int nDs=1;
		float[] temp=new float[width*height];
		temp[0]=D[0];
		Dmap=new int[height][width];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				if((j+i)>0){
					// search through the lookup for matching D values
					boolean found=false;
					float currD=D[j+i*width];
					for(int k=0;k<nDs;k++){
						if(currD==temp[k]){
							Dmap[i][j]=k;
							found=true;
							break;
						}
					}
					if(!found){
						Dmap[i][j]=nDs;
						temp[nDs]=currD;
						nDs++;
					}
				}
			}
		}
		D2=new float[nDs];
		for(int i=0;i<nDs;i++){
			D2[i]=temp[i]*T/(pixsize*pixsize);
		}
	}

	public void init_profiles(){
		gausfunc gf=new gausfunc();
		profiles=new Object[D2.length];
		sizes=new int[D2.length];
		widths=new float[D2.length];
		for(int i=0;i<D2.length;i++){
			if(D2[i]<=0.0f){
				sizes[i]=0;
				profiles[i]=null;
			} else {
				widths[i]=(float)Math.sqrt(2.0f*D2[i]);
				sizes[i]=(int)(widthmult*widths[i]);
				if(sizes[i]<1) sizes[i]=1;
				if(sizes[i]>width) sizes[i]=width;
				if(!oneD && sizes[i]>height) sizes[i]=height;
				if(oneD){
					profiles[i]=gf.get_norm_func(-sizes[i],2*sizes[i]+1,1.0f,widths[i]);
				} else {
					profiles[i]=gf.get_2D_norm_func(-sizes[i],2*sizes[i]+1,1.0f,widths[i]);
				}
			}
		}
	}

	public int getmaxprofilesize(){
		int maxsize=0;
		for(int i=0;i<profiles.length;i++){
			if(sizes[i]>maxsize)
				maxsize=sizes[i];
		}
		return 2*maxsize+1;
	}

	public float[][] output_profiles(){
		int maxsize=getmaxprofilesize();
		float[][] output=new float[profiles.length][maxsize*maxsize];
		for(int i=0;i<profiles.length;i++){
			if(profiles[i]!=null&&sizes[i]>0){
				float[] profile=(float[])profiles[i];
				int tempsize=2*sizes[i]+1;
				for(int j=0;j<tempsize;j++){
					for(int k=0;k<tempsize;k++){
						output[i][k+j*tempsize]=profile[k+j*tempsize];
					}
				}
			}
		}
		return output;
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

	private float[] handle_1D_diffusion(float[] profile){
		float[] temp=new float[width];
		for(int i=0;i<width;i++){
			int index=Dmap[0][i];
			if(sizes[index]>0){
				int start=i-sizes[index];
				int end=i+sizes[index];
				if(boundarycond>1){
					if(start<0)
						start=0;
					if(start>=width)
						start=width-1;
					if(end<0)
						end=0;
					if(end>=width)
						end=width-1;
				}
				for(int j=start;j<=end;j++){
					int index2=j;
					if(index2<0){
						if(boundarycond==0)
							index2+=width;
						if(boundarycond==1)
							index2=-index2;
					}
					if(index2>=width){
						if(boundarycond==0)
							index2-=width;
						if(boundarycond==1)
							index2=width+width-index2;
					}
					temp[index2]+=profile[i]*((float[])profiles[index])[j-i+sizes[index]];
				}
			}
		}
		if(boundarycond==2){
			temp[0]=profile[0];
			temp[width-1]=profile[width-1];
		}
		return temp;
	}

	// need to modify this section to do one d blurring in each dimension
	// separately
	// this will speed up the calculation
	private float[][] handle_2D_diffusion(float[][] profile){
		float[][] temp=new float[height][width];
		for(int i=0;i<height;i++){
			for(int k=0;k<width;k++){
				int index=Dmap[i][k];
				if(sizes[index]>0){
					int startx=k-sizes[index];
					int endx=k+sizes[index];
					int starty=i-sizes[index];
					int endy=i+sizes[index];
					if(boundarycond>1){
						if(startx<0)
							startx=0;
						if(startx>=width)
							startx=width-1;
						if(endx<0)
							endx=0;
						if(endx>=width)
							endx=width-1;
						if(starty<0)
							starty=0;
						if(starty>=height)
							starty=height-1;
						if(endy<0)
							endy=0;
						if(endy>=height)
							endy=height-1;
					}
					for(int j=starty;j<=endy;j++){
						int index2=j;
						if(index2<0){
							if(boundarycond==0)
								index2+=height;
							if(boundarycond==1)
								index2=-index2;
						}
						if(index2>=height){
							if(boundarycond==0)
								index2-=height;
							if(boundarycond==1)
								index2=height+height-index2;
						}
						for(int l=startx;l<=endx;l++){
							int index3=l;
							if(index3<0){
								if(boundarycond==0)
									index3+=width;
								if(boundarycond==1)
									index3=-index3;
							}
							if(index3>=width){
								if(boundarycond==0)
									index3-=width;
								if(boundarycond==1)
									index3=width+width-index3;
							}
							temp[index2][index3]+=profile[i][k]*((float[][])profiles[index])[j-i+sizes[index]][l-k+sizes[index]];
						}
					}
				}else{
					temp[i][k]+=profile[i][k];
				}
			}
		}
		if(boundarycond==2){
			for(int i=0;i<width;i++){
				temp[0][i]=profile[0][i];
				temp[height-1][i]=profile[height-1][i];
			}
			for(int i=0;i<height;i++){
				temp[i][0]=profile[i][0];
				temp[i][width-1]=profile[i][width-1];
			}
		}
		return temp;
	}

	private float[] handle_2D_diffusion(float[] profile){
		float[] temp=new float[height*width];
		for(int i=0;i<height;i++){
			for(int k=0;k<width;k++){
				int pos=i*width+k;
				int index=Dmap[i][k];
				if(sizes[index]>0){
					int startx=k-sizes[index];
					int endx=k+sizes[index];
					int starty=i-sizes[index];
					int endy=i+sizes[index];
					if(boundarycond>2){
						if(startx<0)
							startx=0;
						if(startx>=width)
							startx=width-1;
						if(endx<0)
							endx=0;
						if(endx>=width)
							endx=width-1;
						if(starty<0)
							starty=0;
						if(starty>=height)
							starty=height-1;
						if(endy<0)
							endy=0;
						if(endy>=height)
							endy=height-1;
					}
					for(int j=starty;j<=endy;j++){
						int index2=j;
						if(index2<0){
							if(boundarycond==0)
								index2+=height;
							if(boundarycond==1)
								index2=-index2;
							if(boundarycond==2)
								index2+=height;
						}
						if(index2>=height){
							if(boundarycond==0)
								index2-=height;
							if(boundarycond==1)
								index2=height+height-index2;
							if(boundarycond==2)
								index2-=height;
						}
						for(int l=startx;l<=endx;l++){
							int index3=l;
							if(index3<0){
								if(boundarycond==0)
									index3+=width;
								if(boundarycond==1)
									index3=-index2;
								if(boundarycond==2)
									index3+=width;
							}
							if(index3>=width){
								if(boundarycond==0)
									index3-=width;
								if(boundarycond==1)
									index3=width+width-index3;
								if(boundarycond==2)
									index3-=width;
							}
							temp[index2*width+index3]+=profile[pos]*((float[][])profiles[index])[j-i+sizes[index]][l-k+sizes[index]];
						}
					}
				}
			}
		}
		if(boundarycond==2){
			// reset the edges to the original values
			int temp2=(height-1)*width;
			for(int i=0;i<width;i++){
				temp[i]=profile[i];
				temp[i+temp2]=profile[i+temp2];
			}
			for(int i=0;i<height;i++){
				temp[i*width]=profile[i*width];
				temp[i*width+width-1]=profile[i*width+width-1];
			}
		}
		return temp;
	}
}
