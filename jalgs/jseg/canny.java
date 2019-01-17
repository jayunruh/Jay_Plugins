/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

public class canny{
	public float thresh1,thresh2;
	public int width,height;
	private float quad0,quad1,quad2,quad3;

	public canny(float thresh11,float thresh21,int width1,int height1){
		thresh1=thresh11;
		thresh2=thresh21;
		width=width1;
		height=height1;
		quad0=(float)(Math.PI/8.0);
		quad1=(float)(3.0*Math.PI/8.0);
		quad2=(float)(5.0*Math.PI/8.0);
		quad3=(float)(7.0*Math.PI/8.0);
	}

	public byte[] find_edges(float[] pixels,float[][] sobel){
		float[] gradmag=sobel[0];
		float[] gradangle=sobel[1];
		// now do non-maximum suppression (must be a maximum along gradient)
		byte[] hisedgeimage=new byte[width*height];
		boolean[] maxangle=new boolean[width*height];
		boolean[] notminangle=new boolean[width*height];
		// now that we have the gradient and angle, use maximum detection along
		// the gradient to detect the edges
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				float ga=gradangle[j+i*width];
				maxangle[i*width+j]=ismaximumangle(gradmag,i,j,ga);
				notminangle[i*width+j]=isnotminangle(gradmag,i,j,ga);
				if(gradmag[i*width+j]>=thresh1){
					if(maxangle[i*width+j]){
						hisedgeimage[j+i*width]=(byte)255;
					}
				}
			}
		}
		// now do the hysteresis threshholding
		int nfound;
		int counter=0;
		do{
			nfound=0;
			for(int i=0;i<(height-1);i++){
				for(int j=0;j<(width-1);j++){
					if(hisedgeimage[j+i*width]!=(byte)0){
						// search around each pixel for pixels above the lower
						// threshhold
						int temp=j-1+(i-1)*width;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp++;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp++;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp=j-1+i*width;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp++;
						temp++;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp=j-1+(i+1)*width;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp++;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
						temp++;
						if((gradmag[temp]>=thresh2&&maxangle[temp])&&hisedgeimage[temp]==(byte)0){
							hisedgeimage[temp]=(byte)255;
							nfound++;
						}
					}
				}
			}
			counter++;
			if(counter>100){
				break;
			}
		}while(nfound>0);
		return hisedgeimage;
	}

	public byte[] find_edges(float[] pixels){
		// first calculate the sobel gradients in the x and y direction
		float[][] sobel=(new jsobel(width,height)).do_sobel(pixels);
		return find_edges(pixels,sobel);
	}
	
	public byte[] find_ridges(float[] pixels) {
		//doesn't work!
		float[][] sobel=(new jsobel(width,height)).get_ridges(pixels);
		return find_edges(pixels,sobel);
	}

	private boolean ismaximumangle(float[] gradmag,int i,int j,float ga){
		// 0,1,2
		// 3,4,5
		// 6,7,8

		float[] temp={gradmag[(i-1)*width+j-1],gradmag[(i-1)*width+j],gradmag[(i-1)*width+j+1],gradmag[i*width+j-1],gradmag[i*width+j],gradmag[i*width+j+1],gradmag[(i+1)*width+j-1],
				gradmag[(i+1)*width+j],gradmag[(i+1)*width+j+1]};
		if(ga<quad0||ga>=quad3){
			if(temp[4]>=temp[1]&&temp[4]>=temp[7]){
				return true;
			}
		}else{
			if(ga<quad1){
				if(temp[4]>=temp[0]&&temp[4]>=temp[8]){
					return true;
				}
			}else{
				if(ga<quad2){
					if(temp[4]>=temp[3]&&temp[4]>=temp[5]){
						return true;
					}
				}else{
					if(temp[4]>=temp[6]&&temp[4]>=temp[2]){
						return true;
					}
				}
			}
		}
		return false;
	}

	private boolean isnotminangle(float[] gradmag,int i,int j,float ga){
		// 0,1,2
		// 3,4,5
		// 6,7,8
		float[] temp={gradmag[(i-1)*width+j-1],gradmag[(i-1)*width+j],gradmag[(i-1)*width+j+1],gradmag[i*width+j-1],gradmag[i*width+j],gradmag[i*width+j+1],gradmag[(i+1)*width+j-1],
				gradmag[(i+1)*width+j],gradmag[(i+1)*width+j+1]};
		if(ga<quad0||ga>=quad3){
			if(temp[4]>=temp[1]||temp[4]>=temp[7]){
				return true;
			}
		}else{
			if(ga<quad1){
				if(temp[4]>=temp[0]||temp[4]>=temp[8]){
					return true;
				}
			}else{
				if(ga<quad2){
					if(temp[4]>=temp[3]||temp[4]>=temp[5]){
						return true;
					}
				}else{
					if(temp[4]>=temp[6]||temp[4]>=temp[2]){
						return true;
					}
				}
			}
		}
		return false;
	}
}
