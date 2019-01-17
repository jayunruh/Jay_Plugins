/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

public class jsobel{
	public int width,height;

	public jsobel(int width1,int height1){
		width=width1;
		height=height1;
	}

	public float[][] do_sobel(float[] pixels){
		float[] xgrad=new float[width*height];
		float[] ygrad=new float[width*height];
		float[] gradmag=new float[width*height];
		float[] gradangle=new float[width*height];
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				ygrad[i*width+j]=pixels[(i+1)*width+j-1]+2.0f*pixels[(i+1)*width+j]+pixels[(i+1)*width+j+1]-pixels[(i-1)*width+j-1]-2.0f*pixels[(i-1)*width+j]-pixels[(i-1)*width+j+1];
				xgrad[i*width+j]=pixels[(i+1)*width+j+1]+2.0f*pixels[i*width+j+1]+pixels[(i-1)*width+j+1]-pixels[(i+1)*width+j-1]-2.0f*pixels[i*width+j-1]-pixels[(i-1)*width+j-1];
				int temp=i*width+j;
				float[] mag_angle=get_mag_angle(xgrad[temp],ygrad[temp]);
				gradmag[temp]=mag_angle[0];
				gradangle[temp]=mag_angle[1];
			}
		}
		float[][] retvals={gradmag,gradangle};
		return retvals;
	}
	
	/****************
	 * here we use a second derivative filter to get the ridges instead of edges--this doesn't work at all
	 * @param pixels
	 * @return
	 */
	public float[][] get_ridges(float[] pixels){
		float[] xgrad=new float[width*height];
		float[] ygrad=new float[width*height];
		float[] gradmag=new float[width*height];
		float[] gradangle=new float[width*height];
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				ygrad[i*width+j]=-pixels[(i+1)*width+j-1]+2.0f*pixels[(i+1)*width+j]-pixels[(i+1)*width+j+1]-2.0f*pixels[i*width+j-1]+4.0f*pixels[i*width+j]-2.0f*pixels[i*width+j+1]-pixels[(i-1)*width+j-1]+2.0f*pixels[(i-1)*width+j]-pixels[(i-1)*width+j+1];
				xgrad[i*width+j]=-pixels[(i+1)*width+j-1]-2.0f*pixels[(i+1)*width+j]-pixels[(i+1)*width+j+1]+2.0f*pixels[i*width+j-1]+4.0f*pixels[i*width+j]+2.0f*pixels[i*width+j+1]-pixels[(i-1)*width+j-1]-2.0f*pixels[(i-1)*width+j]-pixels[(i-1)*width+j+1];
				int temp=i*width+j;
				float[] mag_angle=get_mag_angle(xgrad[temp],ygrad[temp]);
				gradmag[temp]=mag_angle[0];
				gradangle[temp]=mag_angle[1];
			}
		}
		float[][] retvals={gradmag,gradangle};
		return retvals;
	}
	
	public float[][] get_gradient(float[] pixels){
		float[] xgrad=new float[width*height];
		float[] ygrad=new float[width*height];
		float[] gradmag=new float[width*height];
		float[] gradangle=new float[width*height];
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				ygrad[i*width+j]=0.5f*(pixels[(i+1)*width+j]-pixels[(i-1)*width+j]);
				xgrad[i*width+j]=0.5f*(pixels[i*width+j+1]-pixels[i*width+j-1]);
				int temp=i*width+j;
				float[] mag_angle=get_mag_angle(xgrad[temp],ygrad[temp]);
				gradmag[temp]=mag_angle[0];
				gradangle[temp]=mag_angle[1];
			}
		}
		float[][] retvals={gradmag,gradangle};
		return retvals;
	}

	public float[][] do_sobel2(float[] pixels){
		float[] xgrad=new float[width*height];
		float[] ygrad=new float[width*height];
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				ygrad[i*width+j]=pixels[(i+1)*width+j-1]+2.0f*pixels[(i+1)*width+j]+pixels[(i+1)*width+j+1]-pixels[(i-1)*width+j-1]-2.0f*pixels[(i-1)*width+j]-pixels[(i-1)*width+j+1];
				xgrad[i*width+j]=pixels[(i+1)*width+j+1]+2.0f*pixels[i*width+j+1]+pixels[(i-1)*width+j+1]-pixels[(i+1)*width+j-1]-2.0f*pixels[i*width+j-1]-pixels[(i-1)*width+j-1];
			}
		}
		float[][] retvals={xgrad,ygrad};
		return retvals;
	}

	public float[][] do_sobel_NaN(float[] pixels){
		float[] gradmag=new float[width*height];
		float[] gradangle=new float[width*height];
		for(int i=1;i<(height-1);i++){
			for(int j=1;j<(width-1);j++){
				float[] temp2=sobel_point(pixels,j,i);
				int temp=i*width+j;
				float[] mag_angle=get_mag_angle(temp2[0],temp2[1]);
				gradmag[temp]=mag_angle[0];
				gradangle[temp]=mag_angle[1];
			}
		}
		float[][] retvals={gradmag,gradangle};
		return retvals;
	}

	public float[] sobel_point(float[] pixels,int x,int y){
		int offset=x+y*width;
		float[] output=new float[2];
		if(x<=0)
			return output;
		if(x>=(width-1))
			return output;
		if(y<=0)
			return output;
		if(y>=(height-1))
			return output;
		if(Float.isNaN(pixels[offset])){
			return output;
		}
		float val=pixels[offset];
		float[][] n=new float[3][3];
		n[0][0]=pixels[offset-width-1];
		n[0][1]=pixels[offset-width];
		n[0][2]=pixels[offset-width+1];
		n[1][0]=pixels[offset-1];
		n[1][1]=pixels[offset];
		n[1][2]=pixels[offset+1];
		n[2][0]=pixels[offset+width-1];
		n[2][1]=pixels[offset+width];
		n[2][2]=pixels[offset+width+1];
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				if(Float.isNaN(n[i][j]))
					n[i][j]=val;
			}
		}
		output[0]=n[2][2]+2.0f*n[1][2]+n[0][2]-n[2][0]-2.0f*n[1][0]-n[0][0];
		output[1]=n[2][0]+2.0f*n[2][1]+n[2][2]-n[0][0]-2.0f*n[0][1]-n[0][2];
		return output;
	}

	public float[] get_mag_angle(float dx,float dy){
		float[] mag_angle=new float[2];
		mag_angle[0]=(float)Math.sqrt(dx*dx+dy*dy);
		if(dx==0.0f){
			if(dy==0.0f){
				mag_angle[1]=0.0f;
			}else{
				mag_angle[1]=(float)(Math.PI/2.0);
			}
		}else{
			mag_angle[1]=(float)Math.atan2(dx,dy);
			if(mag_angle[1]<0.0f){
				mag_angle[1]+=(float)Math.PI;
			}
		}
		return mag_angle;
	}

	public float[][] get_mag_angle(float[] dx,float[] dy){
		float[][] mag_angle=new float[2][width*height];
		for(int i=1;i<height-1;i++){
			for(int j=1;j<width-1;j++){
				int temp=j+i*width;
				float[] temp2=get_mag_angle(dx[temp],dy[temp]);
				mag_angle[0][temp]=temp2[0];
				mag_angle[1][temp]=temp2[1];
			}
		}
		return mag_angle;
	}

	public float[] get_dx_dy(float mag,float angle){
		float[] dxdy=new float[2];
		dxdy[0]=mag*(float)Math.cos(angle);
		dxdy[1]=mag*(float)Math.sin(angle);
		return dxdy;
	}

}
