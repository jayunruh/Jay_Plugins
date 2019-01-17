/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.awt.Polygon;

public class interpolation{
	// this class utilizes simple bilinear and polynomial interpolation methods

	public static float[] shift_image(Object image,int width,int height,float xoff,float yoff){
		float[] retimage=new float[width*height];
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				retimage[j+i*width]=interp2D(image,width,height,j-xoff,i-yoff);
			}
		}
		return retimage;
	}
	
	public static float[] shift_line(Object line,int width,float xoff) {
		float[] retline=new float[width];
		for(int i=0;i<width;i++) {
			retline[i]=interp1D(line,width,(float)i-xoff);
		}
		return retline;
	}
	
	public static float[] affine_transform(Object image,int width,int height,double[][] trans,int newwidth,int newheight){
		//here we transform the coordinates of the new image back to the old image
		float[] retimage=new float[newwidth*newheight];
		float xoff=0.5f*(float)(newwidth-width);
		float yoff=0.5f*(float)(newheight-height);
		for(int i=0;i<newheight;i++){
			float tempy=(float)i-yoff;
			for(int j=0;j<newwidth;j++){
				float tempx=(float)j-xoff;
				float newx=(float)trans[0][0]*tempx+(float)trans[0][1]*tempy+(float)trans[0][2];
				float newy=(float)trans[1][0]*tempx+(float)trans[1][1]*tempy+(float)trans[1][2];
				retimage[j+i*newwidth]=interp2D(image,width,height,newx,newy);
			}
		}
		return retimage;
	}
	
	/**********
	 * here the destination image must be a float image
	 * @param src
	 * @param srcwidth
	 * @param srcheight
	 * @param dest
	 * @param destwidth
	 * @param destheight
	 * @param xoff
	 * @param yoff
	 */
	public static void shift_copy_image(Object src,int srcwidth,int srcheight,float[] dest,int destwidth,int destheight,float xoff,float yoff){
		int ixoff=(int)xoff;
		int iyoff=(int)yoff;
		if(xoff==(float)ixoff && yoff==(float)iyoff){
			int type=algutils.get_array_type(src);
			for(int i=0;i<destheight;i++){
    			int newy=i-(int)yoff;
    			if(newy>=0 && newy<srcheight){
        			for(int j=0;j<destwidth;j++){
        				int newx=j-(int)xoff;
        				if(newx>=0 && newx<srcwidth){
        					dest[j+i*destwidth]=algutils.getArrVal(src,newx+srcwidth*newy,type);
        				}
        			}
    			}
    		}
		} else {
    		for(int i=0;i<destheight;i++){
    			float newy=(float)i-yoff;
    			if(newy>=0.0f && newy<=(float)(srcheight-1)){
        			for(int j=0;j<destwidth;j++){
        				float newx=(float)j-xoff;
        				if(newx>=0.0f && newx<=(float)(srcwidth-1)){
        					dest[j+i*destwidth]=interp2D(src,srcwidth,srcheight,newx,newy);
        				}
        			}
    			}
    		}
		}
	}
	
	/*******************
	 * here the destination can be several array types
	 * @param src
	 * @param srcwidth
	 * @param srcheight
	 * @param dest1
	 * @param destwidth
	 * @param destheight
	 * @param xoff
	 * @param yoff
	 */
	public static void shift_copy_image(Object src,int srcwidth,int srcheight,Object dest1,int destwidth,int destheight,float xoff,float yoff){
		if(src instanceof int[]){
			//this is a color image
			byte[][] rgb=algutils.intval2rgb((int[])src);
			byte[][] temp=new byte[3][destwidth*destheight];
			shift_copy_image(rgb[0],srcwidth,srcheight,temp[0],destwidth,destheight,xoff,yoff);
			shift_copy_image(rgb[1],srcwidth,srcheight,temp[1],destwidth,destheight,xoff,yoff);
			shift_copy_image(rgb[2],srcwidth,srcheight,temp[2],destwidth,destheight,xoff,yoff);
			int[] temp2=algutils.rgb2intval(temp[0],temp[1],temp[2]);
			algutils.copy_subarray(temp2,0,dest1,0,destwidth*destheight);
		} else {
    		float[] dest=new float[destwidth*destheight];
    		shift_copy_image(src,srcwidth,srcheight,dest,destwidth,destheight,xoff,yoff);
    		int type=algutils.get_array_type(src);
    		Object temp=algutils.convert_array(dest,type);
    		algutils.copy_subarray(temp,0,dest1,0,destwidth*destheight);
		}
	}

	public static float[] scale_image(Object image,int width,int height,float s){
		float[] retimage=new float[width*height];
		float xcenter=0.5f*width;
		float ycenter=0.5f*height;
		for(int i=0;i<height;i++){
			float yval=i-ycenter;
			// float yval=(float)i;
			for(int j=0;j<width;j++){
				float xval=j-xcenter;
				// float xval=(float)i;
				float tempx1=xval/s;
				float tempy1=yval/s;
				retimage[j+i*width]=interp2D(image,width,height,tempx1+xcenter,tempy1+ycenter);
			}
		}
		return retimage;
	}

	public static float[] rotate_image(Object image,int width,int height,float angle,float xcenter,float ycenter){
		return rotate_image(image,width,height,angle,xcenter,ycenter,width,height);
	}
	
	public static float[] rotate_image(Object image,int width,int height,float angle,float xcenter,float ycenter,int newwidth,int newheight){
		//ycenter and x center are in the original image coordinates
		float cosval=(float)Math.cos(-angle);
		float sinval=(float)Math.sin(-angle);
		float[] retimage=new float[newwidth*newheight];
		int xoffset=(newwidth-width)/2;
		int yoffset=(newheight-height)/2;
		for(int i=0;i<height;i++){
			float yval=i-ycenter;
			for(int j=0;j<width;j++){
				float xval=j-xcenter;
				float tempx1=xval*cosval-yval*sinval;
				float tempy1=xval*sinval+yval*cosval;
				retimage[j+xoffset+(i+yoffset)*newwidth]=interp2D(image,width,height,tempx1+xcenter,tempy1+ycenter);
			}
		}
		return retimage;
	}

	public static float[][] rotate_image(float[][] image,float angle,float xcenter,float ycenter){
		int width=image.length;
		int height=image[0].length;
		float cosval=(float)Math.cos(-angle);
		float sinval=(float)Math.sin(-angle);
		float[][] retimage=new float[width][height];
		for(int i=0;i<height;i++){
			float yval=i-ycenter;
			for(int j=0;j<width;j++){
				float xval=j-xcenter;
				float tempx1=xval*cosval-yval*sinval;
				float tempy1=xval*sinval+yval*cosval;
				retimage[j][i]=interp2D(image,tempx1+xcenter,tempy1+ycenter);
			}
		}
		return retimage;
	}

	public static boolean[][] rotate_image(boolean[][] image,float angle,float xcenter,float ycenter){
		int width=image.length;
		int height=image[0].length;
		float cosval=(float)Math.cos(-angle);
		float sinval=(float)Math.sin(-angle);
		boolean[][] retimage=new boolean[width][height];
		for(int i=0;i<height;i++){
			float yval=i-ycenter;
			for(int j=0;j<width;j++){
				float xval=j-xcenter;
				float tempx1=xval*cosval-yval*sinval;
				float tempy1=xval*sinval+yval*cosval;
				retimage[j][i]=interp2D(image,tempx1+xcenter,tempy1+ycenter);
			}
		}
		return retimage;
	}
	


	public static float[] rotate_image(Object image,int width,int height,float[][] rmat){
		float[] retimage=new float[width*height];
		float xcenter=0.5f*width;
		float ycenter=0.5f*height;
		for(int i=0;i<height;i++){
			float yval=i-ycenter;
			// float yval=(float)i;
			for(int j=0;j<width;j++){
				float xval=j-xcenter;
				// float xval=(float)i;
				float tempx1=xval*rmat[0][0]+yval*rmat[0][1];
				float tempy1=xval*rmat[1][0]+yval*rmat[1][1];
				retimage[j+i*width]=interp2D(image,width,height,tempx1+xcenter,tempy1+ycenter);
			}
		}
		return retimage;
	}
	
	/********
	 * 
	 * @param stack
	 * @param width
	 * @param height
	 * @param center: this value is in xy pixel units
	 * @param axis: 0,1, and 2 are x, y, and z axes
	 * @param angle
	 * @param zratio
	 * @return
	 */
	public static float[][] rotate_stack(Object[] stack,int[] dims,float[] center,int axis,float angle,float zratio,int[] newdims){
		float[][] newstack=new float[newdims[2]][newdims[0]*newdims[1]];
		if(axis==2){
			//this isn't quite right (bad centerz), have to change it later
			for(int i=0;i<stack.length;i++){
				newstack[i]=interpolation.rotate_image(stack[i],dims[0],dims[1],angle,center[0],center[1],newdims[0],newdims[1]);
			}
			return newstack;
		}
		int[] newcenter={newdims[0]/2,newdims[1]/2,newdims[2]/2};
		float costheta=(float)Math.cos(-angle);
		float sintheta=(float)Math.sin(-angle);
		for(int i=0;i<newdims[2];i++){
			float zpos=zratio*(float)i-newcenter[2];
			for(int j=0;j<newdims[1];j++){
				float ypos=(float)j-newcenter[1];
				for(int k=0;k<newdims[0];k++){
					float xpos=(float)k-newcenter[0];
					float oldzpos=zpos+center[2];
					float oldypos=ypos+center[1];
					float oldxpos=xpos+center[0];
					if(axis==0){ //rotate about x axis
						oldzpos=costheta*zpos+sintheta*ypos+center[2];
						oldypos=costheta*ypos-sintheta*zpos+center[1];
					} else if(axis==1){ //rotate about y axis
						oldzpos=costheta*zpos-sintheta*xpos+center[2];
						oldxpos=costheta*xpos+center[0];
					}
					newstack[i][k+j*newdims[0]]=interpolation.interp3D(stack,dims[0],dims[1],oldxpos,oldypos,oldzpos/zratio);
				}
			}
		}
		return newstack;
	}

	public static Polygon rotate_polygon(Polygon poly,float angle,float xcenter,float ycenter){
		float cosval=(float)Math.cos(-angle);
		float sinval=(float)Math.sin(-angle);
		int[] retxpts=new int[poly.npoints];
		int[] retypts=new int[poly.npoints];
		for(int i=0;i<poly.npoints;i++){
			float yval=poly.ypoints[i]-ycenter;
			float xval=poly.xpoints[i]-xcenter;
			float tempx1=xval*cosval-yval*sinval;
			float tempy1=xval*sinval+yval*cosval;
			retxpts[i]=(int)(tempx1+xcenter);
			retypts[i]=(int)(tempy1+ycenter);
		}
		return new Polygon(retxpts,retypts,poly.npoints);
	}

	public static int[] rotate_color_image(int[] image,int width,int height,float angle,float xcenter,float ycenter){
		byte[][] temp=intval2rgb(image);
		float[][] retimage=new float[3][];
		retimage[0]=rotate_image(temp[0],width,height,angle,xcenter,ycenter);
		retimage[1]=rotate_image(temp[1],width,height,angle,xcenter,ycenter);
		retimage[2]=rotate_image(temp[2],width,height,angle,xcenter,ycenter);
		int[] temp2=rgb2intval(retimage[0],retimage[1],retimage[2]);
		return temp2;
	}

	public static int[] intval2rgb(int value){
		int[] temp=new int[3];
		temp[0]=(value&0xff0000)>>16;
		temp[1]=(value&0xff00)>>8;
		temp[2]=value&0xff;
		return temp;
	}

	public static int[] rgb2intval(float[] r,float[] g,float[] b){
		int[] temp=new int[r.length];
		for(int i=0;i<r.length;i++){
			temp[i]=rgb2intval((byte)((int)r[i]),(byte)((int)g[i]),(byte)((int)b[i]));
		}
		return temp;
	}

	public static byte[][] intval2rgb(int[] values){
		byte[][] temp=new byte[3][values.length];
		for(int i=0;i<values.length;i++){
			int[] temp2=intval2rgb(values[i]);
			temp[0][i]=(byte)temp2[0];
			temp[1][i]=(byte)temp2[1];
			temp[2][i]=(byte)temp2[2];
		}
		return temp;
	}

	public static int rgb2intval(byte r,byte g,byte b){
		int temp=0xff000000|((r&0xff)<<16)|((g&0xff)<<8)|(b&0xff);
		return temp;
	}

	public static float[] rotate_image(Object image,int width,int height,float angle){
		return rotate_image(image,width,height,angle,0.5f*width,0.5f*height);
	}
	
	public static float interp3D(Object[] stack,int width,int height,float x,float y,float z){
		if(z<0.0f) return 0.0f;
		if(z>(stack.length-1)) return 0.0f;
		if(z==0.0f) return interp2D(stack[0],width,height,x,y);
		if(z==stack.length-1) return interp2D(stack[stack.length-1],width,height,x,y);
		int prev=(int)z;
		int next=prev+1;
		float prevval=interp2D(stack[prev],width,height,x,y);
		float nextval=interp2D(stack[next],width,height,x,y);
		float rem=z-(float)prev;
		return rem*(nextval-prevval)+prevval;
	}
	
	public static float[] interpz(Object image1,Object image2,int len,float z){
		if(z<0.0f){
			return null;
		}
		if(z>1.0f){
			return null;
		}
		if(z==0.0f){
			return algutils.convert_arr_float(image1);
		}
		if(z==1.0f){
			return algutils.convert_arr_float(image2);
		}
		float[] retimage=new float[len];
		if(image1 instanceof short[]){
			for(int i=0;i<len;i++){
				float lower=((short[])image1)[i]&0xffff;
				float upper=((short[])image2)[i]&0xffff;
				retimage[i]=lower+z*(upper-lower);
			}
		}else{
			if(image1 instanceof byte[]){
				for(int i=0;i<len;i++){
					float lower=((byte[])image1)[i]&0xff;
					float upper=((byte[])image2)[i]&0xff;
					retimage[i]=lower+z*(upper-lower);
				}
			}else{
				for(int i=0;i<len;i++){
					float lower=((float[])image1)[i];
					float upper=((float[])image2)[i];
					retimage[i]=lower+z*(upper-lower);
				}
			}
		}
		return retimage;
	}

	public static float[] interpz(Object image1,Object image2,int width,int height,float z){
		return interpz(image1,image2,width*height,z);
	}

	public static float interp2D(float[][] image,float x,float y){
		int width=image.length;
		int height=image[0].length;
		if(x<=0.0f){
			return 0.0f;
		}
		if(x>=(width-1)){
			return 0.0f;
		}
		if(y<=0.0f){
			return 0.0f;
		}
		if(y>=(height-1)){
			return 0.0f;
		}
		int intx=(int)x;
		int inty=(int)y;
		float remx=x-intx;
		float remy=y-inty;
		float ul,ur,ll,lr;
		ul=image[intx][inty];
		ur=image[intx+1][inty];
		ll=image[intx][inty+1];
		lr=image[intx+1][inty+1];
		float interpx1=ul+remx*(ur-ul);
		float interpx2=ll+remx*(lr-ll);
		return interpx1+remy*(interpx2-interpx1);
	}

	public static boolean interp2D(boolean[][] image,float x,float y){
		int width=image.length;
		int height=image[0].length;
		if(x<=0.0f){
			return false;
		}
		if(x>=(width-1)){
			return false;
		}
		if(y<=0.0f){
			return false;
		}
		if(y>=(height-1)){
			return false;
		}
		int intx=(int)x;
		int inty=(int)y;
		float remx=x-intx;
		float remy=y-inty;
		float ul,ur,ll,lr;
		ul=image[intx][inty]?1.0f:0.0f;
		ur=image[intx+1][inty]?1.0f:0.0f;
		ll=image[intx][inty+1]?1.0f:0.0f;
		lr=image[intx+1][inty+1]?1.0f:0.0f;
		float interpx1=ul+remx*(ur-ul);
		float interpx2=ll+remx*(lr-ll);
		float newval=interpx1+remy*(interpx2-interpx1);
		return(newval>0.5f);
	}

	public static float interp2D(Object image,int width,int height,float x,float y){
		if(x<0.0f) return 0.0f;
		if(x>(float)(width-1)) return 0.0f;
		if(y<0.0f) return 0.0f;
		if(y>(float)(height-1)) return 0.0f;
		int intx=(int)x;
		int inty=(int)y;
		float remx=x-intx;
		float remy=y-inty;
		float ul,ur,ll,lr;
		if(image instanceof short[]){
			short[] temp=(short[])image;
			ul=temp[intx+width*inty]&0xffff;
			if(intx<(width-1)) ur=temp[intx+1+width*inty]&0xffff;
			else ur=ul;
			if(inty<(height-1)) ll=temp[intx+width*(inty+1)]&0xffff;
			else ll=ul;
			if(intx<(width-1) && inty<(height-1)) lr=temp[intx+1+width*(inty+1)]&0xffff;
			else if(intx>=(width-1) && inty>=(height-1)) lr=ul;
			else if(intx>=(width-1)) lr=ll;
			else lr=ur;
		}else{
			if(image instanceof byte[]){
				byte[] temp=(byte[])image;
				ul=temp[intx+width*inty]&0xff;
				if(intx<(width-1)) ur=temp[intx+1+width*inty]&0xff;
				else ur=ul;
				if(inty<(height-1)) ll=temp[intx+width*(inty+1)]&0xff;
				else ll=ul;
				if(intx<(width-1) && inty<(height-1)) lr=temp[intx+1+width*(inty+1)]&0xff;
				else if(intx>=(width-1) && inty>=(height-1)) lr=ul;
				else if(intx>=(width-1)) lr=ll;
				else lr=ur;
			}else{
				float[] temp=(float[])image;
				ul=temp[intx+width*inty];
				if(intx<(width-1)) ur=temp[intx+1+width*inty];
				else ur=ul;
				if(inty<(height-1)) ll=temp[intx+width*(inty+1)];
				else ll=ul;
				if(intx<(width-1) && inty<(height-1)) lr=temp[intx+1+width*(inty+1)];
				else if(intx>=(width-1) && inty>=(height-1)) lr=ul;
				else if(intx>=(width-1)) lr=ll;
				else lr=ur;
			}
		}
		float interpx1=ul+remx*(ur-ul);
		float interpx2=ll+remx*(lr-ll);
		return interpx1+remy*(interpx2-interpx1);
	}

	public static float interp2D_pad(Object image,int width,int height,float x,float y){
		// here we move out of bounds points orthogonally to be in bounds
		float x2=x;
		float y2=y;
		if(x2<0.0f){
			x2=0.0f;
		}
		if(x2>(width-1)){
			x2=width-1;
		}
		if(y2<0.0f){
			y2=0.0f;
		}
		if(y2>(height-1)){
			y2=height-1;
		}
		int intx=(int)x2;
		int inty=(int)y2;
		float remx=x2-intx;
		float remy=y2-inty;
		if(intx==(width-1)){
			intx-=1;
			remx=1.0f;
		}
		if(inty==(height-1)){
			inty-=1;
			remy=1.0f;
		}
		float ul,ur,ll,lr;
		if(image instanceof short[]){
			ul=((short[])image)[intx+width*inty]&0xffff;
			ur=((short[])image)[intx+width*inty+1]&0xffff;
			ll=((short[])image)[intx+width*(inty+1)]&0xffff;
			lr=((short[])image)[intx+width*(inty+1)+1]&0xffff;
		}else{
			if(image instanceof byte[]){
				ul=((byte[])image)[intx+width*inty]&0xff;
				ur=((byte[])image)[intx+width*inty+1]&0xff;
				ll=((byte[])image)[intx+width*(inty+1)]&0xff;
				lr=((byte[])image)[intx+width*(inty+1)+1]&0xff;
			}else{
				ul=((float[])image)[intx+width*inty];
				ur=((float[])image)[intx+width*inty+1];
				ll=((float[])image)[intx+width*(inty+1)];
				lr=((float[])image)[intx+width*(inty+1)+1];
			}
		}
		float interpx1=ul+remx*(ur-ul);
		float interpx2=ll+remx*(lr-ll);
		return interpx1+remy*(interpx2-interpx1);
	}

	public static float interp1D(Object image,int width,float x){
		if(x<0.0f){
			return 0.0f;
		}
		if(x>(width-1)){
			return 0.0f;
		}
		if(x==(width-1)){
			if(image instanceof short[])
				return ((short[])image)[width-1]&0xffff;
			if(image instanceof byte[])
				return ((byte[])image)[width-1]&0xff;
			if(image instanceof float[])
				return ((float[])image)[width-1];
			return 0.0f;
		}
		int intx=(int)x;
		float remx=x-intx;
		float ul,ur;
		if(image instanceof short[]){
			ul=((short[])image)[intx]&0xffff;
			ur=((short[])image)[intx+1]&0xffff;
		}else{
			if(image instanceof byte[]){
				ul=((byte[])image)[intx]&0xff;
				ur=((byte[])image)[intx+1]&0xff;
			}else{
				ul=((float[])image)[intx];
				ur=((float[])image)[intx+1];
			}
		}
		if(remx==0.0f)
			return ul; // prevent an unneccessary multiplication
		else
			return ul+remx*(ur-ul);
	}

	public static float[] get_local_max(float[] image,int x,int y,int width,int height){
		if(x<=0||x>=(width-1)||y<=0||y>=(height-1)){
			return null;
		}
		float[] temp=algutils.getNeighbors2(image,x,y,width,height);
		float[] temp2=get_neighborhood_max(temp);
		float[] temp3={temp2[0]+x,temp2[1]+y};
		return temp3;
	}

	public static float[] get_local_max3D(Object[] image,int x,int y,int z,int width,int height){
		if(x<=0||x>=(width-1)||y<=0||y>=(height-1)||z<=0||z>=(image.length-1)){
			return null;
		}
		Object temp=algutils.getNeighbors2(image,x,y,z,width,height);
		float[] temp2=get_neighborhood_max3D((float[])temp);
		float[] temp3={temp2[0]+x,temp2[1]+y,temp2[2]+z};
		return temp3;
	}

	public static float[] get_neighborhood_max3D(float[] arr){
		// here we find a neighborhood maximum by fitting the x and y
		// projections to a parabola
		float[] x=cube3Dproj(arr,0);
		float[] y=cube3Dproj(arr,1);
		float[] z=cube3Dproj(arr,2);
		float[] maxima={get_neighborhood_max1D(x),get_neighborhood_max1D(y),get_neighborhood_max1D(z)};
		if(maxima[0]>1.0f)
			maxima[0]=1.0f;
		if(maxima[0]<-1.0f)
			maxima[0]=-1.0f;
		if(maxima[1]>1.0f)
			maxima[1]=1.0f;
		if(maxima[1]<-1.0f)
			maxima[1]=-1.0f;
		if(maxima[2]>1.0f)
			maxima[2]=1.0f;
		if(maxima[2]<-1.0f)
			maxima[2]=-1.0f;
		return maxima;
	}

	public static float[] cube3Dproj(float[] cube,int projaxis){
		if(projaxis==0){
			float[] temp=new float[3];
			for(int z=0;z<3;z++){
				for(int y=0;y<3;y++){
					for(int x=0;x<3;x++){
						temp[x]+=cube[x+y*3+z*9];
					}
				}
			}
			return temp;
		}else{
			if(projaxis==1){
				float[] temp=new float[3];
				for(int z=0;z<3;z++){
					for(int y=0;y<3;y++){
						for(int x=0;x<3;x++){
							temp[y]+=cube[x+y*3+z*9];
						}
					}
				}
				return temp;
			}else{
				float[] temp=new float[3];
				for(int z=0;z<3;z++){
					for(int y=0;y<3;y++){
						for(int x=0;x<3;x++){
							temp[z]+=cube[x+y*3+z*9];
						}
					}
				}
				return temp;
			}
		}
	}

	public static float[] get_neighborhood_max(float[] arr){
		// here we find a neighborhood maximum by fitting the x and y
		// projections to a parabola
		float[] x={arr[0]+arr[3]+arr[6],arr[1]+arr[4]+arr[7],arr[2]+arr[5]+arr[8]};
		float[] y={arr[0]+arr[1]+arr[2],arr[3]+arr[4]+arr[5],arr[6]+arr[7]+arr[8]};
		float[] maxima=new float[2];
		maxima[0]=get_neighborhood_max1D(x);
		maxima[1]=get_neighborhood_max1D(y);
		return maxima;
	}

	public static float get_max1D(float[] data){
		float max=data[0];
		int maxpt=0;
		for(int i=1;i<data.length;i++){
			if(max>data[i]){
				max=data[i];
				maxpt=i;
			}
		}
		if(maxpt<=0)
			maxpt=1;
		if(maxpt>=(data.length-1))
			maxpt=data.length-2;
		return get_local_max1D(data,maxpt);
	}

	public static float get_local_max1D(float[] data,int x){
		int width=data.length;
		if(x<=0||x>=(width-1)){
			return Float.NaN;
		}
		float[] temp={data[x-1],data[x],data[x+1]};
		float temp2=get_neighborhood_max1D(temp);
		return temp2+x;
	}

	public static float get_neighborhood_max1D(float[] arr){
		// here we find a neighborhood maximum by fitting the x and y
		// projections to a parabola
		float ax=0.5f*(arr[2]+arr[0])-arr[1];
		float bx=0.5f*(arr[2]-arr[0]);
		if(ax<0.0f){
			float temp=-0.5f*bx/ax;
			if(temp>1.0f)
				return 1.0f;
			if(temp<-1.0f)
				return -1.0f;
			return temp;
		}else{
			return 0.0f;
		}
	}

	public static float[] get_closest_point(float[][] pl,float[] pt){
		// here we find the closest point on polyline (pl) {X,Y} to point pt
		// find the closest point in the longer polygon
		float mindist=calcdist(pt[0],pt[1],pl[0][0],pl[1][0]);
		int minindex=0;
		for(int j=1;j<pl[0].length;j++){
			float distance=calcdist(pt[0],pt[1],pl[0][j],pl[1][j]);
			if(distance<mindist){
				mindist=distance;
				minindex=j;
			}
		}
		// now find the best neighboring point to form a triangle
		int temp=minindex-1;
		if(temp<0)
			temp=pl[0].length-1;
		float distprev=calcdist(pt[0],pt[1],pl[0][temp],pl[1][temp]);
		float mindist2=distprev;
		int minindex2=temp;
		temp=minindex+1;
		if(temp>=pl[0].length)
			temp=0;
		float distnext=calcdist(pt[0],pt[1],pl[0][temp],pl[1][temp]);
		if(distnext<mindist2){
			mindist2=distnext;
			minindex2=temp;
		}
		// now find the closed point in pl to our current pt point
		float u=(pt[0]-pl[0][minindex])*(pl[0][minindex2]-pl[0][minindex])+(pt[1]-pl[1][minindex])*(pl[1][minindex2]-pl[1][minindex]);
		float dist2=(pl[0][minindex]-pl[0][minindex2])*(pl[0][minindex]-pl[0][minindex2])+(pl[1][minindex]-pl[1][minindex2])*(pl[1][minindex]-pl[1][minindex2]);
		u/=dist2;
		float x=pl[0][minindex]+u*(pl[0][minindex2]-pl[0][minindex]);
		float y=pl[1][minindex]+u*(pl[1][minindex2]-pl[1][minindex]);
		return new float[]{x,y};
	}
	
	public static float[] get_closest_point_3D(float[][] pl,float[] pt){
		// here we find the closest point on polyline (pl) {X,Y} to point pt
		// find the closest point in the longer polygon
		float mindist=calcdist_3D(pt[0],pt[1],pt[2],pl[0][0],pl[1][0],pl[2][0]);
		int minindex=0;
		for(int j=1;j<pl[0].length;j++){
			float distance=calcdist_3D(pt[0],pt[1],pt[2],pl[0][j],pl[1][j],pl[2][j]);
			if(distance<mindist){
				mindist=distance;
				minindex=j;
			}
		}
		// now find the best neighboring point to form a triangle
		int temp=minindex-1;
		if(temp<0)
			temp=pl[0].length-1;
		float distprev=calcdist_3D(pt[0],pt[1],pt[2],pl[0][temp],pl[1][temp],pl[2][temp]);
		float mindist2=distprev;
		int minindex2=temp;
		temp=minindex+1;
		if(temp>=pl[0].length)
			temp=0;
		float distnext=calcdist_3D(pt[0],pt[1],pt[2],pl[0][temp],pl[1][temp],pl[2][temp]);
		if(distnext<mindist2){
			mindist2=distnext;
			minindex2=temp;
		}
		// now find the closed point in pl to our current pt point
		float u=(pt[0]-pl[0][minindex])*(pl[0][minindex2]-pl[0][minindex])+(pt[1]-pl[1][minindex])*(pl[1][minindex2]-pl[1][minindex])+(pt[2]-pl[2][minindex])*(pl[2][minindex2]-pl[2][minindex]);
		float dist2=(pl[0][minindex]-pl[0][minindex2])*(pl[0][minindex]-pl[0][minindex2])+(pl[1][minindex]-pl[1][minindex2])*(pl[1][minindex]-pl[1][minindex2])+(pl[2][minindex]-pl[2][minindex2])*(pl[2][minindex]-pl[2][minindex2]);
		u/=dist2;
		float x=pl[0][minindex]+u*(pl[0][minindex2]-pl[0][minindex]);
		float y=pl[1][minindex]+u*(pl[1][minindex2]-pl[1][minindex]);
		float z=pl[2][minindex]+u*(pl[2][minindex2]-pl[2][minindex]);
		return new float[]{x,y,z};
	}
	
	public static float[] get_closest_index(float[][] pl,float[] pt){
		float mindist=calcdist(pt[0],pt[1],pl[0][0],pl[1][0]);
		int minindex=0;
		for(int j=1;j<pl[0].length;j++){
			float distance=calcdist(pt[0],pt[1],pl[0][j],pl[1][j]);
			if(distance<mindist){
				mindist=distance;
				minindex=j;
			}
		}
		return new float[]{minindex,mindist};
	}
	
	public static float[] get_closest_index_3D(float[][] pl,float[] pt){
		float mindist=calcdist_3D(pt[0],pt[1],pt[2],pl[0][0],pl[1][0],pl[2][0]);
		int minindex=0;
		for(int j=1;j<pl[0].length;j++){
			float distance=calcdist_3D(pt[0],pt[1],pt[2],pl[0][j],pl[1][j],pl[2][j]);
			if(distance<mindist){
				mindist=distance;
				minindex=j;
			}
		}
		return new float[]{minindex,mindist};
	}

	public static float calcdist(float x1,float y1,float x2,float y2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}
	
	public static float calcdist_3D(float x1,float y1,float z1,float x2,float y2,float z2){
		return (float)Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	}
	
	public static float calcdist_ND(float[] x1,float[] x2){
		float dist2=0.0f;
		for(int i=0;i<x1.length;i++){
			dist2+=(x2[i]-x1[i])*(x2[i]-x1[i]);
		}
		return (float)Math.sqrt(dist2);
	}

	public static float get_float_index(float[] arr,float position){
		int length=arr.length;
		if(position<=arr[0])
			return 0.0f;
		if(position>=arr[length-1])
			return length-1.0f;
		int i=1;
		while(arr[i]<position)
			i++;
		return (position-arr[i-1])/(arr[i]-arr[i-1])+(i-1);
	}

	public static float center_of_mass(float[] arr){
		double count=0.0;
		double xsum=0.0;
		for(int i=0;i<arr.length;i++){
			count+=arr[i];
			xsum+=arr[i]*(double)i;
		}
		return (float)(xsum/count+0.5);
	}
	
	public static float[] center_of_mass_2D(Object arr,int width,int height){
		if(arr instanceof float[]) return center_of_mass_2D((float[])arr,width,height);
		else if(arr instanceof byte[]) return center_of_mass_2D((byte[])arr,width,height);
		else return center_of_mass_2D((short[])arr,width,height);
	}
	
	public static float[] center_of_mass_2D(float[] arr,int width,int height){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		int counter=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				count+=arr[counter];
				xsum+=arr[counter]*(double)j;
				ysum+=arr[counter]*(double)i;
				counter++;
			}
		}
		return new float[]{(float)(xsum/count),(float)(ysum/count)};
	}
	
	public static float[] center_of_mass_2D(short[] arr,int width,int height){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		int counter=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				float temp=arr[counter]&0xffff;
				count+=temp;
				xsum+=temp*(double)j;
				ysum+=temp*(double)i;
				counter++;
			}
		}
		return new float[]{(float)(xsum/count),(float)(ysum/count)};
	}
	
	public static float[] center_of_mass_2D(byte[] arr,int width,int height){
		double count=0.0;
		double xsum=0.0;
		double ysum=0.0;
		int counter=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				float temp=arr[counter]&0xff;
				count+=temp;
				xsum+=temp*(double)j;
				ysum+=temp*(double)i;
				counter++;
			}
		}
		return new float[]{(float)(xsum/count),(float)(ysum/count)};
	}

	public static float integrate(Object func,int len,float begin,float end,boolean cumulative){
		int next=1+(int)begin;
		int prev=(int)end;
		float integral=0.0f;
		if(prev<next)
			integral+=(end-begin)*interp1D(func,len,0.5f*(begin+end)); // we are contained within an interval
		else{
			integral+=(next-begin)*interp1D(func,len,0.5f*(next+begin)); // get the interval to the next point
			integral+=(end-prev)*interp1D(func,len,0.5f*(end+prev)); // get the interval to the end point
			for(int i=next;i<prev;i++)
				integral+=interp1D(func,len,i+0.5f); // and all the ones in between
		}
		if(cumulative)
			return integral;
		else
			return integral/(end-begin);
	}
	
	public static float[] makeMirror(float[] avgvals){
		float[] newavgvals=new float[2*avgvals.length-1];
		int halfsize=avgvals.length-1;
		//int newwidth=2*avgvals.length-1;
		//here make sure we don't replicate the center point twice
			newavgvals[halfsize]=avgvals[0];
			for(int j=1;j<avgvals.length;j++){
				newavgvals[j+halfsize]=avgvals[j];
				newavgvals[halfsize-j]=avgvals[j];
			}
		return newavgvals;
	}
	
	public static float[] circavg(Object data,int width,int height,int rsize,float xc,float yc,boolean makeMirror1){
		float[] circavg1=circavg(data,width,height,rsize,xc,yc);
		if(!makeMirror1) return circavg1;
		else return makeMirror(circavg1);
	}
	
	public static float[] circavg(Object data,int width,int height,int rsize,float xc,float yc){
		float twopi=(float)(2.0*Math.PI);
		float[] avgvals=new float[rsize];
		avgvals[0]=interpolation.interp2D(data,width,height,xc,yc);
		for(int j=1;j<rsize;j++){
			float angleincrement=1.0f/(float)j;
			int angles=(int)(twopi/angleincrement);
			angleincrement=twopi/(float)angles;
			for(int k=0;k<angles;k++){
				float angle=angleincrement*(float)k;
				float x=(float)Math.cos(angle)*(float)j+xc;
				float y=(float)Math.sin(angle)*(float)j+yc;
				avgvals[j]+=interp2D(data,width,height,x,y);
			}
			avgvals[j]/=(float)angles;
		}
		return avgvals;
	}
	
	public static float[][] resampleTraj(float[][] traj,float dx){
		//this resamples an n-dimension trajectory with spacing dx
		int dims=traj[0].length;
		float[] dists=new float[traj.length-1];
		float totdist=0.0f;
		for(int i=1;i<traj.length;i++){
			dists[i-1]=calcdist_ND(traj[i],traj[i-1]);
			totdist+=dists[i-1];
		}
		int npts=(int)(totdist/dx)+1;
		float[][] newtraj=new float[npts][traj[0].length];
		float pos=0.0f;
		int ind=0;
		newtraj[0]=traj[0];
		for(int i=1;i<npts;i++){
			float target=dx*(float)i;
			while(pos<target){
				pos+=dists[ind];
				ind++;
			}
			if(pos==target){
				newtraj[i]=traj[ind-1];
			} else {
				float prev=pos-dists[ind-2];
				float frac=(target-prev)/(pos-prev);
				for(int j=0;j<dims;j++) newtraj[i][j]=traj[ind-2][j]+frac*(traj[ind-1][j]-traj[ind-2][j]);
			}
		}
		return newtraj;
	}

}
