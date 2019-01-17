package jalgs.jfft;

public class padding{
	
	public static float[] pad_xy_mirrored(float[] source,int width,int height,int newwidth,int newheight,boolean feather){
		//int newwidth=width+2*padx;
		//int newheight=height+2*padx;
		int padx=(newwidth-width)/2;
		int pady=(newheight-height)/2;
		int xstart=padx;
		int xend=xstart+width;
		int ystart=pady;
		int yend=ystart+height;
		float[] output=new float[newwidth*newheight];
		//first pad in the x direction
		for(int i=0;i<newwidth;i++){
			int index=i-xstart;
			float amp=1.0f;
			if(index<0){
				index=-index;
				if(feather) amp=(float)i/(float)xstart;
			}
			if(index>=width){
				//index=newwidth-index-1;
				index=width-(index-width)-1;
				if(feather) amp=((float)newwidth-(float)i)/padx;
			}
			if(amp<0.0f) amp=0.0f;
			for(int j=0;j<height;j++){
				output[(j+ystart)*newwidth+i]=amp*source[j*width+index];
			}
		}
		//now in the y direction
		//at the top
		for(int i=0;i<ystart;i++){
			int index=ystart+(ystart-i);
			float amp=1.0f;
			if(feather) amp=(float)i/(float)ystart;
			if(amp<0.0f) amp=0.0f;
			for(int j=0;j<newwidth;j++) output[j+i*newwidth]=amp*output[j+index*newwidth];
		}
		//and the bottom
		for(int i=yend;i<newheight;i++){
			int index=yend-(i-yend)-1;
			float amp=1.0f;
			if(feather) amp=((float)newheight-(float)i)/pady;
			if(amp<0.0f) amp=0.0f;
			for(int j=0;j<newwidth;j++) output[j+i*newwidth]=amp*output[j+index*newwidth];
		}
		return output;
	}
	
	public static float[] pad_xy_zeros(float[] source,int width,int height,int newwidth,int newheight){
		//int newwidth=width+2*padx;
		//int newheight=height+2*padx;
		int padx=(newwidth-width)/2;
		int pady=(newheight-height)/2;
		int xstart=padx;
		//int xend=xstart+width;
		int ystart=pady;
		//int yend=ystart+height;
		float[] output=new float[newwidth*newheight];
		//first pad in the x direction
		for(int i=0;i<width;i++){
			for(int j=0;j<height;j++){
				output[(j+ystart)*newwidth+i+xstart]=source[j*width+i];
			}
		}
		return output;
	}
	
	public static Object get_region2_padmirrored(Object source,int x,int y,int rwidth,int rheight,int width,int height){
		if(source instanceof byte[]) return get_region2_padmirrored((byte[])source,x,y,rwidth,rheight,width,height);
		else if(source instanceof short[]) return get_region2_padmirrored((short[])source,x,y,rwidth,rheight,width,height);
		else return get_region2_padmirrored((float[])source,x,y,rwidth,rheight,width,height);
	}
	
	public static byte[] get_region2_padmirrored(byte[] source,int x,int y,int rwidth,int rheight,int width,int height){
		byte[] output=new byte[rwidth*rheight];
		int counter=0;
		for(int i=y;i<(y+rheight);i++){
			int ypos=i;
			if(ypos<0) ypos=-ypos;
			if(ypos>=height) ypos=2*height-ypos-1;
			for(int j=x;j<(x+rwidth);j++){
				int xpos=j;
				if(xpos<0) xpos=-xpos;
				if(xpos>=width) xpos=2*width-xpos-1;
				output[counter]=source[xpos+ypos*width];
				counter++;
			}
		}
		return output;
	}
	
	public static short[] get_region2_padmirrored(short[] source,int x,int y,int rwidth,int rheight,int width,int height){
		short[] output=new short[rwidth*rheight];
		int counter=0;
		for(int i=y;i<(y+rheight);i++){
			int ypos=i;
			if(ypos<0) ypos=-ypos;
			if(ypos>=height) ypos=2*height-ypos-1;
			for(int j=x;j<(x+rwidth);j++){
				int xpos=j;
				if(xpos<0) xpos=-xpos;
				if(xpos>=width) xpos=2*width-xpos-1;
				output[counter]=source[xpos+ypos*width]; //getting an array out of bounds here on the source side
				counter++;
			}
		}
		return output;
	}
	
	public static float[] get_region2_padmirrored(float[] source,int x,int y,int rwidth,int rheight,int width,int height){
		float[] output=new float[rwidth*rheight];
		int counter=0;
		for(int i=y;i<(y+rheight);i++){
			int ypos=i;
			if(ypos<0) ypos=-ypos;
			if(ypos>=height) ypos=2*height-ypos-1;
			for(int j=x;j<(x+rwidth);j++){
				int xpos=j;
				if(xpos<0) xpos=-xpos;
				if(xpos>=width) xpos=2*width-xpos-1;
				output[counter]=source[xpos+ypos*width];
				counter++;
			}
		}
		return output;
	}
	
	public static float[] get_region2_padzeros(float[] source,int x,int y,int rwidth,int rheight,int width,int height){
		float[] output=new float[rwidth*rheight];
		for(int i=y;i<(y+rheight);i++){
			if(i<0) continue;
			if(i>=height) continue;
			for(int j=x;j<(x+rwidth);j++){
				if(j<0) continue;
				if(j>=width) continue;
				output[j-x+(i-y)*rwidth]=source[j+i*width];
			}
		}
		return output;
	}

}
