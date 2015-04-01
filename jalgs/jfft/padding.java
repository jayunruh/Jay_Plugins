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
	
	

}
