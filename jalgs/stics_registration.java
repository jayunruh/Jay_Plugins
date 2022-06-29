package jalgs;

import jalgs.jfft.crosscorr2D;

public class stics_registration{
	//this class uses spatio-temporal image correlation to generate a displacement map from a first image to a second
	//that map is then interpolated bicubically and used to align the first image to the second
	//NaN values are omitted in the calculation
	
	public int width,height,typeindex,newwidth,newheight;
	public float xmin,xmax,ymin,ymax;
	crosscorr2D ccclass;
	public int fftdim;
	public float[] xdisplacements,ydisplacements;
	public gui_interface gui;

	public stics_registration(int width,int height,int typeindex,gui_interface gui){
		// TODO Auto-generated constructor stub
	}

}
