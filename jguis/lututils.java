package jguis;

import jalgs.interpolation;
import jalgs.jsim.rngs;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class lututils{
	// this class contains utility methods creating various luts
	public static byte[][] vis_spectrum(float wlstart,float wlend){
		// here we interpolate from the 380-780 spectrum below
		float[][] spectrum=vis_spectrum();
		byte[][] lut=new byte[3][255];
		float start=wlstart-380f;
		float end=wlend-380f;
		float inc=(end-start)/255f;
		for(int i=0;i<256;i++){
			float pos=start+inc*i;
			float R=interpolation.interp1D(spectrum[0],spectrum[0].length,pos);
			float G=interpolation.interp1D(spectrum[1],spectrum[1].length,pos);
			float B=interpolation.interp1D(spectrum[2],spectrum[2].length,pos);
			lut[0][i]=(byte)((int)(R*255.999f));
			lut[1][i]=(byte)((int)(G*255.999f));
			lut[2][i]=(byte)((int)(B*255.999f));
		}
		return lut;
	}

	public static float[][] vis_spectrum(){
		// this version goes from 380 to 780 in 1 nm increments
		// it is an approximation using linear interpolation
		// the code is adapted from
		// http://www.physics.sfasu.edu/astro/color.html
		// originally written by Dan Bruton
		float[][] lut=new float[3][401];
		for(int i=0;i<=400;i++){
			float R=0f;
			float G=0f;
			float B=0f;
			float wl=380f+i;
			if(wl<440f){
				R=(440f-wl)/(440f-380f);
				G=0f;
				B=1f;
			}else if(wl<490f){
				R=0f;
				G=(wl-440f)/(490f-440f);
				B=1f;
			}else if(wl<510f){
				R=0f;
				G=1f;
				B=(510f-wl)/(510f-490f);
			}else if(wl<580){
				R=(wl-510f)/(580f-510f);
				G=1f;
				B=0f;
			}else if(wl<645f){
				R=1f;
				G=(wl-645f)/(645f-580f);
				B=0f;
			}else{
				R=1f;
				G=0f;
				B=0f;
			}
			// now the attenuation factor
			float SSS=1f;
			if(wl>=700f){
				SSS=0.3f+0.7f*(780f-wl)/(780f-700f);
			}else if(wl<420f){
				SSS=0.3f+0.7f*(wl-380f)/(420f-380f);
			}
			lut[0][i]=SSS*R;
			lut[1][i]=SSS*G;
			lut[2][i]=SSS*B;
		}
		return lut;
	}

	public static byte[][] angle_lut(){
		// this lookup table varies slowly across angles including 360-0
		byte[][] ctable=new byte[3][256];
		// start with the transition from blue to cyan
		int off=0;
		for(int i=0;i<43;i++){
			ctable[2][off]=(byte)255;
			int temp=(int)(i*255.0f/43.0f);
			ctable[1][off]=(byte)temp;
			off++;
		}
		// now from cyan to green
		for(int i=0;i<43;i++){
			ctable[1][off]=(byte)255;
			int temp=255-(int)(i*255.0f/43.0f);
			ctable[2][off]=(byte)temp;
			off++;
		}
		// now green to yellow
		for(int i=0;i<42;i++){
			ctable[1][off]=(byte)255;
			int temp=(int)(i*255.0f/42.0f);
			ctable[0][off]=(byte)temp;
			off++;
		}
		// now yellow to red
		for(int i=0;i<43;i++){
			ctable[0][off]=(byte)255;
			int temp=255-(int)(i*255.0f/43.0f);
			ctable[1][off]=(byte)temp;
			off++;
		}
		// now red to magenta
		for(int i=0;i<43;i++){
			ctable[0][off]=(byte)255;
			int temp=(int)(i*255.0f/43.0f);
			ctable[2][off]=(byte)temp;
			off++;
		}
		// and finally magenta back to blue
		for(int i=0;i<42;i++){
			ctable[2][off]=(byte)255;
			int temp=255-(int)(i*255.0f/42.0f);
			ctable[0][off]=(byte)temp;
			off++;
		}
		return ctable;
	}

	public static byte[][] nice_lut(boolean whiteback){
		int[] r={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,
				200,204,208,212,216,220,224,228,232,236,240,244,248,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255};
		int[] g={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,
				144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,127,123,119,115,111,107,103,99,95,91,87,83,79,
				75,71,67,63,59,55,51,47,43,39,35,31,27,23,19,15,11,7,3,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,255};
		int[] b={0,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
				255,255,255,255,255,255,251,247,243,239,235,231,227,223,219,215,211,207,203,199,195,191,187,183,179,175,171,167,163,159,155,151,147,143,139,135,131,127,123,119,115,111,107,103,99,95,
				91,87,83,79,75,71,67,63,59,55,51,47,43,39,35,31,27,23,19,15,11,7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				0,0,0,0,0,0,0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,255};
		byte[][] lut=new byte[3][256];
		for(int i=0;i<256;i++){
			lut[0][i]=(byte)r[i];
			lut[1][i]=(byte)g[i];
			lut[2][i]=(byte)b[i];
		}
		if(whiteback){
			lut[0][0]=(byte)255; lut[1][0]=(byte)255; lut[2][0]=(byte)255;
			lut[0][255]=(byte)255; lut[1][255]=(byte)120; lut[2][255]=(byte)120;
		}
		return lut;
	}

	public static byte[][] ryan_lut(){
		// this design came from ryan colyer
		byte[][] ctable=new byte[3][256];
		for(int i=1;i<256;i++){
			float val=(i)/255.0f;
			float bright=128.0f;
			float sat=124.0f;
			float hue=(1.0f-0.95f*val-0.2f)*6.28318530717958f;
			if(val<0.125){
				bright-=(0.35355339059f-(float)Math.sqrt(val))*2.3f*sat;
				sat-=(0.35355339059f-(float)Math.sqrt(val))*2.3f*sat;
				hue=(1.0f-0.95f*0.125f-0.2f)*6.28318530717958f;
			}else if(val>0.875){
				bright+=(0.35355339059f-(float)Math.sqrt(1.0f-val))*2.3f*sat;
				sat-=(0.35355339059f-(float)Math.sqrt(1.0f-val))*2.3f*sat;
				hue=(1.0f-0.95f*0.875f-0.2f)*6.28318530717958f;
			}
			float r=bright+sat*(float)Math.cos(hue);
			float g=bright+sat*(float)Math.cos(hue+4.188790204786387f);
			float b=bright+sat*(float)Math.cos(hue+2.094395102393193f);
			ctable[0][i]=(byte)r;
			ctable[1][i]=(byte)g;
			ctable[2][i]=(byte)b;
		}
		return ctable;
	}

	public static byte[][] rand_lut(){
		// this is good for distinguishing closely colored objects or intensity
		// coded objects
		byte[][] ctable=new byte[3][256];
		rngs random=new rngs();
		for(int i=1;i<256;i++){
			double r=random.unidev(255.0,0.0);
			double g=random.unidev(255.0,0.0);
			double b=random.unidev(255.0,0.0);
			ctable[0][i]=(byte)r;
			ctable[1][i]=(byte)g;
			ctable[2][i]=(byte)b;
		}
		return ctable;
	}

	public static byte[][] log_lut(int color){
		// colors are gray,red,green,blue
		byte[][] ctable=new byte[3][256];
		float logmax=(float)Math.log(255.0);
		for(int i=1;i<256;i++){
			float logval=(float)Math.log(i);
			logval/=logmax;
			logval*=255.0f;
			logval=((int)logval);
			byte logvalb=(byte)logval;
			if(color==0){
				ctable[0][i]=logvalb;
				ctable[1][i]=logvalb;
				ctable[2][i]=logvalb;
			}else if(color==1){
				ctable[0][i]=logvalb;
			}else if(color==2){
				ctable[1][i]=logvalb;
			}else{
				ctable[2][i]=logvalb;
			}
		}
		return ctable;
	}

	public static byte[][] greenyellowred_lut(){
		byte[][] ctable=new byte[3][256];
		for(int i=1;i<128;i++){
			ctable[0][i]=(byte)(2*i);
			ctable[1][i]=(byte)255;
			ctable[2][i]=(byte)0;
		}
		for(int i=128;i<256;i++){
			ctable[0][i]=(byte)255;
			ctable[1][i]=(byte)(255-2*(i-128));
			ctable[2][i]=(byte)0;
		}
		return ctable;
	}

	public static byte[][] bluegreenred_lut(){
		byte[][] ctable=new byte[3][256];
		for(int i=1;i<64;i++){
			ctable[1][i]=(byte)(4*i);
			ctable[2][i]=(byte)255;
			ctable[0][i]=(byte)0;
		}
		for(int i=64;i<128;i++){
			ctable[1][i]=(byte)255;
			ctable[2][i]=(byte)(255-4*(i-64));
			ctable[0][i]=(byte)0;
		}
		for(int i=128;i<192;i++){
			ctable[0][i]=(byte)(4*(i-128));
			ctable[1][i]=(byte)255;
			ctable[2][i]=(byte)0;
		}
		for(int i=192;i<256;i++){
			ctable[0][i]=(byte)255;
			ctable[1][i]=(byte)(255-4*(i-192));
			ctable[2][i]=(byte)0;
		}
		return ctable;
	}

	public static boolean write_lut(byte[][] lut,String directory,String name){
		try{
			OutputStream outstream=new BufferedOutputStream(new FileOutputStream(directory+name));
			for(int i=0;i<256;i++){
				outstream.write(lut[0][i]);
			}
			for(int i=0;i<256;i++){
				outstream.write(lut[1][i]);
			}
			for(int i=0;i<256;i++){
				outstream.write(lut[2][i]);
			}
			outstream.close();
			return true;
		}catch(IOException e){
			return false;
		}
	}

}
