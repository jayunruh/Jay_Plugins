package jalgs.jfft;

import jalgs.jstatistics;

public class jwavelets{
	
	public static String[] coefs={"2","4","6","8","10","12","14","16","18","20"};
	
	public static double[][] coef={{0.707106781,0.707106781},
		{0.482962913,0.836516304,0.22414387,-0.129409521},
		{0.332670553,0.80689151,0.4598775,-0.135011023,-0.085441275,0.035226292},
		{0.230377684,0.714846169,0.630881303,-0.027983756,-0.187034708,0.030841364,0.032882992,-0.010597396},
		{0.111540742,0.494623889,0.75113391,0.31525035,-0.226264695,-0.129766865,0.097501604,0.027522866,-0.031582041,0.000553842,0.004777258,-0.001077301},
		{0.077852053,0.396539313,0.729132083,0.469782282,-0.143906001,-0.224036182,0.071309385,0.080612612,-0.038029935,-0.01657454,0.012550997,0.000429578,-0.001801641,0.000353714},
		{0.054415841,0.312871593,0.675630737,0.585354682,-0.015829109,-0.284015541,0.000472485,0.128747429,-0.0173693,-0.044088256,0.013981028,0.008746091,-0.004870353,-0.00039174,0.000675449,-0.000117477},
		{0.038077948,0.243834673,0.604823124,0.657288078,0.133197387,-0.293273784,-0.096840784,0.148540748,0.030725681,-0.067632826,0.000250947,0.022361665,-0.004723205,-0.004281504,0.001847647,0.000230386,-0.000251963,3.93473E-05},
		{0.02667006,0.188176798,0.527201189,0.68845904,0.281172343,-0.249846423,-0.195946276,0.127369342,0.093057367,-0.071394146,-0.029457538,0.033212671,0.003606554,-0.010733174,0.001395352,0.001992405,-0.000685857,-0.000116467,9.35887E-05,-1.32642E-05}
	};
	
	public static float[][] coeff={{0.707106781f,0.707106781f},
		{0.482962913f,0.836516304f,0.22414387f,-0.129409521f},
		{0.332670553f,0.80689151f,0.4598775f,-0.135011023f,-0.085441275f,0.035226292f},
		{0.230377684f,0.714846169f,0.630881303f,-0.027983756f,-0.187034708f,0.030841364f,0.032882992f,-0.010597396f},
		{0.111540742f,0.494623889f,0.75113391f,0.31525035f,-0.226264695f,-0.129766865f,0.097501604f,0.027522866f,-0.031582041f,0.000553842f,0.004777258f,-0.001077301f},
		{0.077852053f,0.396539313f,0.729132083f,0.469782282f,-0.143906001f,-0.224036182f,0.071309385f,0.080612612f,-0.038029935f,-0.01657454f,0.012550997f,0.000429578f,-0.001801641f,0.000353714f},
		{0.054415841f,0.312871593f,0.675630737f,0.585354682f,-0.015829109f,-0.284015541f,0.000472485f,0.128747429f,-0.0173693f,-0.044088256f,0.013981028f,0.008746091f,-0.004870353f,-0.00039174f,0.000675449f,-0.000117477f},
		{0.038077948f,0.243834673f,0.604823124f,0.657288078f,0.133197387f,-0.293273784f,-0.096840784f,0.148540748f,0.030725681f,-0.067632826f,0.000250947f,0.022361665f,-0.004723205f,-0.004281504f,0.001847647f,0.000230386f,-0.000251963f,3.93473E-05f},
		{0.02667006f,0.188176798f,0.527201189f,0.68845904f,0.281172343f,-0.249846423f,-0.195946276f,0.127369342f,0.093057367f,-0.071394146f,-0.029457538f,0.033212671f,0.003606554f,-0.010733174f,0.001395352f,0.001992405f,-0.000685857f,-0.000116467f,9.35887E-05f,-1.32642E-05f}
	};
	
	public static float[] waveletFilter(float[] data,float startpos,float endpos,float thresh,float[] coef){
		float[] trans=data.clone();
		waveletTransform(trans,coef,false);
		for(int i=(int)(startpos*data.length);i<(int)(endpos*data.length);i++){
			trans[i]=0.0f;
		}
		float max=jstatistics.getstatistic("Max",trans,null);
		int countfilt=0;
		for(int i=0;i<trans.length;i++) if(Math.abs(trans[i])<thresh*max){trans[i]=0.0f; countfilt++;}
		float ffilt=(float)countfilt/(float)data.length;
		waveletTransform(trans,coef,true);
		return trans;
	}
	
	public static float[] waveletFilter2D(float[] data,int width,int height,float startpos,float endpos,float thresh,float[] coef){
		float[] trans=data.clone();
		waveletTransform2D(trans,width,height,coef,false);
		for(int i=(int)(startpos*height);i<(int)(endpos*height);i++){
			for(int j=(int)(startpos*width);j<(int)(endpos*width);j++){
				trans[j+i*width]=0.0f;
			}
		}
		float max=jstatistics.getstatistic("Max",trans,null);
		int countfilt=0;
		for(int i=0;i<trans.length;i++) if(Math.abs(trans[i])<thresh*max){trans[i]=0.0f; countfilt++;}
		float ffilt=(float)countfilt/(float)data.length;
		waveletTransform2D(trans,width,height,coef,true);
		return trans;
	}
	
	public static void waveletTransform2D(float[] data,int width,int height,float[] coef,boolean inverse){
		float[][] columns=new float[width][height];
		int pos=0;
		for(int i=0;i<height;i++){
			float[] line=new float[width];
			System.arraycopy(data,pos,line,0,width);
			waveletTransform(line,coef,inverse);
			for(int j=0;j<width;j++) columns[j][i]=line[j];
			pos+=width;
		}
		for(int i=0;i<width;i++){
			waveletTransform(columns[i],coef,inverse);
			pos=i;
			for(int j=0;j<height;j++){data[pos]=columns[i][j]; pos+=width;}
		}
	}

	public static void waveletTransform(float[] data,float[] coef,boolean inverse){
		if(!inverse){
			for(int i=data.length;i>=coef.length;i>>=1){
				waveletFilter(data,coef,i,false);
			}
		} else {
			for(int i=coef.length;i<=data.length;i<<=1){
				waveletFilter(data,coef,i,true);
			}
		}
	}


	public static void waveletFilter(float[] data,float[] coef,int len,boolean inverse){
		float[] filtered=new float[len];
		int halflen=len/2;
		float[] rcoef=revCoef(coef);
		int len1=len-1;
		int lenmod=coef.length*len;
		//float[] filtered=new float[len];
		if(!inverse){
			int index=0;
			for(int j=0;j<len;j+=2){
				for(int i=0;i<coef.length;i++){
					int shift=len1&(j+1+lenmod+i+1);
					//int shift=wrap(j+1+lenmod+i+1,len);
					filtered[index]+=coef[i]*data[shift];
					filtered[index+halflen]+=rcoef[i]*data[shift];
				}
				index++;
			}
		} else {
			int index=0;
			for(int j=0;j<len;j+=2){
				float t=data[index];
				float t2=data[index+halflen];
				for(int i=0;i<coef.length;i++){
					int shift=len1&(j+1+lenmod+i+1);
					//int shift=wrap(j+1+lenmod+i+1,len);
					filtered[shift]+=coef[i]*t;
					filtered[shift]+=rcoef[i]*t2;
				}
				index++;
			}
		}
		System.arraycopy(filtered,0,data,0,len);
	}

	public static double[] revCoef(double[] coef){
		double[] coef2=new double[coef.length];
		double mult=-1.0;
		for(int i=0;i<coef.length;i++){
			if(mult<0) coef2[coef.length-i-1]=-coef[i];
			else coef2[coef.length-i-1]=coef[i];
			mult=-mult;
		}
		return coef2;
	}

	public static float[] revCoef(float[] coef){
		float[] coef2=new float[coef.length];
		double mult=-1.0;
		for(int i=0;i<coef.length;i++){
			if(mult<0) coef2[coef.length-i-1]=-coef[i];
			else coef2[coef.length-i-1]=coef[i];
			mult=-mult;
		}
		return coef2;
	}

	public static int wrap(int var,int length){
		int var2=var;
		if(var2<0){
			int rem=Math.abs(var)%length;
			var2=length-rem;
		} else if(var2>(length-1)){
			var2%=length;
		}
		return var2;
	}

}
