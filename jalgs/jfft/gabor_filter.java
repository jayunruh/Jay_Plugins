package jalgs.jfft;

import jalgs.algutils;

public class gabor_filter{
	
	public float[][] kernel;
	public int N,nthetas,nphis;
	public int[] fftindex1,fftindex2;

	public gabor_filter(float lambda,float sigma,float[] theta,float[] phi,float gamma){
		N=1+2*(int)Math.ceil(2.5*sigma/gamma);
		nthetas=theta.length; nphis=phi.length;
		kernel=new float[theta.length*phi.length][];
		//the function is exp((-xp*xp-g*g*yp*yp)/(2*sig*sig))*cos(2*pi*xp/lambda+phi)
		//xp=x*cos(theta)+y*sin(theta)
		//yp=-x*sin(theta)+y*cos(theta)
		int counter=0;
		for(int k=0;k<theta.length;k++){
			for(int l=0;l<phi.length;l++){
				kernel[counter]=get_kernel(N,lambda,sigma,theta[k],phi[l],gamma);
        		counter++;
			}
		}
	}
	
	public gabor_filter(float lambda,float sigma,float[] theta,float[] phi,float gamma,int N){
		this.N=N;
		nthetas=theta.length; nphis=phi.length;
		kernel=new float[theta.length*phi.length][];
		//the function is exp((-xp*xp-g*g*yp*yp)/(2*sig*sig))*cos(2*pi*xp/lambda+phi)
		//xp=x*cos(theta)+y*sin(theta)
		//yp=-x*sin(theta)+y*cos(theta)
		int counter=0;
		for(int k=0;k<theta.length;k++){
			for(int l=0;l<phi.length;l++){
				kernel[counter]=get_kernel(N,lambda,sigma,theta[k],phi[l],gamma);
        		counter++;
			}
		}
	}
	
	public static float[] get_kernel(int size,float lambda,float sigma,float theta,float phi,float gamma){
		float[] kernel2=new float[size*size];
		float xc=0.5f*size;
		float yc=xc;
		double avg=0.0;
		//float sumneg=0.0f;
		float sintheta=(float)Math.sin(theta);
		float costheta=(float)Math.cos(theta);
		for(int i=0;i<size;i++){
			for(int j=0;j<size;j++){
				float x=j-xc;
				float y=i-yc;
				float temp=get_kernel(x,y,lambda,sigma,costheta,sintheta,phi,gamma);
				kernel2[j+i*size]=temp;
				//if(temp>0) sumpos+=temp;
				//else sumneg+=Math.abs(temp);
				avg+=(double)temp;
			}
		}
		avg/=(double)size*(double)size;
		/*float meansum=(sumpos+sumneg)/2.0f;
		if(meansum>0){
			sumpos/=meansum;
			sumneg/=meansum;
		}*/
		for(int i=0;i<size*size;i++){
			//if(kernel2[i]>0.0f) kernel2[i]*=sumpos;
			//if(kernel2[i]<0.0f) kernel2[i]*=sumneg;
			kernel2[i]-=(float)avg;
		}
		return kernel2;
	}
	
	public static float get_kernel(float x,float y,float lambda,float sigma,float costheta,float sintheta,float phi,float gamma){
		/*float xp=x*costheta-y*sintheta;
		float yp=y*costheta+x*sintheta;
		return (float)Math.exp(-(xp*xp+gamma*gamma*yp*yp)/(2.0*sigma*sigma))*(float)Math.cos(-2.0*Math.PI*yp/lambda-phi);*/
		float xp=x*costheta+y*sintheta;
		float yp=y*costheta-x*sintheta;
		return (float)Math.exp(-(xp*xp+gamma*gamma*yp*yp)/(2.0*sigma*sigma))*(float)Math.cos(2.0*Math.PI*xp/lambda+phi);
	}
	
	public static float get_kernel(float x,float y,float lambda,float sigma,float theta,float phi,float gamma){
		return get_kernel(x,y,lambda,sigma,(float)Math.cos(theta),(float)Math.sin(theta),phi,gamma);
	}
	
	public static float calc_sigma(float lambda,float bandwidth){
		double ratio=Math.sqrt(0.5*Math.log(2))*(1.0+Math.pow(2.0,bandwidth))/(Math.PI*(Math.pow(2.0,bandwidth)-1.0));
		return (float)ratio*lambda;
	}
	
	public static float calc_lambda(float sigma,float bandwidth){
		double ratio=Math.sqrt(0.5*Math.log(2))*(1.0+Math.pow(2.0,bandwidth))/(Math.PI*(Math.pow(2.0,bandwidth)-1.0));
		return sigma/(float)ratio;
	}
	
	public float[] applyFilter(float[] image,int width,int height,float mult,String stat,boolean trunc){
		float[][] filtered=getFilterStack(image,width,height);
		float[] proj=algutils.get_stack_proj_stat(stat,filtered,width,height,filtered.length,null);
		for(int i=0;i<proj.length;i++){
			proj[i]*=mult;
			if(trunc){
				if(proj[i]<0.0f) proj[i]=0.0f;
			}
		}
		return proj;
	}
	
	public float[][] getFilterStack(float[] image,int width,int height){
		//convolution2D cc=new convolution2D();
		fftindex1=fftutils.get_best_index(width,true,19);
		fftindex2=fftutils.get_best_index(height,true,19);
		convolution2D cc=new convolution2D(fftindex1[1],fftindex2[1],fftindex1[0],fftindex2[0]);
		float[] image2=image.clone();
		//pad if necessary
		if(fftindex1[1]!=width || fftindex2[1]!=height){
			image2=padding.pad_xy_mirrored(image2,width,height,fftindex1[1],fftindex2[1],false);
		}
		float[] im=new float[fftindex1[1]*fftindex2[1]];
		cc.fft.dorealfft2D(image2,im,false);
		//we will almost certainly need to pad the kernels as we go
		float[][] filtered=new float[nthetas*nphis][];
		for(int i=0;i<nthetas;i++){
			for(int j=0;j<nphis;j++){
				float[] kernel2=padding.pad_xy_zeros(kernel[j+i*nphis],N,N,fftindex1[1],fftindex2[1]);
				kernel2=(new manipulate_quads()).shiftxycenter(kernel2,fftindex1[1],fftindex2[1]);
				float[] temp=cc.convolve2D(kernel2,image2,im);
				//temp=(new manipulate_quads()).shiftxycenter(temp,fftindex1[1],fftindex2[1]);
				filtered[j+i*nphis]=algutils.get_region(temp,fftindex1[1]/2,fftindex2[1]/2,width,height,fftindex1[1],fftindex2[1]);
			}
		}
		return filtered;
	}
	
}
