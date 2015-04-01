package jalgs.jfit;

import jalgs.algutils;

public class linear_unmix{
	//here we have some simple utilities for linear unmixing 
	//(basically an automation class for linleastsquares
	public float[][] spectra;
	public int startch,endch;
	public linleastsquares lls;
	
	public linear_unmix(float[][] spectra,int startch,int endch){
		this.spectra=spectra;
		this.startch=startch;
		this.endch=endch;
		lls=new linleastsquares(spectra,false,startch,endch);
	}
	
	public void update_spectra(float[][] spectra,int startch,int endch){
		this.spectra=spectra;
		this.startch=startch;
		this.endch=endch;
		lls=new linleastsquares(spectra,false,startch,endch);
	}
	
	public float[][] unmix(Object[] stack,int npix,boolean truncneg){
		float[][] output=new float[spectra.length][npix];
		for(int i=0;i<npix;i++){
			float[] col=algutils.convert_arr_float(algutils.get_stack_col(stack,npix,1,i,stack.length));
			double[] contributions=lls.fitdata(col,null);
			for(int j=0;j<contributions.length;j++){
				if(truncneg && contributions[j]<0.0f) contributions[j]=0.0f;
				output[j][i]=(float)contributions[j];
			}
		}
		return output;
	}
	
	public float[][] get_resid(Object[] stack,float[][] contr,int npix){
		float[][] output=new float[stack.length][npix];
		for(int i=0;i<npix;i++){
			float[] col=algutils.convert_arr_float(algutils.get_stack_col(stack,npix,1,i,stack.length));
			double[] contrcol=algutils.convert_arr_double(algutils.get_stack_col(contr,npix,1,i,stack.length));
			float[] resid=lls.get_fresid(contrcol,col,null);
			for(int j=0;j<resid.length;j++){
				output[j][i]=resid[j];
			}
		}
		return output;
	}
	
	public float[] get_c2(Object[] stack,float[][] contr,int npix){
		float[] output=new float[npix];
		for(int i=0;i<npix;i++){
			float[] col=algutils.convert_arr_float(algutils.get_stack_col(stack,npix,1,i,stack.length));
			double[] contrcol=algutils.convert_arr_double(algutils.get_stack_col(contr,npix,1,i,stack.length));
			output[i]=(float)lls.get_c2(contrcol,col,null);
		}
		return output;
	}

}
