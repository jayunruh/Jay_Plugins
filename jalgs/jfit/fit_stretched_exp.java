package jalgs.jfit;

import jalgs.algutils;
import jalgs.gui_interface;
import jalgs.jstatistics;

public class fit_stretched_exp implements NLLSfitinterface_v2{
	public float[] xvals,yvals;
	public double[][] constraints;
	public double[] sparams;
	public boolean redirect;
	gui_interface gi;
	public String[] paramsnames= {"baseline","amp","EC50","alpha","xshift"};
	
	public fit_stretched_exp(float[] xvals,float[] yvals,float minec50,float maxec50,gui_interface gi) {
		this.xvals=xvals;
		this.yvals=yvals;
		this.gi=gi;
		redirect=false;
		double[] guesses=initEC50(minec50,maxec50);
		constraints=new double[2][5];
		sparams=new double[]{0.0,1.0,guesses[0],guesses[1],0.0};
		constraints[0][0]=sparams[0]-sparams[1]; constraints[1][0]=sparams[0]+sparams[1];
		constraints[0][1]=0.1*sparams[1]; constraints[1][1]=10.0*sparams[1];
		constraints[0][2]=minec50; constraints[1][2]=maxec50;
		constraints[0][3]=0.1; constraints[1][3]=10.0;
		constraints[0][4]=0.0; constraints[1][4]=maxec50;
	}
	
	public Object[] autoFit(int[] fixes,double[] params,float[] weights,boolean errors,boolean output) {
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0.0001,50,0.1);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,constraints,yvals,weights,stats,output);
		double[] bestparams=params.clone();
		float[] errstdevs=null;
		if(errors) {
    		monte_carlo_errors_v2 erclass=new monte_carlo_errors_v2(this,0.0001,50,false,0.1);
    		redirect=true;
    		double[][] errs=erclass.geterrors(params,fixes,constraints,yvals,weights,100);
    		redirect=false;
    		errstdevs=new float[errs.length];
    		for(int i=0;i<errs.length;i++) {
    			errstdevs[i]=jstatistics.getstatistic("StDev",algutils.convert_arr_float(errs[i]),null);
    		}
		}
		return new Object[] {bestparams,stats,fit,errstdevs};
	}
	
	public double getc2(double[] fit){
		double c2=0.0;
		for(int i=0;i<fit.length;i++){
			c2+=(fit[i]-(double)yvals[i])*(fit[i]-(double)yvals[i]);
		}
		c2/=(double)(fit.length-3);
		return c2;
	}

	public double[] initEC50(double min,double max){
		double minc2=-1.0;
		double minec50=min;
		double minamp=0.0;
		double minoff=0.0;
		double minalpha=2.0;
		//added to check alpha values 2, 4, and 6
		for(double ec50=min;ec50<=max;ec50*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,ec50,2.0,0.0});
			//double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double[] ampoff={1.0,0.0};
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minec50=ec50; minamp=ampoff[0]; minoff=ampoff[1]; minalpha=2.0;}
		}
		for(double ec50=min;ec50<=max;ec50*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,ec50,4.0,0.0});
			//double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double[] ampoff={1.0,0.0};
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minec50=ec50; minamp=ampoff[0]; minoff=ampoff[1]; minalpha=4.0;}
		}
		for(double ec50=min;ec50<=max;ec50*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,ec50,6.0,0.0});
			//double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double[] ampoff={1.0,0.0};
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minec50=ec50; minamp=ampoff[0]; minoff=ampoff[1]; minalpha=6.0;}
		}
		return new double[]{minec50,minalpha,minc2,minamp,minoff};
	}

	public double[] fitfunc(double[] params){
		//params are 0baseline,1amp,2ec50,3alpha,4xshift
		double[] func=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			double x=(double)xvals[i];
			if(x<params[4]) func[i]=params[0];
			else {
				double temp=(x-params[4])/(params[2]-params[4]);
				func[i]=params[0]+params[1]*(1.0-Math.exp(-0.693*Math.pow(temp,params[3])));
			}
		}
		return func;
	}

	public void showresults(String results){
		if(!redirect)
			gi.showMessage(results);
	}
}
