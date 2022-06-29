package jalgs.jfit;

import java.util.List;

import jalgs.algutils;
import jalgs.jstatistics;

public class fit_binding implements NLLSfitinterface_v2{
	//this class has fitting utilities for binding functions--implement logistic first
	public float[] yvals,xvals,weights;
	double[][] constraints;
	double[] params;
	public static String[] paramsnames={"baseline","amp","Kd","n"};

	public fit_binding(float[] xvals,float[] yvals,double minkd,double maxkd,double guessn,float[] errs){
		this.xvals=xvals;
		this.yvals=yvals;
		weights=null;
		if(errs!=null){
			weights=new float[errs.length];
			for(int i=0;i<errs.length;i++) weights[i]=1.0f/(errs[i]*errs[i]);
		}
		params=null;
		double[] guesses=initKd(this.yvals,minkd,maxkd,guessn);
		params=new double[]{guesses[3],guesses[2],guesses[0],guessn};
		if(params[1]<0.0){
			params[1]=jstatistics.getstatistic("Max",yvals,null)-jstatistics.getstatistic("Min",yvals,null);
		}
		constraints=new double[2][4];
		constraints[0][0]=params[0]-params[1]; constraints[1][0]=params[0]+params[1];
		constraints[0][1]=0.2*params[1]; constraints[1][1]=5.0*params[1];
		constraints[0][2]=minkd; constraints[1][2]=maxkd;
		constraints[0][3]=0.1; constraints[1][3]=10.0;
	}
	
	public fit_binding(List<Float> xvals,List<Float> yvals,double minkd,double maxkd,double guessn,List<Float> errs1){
		this.xvals=algutils.convert_arr_float(xvals);
		this.yvals=algutils.convert_arr_float(yvals);
		weights=null;
		if(errs1!=null){
			float[] errs=algutils.convert_arr_float(errs1);
			weights=new float[errs.length];
			for(int i=0;i<errs.length;i++) weights[i]=1.0f/(errs[i]*errs[i]);
		}
		params=null;
		double[] guesses=initKd(this.yvals,minkd,maxkd,guessn);
		params=new double[]{guesses[3],guesses[2],guesses[0],guessn};
		if(params[1]<0.0){
			params[1]=jstatistics.getstatistic("Max",yvals,null)-jstatistics.getstatistic("Min",yvals,null);
		}
		constraints=new double[2][4];
		constraints[0][0]=params[0]-params[1]; constraints[1][0]=params[0]+params[1];
		constraints[0][1]=0.2*params[1]; constraints[1][1]=5.0*params[1];
		constraints[0][2]=minkd; constraints[1][2]=maxkd;
		constraints[0][3]=0.1; constraints[1][3]=10.0;
	}
	
	public double[][] runFit(double[] params,int[] fixes,boolean output){
		NLLSfit_v2 fitclass=new NLLSfit_v2(this,0.0001,20,0.1);
		double[] stats=new double[2];
		if(params!=null) this.params=params.clone();
		int[] fixes1=new int[paramsnames.length];
		if(fixes!=null) fixes1=fixes.clone();
		float[] fit=fitclass.fitdata(this.params,fixes1,constraints,yvals,weights,stats,output);
		return new double[][]{this.params,stats,algutils.convert_arr_double(fit)};
	}
	
	public double[] initKd(float[] data,double min,double max,double guessn){
		double minc2=-1.0;
		double minkd=min;
		double minamp=0.0;
		double minoff=0.0;
		for(double kd=min;kd<=max;kd*=1.01){
			double[] func=fitfunc(new double[]{0.0,1.0,kd,guessn});
			double[] ampoff=(new linleastsquares()).get_amp_offset(func,yvals,true);
			double c2=(new linleastsquares()).get_amp_offset_c2(func,yvals,ampoff);
			if(minc2<0 || c2<minc2){minc2=c2; minkd=kd; minamp=ampoff[0]; minoff=ampoff[1];}
		}
		return new double[]{minkd,minc2,minamp,minoff};
	}
	
	public double[] fitfunc(double[] params){
		//params are 0baseline,1amp,2kd,3n
		double[] func=new double[xvals.length];
		for(int i=0;i<xvals.length;i++){
			double x=(double)xvals[i];
			func[i]=params[0]+params[1]/(1.0+Math.exp(-params[3]*(x-params[2])));
		}
		return func;
	}
	
	public void showresults(String results){
		System.out.println(results);
	}

}
