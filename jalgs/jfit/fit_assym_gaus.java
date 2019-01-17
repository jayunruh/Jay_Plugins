package jalgs.jfit;

import jalgs.gui_interface;
import jalgs.jstatistics;


public class fit_assym_gaus implements NLLSfitinterface_v3{
	//params are 0baseline, 1amp, 2avgstdev, 3stdevratio, 4xc, 5yc, 6angle
	public int xpts,ypts;
	public int nfit;
	public double[][] constraints;
	public double[] params;
	public float[] angles,c2plot;
	public gui_interface gui;
	public static String[] paramsnames={"baseline","amp","avgstdev","stdevratio","xc","yc","angle"};
	
	public fit_assym_gaus(int xpts,int ypts,gui_interface gui){
		this.xpts=xpts;
		this.ypts=ypts;
		this.gui=gui;
	}
	
	public float[] fitdata(float[] data,double[] stats,double stdevguess,float dangle,float[] weights,float guessx,float guessy){
		constraints=new double[2][7];
		nfit=0;
		float max=jstatistics.getstatistic("Max",data,null);
		float min=jstatistics.getstatistic("Min",data,null);
		float amp=max-min;
		//float maxpos=jstatistics.getstatistic("Maxpos",data,null);
		//int yc=(int)(maxpos/(float)xpts);
		//int xc=(int)(maxpos-yc*xpts);
		int yc=(int)guessy;
		int xc=(int)guessx;
		params=new double[]{(double)min,(double)amp,stdevguess,2.0,(double)xc,(double)yc,0.0};
		constraints[0][0]=min-0.1*amp; constraints[1][0]=min+0.1*amp;
		constraints[0][1]=0; constraints[1][1]=5.0*amp;
		constraints[0][2]=0.5; constraints[1][2]=2.0*stdevguess;
		constraints[0][3]=1.5; constraints[1][3]=4.0;
		constraints[0][4]=(double)(xc-2); constraints[1][4]=(double)(xc+2);
		constraints[0][5]=(double)(yc-2); constraints[1][5]=(double)(yc+2);
		constraints[0][6]=-0.1; constraints[1][6]=Math.PI+0.1;
		//float[] weights=null;
		//float[] weights=new float[data.length];
		//for(int i=0;i<weights.length;i++) weights[i]=data[i];
		NLLSfit_v3 fitclass=new NLLSfit_v3(this,0.0001,10,0.0);
		//search over the angles
		int[] fixes={0,0,0,0,0,0,1};
		float dx=0.5f;
		double[] minparams=params.clone();
		double minc2=Float.MAX_VALUE;
		c2plot=new float[(int)(180.0/dangle)];
		angles=new float[(int)(180.0/dangle)];
		int counter=0;
		for(float i=0.0f;i<180.0f;i+=dangle){
			double minc22=Float.MAX_VALUE;
			double[] minparams2=new double[params.length];
			for(float j=(yc-1.5f);j<=(yc+1.5f);j+=dx){
				for(float k=(xc-1.5f);k<=(xc+1.5f);k+=dx){
    				params[6]=Math.toRadians((double)i);
    				params[5]=(double)j;
    				params[4]=(double)k;
    				double[] tempparams=params.clone();
    				fitclass.fitdata(tempparams,fixes,data,weights,stats,false);
    				if(stats[1]<minc22){
    					minc22=stats[1];
    					minparams2=tempparams.clone();
    				}
				}
			}
			if(minc22<minc2){
				minc2=minc22;
				minparams=minparams2.clone();
			}
			angles[counter]=i;
			c2plot[counter]=(float)minc22;
			counter++;
			showresults("angle = \t"+i+"\t c2 = \t"+minc22);
		}
		params=minparams.clone();
		constraints[0][4]=params[4]-2.0; constraints[1][4]=params[4]+2.0;
		constraints[0][5]=params[5]-2.0; constraints[1][5]=params[5]+2.0;
		constraints[0][6]=params[6]-0.1; constraints[1][6]=params[6]+0.1;
		//now minimize the angle, fixing everything but the amplitudes
		fixes=new int[]{0,0,1,1,1,1,0};
		//fitclass.fitdata(params,fixes,data,weights,stats,false);
		//now minimize everything at once
		//fixes=new int[]{0,0,0,0,0,0,0};
		return fitclass.fitdata(params,fixes,data,weights,stats,true);
	}
	
	public double[] guessParams(float[] data,double stdevguess){
		float max=jstatistics.getstatistic("Max",data,null);
		float min=jstatistics.getstatistic("Min",data,null);
		float amp=max-min;
		float maxpos=jstatistics.getstatistic("Maxpos",data,null);
		int yc=(int)(maxpos/(float)xpts);
		int xc=(int)(maxpos-yc*xpts);
		double[] temp=new double[]{(double)min,(double)(max-min),stdevguess,2.0,(double)xc,(double)yc,0.0};
		return temp;
	}

	public double[] fitfunc(double[] params){
		double cost=Math.cos(params[6]);
		double sint=Math.sin(params[6]);
		double sin2t=Math.sin(2.0*params[6]);
		double stdev2=params[2]*2.0/(params[3]+1.0);
		double stdev1=params[2]*2.0-stdev2;
		double a=cost*cost/(2.0*stdev1*stdev1)+sint*sint/(2.0*stdev2*stdev2);
		double b=sin2t/(4.0*stdev2*stdev2)-sin2t/(4.0*stdev1*stdev1);
		double c=sint*sint/(2.0*stdev1*stdev1)+cost*cost/(2.0*stdev2*stdev2);
		double[] func=new double[xpts*ypts];
		for(int i=0;i<ypts;i++){
			for(int j=0;j<xpts;j++){
				func[j+i*xpts]=params[0]+params[1]*Math.exp(-a*(j-params[4])*(j-params[4])+2.0*b*(j-params[4])*(i-params[5])-c*(i-params[5])*(i-params[5]));
			}
		}
		return func;
	}

	public void applyconstraints(double[] params,int[] fixes){
		if(constraints!=null){
    		for(int i=0;i<constraints[0].length;i++){
    			if(fixes[i]==0){
    				if(Double.isNaN(params[i])) params[i]=constraints[0][i];
    				if(Double.isInfinite(params[i])) params[i]=constraints[1][i];
    				if(params[i]<constraints[0][i]) params[i]=constraints[0][i];
    				if(params[i]>constraints[1][i]) params[i]=constraints[1][i];
    			}
    		}
		}
	}

	public double[][] derivfunc(double[] params,int[] fixes,double[] fit){
		//params are 0baseline, 1amp, 2avgstdev, 3stdevratio, 4xc, 5yc, 6angle
		double[][] deriv=new double[xpts*ypts][7];
		double t=params[6];
		double sr=params[3];
		double sa=params[2];
		double xc=params[4];
		double yc=params[5];
		double cost=Math.cos(t);
		double sint=Math.sin(t);
		double cost2=cost*cost;
		double sint2=sint*sint;
		double cos2t=Math.cos(2.0*t);
		double sin2t=Math.sin(2.0*t);
		//double sin4t=Math.sin(4.0*t);
		double sradd=(1.0+sr);
		double sradd2=sradd*sradd;
		double temp=1.0/(sa*sa*sr*sr);
		double sr3=sr*sr*sr;

		for(int i=0;i<ypts;i++){
			double yd=(double)i-yc;
			double yd2=yd*yd;
			int off=i*xpts;
			for(int j=0;j<xpts;j++){
				double xd=(double)j-xc;
				double xd2=xd*xd;
				int ind=j+off;
				double nobase=fit[ind]-params[0];
    			deriv[ind][0]=1.0f;
    			deriv[ind][1]=nobase/params[1];
    			deriv[ind][2]=0.25*sradd2*(temp/sa)*nobase*((xd2+yd2*sr*sr)*cos2t+(sr*xd2+yd2)*sin2t+(1.0-sr)*sradd*xd*yd*sin2t);
    			deriv[ind][3]=nobase*0.25*sradd*(temp/sr)*((xd2-yd2*sr3)*cost2+(yd2-xd2*sr3)*sint2+(1.0+sr3)*xd*yd*sin2t);
    			deriv[ind][4]=nobase*sradd2*0.125*temp*(2.0*xd*cost2+2.0*xd*sint2*sr*sr-sradd*(sr-1.0)*yd*sin2t);
    			deriv[ind][5]=nobase*sradd2*0.125*temp*(2.0*yd*cost2*sr*sr+2.0*yd*sint2-sradd*(sr-1.0)*xd*sin2t);
    			deriv[ind][6]=0.25*temp*nobase*sradd2*((xd2-sr*sr*yd2)*cost*sint-(yd2+xd2*sr*sr)*cost*sint+(sr-1.0)*sradd*xd*yd*cos2t);
			}
		}
		return deriv;
	}

	public void showresults(String results){
		gui.showMessage(results);
	}

}
