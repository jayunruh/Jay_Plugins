package jalgs.jfit;


public class fit_assym_gaus implements NLLSfitinterface_v3{
	//params are 0baseline, 1amp, 2avgstdev, 3stdevratio, 4xc, 5yc, 6angle
	public int xpts,ypts;
	//public NLLSfitinterface_v2 fitclass;
	
	/*public fit_assym_gaus(NLLSfitinterface_v2 fitclass){
		this.fitclass=fitclass;
	}*/

	public double[] fitfunc(double[] params){
		double cost=Math.cos(params[6]);
		double sint=Math.sin(params[6]);
		double sin2t=Math.sin(2.0*params[6]);
		double stdev2=params[2]*2.0/(params[3]+1.0);
		double stdev1=params[2]*2.0-stdev2;
		double a=cost*cost/(2.0*stdev1*stdev1)+sint*sint/(2.0*stdev2*stdev2);
		double b=sin2t/(4.0*stdev2*stdev2)-sin2t/(4.0*stdev2*stdev2);
		double c=sint*sint/(2.0*stdev1*stdev1)+cost*cost/(2.0*stdev2*stdev2);
		double[] func=new double[xpts*ypts];
		for(int i=0;i<ypts;i++){
			for(int j=0;j<xpts;j++){
				func[j+i*xpts]=params[0]+params[1]*Math.exp(-a*(j-params[4])*(j-params[4])+2.0*b*(j-params[4])*(i-params[5])+c*(i-params[5])*(i-params[5]));
			}
		}
		return func;
	}

	public void applyconstraints(double[] params,int[] fixes){
		// TODO Auto-generated method stub

	}

	public double[][] derivfunc(double[] params,int[] fixes,double[] fit){
		// TODO Auto-generated method stub
		return null;
	}

	public void showresults(String results){
		// TODO Auto-generated method stub
		//IJ.log(results);
		//fitclass.showresults(results);
		System.out.println(results);
	}

}
