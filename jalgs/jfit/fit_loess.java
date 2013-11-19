package jalgs.jfit;

public class fit_loess{
	//this class performs local regression parabolic (loess) or linear (lowess) smoothing
	//assume equally spaced x values
	//weighting function is tricube
	public float[] centerweights;
	public int span,length,halfspan;
	public fitpoly fp;
	
	public fit_loess(int length,int span){
		this(length,span,2);
	}

	public fit_loess(int length,int span,int order){
		//span is truncated to an odd number
		this.length=length;
		this.span=span;
		if(this.span>=length) this.span=length;
		if(this.span%2==0) this.span--;
		float[] xvals=new float[span];
		for(int i=0;i<this.span;i++) xvals[i]=(float)i;
		int torder=order;
		if(torder<1) torder=1;
		//if(torder>3) torder=3;
		fp=new fitpoly(torder,xvals,false);
		halfspan=(int)(0.5f*(float)this.span);
		initweights();
	}
	
	public float[] getfit(float[] yvals){
		float[] fit=new float[length];
		for(int i=0;i<length;i++){
			float[][] weightsdata=getweightsdata(i,yvals);
			double[] coefs=fp.fitdatapolyd(weightsdata[1],weightsdata[0]);
			int offset=(int)weightsdata[2][0];
			fit[i]=fp.getfitpt(coefs,offset);
		}
		return fit;
	}
	
	public void initweights(){
		centerweights=new float[span];
		int dist=halfspan;
		int center=dist+1;
		for(int i=0;i<span;i++){
			int dist2=i-center;
			if(dist2<0) dist2=-dist2;
			float ratio=(float)dist2/(float)dist;
			float temp=1.0f-ratio*ratio*ratio;
			centerweights[i]=temp*temp*temp;
		}
	}
	
	public float[][] getweightsdata(int center,float[] yvals){
		float[] weights=new float[span];
		float[] data=new float[span];
		int dist=halfspan;
		int start=center-halfspan;
		int end=center+halfspan;
		if(start<0){
			dist=span-center;
			start=0;
			end=span-1;
			System.arraycopy(yvals,start,data,0,span);
		}
		else if(end>=length){
			dist=center-(length-span-1);
			end=length-1;
			start=end-span+1;
			System.arraycopy(yvals,start,data,0,span);
		} else {
			System.arraycopy(yvals,start,data,0,span);
			return new float[][]{centerweights,data,new float[]{center-start}};
		}
		for(int i=start;i<=end;i++){
			int dist2=i-center;
			if(dist2<0) dist2=-dist2;
			float ratio=(float)dist2/(float)dist;
			float temp=1.0f-ratio*ratio*ratio;
			weights[i-start]=temp*temp*temp;
		}
		return new float[][]{weights,data,new float[]{center-start}};
	}
	
}
