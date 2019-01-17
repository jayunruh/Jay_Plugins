package jalgs.jfit;

import jalgs.jstatistics;

public class fit_kde_gaus{
	
	public static float get_kde_stdev(float[] data){
		//this uses the scott (1992) method
		//silverman used a leading factor of 0.9--narrower bandwidth
		//assumes the data is normal and univariate
		//will tend to over estimate the kernel bandwidth (thus oversmoothing the data)
		float[] quartiles={75.0f,25.0f};
		jstatistics.getstatistic("Percentile",data,quartiles);
		float iqr=(float)Math.abs(quartiles[0]-quartiles[1]);
		float stdev=jstatistics.getstatistic("stdev",data,null);
		float h=iqr/1.34f;
		return 1.06f*(float)Math.min(stdev,h)*(float)Math.pow((float)data.length,-0.2f);
	}
	
	public static float[] get_kde(float[] data,float bandwidth,float start,float dx,int npts){
		float[] dist=new float[npts];
		gausfunc gf=new gausfunc();
		for(int i=0;i<data.length;i++){
			float[] temp=gf.get_norm_func(start-data[i],npts,dx,bandwidth);
			for(int j=0;j<dist.length;j++){
				dist[j]+=temp[j];
			}
		}
		return dist;
	}
	
	public static float[] get_kde(float[] data,float start,float dx,int npts){
		return get_kde(data,get_kde_stdev(data),start,dx,npts);
	}
	
	public static float[][] get_kde(float[] data,float bandwidth){
		float min=jstatistics.getstatistic("min",data,null);
		float max=jstatistics.getstatistic("max",data,null);
		float minkde=min-2.0f*bandwidth;
		float maxkde=max+2.0f*bandwidth;
		float dx=0.5f*bandwidth;
		int npts=1+(int)((maxkde-minkde)/dx);
		float[] xvals=new float[npts];
		for(int i=0;i<xvals.length;i++) xvals[i]=minkde+dx*(float)i;
		return new float[][]{xvals,get_kde(data,bandwidth,minkde,dx,npts)};
	}
	
	public static float[][] get_kde(float[] data){
		return get_kde(data,get_kde_stdev(data));
	}

}
