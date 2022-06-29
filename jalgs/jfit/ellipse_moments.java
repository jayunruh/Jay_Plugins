package jalgs.jfit;


public class ellipse_moments{
	
	public static float[] get_ellipse_moments(float[] image,int width,int height){
		// this is a adapted from Bod Rodieck's code in ImageJ (EllipseFitter
		// class) that finds the ellipse
		// parameters that best match the spatial moments of the Roi
		double sum=0;
		double xsum=0;
		double ysum=0;
		double x2sum=0;
		double y2sum=0;
		double xysum=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				double val=image[j+i*width];
				sum+=val;
				xsum+=val*j;
				ysum+=val*i;
				x2sum+=val*j*j;
				y2sum+=val*i*i;
				xysum+=val*j*i;
			}
		}
		xsum/=sum;
		ysum/=sum;
		x2sum/=sum;
		y2sum/=sum;
		xysum/=sum;
		double xvar=x2sum-xsum*xsum;
		double yvar=y2sum-ysum*ysum;
		double xyvar=xysum-xsum*ysum;
		//double m4=4.0*Math.abs(yvar*xvar-xyvar*xyvar);
		//if(m4<0.000001)
		//	m4=0.000001;
		//a is the covariance matrix
		double m4=1.0f;
		double a11=yvar/m4;
		double a12=xyvar/m4;
		double a22=xvar/m4;
		double tmp=a11-a22;
		if(tmp==0.0)
			tmp=0.000001;
		double theta=0.5*Math.atan(2.0*a12/tmp);
		if(theta<0.0)
			theta+=0.5*Math.PI;
		if(a12>0.0)
			theta+=0.5*Math.PI;
		else if(a12==0.0){
			if(a22>a11){
				theta=0.0;
				tmp=a22;
				a22=a11;
				a11=tmp;
			}else if(a11!=a22)
				theta=0.5*Math.PI;
		}
		//tmp=Math.sin(theta);
		//if(tmp==0.0)
		//	tmp=0.000001;
		//double z=a12*Math.cos(theta)/tmp;
		//double major=Math.sqrt(1.0/Math.abs(a22+z));
		//double minor=Math.sqrt(1.0/Math.abs(a11-z));
		//double scale=Math.sqrt(sum/(Math.PI*major*minor)); // equalize
		// areas
		//major=major*scale*2.0;
		//minor=minor*scale*2.0;
		tmp=a11-a22;
		if(tmp==0.0)
			tmp=0.000001;
		double tempval=0.5f*Math.sqrt(4.0*a12*a12+tmp*tmp);
		double major=Math.sqrt(0.5*(a11+a22)+tempval);
		double minor=Math.sqrt(0.5*(a11+a22)-tempval);
		/*
		 * double angle = 180.0 * theta / Math.PI; if (angle == 180.0) angle =
		 * 0.0;
		 */
		double angle=theta;
		if(angle==Math.PI)
			angle=0.0;
		if(major<minor){
			tmp=major;
			major=minor;
			minor=tmp;
		}
		float[] output={(float)xsum,(float)ysum,(float)angle,(float)major,(float)minor};
		//float[] output={(float)xsum,(float)ysum,(float)xvar,(float)yvar,(float)xyvar};
		return output;
	}
	
	public static float[] get_spatial_moments(float[] image,int width,int height){
		// this is a adapted from Bod Rodieck's code in ImageJ (EllipseFitter
		// class) that finds the ellipse
		// parameters that best match the spatial moments of the Roi
		double sum=0;
		double xsum=0;
		double ysum=0;
		double x2sum=0;
		double y2sum=0;
		double xysum=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				double val=image[j+i*width];
				sum+=val;
				xsum+=val*j;
				ysum+=val*i;
				x2sum+=val*j*j;
				y2sum+=val*i*i;
				xysum+=val*j*i;
			}
		}
		xsum/=sum;
		ysum/=sum;
		x2sum/=sum;
		y2sum/=sum;
		xysum/=sum;
		double xvar=x2sum-xsum*xsum;
		double yvar=y2sum-ysum*ysum;
		double xyvar=xysum-xsum*ysum;
		float[] output={(float)xsum,(float)ysum,(float)xvar,(float)yvar,(float)xyvar};
		return output;
	}
	
	public static float[] get_ellipse_moments(byte[] mask,int width,int height){
		// this is a adapted from Bod Rodieck's code in ImageJ (EllipseFitter
		// class) that finds the ellipse
		// parameters that best match the spatial moments of the Roi
		int npixels=0;
		double xsum=0;
		double ysum=0;
		double x2sum=0;
		double y2sum=0;
		double xysum=0;
		for(int i=0;i<height;i++){
			for(int j=0;j<width;j++){
				int val=mask[j+i*width]&0xff;
				if(val>0){
					xsum+=j;
					ysum+=i;
					x2sum+=j*j;
					y2sum+=i*i;
					xysum+=j*i;
					npixels++;
				}
			}
		}
		/*
		 * x2sum+=0.08333333*npixels; y2sum+=0.08333333*npixels;
		 */
		xsum/=npixels;
		ysum/=npixels;
		x2sum/=npixels;
		y2sum/=npixels;
		xysum/=npixels;
		double xvar=x2sum-xsum*xsum;
		double yvar=y2sum-ysum*ysum;
		double xyvar=xysum-xsum*ysum;
		double m4=4.0*Math.abs(yvar*xvar-xyvar*xyvar);
		if(m4<0.000001)
			m4=0.000001;
		double a11=yvar/m4;
		double a12=xyvar/m4;
		double a22=xvar/m4;
		double tmp=a11-a22;
		if(tmp==0.0)
			tmp=0.000001;
		double theta=0.5*Math.atan(2.0*a12/tmp);
		if(theta<0.0)
			theta+=0.5*Math.PI;
		if(a12>0.0)
			theta+=0.5*Math.PI;
		else if(a12==0.0){
			if(a22>a11){
				theta=0.0;
				tmp=a22;
				a22=a11;
				a11=tmp;
			}else if(a11!=a22)
				theta=0.5*Math.PI;
		}
		tmp=Math.sin(theta);
		if(tmp==0.0)
			tmp=0.000001;
		double z=a12*Math.cos(theta)/tmp;
		double major=Math.sqrt(1.0/Math.abs(a22+z));
		double minor=Math.sqrt(1.0/Math.abs(a11-z));
		double scale=Math.sqrt(npixels/(Math.PI*major*minor)); // equalize
		// areas
		major=major*scale*2.0;
		minor=minor*scale*2.0;
		/*
		 * double angle = 180.0 * theta / Math.PI; if (angle == 180.0) angle =
		 * 0.0;
		 */
		double angle=theta;
		if(angle==Math.PI)
			angle=0.0;
		if(major<minor){
			tmp=major;
			major=minor;
			minor=tmp;
		}
		float[] output={(float)xsum,(float)ysum,(float)angle,(float)major,(float)minor};
		return output;
	}

}
