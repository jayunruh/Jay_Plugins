package jalgs.jfit;

import jalgs.*;

public class fit_mixed_gaus implements NLLSfitinterface_v3{
	public int xpts,ypts;
	public gui_interface gui;
	gausfunc gf;
	double[][] constraints;
	double oldtheta,costheta,sintheta;
	
	public fit_mixed_gaus(int xpts,int ypts,gui_interface gui){
		this.xpts=xpts;
		this.ypts=ypts;
		this.gui=gui;
		gf=new gausfunc();
	}

	public float[] fit_data(float[] data) {
		//this plugin takes an roi or small 2D image and fits it to a mixture of gaussians
		//rather than fitting to two independent gaussians, we fit to a pair of gaussians
		//with a center and an angle between them and an amplitude ratio
		//in that way, we can constrain the fit to prevent too small a distance and too large
		//an intensity ratio.
		//params are baseline,ampavg,ampratio,stdevavg,stdevratio,xc,yc,angle,distance
		String[] paramsnames={"baseline","ampavg","ampratio","stdevavg","stdevratio","xc","yc","angle","distance"};

		//make vague guesses as to the parameters
		//the stdev is typically 1 pixel for yeast spb with sim
		float stdevguess=1.0f;
		//first find the maximum point
		int maxx=0;
		int maxy=0;
		float maxval=data[0];
		float minval=data[0];
		for(int i=0;i<ypts;i++){
			for(int j=0;j<xpts;j++){
				if(data[j+i*xpts]>maxval){
					maxx=j;
					maxy=i;
					maxval=data[j+i*xpts];
				}
				if(data[j+i*xpts]<minval) minval=data[j+i*xpts];
			}
		}
		//now find the max point more than 2 stdevs away from the first (could be a red herring)
		int maxx2=0; int maxy2=0; float maxval2=data[0];
		for(int i=0;i<ypts;i++){
			for(int j=0;j<xpts;j++){
				if((float)((i-maxy)*(i-maxy)+(j-maxx)*(j-maxx))>(4.0f*stdevguess*stdevguess)){
					if(data[j+i*xpts]>maxval2){
						maxval2=data[j+i*xpts];
						maxx2=j;
						maxy2=i;
					}
				}
			}
		}
		//use the middle point as our guess center
		//use the point-to-point vector as our guess angle
		double dist=Math.sqrt((maxx2-maxx)*(maxx2-maxx)+(maxy2-maxy)*(maxy2-maxy));
		double angle=Math.acos((double)(maxx-maxx2)/dist);
		double[] gparams={(double)minval,(double)(0.5f*(maxval+maxval2)-minval),(double)(maxval2-minval)/(maxval-minval),(double)stdevguess,1.0,(double)(0.5f*(maxx+maxx2)),(double)(0.5f*(maxy+maxy2)),angle,dist};
		showresults("init params");
		for(int i=0;i<9;i++){
			showresults(paramsnames[i]+" , \t"+gparams[i]);
		}

		float max=jstatistics.getstatistic("Max",data,null);
		float min=jstatistics.getstatistic("Min",data,null);
		/*double stdev=2.0*(double)eparams[4]/2.35;
		double mindist=2.0*(double)eparams[3]-2.0*(double)eparams[4];
		double maxdist=2.0*(double)eparams[4];*/
		//params are 0baseline,1ampavg,2ampratio,3stdevavg,4stdevratio,5xc,6yc,7angle,8dist
		double[] params=gparams.clone();
		//double[] params={min,(double)(max-min),1.0,stdev,1.0,(double)eparams[0],(double)eparams[1],(double)eparams[2],0.5*(mindist+maxdist)};
		double[] paramscopy=params.clone();
		int[] fixes={0,0,0,0,1,1,1,1,0}; //set temporarily for troubleshooting
		constraints=new double[2][9];
		constraints[0][0]=-max; constraints[1][0]=max;
		constraints[0][1]=0.05*(max-min); constraints[1][1]=5.0*(max-min);
		constraints[0][2]=0.25; constraints[1][2]=4.0;
		constraints[0][3]=0.5; constraints[1][3]=2.0*gparams[3];
		constraints[0][4]=0.5; constraints[1][4]=2.0;
		constraints[0][5]=params[5]-3.0; constraints[1][5]=params[5]+3.0; //xc (was 3)
		constraints[0][6]=params[6]-3.0; constraints[1][6]=params[6]+3.0; //yc
		constraints[0][7]=params[7]-1.6; constraints[1][7]=params[7]+1.6; //angle
		constraints[0][8]=gparams[3]; constraints[1][8]=params[8]*4.0; //dist

		NLLSfit_v3 fitclass=new NLLSfit_v3(this,0.0001,10,0.0);
		double[] stats=new double[2];
		float[] fit=fitclass.fitdata(params,fixes,data,null,stats,true);

		//search for the best angle +/- 45 degrees
		/*double minangle=(double)eparams[2]-0.5*Math.PI;
		double maxangle=(double)minangle+Math.PI;
		double bestangle=(double)params[7];
		double minx=(double)eparams[0]-4.0;
		double maxx=(double)eparams[0]+4.0;
		double bestx=(double)params[5];
		double miny=(double)eparams[1]-4.0;
		double maxy=(double)eparams[1]+4.0;
		double besty=(double)params[6];
		double bestc2=-0.1;
		ImageStack tempstack=new ImageStack(xpts,ypts);
		searchlabel:
		for(double angle=minangle;angle<=maxangle;angle+=0.2){
			for(double xpos=minx;xpos<=maxx;xpos+=0.2){
				for(double ypos=miny;ypos<=maxy;ypos+=0.2){
					//double angle=findBestAngle(data,xpos,ypos,1.0,5.0,0.5,100);
					params=paramscopy.clone();
					params[7]=angle;
					params[6]=ypos;
					params[5]=xpos;
					float[] temp=fitclass.fitdata(params,fixes,data,null,stats,false);
					String label=""+angle+" , "+xpos+" , "+ypos;
					tempstack.addSlice(label,temp);
					if(bestc2<0.0 || stats[1]<bestc2){bestangle=params[7]; besty=params[6]; bestx=params[5]; bestc2=stats[1];}
					IJ.log(label+" , "+stats[1]);
					if(IJ.escapePressed()) break searchlabel;
				}
			}
		}
		String label2="best = "+bestangle+" , "+bestx+" , "+besty+" , "+bestc2;
		IJ.log(label2);
		new ImagePlus("angles",tempstack).show();
		params=paramscopy.clone();
		params[7]=bestangle;
		params[6]=besty;
		params[5]=bestx;
		fixes=new int[]{0,0,0,0,1,0,0,0,0};
		float[] fit=fitclass.fitdata(params,fixes,data,null,stats,false);*/

		return fit;
	}

	public double[] fitfunc(double[] params){
		//params are 0baseline,1ampavg,2ampratio,3stdevavg,4stdevratio,5xc,6yc,7angle,8dist
		//ampavg=0.5*(amp1+amp2)
		//ampratio=amp1/amp2
		if(params[7]!=oldtheta){
			oldtheta=params[7];
			costheta=Math.cos(oldtheta);
			sintheta=Math.sin(oldtheta);
		}
		double amp2=params[1]*2.0/(params[2]+1.0);
		double amp1=params[1]*2.0-amp2;
		double stdev2=params[3]*2.0/(params[4]+1.0);
		double stdev1=params[3]*2.0-stdev2;
		double xinc=0.5*params[8]*costheta;
		double yinc=-0.5*params[8]*sintheta;
		double x1=params[5]-xinc;
		double x2=params[5]+xinc;
		double y1=params[6]-yinc;
		double y2=params[6]+yinc;
		double[] temp=new double[xpts*ypts];
		gf.draw_2D_func(temp,x1,y1,xpts,ypts,stdev1,amp1);
		gf.draw_2D_func(temp,x2,y2,xpts,ypts,stdev2,amp2);
		for(int i=0;i<temp.length;i++) temp[i]+=params[0];
		return temp;
	}

	public void applyconstraints(double[] params,int[] fixes){
		//now implement the constraints
		//parameters are truncated at the constraint boundary
		for(int i=0;i<params.length;i++){
			if(fixes[i]==0){
				if(constraints!=null){
					if(params[i]<constraints[0][i]){
						params[i]=constraints[0][i];
					}
					if(params[i]>constraints[1][i]){
						params[i]=constraints[1][i];
					}
				}
			}
		}
	}

	public double[][] derivfunc(double[] params,int[] fixes,double[] fit){
		//params are 0baseline,1aa,2ar,3sa,4sr,5xc,6yc,7angle,8dist
		//this function should return the derivatives for each parameter (except the fixed ones)
		if(params[7]!=oldtheta){
			oldtheta=params[7];
			costheta=Math.cos(oldtheta);
			sintheta=Math.sin(oldtheta);
		}
		double xinc=0.5*params[8]*costheta;
		double yinc=-0.5*params[8]*sintheta;
		double x1=params[5]-xinc;
		double x2=params[5]+xinc;
		double y1=params[6]-yinc;
		double y2=params[6]+yinc;
		double stdev2=params[3]*2.0/(params[4]+1.0);
		double stdev1=params[3]*2.0-stdev2;
		double[] g1f=new double[xpts*ypts]; gf.draw_2D_func(g1f,x1,y1,xpts,ypts,stdev1,1.0f);
		double[] g2f=new double[xpts*ypts]; gf.draw_2D_func(g2f,x2,y2,xpts,ypts,stdev2,1.0f);
		double[][] ders=new double[xpts*ypts][params.length];
		int index=0;
		for(int y=0;y<ypts;y++){
			for(int x=0;x<xpts;x++){
				double tempx1=(double)x-x1;
				double tempx2=(double)x-x2;
				double tempy1=(double)y-y1;
				double tempy2=(double)y-y2;
				double temp1sr2=(1.0+params[4])*(1.0+params[4]);
				double a1=temp1sr2*(-tempx1*tempx1-tempy1*tempy1)/(8.0*params[3]*params[3]*params[4]*params[4]);
				double a2=temp1sr2*(-tempx2*tempx2-tempy2*tempy2)/(8.0*params[3]*params[3]);
				//double a1=Math.log(g1f[index]);
				//double a2=Math.log(g1f[index]);
				double g1=g1f[index];
				double g2=g2f[index];
				//double g1=Math.exp(a1);
				//double g2=Math.exp(a2);
				double temp5=params[2]*g1/(params[4]*params[4]);
				double temp6=params[1]*temp1sr2/(2.0*params[3]*params[3]*(1.0+params[2]));
				double[] der=new double[params.length];
				der[0]=1.0;
				//if(fixes[1]==0) der[1]=(fit[index]-params[0])/params[1];
				if(fixes[1]==0) der[1]=2.0*(params[2]*g1+g2)/(1.0+params[2]);
				if(fixes[2]==0) der[2]=(2.0*params[1]*g1-fit[index]+params[0])/(1.0+params[2]);
				if(fixes[3]==0) der[3]=-4.0*params[1]*(params[2]*g1*a1+g2*a2)/(params[3]*(1.0+params[2]));
				if(fixes[4]==0) der[4]=4.0*params[1]*(g2*a2+params[2]*g1*a1*(1.0-1.0/params[4]))/(1.0+params[2]);
				if(fixes[5]==0) der[5]=temp6*(temp5*tempx1+g2*tempx2);
				if(fixes[6]==0) der[6]=temp6*(temp5*tempy1+g2*tempy2);
				if(fixes[7]==0) der[7]=0.5*temp6*(temp5*(params[8]*sintheta*tempx1+params[8]*costheta*tempy1)-g2*(params[8]*sintheta*tempx2+params[8]*costheta*tempy2));
				if(fixes[8]==0) der[8]=0.5*temp6*(temp5*(sintheta*tempy1-costheta*tempx1)+g2*(costheta*tempx2-sintheta*tempy2));
				ders[index]=der;
				index++;
			}
		}
		//ders=NLLSfit_v3.derivfunc2(params,fixes,fit,this,0.0001);
		/*float[][] fders=new float[params.length][xpts*ypts];
		for(int i=0;i<params.length;i++){
			for(int j=0;j<xpts*ypts;j++){
				fders[i][j]=(float)ders[j][i];
			}
		}
		new ImagePlus("derivatives",jutils.array2stack(fders,xpts,ypts)).show();*/
		return ders;
	}

	public void showresults(String results){
		gui.showMessage(results);
	}

}

