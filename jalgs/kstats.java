/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class kstats{
	// this class contains methods for calculating k-statistics and variances of
	// k-statistics from raw moments
	public double[] conf={1.0,1.15,1.28,1.44,1.64,1.96,2.58};
	// these are the confidence intervals (in units of the stdev) that
	// correspond to 0.68, 0.75, 0.8, 0.85, 0.9, 0.95, and 0.99 probabilities
	public binning_function bf;
	public boolean binning;

	public kstats(boolean binning1){
		if(binning1){
			bf=new binning_function(5.0);
			binning=true;
		}
	}

	public kstats(){
		binning=false;
	}

	public double[] rawmoments(float[] data){
		double[] moms=new double[9];
		double dlength=data.length;
		for(int i=0;i<data.length;i++){
			double temp=data[i];
			for(int j=1;j<=8;j++){
				moms[j]+=intpow(temp,j)/dlength;
			}
		}
		return moms;
	}

	public double[] rawmomentsshort(float[] data){
		double[] moms=new double[5];
		double dlength=data.length;
		for(int i=0;i<data.length;i++){
			double temp=data[i];
			for(int j=1;j<=4;j++){
				moms[j]+=intpow(temp,j)/dlength;
			}
		}
		return moms;
	}

	public double[][] rawmoments(float[] data1,float[] data2){
		double[][] moms=new double[9][9];
		int length=data1.length;
		if(data2.length<data1.length){
			length=data2.length;
		}
		double dlength=length;

		double[] moms1=rawmoments(data1);
		double[] moms2=rawmoments(data2);
		for(int i=1;i<=8;i++){
			moms[i][0]=moms1[i];
			moms[0][i]=moms2[i];
		}
		for(int i=0;i<length;i++){
			double temp1=data1[i];
			double temp2=data2[i];
			for(int a=1;a<=6;a++){
				for(int b=1;b<=6;b++){
					moms[a][b]+=(intpow(temp1,a)*intpow(temp2,b))/dlength;
				}
			}
		}
		return moms;
	}

	public double[][] rawmomentsshort(float[] data1,float[] data2){
		double[][] moms=new double[5][5];
		int length=data1.length;
		if(data2.length<data1.length){
			length=data2.length;
		}
		double dlength=length;

		double[] moms1=rawmomentsshort(data1);
		double[] moms2=rawmomentsshort(data2);
		for(int i=1;i<=4;i++){
			moms[i][0]=moms1[i];
			moms[0][i]=moms2[i];
		}
		for(int i=0;i<length;i++){
			double temp1=data1[i];
			double temp2=data2[i];
			for(int a=1;a<=3;a++){
				for(int b=1;b<=3;b++){
					moms[a][b]+=(intpow(temp1,a)*intpow(temp2,b))/dlength;
				}
			}
		}
		return moms;
	}

	public double[][] kstatistics(float[] data1,float[] data2){
		double[][] u=rawmoments(data1,data2);
		if(data1.length<data2.length){
			return kstatistics(u,data1.length);
		}else{
			return kstatistics(u,data2.length);
		}
	}

	public double[] kstatistics(float[] data1){
		double[] u=rawmoments(data1);
		return kstatistics(u,data1.length);
	}

	public double[] kstatisticsshort(float[] data1){
		double[] u=rawmomentsshort(data1);
		return kstatisticsshort(u,data1.length);
	}

	public double[] kstatisticsshort(float[] data1,int length){
		if(length<data1.length){
			float[] temp=new float[length];
			System.arraycopy(data1,0,temp,0,length);
			return kstatisticsshort(temp);
		}else{
			return kstatisticsshort(data1);
		}
	}

	public double[][] kstatisticsshort(float[] data1,float[] data2){
		double[][] u=rawmomentsshort(data1,data2);
		if(data1.length<data2.length){
			return kstatisticsshort(u,data1.length);
		}else{
			return kstatisticsshort(u,data2.length);
		}
	}

	public double[][] kstatisticsshort(float[] data1,float[] data2,int length){
		if(length<data1.length&&length<data2.length){
			float[] temp1=new float[length];
			float[] temp2=new float[length];
			System.arraycopy(data1,0,temp1,0,length);
			System.arraycopy(data2,0,temp2,0,length);
			return kstatisticsshort(temp1,temp2);
		}else{
			return kstatisticsshort(data1,data2);
		}
	}

	public double[][] kstatistics(double[][] u,int n1){
		double n=n1;
		double[][] k=new double[9][9];
		// start with the univariate cases
		double[] u1=new double[9];
		double[] u2=new double[9];
		for(int i=0;i<9;i++){
			u1[i]=u[i][0];
			u2[i]=u[0][i];
		}
		double[] k1=kstatistics(u1,n1);
		double[] k2=kstatistics(u2,n1);

		for(int i=0;i<9;i++){
			k[i][0]=k1[i];
			k[0][i]=k2[i];
		}

		k[2][2]=-((n*n*(-2.0*u[1][1]*u[1][1]+2.0*n*u[1][1]*u[1][1]+2.0*u[1][0]*u[1][2]+2.0*n*u[1][0]*u[1][2]-u[0][2]*(2.0*n*u[1][0]*u[1][0]-(-1.0+n)*u[2][0])+u[0][1]*u[0][1]
				*(6.0*n*u[1][0]*u[1][0]-2.0*n*u[2][0])+2.0*u[0][1]*(-4.0*n*u[1][0]*u[1][1]+(1.0+n)*u[2][1])-u[2][2]-n*u[2][2]))/(-6.0+11.0*n-6.0*n*n+n*n*n));

		k[3][1]=(-6.0*u[0][1]*u[1][0]*u[1][0]*u[1][0]-3.0*(-1.0+n)*n*u[1][1]*u[2][0]+2.0*n*(3.0*u[1][0]*u[1][0]*u[1][1]+3.0*u[0][1]*u[1][0]*u[2][0])-n*(1.0+n)*(3.0*u[1][0]*u[2][1]+u[0][1]*u[3][0])+n
				*n*(1.0+n)*u[3][1])
				/((-3.0+n)*(-2.0+n)*(-1.0+n)*n);

		k[1][3]=(-6.0*u[1][0]*u[0][1]*u[0][1]*u[0][1]-3.0*(-1.0+n)*n*u[1][1]*u[0][2]+2.0*n*(3.0*u[0][1]*u[0][1]*u[1][1]+3.0*u[1][0]*u[0][1]*u[0][2])-n*(1.0+n)*(3.0*u[0][1]*u[1][2]+u[1][0]*u[0][3])+n
				*n*(1.0+n)*u[1][3])
				/((-3.0+n)*(-2.0+n)*(-1.0+n)*n);

		k[2][1]=(n*n*2.0*u[0][1]*u[1][0]*u[1][0]-n*n*(2.0*u[1][0]*u[1][1]+u[0][1]*u[2][0])+n*n*u[2][1])/((-2.0+n)*(-1.0+n));

		k[1][2]=(n*n*2.0*u[1][0]*u[0][1]*u[0][1]-n*n*(2.0*u[0][1]*u[1][1]+u[1][0]*u[0][2])+n*n*u[1][2])/((-2.0+n)*(-1.0+n));

		k[1][1]=(-n*u[0][1]*u[1][0]+n*u[1][1])/(-1.0+n);

		return k;
	}

	public double[][] kstatisticsshort(double[][] u,int n1){
		double n=n1;
		double[][] k=new double[5][5];
		// start with the univariate cases
		double[] u1=new double[5];
		double[] u2=new double[5];
		for(int i=0;i<5;i++){
			u1[i]=u[i][0];
			u2[i]=u[0][i];
		}
		double[] k1=kstatisticsshort(u1,n1);
		double[] k2=kstatisticsshort(u2,n1);

		for(int i=0;i<5;i++){
			k[i][0]=k1[i];
			k[0][i]=k2[i];
		}

		k[2][1]=(n*n*2.0*u[0][1]*u[1][0]*u[1][0]-n*n*(2.0*u[1][0]*u[1][1]+u[0][1]*u[2][0])+n*n*u[2][1])/((-2.0+n)*(-1.0+n));

		k[1][2]=(n*n*2.0*u[1][0]*u[0][1]*u[0][1]-n*n*(2.0*u[0][1]*u[1][1]+u[1][0]*u[0][2])+n*n*u[1][2])/((-2.0+n)*(-1.0+n));

		k[1][1]=(-n*u[0][1]*u[1][0]+n*u[1][1])/(-1.0+n);

		return k;
	}

	public double[] kstatistics(double[] u,int nn){
		double[] k=new double[9];
		double n=nn;

		k[8]=-((n*n*(n
				*n
				*(-3780.0*u[2]*u[2]*u[2]*u[2]+840.0*u[2]*u[2]*(12.0*u[1]*u[3]-23.0*u[4])-4235.0*u[4]*u[4]-6048.0*u[1]*u[1]*u[1]*u[5]+6216.0*u[3]*u[5]+28.0*u[2]
						*(400.0*u[3]*u[3]+540.0*u[1]*u[1]*u[4]+276.0*u[1]*u[5]-157.0*u[6])-224.0*u[1]*u[1]*(75.0*u[3]*u[3]-2.0*u[6])+8.0*u[1]*(350.0*u[3]*u[4]+141.0*u[7])-141.0*u[8])
				+n
				*n
				*n
				*n
				*(-3780.0*u[2]*u[2]*u[2]*u[2]+33600.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[3]-840.0*u[2]*u[2]*u[4]+35.0*u[4]*u[4]-1680.0*u[1]*u[1]*u[1]*u[1]*(15.0*u[2]*u[2]+17.0*u[4])-672.0*u[1]*u[1]*u[1]
						*(40.0*u[2]*u[3]-19.0*u[5])+504.0*u[3]*u[5]+112.0*u[1]*u[1]*(270.0*u[2]*u[2]*u[2]+30.0*u[3]*u[3]+225.0*u[2]*u[4]-34.0*u[6])+28.0*u[2]*(80.0*u[3]*u[3]+37.0*u[6])-8.0*u[1]
						*(1260.0*u[2]*u[2]*u[3]+350.0*u[3]*u[4]+966.0*u[2]*u[5]-99.0*u[7])-99.0*u[8])
				+n
				*n
				*n
				*n
				*n
				*(5040.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]-20160.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]*u[2]+630.0*u[2]*u[2]*u[2]*u[2]+6720.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[3]+1680.0*u[1]*u[1]*u[1]*u[1]
						*(15.0*u[2]*u[2]-u[4])-420.0*u[2]*u[2]*u[4]+35.0*u[4]*u[4]-336.0*u[1]*u[1]*u[1]*(40.0*u[2]*u[3]-u[5])+56.0*u[3]*u[5]+28.0*u[2]*(-20.0*u[3]*u[3]+u[6])-56.0*u[1]*u[1]
						*(180.0*u[2]*u[2]*u[2]-30.0*u[3]*u[3]-45.0*u[2]*u[4]+u[6])+8.0*u[1]*(630.0*u[2]*u[2]*u[3]-70.0*u[3]*u[4]-42.0*u[2]*u[5]+u[7])-u[8])
				-120.0
				*(35.0*u[4]*u[4]-56.0*u[3]*u[5]+28.0*u[2]*u[6]-8.0*u[1]*u[7]+u[8])
				+2.0
				*n
				*(5040.0*u[2]*u[2]*u[4]+3605.0*u[4]*u[4]-6104.0*u[3]*u[5]-28.0*u[2]*(120.0*u[3]*u[3]+180.0*u[1]*u[5]-139.0*u[6])+1680.0*u[1]*u[1]*u[6]+8.0*u[1]*(210.0*u[3]*u[4]-199.0*u[7])+199.0*u[8])-n
				*n
				*n
				*(-6930.0*u[2]*u[2]*u[2]*u[2]+20160.0*u[1]*u[1]*u[1]*u[1]*u[4]-10500.0*u[2]*u[2]*u[4]-1155.0*u[4]*u[4]+1288.0*u[3]*u[5]-1008.0*u[1]*u[1]*u[1]*(40.0*u[2]*u[3]+33.0*u[5])+28.0*u[2]
						*(220.0*u[3]*u[3]+39.0*u[6])+56.0*u[1]*u[1]*(360.0*u[2]*u[2]*u[2]-210.0*u[3]*u[3]+765.0*u[2]*u[4]+359.0*u[6])+8.0*u[1]
						*(630.0*u[2]*u[2]*u[3]+350.0*u[3]*u[4]-1302.0*u[2]*u[5]-757.0*u[7])+757.0*u[8])))/(-5040.0+13068.0*n-13132.0*n*n+6769.0*n*n*n-1960.0*n*n*n*n+322.0*n*n*n*n*n-28.0*n*n*n*n*n*n+n
				*n*n*n*n*n*n));

		k[7]=(n*n*n*(-7.0
				*n
				*(-60.0*u[2]*u[2]*u[3]+180.0*u[1]*u[1]*u[1]*u[4]+25.0*u[3]*u[4]-93.0*u[2]*u[5]-24.0*u[1]*u[1]*(15.0*u[2]*u[3]+11.0*u[5])+u[1]
						*(180.0*u[2]*u[2]*u[2]-100.0*u[3]*u[3]+390.0*u[2]*u[4]+119.0*u[6])-17.0*u[7])
				-42.0
				*(-5.0*u[3]*u[4]+12.0*u[1]*u[1]*u[5]+9.0*u[2]*u[5]+u[1]*(20.0*u[3]*u[3]-30.0*u[2]*u[4]-7.0*u[6])+u[7])
				+42.0
				*n
				*n
				*(80.0*u[1]*u[1]*u[1]*u[1]*u[3]-15.0*u[2]*u[2]*u[3]-5.0*u[1]*u[1]*u[1]*(12.0*u[2]*u[2]+13.0*u[4])-6.0*u[2]*u[5]+u[1]*u[1]*(-30.0*u[2]*u[3]+27.0*u[5])+u[1]
						*(45.0*u[2]*u[2]*u[2]+30.0*u[2]*u[4]-7.0*u[6])+u[7])+n
				*n
				*n
				*(720.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]-2520.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[2]+840.0*u[1]*u[1]*u[1]*u[1]*u[3]+210.0*u[2]*u[2]*u[3]+210.0*u[1]*u[1]*u[1]*(12.0*u[2]*u[2]-u[4])-35.0
						*u[3]*u[4]-42.0*u[1]*u[1]*(30.0*u[2]*u[3]-u[5])-21.0*u[2]*u[5]-7.0*u[1]*(90.0*u[2]*u[2]*u[2]-20.0*u[3]*u[3]-30.0*u[2]*u[4]+u[6])+u[7])))
				/(720.0-1764.0*n+1624.0*n*n-735.0*n*n*n+175.0*n*n*n*n-21.0*n*n*n*n*n+n*n*n*n*n*n);

		k[6]=-((n*n*(n
				*(-60.0*u[2]*u[2]*u[2]+50.0*u[3]*u[3]+15.0*u[2]*(8.0*u[1]*u[3]-7.0*u[4])-60.0*u[1]*u[1]*u[4]+66.0*u[1]*u[5]-11.0*u[6])
				+2.0
				*n
				*n
				*(45.0*u[2]*u[2]*u[2]+180.0*u[1]*u[1]*u[1]*u[3]-10.0*u[3]*u[3]+15.0*u[2]*u[4]-135.0*u[1]*u[1]*(u[2]*u[2]+u[4])+48.0*u[1]*u[5]-8.0*u[6])
				+n
				*n
				*n
				*(120.0*u[1]*u[1]*u[1]*u[1]*u[1]*u[1]-360.0*u[1]*u[1]*u[1]*u[1]*u[2]-30.0*u[2]*u[2]*u[2]+120.0*u[1]*u[1]*u[1]*u[3]+10.0*u[3]*u[3]+30.0*u[1]*u[1]*(9.0*u[2]*u[2]-u[4])+15.0*u[2]*u[4]
						+6.0*u[1]*(-20.0*u[2]*u[3]+u[5])-u[6])+4.0*(-10.0*u[3]*u[3]+15.0*u[2]*u[4]-6.0*u[1]*u[5]+u[6])))/(-120.0+274.0*n-225.0*n*n+85.0*n*n*n-15.0*n*n*n*n+n*n*n*n*n));

		k[5]=(n*n*n*(n*(24.0*u[1]*u[1]*u[1]*u[1]*u[1]-60.0*u[1]*u[1]*u[1]*u[2]+20.0*u[1]*u[1]*u[3]-10.0*u[2]*u[3]+5.0*u[1]*(6.0*u[2]*u[2]-u[4])+u[5])+5.0*(8.0*u[1]*u[1]*u[3]+2.0*u[2]*u[3]-u[1]
				*(6.0*u[2]*u[2]+5.0*u[4])+u[5])))
				/(24.0-50.0*n+35.0*n*n-10.0*n*n*n+n*n*n*n);

		k[4]=-(-3.0*u[2]*u[2]+4.0*u[1]*u[3]+n*(6.0*u[1]*u[1]*u[1]*u[1]-12.0*u[1]*u[1]*u[2]+3.0*u[2]*u[2]+4.0*u[1]*u[3]-u[4])-u[4])*((n*n)/(-6.0+11.0*n-6.0*n*n+n*n*n));

		k[3]=(2.0*u[1]*u[1]*u[1]-3.0*u[1]*u[2]+u[3])*((n*n)/(2.0-3.0*n+n*n));

		k[2]=(n*(-u[1]*u[1]+u[2]))/(-1.0+n);

		k[1]=u[1];

		return k;
	}

	public double[] kstatisticsshort(double[] u,int nn){
		double[] k=new double[5];
		double n=nn;

		k[4]=-(-3.0*u[2]*u[2]+4.0*u[1]*u[3]+n*(6.0*u[1]*u[1]*u[1]*u[1]-12.0*u[1]*u[1]*u[2]+3.0*u[2]*u[2]+4.0*u[1]*u[3]-u[4])-u[4])*((n*n)/(-6.0+11.0*n-6.0*n*n+n*n*n));

		k[3]=(2.0*u[1]*u[1]*u[1]-3.0*u[1]*u[2]+u[3])*((n*n)/(2.0-3.0*n+n*n));

		k[2]=(n*(-u[1]*u[1]+u[2]))/(-1.0+n);

		k[1]=u[1];

		return k;
	}

	public double[] kstatsvars(double[] K,int nn){
		double[] kvars=new double[5];
		double n=nn;

		kvars[1]=K[2]/n;

		kvars[2]=(2.0*K[2]*K[2])/(-1.0+n)+K[4]/n;

		kvars[3]=(6.0*n*K[2]*K[2]*K[2])/(2.0-3.0*n+n*n)+(9.0*(K[3]*K[3]+K[2]*K[4]))/(-1.0+n)+K[6]/n;

		kvars[4]=(24.0*n*n*K[2]*K[2]*K[2]*K[2])/(-6.0+11.0*n-6.0*n*n+n*n*n)+(24.0*n*K[2]*(K[2]*K[2]*K[2]+6.0*(-3.0+n)*K[3]*K[3]+3.0*(-3.0+n)*K[2]*K[4]))/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(2.0*(17.0*K[4]*K[4]+24.0*K[3]*K[5]+8.0*K[2]*K[6]))/(-1.0+n)+K[8]/n;

		return kvars;
	}

	public double[][] kstatsvars(double[][] K,int n){
		double[][] kvars=new double[5][5];

		kvars[1][0]=K[2][0]/n;

		kvars[0][1]=K[2][0]/n;

		kvars[2][0]=(2.0*K[2][0]*K[2][0])/(-1.0+n)+K[4][0]/n;

		kvars[0][2]=(2.0*K[0][2]*K[0][2])/(-1.0+n)+K[0][4]/n;

		kvars[3][0]=(6.0*n*K[2][0]*K[2][0]*K[2][0])/(2.0-3.0*n+n*n)+(9.0*(K[3][0]*K[3][0]+K[2][0]*K[4][0]))/(-1.0+n)+K[6][0]/n;

		kvars[0][3]=(6.0*n*K[0][2]*K[0][2]*K[0][2])/(2.0-3.0*n+n*n)+(9.0*(K[0][3]*K[0][3]+K[0][2]*K[0][4]))/(-1.0+n)+K[0][6]/n;

		kvars[4][0]=(24.0*n*n*K[2][0]*K[2][0]*K[2][0]*K[2][0])/(-6.0+11.0*n-6.0*n*n+n*n*n)+(24.0*n*K[2][0]*(K[2][0]*K[2][0]*K[2][0]+6.0*(-3.0+n)*K[3][0]*K[3][0]+3.0*(-3.0+n)*K[2][0]*K[4][0]))
				/(-6.0+11.0*n-6.0*n*n+n*n*n)+(2.0*(17.0*K[4][0]*K[4][0]+24.0*K[3][0]*K[5][0]+8.0*K[2][0]*K[6][0]))/(-1.0+n)+K[8][0]/n;

		kvars[0][4]=(24.0*n*n*K[0][2]*K[0][2]*K[0][2]*K[0][2])/(-6.0+11.0*n-6.0*n*n+n*n*n)+(24.0*n*K[0][2]*(K[0][2]*K[0][2]*K[0][2]+6.0*(-3.0+n)*K[0][3]*K[0][3]+3.0*(-3.0+n)*K[0][2]*K[0][4]))
				/(-6.0+11.0*n-6.0*n*n+n*n*n)+(2.0*(17.0*K[0][4]*K[0][4]+24.0*K[0][3]*K[0][5]+8.0*K[0][2]*K[0][6]))/(-1.0+n)+K[0][8]/n;

		kvars[1][1]=K[1][1]*K[1][1]/(-1.0+n)+(K[0][2]*K[2][0])/(-1.0+n)+K[2][2]/n;

		kvars[2][1]=(2.0*n*K[2][0]*(2.0*K[1][1]*K[1][1]+K[0][2]*K[2][0]))/(2.0-3.0*n+n*n)+(5.0*K[2][1]*K[2][1]+4.0*K[2][0]*K[2][2]+4.0*K[1][2]*K[3][0]+4.0*K[1][1]*K[3][1]+K[0][2]*K[4][0])/(-1.0+n)
				+K[4][2]/n;

		kvars[1][2]=(2.0*n*K[0][2]*(2.0*K[1][1]*K[1][1]+K[0][2]*K[2][0]))/(2.0-3.0*n+n*n)+(5.0*K[1][2]*K[1][2]+4.0*K[1][1]*K[1][3]+K[0][4]*K[2][0]+4.0*K[0][3]*K[2][1]+4.0*K[0][2]*K[2][2])/(-1.0+n)
				+K[2][4]/n;

		kvars[3][1]=(6.0*n*n*K[2][0]*K[2][0]*(3.0*K[1][1]*K[1][1]+K[0][2]*K[2][0]))
				/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(3.0*n*(3.0*(-3.0+n)*K[2][0]*(5.0*K[2][1]*K[2][1]+2.0*K[2][0]*K[2][2]+4.0*K[1][2]*K[3][0])+6.0*(-3.0+n)*K[1][1]*(3.0*K[2][1]*K[3][0]+2.0*K[2][0]*K[3][1])+3.0*K[1][1]*K[1][1]
						*(2.0*K[2][0]*K[2][0]+(-3.0+n)*K[4][0])+K[0][2]*(2.0*K[2][0]*K[2][0]*K[2][0]+3.0*(-3.0+n)*K[3][0]*K[3][0]+3.0*(-3.0+n)*K[2][0]*K[4][0])))/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(19.0*K[3][1]*K[3][1]+18.0*K[3][0]*K[3][2]+15.0*K[2][2]*K[4][0]+24.0*K[2][1]*K[4][1]+9.0*K[2][0]*K[4][2]+6.0*K[1][2]*K[5][0]+6.0*K[1][1]*K[5][1]+K[0][2]*K[6][0])/(-1.0+n)+K[6][2]/n;

		kvars[1][3]=(6.0*n*n*K[0][2]*K[0][2]*(3.0*K[1][1]*K[1][1]+K[0][2]*K[2][0]))
				/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(3.0*n*(2.0*K[0][2]*K[0][2]*K[0][2]*K[2][0]+3.0*(-3.0+n)*(K[0][4]*K[1][1]*K[1][1]+K[0][3]*(6.0*K[1][1]*K[1][2]+K[0][3]*K[2][0]))+3.0*(-3.0+n)*K[0][2]
						*(5.0*K[1][2]*K[1][2]+4.0*K[1][1]*K[1][3]+K[0][4]*K[2][0]+4.0*K[0][3]*K[2][1])+6.0*K[0][2]*K[0][2]*(K[1][1]*K[1][1]+(-3.0+n)*K[2][2])))/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(19.0*K[1][3]*K[1][3]+24.0*K[1][2]*K[1][4]+6.0*K[1][1]*K[1][5]+K[0][6]*K[2][0]+6.0*K[0][5]*K[2][1]+15.0*K[0][4]*K[2][2]+18.0*K[0][3]*K[2][3]+9.0*K[0][2]*K[2][4])/(-1.0+n)+K[2][6]/n;

		kvars[2][2]=(4.0*n*n*(K[1][1]*K[1][1]*K[1][1]*K[1][1]+4.0*K[0][2]*K[1][1]*K[1][1]*K[2][0]+K[0][2]*K[0][2]*K[2][0]*K[2][0]))
				/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(2.0*n*(2.0*K[1][1]*K[1][1]*K[1][1]*K[1][1]+10.0*(-3.0+n)*K[1][2]*K[1][2]*K[2][0]+2.0*K[0][2]*K[0][2]*K[2][0]*K[2][0]-3.0*K[0][4]*K[2][0]*K[2][0]+n*K[0][4]*K[2][0]*K[2][0]-24.0
						*K[0][3]*K[2][0]*K[2][1]+8.0*n*K[0][3]*K[2][0]*K[2][1]-30.0*K[0][2]*K[2][1]*K[2][1]+10.0*n*K[0][2]*K[2][1]*K[2][1]-24.0*K[0][2]*K[2][0]*K[2][2]+8.0*n*K[0][2]*K[2][0]*K[2][2]
						+2.0*K[1][1]*K[1][1]*(4.0*K[0][2]*K[2][0]+5.0*(-3.0+n)*K[2][2])+8.0*(-3.0+n)*K[0][2]*K[1][2]*K[3][0]+4.0*(-3.0+n)*K[1][1]
						*(2.0*K[1][3]*K[2][0]+8.0*K[1][2]*K[2][1]+K[0][3]*K[3][0]+2.0*K[0][2]*K[3][1])-3.0*K[0][2]*K[0][2]*K[4][0]+n*K[0][2]*K[0][2]*K[4][0]))
				/(-6.0+11.0*n-6.0*n*n+n*n*n)
				+(17.0*K[2][2]*K[2][2]+20.0*K[2][1]*K[2][3]+4.0*K[2][0]*K[2][4]+4.0*K[1][4]*K[3][0]+16.0*K[1][3]*K[3][1]+20.0*K[1][2]*K[3][2]+8.0*K[1][1]*K[3][3]+K[0][4]*K[4][0]+4.0*K[0][3]*K[4][1]+4.0
						*K[0][2]*K[4][2])/(-1.0+n)+K[4][4]/n;

		return kvars;
	}

	public double[] cumslopes(double tau1,double amp1,double tau2){
		double[] rawmoments=new double[9];
		for(int i=1;i<=8;i++){
			rawmoments[i]=amp1*intpow(tau1,i)+(1.0-amp1)*intpow(tau2,i);
			rawmoments[i]*=factorial(i);
		}
		double avgint=rawmoments[1];
		for(int i=1;i<8;i++){
			rawmoments[i]/=avgint;
		}
		return rawmoments;
	}

	public double[] cumslopes(double tau){
		double[] rawmoments=new double[9];
		for(int i=1;i<=8;i++){
			rawmoments[i]=intpow(tau,i);
			rawmoments[i]*=factorial(i);
		}
		double avgint=rawmoments[1];
		for(int i=1;i<8;i++){
			rawmoments[i]/=avgint;
		}
		return rawmoments;
	}

	public double[] cumslopes(double[] spcentmoms){
		double[] slopes=new double[9];
		slopes[1]=1.0f;
		slopes[2]=spcentmoms[1]+spcentmoms[2]/spcentmoms[1];
		slopes[3]=spcentmoms[1]*spcentmoms[1]+3.0*spcentmoms[2]+spcentmoms[3]/spcentmoms[1];
		slopes[4]=intpow(spcentmoms[1],3)+6.0*spcentmoms[1]*spcentmoms[2]+4.0*spcentmoms[3]+spcentmoms[4]/spcentmoms[1];
		slopes[5]=intpow(spcentmoms[1],4)+10.0*intpow(spcentmoms[1],2)*spcentmoms[2]+10.0*spcentmoms[1]*spcentmoms[3]+5.0*spcentmoms[4]+spcentmoms[5]/spcentmoms[1];
		slopes[6]=intpow(spcentmoms[1],5)+15.0*intpow(spcentmoms[1],3)*spcentmoms[2]+20.0*intpow(spcentmoms[1],2)*spcentmoms[3]+15.0*spcentmoms[1]*spcentmoms[4]+6.0*spcentmoms[5]+spcentmoms[6]
				/spcentmoms[1];
		slopes[7]=intpow(spcentmoms[1],6)+21.0*intpow(spcentmoms[1],4)*spcentmoms[2]+35.0*intpow(spcentmoms[1],3)*spcentmoms[3]+35.0*intpow(spcentmoms[1],2)*spcentmoms[4]+21.0*spcentmoms[1]
				*spcentmoms[5]+7.0*spcentmoms[6]+spcentmoms[7]/spcentmoms[1];
		slopes[8]=intpow(spcentmoms[1],7)+28.0*intpow(spcentmoms[1],5)*spcentmoms[2]+56.0*intpow(spcentmoms[1],4)*spcentmoms[3]+70.0*intpow(spcentmoms[1],3)*spcentmoms[4]+56.0*intpow(spcentmoms[1],2)
				*spcentmoms[5]+28.0*spcentmoms[1]*spcentmoms[6]+8.0*spcentmoms[7]+spcentmoms[8]/spcentmoms[1];
		return slopes;
	}

	public double[] spcentmoms(double tau1,double amp1,double tau2){
		double[] rawmoments=new double[9];
		for(int i=1;i<=8;i++){
			rawmoments[i]=amp1*intpow(tau1,i)+(1.0-amp1)*intpow(tau2,i);
			rawmoments[i]*=factorial(i);
		}
		double[] centmoms=rawmoms2centmoms(rawmoments);
		return centmoms;
	}

	public double[] analog_cumulants(float[] brightnesses,float[] numbers,double[] gammas,double[] slopes,double[] offsets,double gain){
		double[] afcums=new double[9];
		double[] acums=new double[9];
		// first calculate the analog factorial cumulants
		for(int i=1;i<=4;i++){
			for(int j=0;j<brightnesses.length;j++){
				afcums[i]+=intpow(brightnesses[j],i)*numbers[j];
			}
		}
		for(int i=1;i<=4;i++){
			afcums[i]*=gammas[i]*intpow(gain,i);
		}
		// now calculate the analog cumulants from these
		acums[1]=offsets[1]+afcums[1];
		acums[2]=offsets[2]+afcums[2]+afcums[1]*slopes[2];
		acums[3]=offsets[3]+afcums[3]+3.0*slopes[2]*afcums[2]+afcums[1]*slopes[3];
		acums[4]=offsets[4]+afcums[4]+6.0*afcums[3]*slopes[2]+3.0*afcums[2]*intpow(slopes[2],2)+4.0*afcums[2]*slopes[3]+afcums[1]*slopes[4];
		acums[5]=offsets[5]+afcums[5]+10.0*afcums[4]*slopes[2]+15.0*afcums[3]*intpow(slopes[2],2)+10.0*afcums[3]*slopes[3]+10.0*afcums[2]*slopes[2]*slopes[3]+5.0*afcums[2]*slopes[4]+afcums[1]
				*slopes[5];
		acums[6]=offsets[6]+afcums[6]+15.0*afcums[5]*slopes[2]+45.0*afcums[4]*intpow(slopes[2],2)+15.0*afcums[3]*intpow(slopes[2],3)+20.0*afcums[4]*slopes[3]+60.0*afcums[3]*slopes[2]*slopes[3]+10.0
				*afcums[2]*intpow(slopes[3],2)+15.0*afcums[3]*slopes[4]+15.0*afcums[2]*slopes[2]*slopes[4]+6.0*afcums[2]*slopes[5]+afcums[1]*slopes[6];
		acums[7]=offsets[7]+afcums[7]+21.0*afcums[6]*slopes[2]+105.0*afcums[5]*intpow(slopes[2],2)+105.0*afcums[4]*intpow(slopes[2],3)+35.0*afcums[5]*slopes[3]+210.0*afcums[4]*slopes[2]*slopes[3]
				+105.0*afcums[3]*intpow(slopes[2],2)*slopes[3]+70.0*afcums[3]*intpow(slopes[3],2)+35.0*afcums[4]*slopes[4]+105.0*afcums[3]*slopes[2]*slopes[4]+35.0*afcums[2]*slopes[3]*slopes[4]+21.0
				*afcums[3]*slopes[5]+21.0*afcums[2]*slopes[2]*slopes[5]+7.0*afcums[2]*slopes[6]+afcums[1]*slopes[7];
		acums[8]=offsets[8]+afcums[8]+28.0*afcums[7]*slopes[2]+210.0*afcums[6]*intpow(slopes[2],2)+420.0*afcums[5]*intpow(slopes[2],3)+105.0*afcums[4]*intpow(slopes[2],4)+56.0*afcums[6]*slopes[3]
				+560.0*afcums[5]*slopes[2]*slopes[3]+840.0*afcums[4]*intpow(slopes[2],2)*slopes[3]+280.0*afcums[4]*intpow(slopes[3],2)+280.0*afcums[3]*slopes[2]*intpow(slopes[3],2)
				+70.0*afcums[5]*slopes[4]+420.0*afcums[4]*slopes[2]*slopes[4]+210.0*afcums[3]*intpow(slopes[2],2)*slopes[4]+280.0*afcums[3]*slopes[3]*slopes[4]+35.0*afcums[2]*intpow(slopes[4],2)
				+56.0*afcums[4]*slopes[5]+168.0*afcums[3]*slopes[2]*slopes[5]+56.0*afcums[2]*slopes[3]*slopes[5]+28.0*afcums[3]*slopes[6]+28.0*afcums[2]*slopes[2]*slopes[6]+8.0*afcums[2]*slopes[7]
				+afcums[1]*slopes[8];
		return acums;
	}

	public double[] cums2afcums(double[] k,double[] slopes,double[] offsets){
		double[] afcums=new double[9];
		// first calculate the analog factorial cumulants
		// now calculate the analog cumulants from these
		afcums[1]=k[1]-offsets[1];
		afcums[2]=k[2]-offsets[2]-afcums[1]*slopes[2];
		afcums[3]=k[3]-offsets[3]-3.0*slopes[2]*afcums[2]-afcums[1]*slopes[3];
		afcums[4]=k[4]-offsets[4]-6.0*afcums[3]*slopes[2]-3.0*afcums[2]*intpow(slopes[2],2)-4.0*afcums[2]*slopes[3]-afcums[1]*slopes[4];
		afcums[5]=k[5]-offsets[5]-10.0*afcums[4]*slopes[2]-15.0*afcums[3]*intpow(slopes[2],2)-10.0*afcums[3]*slopes[3]-10.0*afcums[2]*slopes[2]*slopes[3]-5.0*afcums[2]*slopes[4]-afcums[1]*slopes[5];
		afcums[6]=k[6]-offsets[6]-15.0*afcums[5]*slopes[2]-45.0*afcums[4]*intpow(slopes[2],2)-15.0*afcums[3]*intpow(slopes[2],3)-20.0*afcums[4]*slopes[3]-60.0*afcums[3]*slopes[2]*slopes[3]-10.0
				*afcums[2]*intpow(slopes[3],2)-15.0*afcums[3]*slopes[4]-15.0*afcums[2]*slopes[2]*slopes[4]-6.0*afcums[2]*slopes[5]-afcums[1]*slopes[6];
		afcums[7]=k[7]-offsets[7]-21.0*afcums[6]*slopes[2]-105.0*afcums[5]*intpow(slopes[2],2)-105.0*afcums[4]*intpow(slopes[2],3)-35.0*afcums[5]*slopes[3]-210.0*afcums[4]*slopes[2]*slopes[3]-105.0
				*afcums[3]*intpow(slopes[2],2)*slopes[3]-70.0*afcums[3]*intpow(slopes[3],2)-35.0*afcums[4]*slopes[4]-105.0*afcums[3]*slopes[2]*slopes[4]-35.0*afcums[2]*slopes[3]*slopes[4]-21.0
				*afcums[3]*slopes[5]-21.0*afcums[2]*slopes[2]*slopes[5]-7.0*afcums[2]*slopes[6]-afcums[1]*slopes[7];
		afcums[8]=k[8]-offsets[8]-28.0*afcums[7]*slopes[2]-210.0*afcums[6]*intpow(slopes[2],2)-420.0*afcums[5]*intpow(slopes[2],3)-105.0*afcums[4]*intpow(slopes[2],4)-56.0*afcums[6]*slopes[3]-560.0
				*afcums[5]*slopes[2]*slopes[3]-840.0*afcums[4]*intpow(slopes[2],2)*slopes[3]-280.0*afcums[4]*intpow(slopes[3],2)-280.0*afcums[3]*slopes[2]*intpow(slopes[3],2)-70.0*afcums[5]*slopes[4]
				-420.0*afcums[4]*slopes[2]*slopes[4]-210.0*afcums[3]*intpow(slopes[2],2)*slopes[4]-280.0*afcums[3]*slopes[3]*slopes[4]-35.0*afcums[2]*intpow(slopes[4],2)-56.0*afcums[4]*slopes[5]
				-168.0*afcums[3]*slopes[2]*slopes[5]-56.0*afcums[2]*slopes[3]*slopes[5]-28.0*afcums[3]*slopes[6]-28.0*afcums[2]*slopes[2]*slopes[6]-8.0*afcums[2]*slopes[7]-afcums[1]*slopes[8];
		return afcums;
	}

	public double[][] cums2afcums(double[][] k,double[][] slopes,double[][] offsets){
		double[][] afcums=new double[9][9];
		double[] tempu1=new double[9];
		double[] off1=new double[9];
		double[] tempu2=new double[9];
		double[] off2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=k[i][0];
			tempu2[i]=k[0][i];
			off1[i]=offsets[i][0];
			off2[i]=offsets[0][i];
		}
		tempu1=cums2afcums(tempu1,slopes[0],off1);
		tempu2=cums2afcums(tempu2,slopes[1],off2);
		for(int i=0;i<9;i++){
			afcums[i][0]=tempu1[i];
			afcums[0][i]=tempu2[i];
		}
		afcums[1][1]=k[1][1]-offsets[1][1];
		afcums[2][1]=k[2][1]-offsets[2][1]-afcums[1][1]*slopes[0][2];
		afcums[1][2]=k[1][2]-offsets[1][2]-afcums[1][1]*slopes[1][2];
		afcums[3][1]=k[3][1]-offsets[3][1]-3.0*slopes[0][2]*afcums[2][1]-slopes[0][3]*afcums[1][1];
		afcums[1][3]=k[1][3]-offsets[1][3]-3.0*slopes[1][2]*afcums[1][2]-slopes[1][3]*afcums[1][1];
		afcums[2][2]=k[2][2]-offsets[2][2]-slopes[0][2]*afcums[1][2]-slopes[1][2]*afcums[2][1]-slopes[0][2]*slopes[1][2]*afcums[1][1];
		return afcums;
	}

	public double[] rawmoms2centmoms(double[] u){
		double[] centmoms=new double[9];
		centmoms[1]=u[1];
		centmoms[2]=-intpow(u[1],2)+u[2];
		centmoms[3]=2.0*intpow(u[1],3)-3.0*u[1]*u[2]+u[3];
		centmoms[4]=-3.0*intpow(u[1],4)+6.0*intpow(u[1],2)*u[2]-4.0*u[1]*u[3]+u[4];
		centmoms[5]=4.0*intpow(u[1],5)-10.0*intpow(u[1],3)*u[2]+10.0*intpow(u[1],2)*u[3]-5.0*u[1]*u[4]+u[5];
		centmoms[6]=-5.0*intpow(u[1],6)+15.0*intpow(u[1],4)*u[2]-20.0*intpow(u[1],3)*u[3]+15.0*intpow(u[1],2)*u[4]-6.0*u[1]*u[5]+u[6];
		centmoms[7]=6.0*intpow(u[1],7)-21.0*intpow(u[1],5)*u[2]+35.0*intpow(u[1],4)*u[3]-35.0*intpow(u[1],3)*u[4]+21.0*intpow(u[1],2)*u[5]-7.0*u[1]*u[6]+u[7];
		centmoms[8]=-7.0*intpow(u[1],8)+28.0*intpow(u[1],6)*u[2]-56.0*intpow(u[1],5)*u[3]+70.0*intpow(u[1],4)*u[4]-56.0*intpow(u[1],3)*u[5]+28.0*intpow(u[1],2)*u[6]-8.0*u[1]*u[7]+u[8];
		return centmoms;
	}

	public double[][] rawmoms2centmoms(double[][] u){
		double[][] centmoms=new double[9][9];
		double[] tempu1=new double[9];
		double[] tempu2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=u[i][0];
			tempu2[i]=u[0][i];
		}
		tempu1=rawmoms2centmoms(tempu1);
		tempu2=rawmoms2centmoms(tempu2);
		for(int i=0;i<9;i++){
			centmoms[i][0]=tempu1[i];
			centmoms[0][i]=tempu2[i];
		}
		centmoms[1][1]=u[1][1]-u[1][0]*u[0][1];
		centmoms[2][1]=2.0*u[0][1]*u[1][0]*u[1][0]-2.0*u[1][0]*u[1][1]-u[0][1]*u[2][0]+u[2][1];
		centmoms[1][2]=2.0*u[0][1]*u[0][1]*u[1][0]-u[0][2]*u[1][0]-2.0*u[0][1]*u[1][1]+u[1][2];
		centmoms[3][1]=-3.0*u[0][1]*u[1][0]*u[1][0]*u[1][0]+3.0*u[1][0]*u[1][0]*u[1][1]+3.0*u[0][1]*u[1][0]*u[2][0]-3.0*u[1][0]*u[2][1]-u[0][1]*u[3][0]+u[3][1];
		centmoms[1][3]=-3.0*u[0][1]*u[0][1]*u[0][1]*u[1][0]+3.0*u[0][1]*u[0][2]*u[1][0]-u[0][3]*u[1][0]+3.0*u[0][1]*u[0][1]*u[1][1]-3.0*u[0][1]*u[1][2]+u[1][3];
		centmoms[2][2]=-3.0*u[0][1]*u[0][1]*u[1][0]*u[1][0]+u[0][2]*u[1][0]*u[1][0]+4.0*u[0][1]*u[1][0]*u[1][1]-2.0*u[1][0]*u[1][2]+u[0][1]*u[0][1]*u[2][0]-2.0*u[0][1]*u[2][1]+u[2][2];
		return centmoms;
	}

	double[] centmoms2rawmoms(double[] uc){
		double[] rawmoms=new double[9];
		rawmoms[1]=uc[1];
		rawmoms[2]=intpow(uc[1],2)+uc[2];
		rawmoms[3]=intpow(uc[1],3)+3.0*uc[1]*uc[2]+uc[3];
		rawmoms[4]=intpow(uc[1],4)+6.0*intpow(uc[1],2)*uc[2]+4.0*uc[1]*uc[3]+uc[4];
		rawmoms[5]=intpow(uc[1],5)+10.0*intpow(uc[1],3)*uc[2]+10.0*intpow(uc[1],2)*uc[3]+5.0*uc[1]*uc[4]+uc[5];
		rawmoms[6]=intpow(uc[1],6)+15.0*intpow(uc[1],4)*uc[2]+20.0*intpow(uc[1],3)*uc[3]+15.0*intpow(uc[1],2)*uc[4]+6.0*uc[1]*uc[5]+uc[6];
		rawmoms[7]=intpow(uc[1],7)+21.0*intpow(uc[1],5)*uc[2]+35.0*intpow(uc[1],4)*uc[3]+35.0*intpow(uc[1],3)*uc[4]+21.0*intpow(uc[1],2)*uc[5]+7.0*uc[1]*uc[6]+uc[7];
		rawmoms[8]=intpow(uc[1],8)+28.0*intpow(uc[1],6)*uc[2]+56.0*intpow(uc[1],5)*uc[3]+70.0*intpow(uc[1],4)*uc[4]+56.0*intpow(uc[1],3)*uc[5]+28.0*intpow(uc[1],2)*uc[6]+8.0*uc[1]*uc[7]+uc[8];
		return rawmoms;
	}

	public double[][] centmoms2rawmoms(double[][] uc){
		double[][] rawmoms=new double[9][9];
		double[] tempu1=new double[9];
		double[] tempu2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=uc[i][0];
			tempu2[i]=uc[0][i];
		}
		tempu1=centmoms2rawmoms(tempu1);
		tempu2=centmoms2rawmoms(tempu2);
		for(int i=0;i<9;i++){
			rawmoms[i][0]=tempu1[i];
			rawmoms[0][i]=tempu2[i];
		}
		rawmoms[1][1]=uc[0][1]*uc[1][0]+uc[1][1];
		rawmoms[2][1]=uc[0][1]*uc[1][0]*uc[1][0]+2.0*uc[1][0]*uc[1][1]+uc[0][1]*uc[2][0]+uc[2][1];
		rawmoms[1][2]=uc[0][1]*uc[0][1]*uc[1][0]+uc[0][2]*uc[1][0]+2.0*uc[0][1]*uc[1][1]+uc[1][2];
		rawmoms[3][1]=uc[0][1]*uc[1][0]*uc[1][0]*uc[1][0]+3.0*uc[1][0]*uc[1][0]*uc[1][1]+3.0*uc[0][1]*uc[1][0]*uc[2][0]+3.0*uc[1][0]*uc[2][1]+uc[0][1]*uc[3][0]+uc[3][1];
		rawmoms[1][3]=uc[0][1]*uc[0][1]*uc[0][1]*uc[1][0]+3.0*uc[0][1]*uc[0][2]*uc[1][0]+uc[0][3]*uc[1][0]+3.0*uc[0][1]*uc[0][1]*uc[1][1]+3.0*uc[0][1]*uc[1][2]+uc[1][3];
		rawmoms[2][2]=uc[0][1]*uc[0][1]*uc[1][0]*uc[1][0]+uc[0][2]*uc[1][0]*uc[1][0]+4.0*uc[0][1]*uc[1][0]*uc[1][1]+2.0*uc[1][0]*uc[1][2]+uc[0][1]*uc[0][1]*uc[2][0]+2.0*uc[0][1]*uc[2][1]+uc[2][2];
		return rawmoms;
	}

	public double[] centmoms2cums(double[] uc){
		double[] cums=new double[9];
		cums[1]=uc[1];
		cums[2]=uc[2];
		cums[3]=uc[3];
		cums[4]=uc[4]-3.0*intpow(uc[2],2);
		cums[5]=uc[5]-10.0*uc[2]*uc[3];
		cums[6]=uc[6]-10.0*intpow(uc[3],2)-15.0*uc[2]*uc[4]+30.0*intpow(uc[2],3);
		cums[7]=210.0*intpow(uc[2],2)*uc[3]-35.0*uc[3]*uc[4]-21.0*uc[2]*uc[5]+uc[7];
		cums[8]=-630.0*intpow(uc[2],4)+560.0*uc[2]*intpow(uc[3],2)+420.0*intpow(uc[2],2)*uc[4]-35.0*intpow(uc[4],2)-56.0*uc[3]*uc[5]-28.0*uc[2]*uc[6]+uc[8];
		return cums;
	}

	public double[][] centmoms2cums(double[][] uc){
		double[][] cums=new double[9][9];
		double[] tempu1=new double[9];
		double[] tempu2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=uc[i][0];
			tempu2[i]=uc[0][i];
		}
		tempu1=centmoms2cums(tempu1);
		tempu2=centmoms2cums(tempu2);
		for(int i=0;i<9;i++){
			cums[i][0]=tempu1[i];
			cums[0][i]=tempu2[i];
		}
		cums[1][1]=uc[1][1];
		cums[2][1]=uc[2][1];
		cums[1][2]=uc[1][2];
		cums[3][1]=uc[3][1]-3.0*uc[1][1]*uc[2][0];
		cums[1][3]=uc[1][3]-3.0*uc[1][1]*uc[0][2];
		cums[2][2]=uc[2][2]-2.0*uc[1][1]*uc[1][1]-uc[2][0]*uc[0][2];
		return cums;
	}

	public double[] cums2centmoms(double[] K){
		double[] centmoms=new double[9];
		centmoms[1]=K[1];
		centmoms[2]=K[2];
		centmoms[3]=K[3];
		centmoms[4]=K[4]+3.0*intpow(K[2],2);
		centmoms[5]=K[5]+10.0*K[2]*K[3];
		centmoms[6]=K[6]+10.0*intpow(K[3],2)+15.0*K[2]*K[4]+15.0*intpow(K[2],3);
		centmoms[7]=105.0*intpow(K[2],2)*K[3]+35.0*K[3]*K[4]+21.0*K[2]*K[5]+K[7];
		centmoms[8]=105.0*intpow(K[2],4)+280.0*K[2]*intpow(K[3],2)+210.0*intpow(K[2],2)*K[4]+35.0*intpow(K[4],2)+56.0*K[3]*K[5]+28.0*K[2]*K[6]+K[8];
		return centmoms;
	}

	public double[][] cums2centmoms(double[][] K){
		double[][] centmoms=new double[9][9];
		double[] tempu1=new double[9];
		double[] tempu2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=K[i][0];
			tempu2[i]=K[0][i];
		}
		tempu1=cums2centmoms(tempu1);
		tempu2=cums2centmoms(tempu2);
		for(int i=0;i<9;i++){
			centmoms[i][0]=tempu1[i];
			centmoms[0][i]=tempu2[i];
		}
		centmoms[1][1]=K[1][1];
		centmoms[2][1]=K[2][1];
		centmoms[1][2]=K[1][2];
		centmoms[3][1]=K[3][1]+3.0*K[1][1]*K[2][0];
		centmoms[1][3]=K[1][3]-3.0*K[1][1]*K[0][2];
		centmoms[2][2]=K[2][2]+2.0*K[1][1]*K[1][1]+K[2][0]*K[0][2];
		return centmoms;
	}

	public double[] factcums2cums(double[] factcums){
		double[] cums=new double[9];
		cums[1]=factcums[1];
		cums[2]=factcums[1]+factcums[2];
		cums[3]=factcums[1]+3.0*factcums[2]+factcums[3];
		cums[4]=factcums[1]+7.0*factcums[2]+6.0*factcums[3]+factcums[4];
		cums[5]=factcums[1]+15.0*factcums[2]+25.0*factcums[3]+10.0*factcums[4]+factcums[5];
		cums[6]=factcums[1]+31.0*factcums[2]+90.0*factcums[3]+65.0*factcums[4]+15.0*factcums[5]+factcums[6];
		cums[7]=factcums[1]+63.0*factcums[2]+301.0*factcums[3]+350.0*factcums[4]+140.0*factcums[5]+21.0*factcums[6]+factcums[7];
		cums[8]=factcums[1]+127.0*factcums[2]+966.0*factcums[3]+1701.0*factcums[4]+1050.0*factcums[5]+266.0*factcums[6]+28.0*factcums[7]+factcums[8];
		return cums;
	}

	public double[][] factcums2cums(double[][] factcums){
		double[][] cums=new double[9][9];
		double[] tempu1=new double[9];
		double[] tempu2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=factcums[i][0];
			tempu2[i]=factcums[0][i];
		}
		tempu1=factcums2cums(tempu1);
		tempu2=factcums2cums(tempu2);
		for(int i=0;i<9;i++){
			cums[i][0]=tempu1[i];
			cums[0][i]=tempu2[i];
		}
		cums[1][1]=factcums[1][1];
		cums[2][1]=factcums[2][1]+factcums[1][1];
		cums[1][2]=factcums[1][2]+factcums[1][1];
		cums[3][1]=factcums[3][1]+3.0*factcums[2][1]+factcums[1][1];
		cums[1][3]=factcums[1][3]+3.0*factcums[1][2]+factcums[1][1];
		cums[2][2]=factcums[2][2]+factcums[2][1]+factcums[1][2]+factcums[1][1];
		return cums;
	}

	public double[] cums2factcums(double[] cums){
		double[] factcums=new double[9];
		factcums[1]=cums[1];
		factcums[2]=cums[2]-cums[1];
		factcums[3]=cums[3]-3.0*cums[2]+2.0*cums[1];
		factcums[4]=cums[4]-6.0*cums[3]+11.0*cums[2]-6.0*cums[1];
		factcums[5]=cums[5]-10.0*cums[4]+35.0*cums[3]-50.0*cums[2]+24.0*cums[1];
		factcums[6]=cums[6]-15.0*cums[5]+85.0*cums[4]-225.0*cums[3]+274.0*cums[2]-120.0*cums[1];
		factcums[7]=cums[7]-21.0*cums[6]+175.0*cums[5]-735.0*cums[4]+1624.0*cums[3]-1764.0*cums[2]+720.0*cums[1];
		factcums[8]=cums[8]-28.0*cums[7]+322.0*cums[6]-1960.0*cums[5]+6769.0*cums[4]-13132.0*cums[3]+13068.0*cums[2]-5040.0*cums[1];
		return factcums;
	}

	public double[][] cums2factcums(double[][] cums){
		double[][] factcums=new double[9][9];
		double[] tempu1=new double[9];
		double[] tempu2=new double[9];
		for(int i=0;i<9;i++){
			tempu1[i]=cums[i][0];
			tempu2[i]=cums[0][i];
		}
		tempu1=cums2factcums(tempu1);
		tempu2=cums2factcums(tempu2);
		for(int i=0;i<9;i++){
			factcums[i][0]=tempu1[i];
			factcums[0][i]=tempu2[i];
		}
		factcums[1][1]=cums[1][1];
		factcums[2][1]=cums[2][1]-cums[1][1];
		factcums[1][2]=cums[1][2]-cums[1][1];
		factcums[3][1]=cums[3][1]-3.0*cums[2][1]+cums[1][1];
		factcums[1][3]=cums[1][3]-3.0*cums[1][2]+cums[1][1];
		factcums[2][2]=cums[2][2]-cums[2][1]-cums[1][2]+cums[1][1];
		return cums;
	}

	public double[] gamma2GL(){
		double[] gammas=new double[9];
		for(int i=1;i<=8;i++){
			gammas[i]=factorial(4*i-4)/(intpow(4.0*Math.PI*Math.PI,i-1)*i*intpow(factorial(2*i-2),2));
		}
		return gammas;
	}

	public double[] gamma2GLxz(){
		double[] gammas=new double[9];
		for(int i=1;i<=8;i++){
			gammas[i]=(intpow(2.0,6*i-5)*intpow(factorial(2*i-2),2)*(2*i-1))/(Math.sqrt(i)*intpow(Math.PI,2*i-2)*factorial(4*i-2));
		}
		return gammas;
	}

	public double[] gamma3DG(){
		double[] gammas=new double[9];
		for(int i=1;i<=8;i++){
			gammas[i]=Math.pow(i,-1.5);
		}
		return gammas;
	}

	public double[] gamma2DG(){
		double[] gammas=new double[9];
		for(int i=1;i<=8;i++){
			gammas[i]=1.0/i;
		}
		return gammas;
	}

	public double[] pch_cumulants(float[] brightnesses,float[] numbers,double[] gammas){
		// first calculate the factorial cumulants
		double[] factcums=new double[9];
		for(int i=1;i<=8;i++){
			for(int j=0;j<brightnesses.length;j++){
				factcums[i]+=gammas[i]*intpow(brightnesses[j],i)*numbers[j];
			}
		}
		// now calculate the cumulants from these
		return factcums2cums(factcums);
	}

	public double pch_cumulant(float[] brightnesses,float[] numbers,double[] gammas,int order){
		double[] temp=pch_cumulants(brightnesses,numbers,gammas);
		return temp[order];
	}

	double intpow(double value,int power){
		double retval=value;
		for(int i=1;i<power;i++){
			retval*=value;
		}
		if(power==0){
			retval=1.0;
		}
		return retval;
	}

	int factorial(int n){
		int retval=1;
		for(int i=n;i>1;i--){
			retval*=i;
		}
		return retval;
	}

	public double[] two_comp_n_b(double[] k,double gamma2,double gamma3,double e1,double mine2,double maxe2,double delta){
		double[] params=new double[3];
		double minc2=1000.0;
		double[] gammas=new double[9];
		gammas[1]=1.0;
		gammas[2]=gamma2;
		gammas[3]=gamma3;
		int pts=1+(int)((maxe2-mine2)/delta);
		boolean underzero=true;
		for(int i=0;i<pts;i++){
			double e2=mine2+delta*i;
			double N2=(k[2]-k[1]-gamma2*k[1]*e1)/(gamma2*e2*(e2-e1));
			double N1=(k[1]-e2*N2)/e1;
			float[] brightnesses={(float)e1,(float)e2};
			float[] numbers={(float)N1,(float)N2};
			double k3fit=pch_cumulant(brightnesses,numbers,gammas,3);
			double c2=Math.abs(k3fit-k[3]);
			if(underzero){
				if(N1>0.001&&N2>0.001){
					underzero=false;
					minc2=c2;
					params[0]=N1;
					params[1]=N2;
					params[2]=e2;
				}
			}else{
				if(c2<minc2&&(N1>0.001&&N2>0.001)){
					minc2=c2;
					params[0]=N1;
					params[1]=N2;
					params[2]=e2;
				}
			}
		}
		if(underzero){
			return null;
		}else{
			return params;
		}
	}

	public double[] gammaval(double[] k,double gamma2){
		double[] factcums=cums2factcums(k);
		double[] gvals=new double[9];
		gvals[1]=1.0;
		gvals[2]=gamma2;
		double e=factcums[2]/(gamma2*factcums[1]);
		double N=factcums[1]/e;
		gvals[3]=factcums[3]/(e*e*e*N);
		gvals[4]=factcums[4]/(e*e*e*e*N);
		gvals[5]=factcums[5]/(e*e*e*e*e*N);
		gvals[6]=factcums[6]/(e*e*e*e*e*e*N);
		gvals[7]=factcums[7]/(e*e*e*e*e*e*e*N);
		gvals[8]=factcums[8]/(e*e*e*e*e*e*e*e*N);
		return gvals;
	}

	public double[] analog_gammaval(double[] k,double gamma2,double[] offsets,double[] slopes,double gain){
		double[] afcums=cums2afcums(k,slopes,offsets);
		double[] gvals=new double[9];
		gvals[1]=1.0;
		gvals[2]=gamma2;
		double eg=afcums[2]/(gamma2*afcums[1]);
		double N=afcums[1]/eg;
		gvals[3]=afcums[3]/(eg*eg*eg*N);
		gvals[4]=afcums[4]/(eg*eg*eg*eg*N);
		gvals[5]=afcums[5]/(eg*eg*eg*eg*eg*N);
		gvals[6]=afcums[6]/(eg*eg*eg*eg*eg*eg*N);
		gvals[7]=afcums[7]/(eg*eg*eg*eg*eg*eg*eg*N);
		gvals[8]=afcums[8]/(eg*eg*eg*eg*eg*eg*eg*eg*N);
		return gvals;
	}

	public double[] concs_n_b(double[] brightnesses,double[] gammas,double[] k,int order){
		double[] factk=cums2factcums(k);
		double[][] A=new double[order][order];
		double[] b=new double[order];
		for(int i=0;i<order;i++){
			for(int j=0;j<order;j++){
				A[i][j]=gammas[i+1]*intpow(brightnesses[j],i+1);
			}
			b[i]=factk[i+1];
		}
		double[] x=new double[order];
		new matrixsolve().gjsolve(A,b,x,order);
		return x;
	}

	public double[] oligomers_n_b(double monbright,double[] gammas,double[] k,int order){
		double[] temp={monbright,2.0*monbright,3.0*monbright,4.0*monbright};
		return concs_n_b(temp,gammas,k,order);
	}

	public float[] is_single_species(double[] k,int psftype,double gamma3,int npts){
		// this method tests whether or not a data set has a single species
		double[] gammas;
		if(psftype==0){
			gammas=gamma3DG();
		}else{
			if(psftype==1){
				gammas=gamma2DG();
			}else{
				if(psftype==2){
					gammas=gamma2GL();
				}else{
					gammas=gamma2GLxz();
				}
			}
		}
		gammas[3]=gamma3;
		// first calculate the brightness
		float[] brightness=new float[1];
		brightness[0]=(float)(((k[2]/k[1])-1.0)/gammas[2]);
		float[] number=new float[1];
		number[0]=(float)k[1]/brightness[0];

		// now calculate what the k[3] value would be for single species and its
		// error
		double[] singlekvals=pch_cumulants(brightness,number,gammas);
		double[] singlekvars=kstatsvars(singlekvals,npts);
		double stdev=Math.sqrt(singlekvars[3]);
		double diff=Math.abs((singlekvals[3]-k[3])/stdev);
		float[] out=new float[5];
		out[0]=brightness[0];
		out[1]=(float)diff;
		out[2]=(float)singlekvals[3];
		out[3]=(float)stdev;
		out[4]=(float)k[3];
		return out;
	}

	public float[] td_brightnesses(double td1,double td2,double[] kT1,double[] kT2,double minl2,double maxl2,double step,double gamma2,double gamma3){
		double r=5.0;
		binning_function bf=new binning_function(r);
		double[] kfactT1=cums2factcums(kT1);
		double[] kfactT2=cums2factcums(kT2);
		boolean firstc2=false;
		double B1T1=bf.b2(td1,td1);
		double B1T2=bf.b2(td1,td2);
		double B2T1=bf.b2(td2,td1);
		double B2T2=bf.b2(td2,td2);
		double B31T1=bf.b3(td1,td1);
		// double B31T2=bf.b3(td1, td2);
		double B32T1=bf.b3(td2,td1);
		// double B32T2=bf.b3(td2, td2);
		double c2min=0.0;
		double l1min=0.0;
		double N1min=0.0;
		double l2min=0.0;
		double N2min=0.0;
		for(double l2=minl2;l2<=maxl2;l2+=step){
			double l1=(B2T1*l2*kfactT2[2]*td1-B2T2*l2*kfactT1[2]*td1)/(B1T2*B2T1*gamma2*kfactT1[1]*l2-B1T1*B2T2*gamma2*kfactT1[1]*l2-B1T2*kfactT1[2]*td1+B1T1*kfactT2[2]*td1);
			double N2=(kfactT1[2]*td1-B1T1*gamma2*kfactT1[1]*l1)/(gamma2*l2*td1*(B2T1*l2-B1T1*l1));
			double N1=(kfactT1[1]-l2*N2*td1)/(l1*td1);
			double K3T1=gamma3*(l1*l1*l1*N1*B31T1+l2*l2*l2*N2*B32T1);
			double K3T2=gamma3*(l1*l1*l1*N1*B31T1+l2*l2*l2*N2*B32T1);
			double c2=0.0;
			if(td1<td2){
				c2=Math.abs(K3T1-kfactT1[3]);
			}else{
				c2=Math.abs(K3T2-kfactT2[3]);
			}
			if(!firstc2){
				if(l1>0.0&&N1>0.0&&N2>0.0){
					c2min=c2;
					l1min=l1;
					N1min=N1;
					N2min=N2;
					l2min=l2;
					firstc2=true;
				}
			}else{
				if(c2<c2min){
					c2min=c2;
					l1min=l1;
					N1min=N1;
					N2min=N2;
					l2min=l2;
				}
			}
		}
		if(!firstc2){
			return null;
		}else{
			float[] temp={(float)l1min,(float)N1min,(float)l2min,(float)N2min,(float)c2min};
			return temp;
		}
	}

}
