/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfft;

public class binmultilog{
	public int tauwidth;

	/*
	 * this class quasilogarithmically bins a data set the bin size is doubled
	 * every tauwidth bins Copyright Jay Unruh Stowers Institute for Medical
	 * Research 4/25/08
	 */
	public binmultilog(){
		tauwidth=8;
	}

	public binmultilog(int tauwidth1){
		tauwidth=tauwidth1;
	}

	public float[] dobinmultilog(float[] ac,int length){
		// length here is half the length of the original autocorrelated data
		// set
		float[] tempac=new float[ac.length];
		System.arraycopy(ac,0,tempac,0,ac.length);
		int binsize=2;
		int counter1=tauwidth+1;
		int counter2=tauwidth+2;
		int counter3=1;
		tempac[tauwidth]=(tempac[tauwidth]+tempac[tauwidth+1])/2.0f;
		do{
			tempac[counter1]=0.0f;
			for(int i=0;i<binsize;i++){
				tempac[counter1]+=tempac[counter2+i]/binsize;
			}
			counter1++;
			counter2+=binsize;
			counter3++;
			if(counter3>=tauwidth){
				binsize+=binsize;
				counter3=0;
			}
		}while((counter2+binsize)<length);
		int size=counter1-1;
		float[] newresult=new float[size];
		for(int i=0;i<size;i++){
			newresult[i]=tempac[i];
		}
		return newresult;
	}

	public float[] unbinmultilog(float[] binned){
		// length here is half the length of the original autocorrelated data
		// set
		// here we set the points in every bin equal to the average of the bin
		// this allows for arbitrary rebinning of correlation data
		int length=getoriglength(binned.length);
		float[] tempac=new float[length];
		System.arraycopy(binned,0,tempac,0,tauwidth);
		int binsize=2;
		int counter1=tauwidth+1;
		int counter2=tauwidth+2;
		int counter3=1;
		tempac[tauwidth]=(tempac[tauwidth]+tempac[tauwidth+1])/2.0f;
		tempac[tauwidth]=binned[tauwidth];
		tempac[tauwidth+1]=binned[tauwidth];
		do{
			for(int i=0;i<binsize;i++){
				tempac[counter2+i]=binned[counter1];
			}
			counter1++;
			counter2+=binsize;
			counter3++;
			if(counter3>=tauwidth){
				binsize+=binsize;
				counter3=0;
			}
			if(counter1>=binned.length)
				counter1=binned.length-1;
		}while((counter2+binsize)<length);
		return tempac;
	}

	public float[] rebin_data(float[] prebinned,int binsize){
		float[] orig=unbinmultilog(prebinned);
		int binlength=(int)((float)orig.length/(float)binsize);
		float[] newbinned=new float[binlength];
		for(int i=0;i<binlength;i++){
			for(int j=0;j<binsize;j++){
				newbinned[i]+=orig[i*binsize+j]/binsize;
			}
		}
		return dobinmultilog(newbinned,binlength);
	}

	public float[] getxvals(int length){
		// length here is half the length of the original autocorrelated data
		// set
		// this returns the x axis of the binned autocorrelation
		float[] xvals=new float[length];
		for(int i=0;i<length;i++){
			xvals[i]=i;
		}
		int binsize=2;
		int counter1=tauwidth+1;
		int counter2=tauwidth+2;
		int counter3=1;
		xvals[tauwidth]=(xvals[tauwidth]+xvals[tauwidth+1])/2.0f;
		do{
			xvals[counter1]=0.0f;
			for(int i=0;i<binsize;i++){
				xvals[counter1]+=xvals[counter2+i]/binsize;
			}
			counter1++;
			counter2+=binsize;
			counter3++;
			if(counter3>=tauwidth){
				binsize+=binsize;
				counter3=0;
			}
		}while((counter2+binsize)<length);
		int size=counter1-1;
		float[] newxvals=new float[size];
		for(int i=0;i<size;i++){
			newxvals[i]=xvals[i];
		}
		return newxvals;
	}

	public float[] getxvals2(int binnedlength){
		// here we get the xvals array from the length of the binned
		// autocorrelation
		return getxvals(getoriglength(binnedlength));
	}

	public int getoriglength(int binnedlength){
		// here we get the length of the original autocorrelation from the
		// length of the binned autocorrelation
		int binsize=2;
		int counter1=tauwidth+1;
		int counter2=tauwidth+2;
		int counter3=1;
		do{
			counter1++;
			counter2+=binsize;
			counter3++;
			if(counter3>=tauwidth){
				binsize+=binsize;
				counter3=0;
			}
		}while((counter1-1)<binnedlength);
		int length=counter2+binsize;
		return length;
	}
}
