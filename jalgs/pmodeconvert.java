/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class pmodeconvert{

	public float[] pm2tm(int[] pmdata,double sfreq,int pmfreq){
		// this method converts photon mode data to time mode
		int total_time=0;
		for(int i=0;i<pmdata.length;i++){
			total_time+=pmdata[i];
		}
		double chan_clocks=pmfreq/sfreq;
		int numpts=(int)(total_time/chan_clocks)-1;
		float[] data=new float[numpts];

		int j=0;
		double temptime=0.0;
		int photons=0;
		for(int i=0;i<(numpts-1);i++){
			while(temptime<chan_clocks){
				temptime+=pmdata[j];
				photons++;
				j++;
			}
			data[i]=photons-1;
			temptime-=chan_clocks;
			photons=1;
		}
		while(temptime<chan_clocks&&j<pmdata.length){
			temptime+=pmdata[j];
			photons++;
			j++;
		}
		if(j>=pmdata.length){
			data[numpts-1]=photons;
		}else{
			data[numpts-1]=photons-1;
		}
		return data;
	}

	public float[] pm2tm(long[] pmdata,double sfreq,long pmfreq){
		// this method converts photon mode data to time mode
		long total_time=0L;
		for(int i=0;i<pmdata.length;i++){
			total_time+=pmdata[i];
		}
		double chan_clocks=pmfreq/sfreq;
		int numpts=(int)(total_time/chan_clocks)-1;
		float[] data=new float[numpts];

		int j=0;
		double temptime=0.0;
		int photons=0;
		for(int i=0;i<(numpts-1);i++){
			while(temptime<chan_clocks){
				temptime+=pmdata[j];
				photons++;
				j++;
			}
			data[i]=photons-1;
			temptime-=chan_clocks;
			photons=1;
		}
		while(temptime<chan_clocks&&j<pmdata.length){
			temptime+=pmdata[j];
			photons++;
			j++;
		}
		if(j>=pmdata.length){
			data[numpts-1]=photons;
		}else{
			data[numpts-1]=photons-1;
		}
		return data;
	}

	public float[] pm2tm(int[] pmdata,double sfreq,int pmfreq,double ill_delay){
		// this method converts photon mode data to time mode with an
		// illumination delay for laser switching
		// here we essentially have 3 time bins for each actual time bin
		// one at the beginning and end for the delay and the other for the
		// actual data
		int total_time=0;
		for(int i=0;i<pmdata.length;i++){
			total_time+=pmdata[i];
		}
		double chan_clocks=pmfreq/sfreq;
		int numpts=(int)(total_time/chan_clocks)-1;
		float[] data=new float[numpts];
		double frac_off=ill_delay*sfreq;
		double clocks_off=chan_clocks*frac_off;
		double clocks_on=chan_clocks-2.0*clocks_off;

		int j=0;
		double temptime=0.0;
		int photons=0;
		for(int i=0;i<(numpts-1);i++){
			while(temptime<clocks_off){
				temptime+=pmdata[j];
				j++;
			}
			temptime-=clocks_off;
			while(temptime<clocks_on){
				photons++;
				temptime+=pmdata[j];
				j++;
			}
			data[i]=photons-1;
			temptime-=clocks_on;
			photons=1;
			while(temptime<clocks_off){
				temptime+=pmdata[j];
				j++;
			}
			temptime-=clocks_off;
		}
		while(temptime<clocks_off&&j<pmdata.length){
			temptime+=pmdata[j];
			j++;
		}
		temptime-=clocks_off;
		while(temptime<chan_clocks&&j<pmdata.length){
			temptime+=pmdata[j];
			photons++;
			j++;
		}
		if(j>=pmdata.length){
			data[numpts-1]=photons;
		}else{
			data[numpts-1]=photons-1;
		}
		return data;
	}

	public float[] pm2tm(int[] pmdata,double sfreq,int pmfreq,int offset){
		// this method converts photon mode data to time mode
		// the offset is in inverse pmfreq units and simply gets removed from
		// the data set
		return pm2tm(subpmoffset(pmdata,offset),sfreq,pmfreq);
	}

	public float[] pm2tm(int[] pmdata,double sfreq,int pmfreq,int offset,double ill_delay){
		// this method converts photon mode data to time mode
		// the offset is in inverse pmfreq units and simply gets removed from
		// the data set
		return pm2tm(subpmoffset(pmdata,offset),sfreq,pmfreq,ill_delay);
	}

	public int[] subpmoffset(int[] pmdata,int offset){
		int tempoff=offset;
		int counter=0;
		while(tempoff>=0){
			tempoff-=pmdata[counter];
			counter++;
		}

		int[] newpmdata=new int[pmdata.length-counter+1];
		newpmdata[0]=-tempoff;
		System.arraycopy(pmdata,counter,newpmdata,1,pmdata.length-counter);
		return newpmdata;
	}

	public Object[] subdivide_pmdata(int[] pmdata,double sfreq,double sublength){
		// sublength here is in seconds
		int sublength2=(int)(sublength*sfreq); // sublength 2 is in channel
		// clocks
		int total_time=0;
		for(int i=0;i<pmdata.length;i++){
			total_time+=pmdata[i];
		}
		int sections=(int)((double)total_time/(double)sublength2);
		if((sublength2*sections)>=total_time){
			sections--;
		}
		int secttemplength=2*pmdata.length/sections;
		Object[] outdata=new Object[sections];
		int counter1=1;
		int remainder=pmdata[0];
		for(int i=0;i<sections;i++){
			int[] subdata=new int[secttemplength];
			subdata[0]=remainder;
			int tempoff=sublength2-remainder;
			int counter=1;
			while(tempoff>=0){
				tempoff-=pmdata[counter1];
				if(tempoff>=0){
					subdata[counter]=pmdata[counter1];
				}
				counter1++;
				counter++;
			}
			remainder=(-tempoff);
			int[] subdata2=new int[counter];
			System.arraycopy(subdata,0,subdata2,0,counter);
			outdata[i]=subdata2;
		}
		return outdata;
	}

	public float[][] pm2tm_alex(int[] pmdata,int[] pmdata2,double swfreq,int pmfreq,int offset,int binby){
		// this method converts photon mode data to time mode
		// for wavelength switching data. The switching frequency must be known
		// here the phase offset must be known as well
		double swfreq2=2.0*swfreq;
		int chan_clocks=(int)((pmfreq)/swfreq2);
		float[] tmdata1=pm2tm(pmdata,swfreq2,pmfreq,offset);
		float[] tmdata2=pm2tm(pmdata2,swfreq2,pmfreq,offset+chan_clocks);
		int length=tmdata1.length;
		if(tmdata2.length<length){
			length=tmdata2.length;
		}
		int newlength=length/(2*binby);
		float[][] retvals=new float[2][newlength];
		for(int i=0;i<newlength;i++){
			for(int j=0;j<binby;j++){
				retvals[0][i]+=tmdata1[2*(i*binby+j)];
				retvals[1][i]+=tmdata2[2*(i*binby+j)];
			}
		}
		return retvals;
	}

	public float[][] pm2tm_alex(int[] pmdata,int[] pmdata2,double swfreq,int pmfreq,int offset,int binby,double ill_delay){
		// this method converts photon mode data to time mode
		// for wavelength switching data. The switching frequency must be known
		// here the phase offset must be known as well
		double swfreq2=2.0*swfreq;
		int chan_clocks=(int)((pmfreq)/swfreq2);
		float[] tmdata1=pm2tm(pmdata,swfreq2,pmfreq,offset,ill_delay);
		float[] tmdata2=pm2tm(pmdata2,swfreq2,pmfreq,offset+chan_clocks,ill_delay);
		int length=tmdata1.length;
		if(tmdata2.length<length){
			length=tmdata2.length;
		}
		int newlength=length/(2*binby);
		float[][] retvals=new float[2][newlength];
		for(int i=0;i<newlength;i++){
			for(int j=0;j<binby;j++){
				retvals[0][i]+=tmdata1[2*(i*binby+j)];
				retvals[1][i]+=tmdata2[2*(i*binby+j)];
			}
		}
		return retvals;
	}

	public float[][] pm2tm_multiposition(int[] pmdata2,int[] markers,int pmfreq){
		// this method converts photon mode to time mode for position switching
		// data
		// here we use the markers to delimit the data (will only work for slow
		// switching)
		// first find out what the offset is to the first position zero
		int counter=1;
		int offset=markers[0];
		while(markers[counter]!=0x0000000000){
			offset+=markers[counter+1];
			counter+=2;
		}
		offset+=markers[counter+1];
		int startpt=counter+3;
		int numpts=(markers.length-startpt)/4;
		float[] data=new float[numpts];
		int j=0;
		double temptime=0.0;
		int photons=0;
		int[] pmdata=subpmoffset(pmdata2,offset);
		double avgperiod=0.0;
		for(int i=0;i<(numpts-1);i++){
			// read in the photons
			while(temptime<markers[startpt+i*4]){
				temptime+=pmdata[j];
				photons++;
				j++;
			}
			data[i]=photons-1;
			temptime-=markers[startpt+i*4];
			// read in the scanner motion period
			while(temptime<markers[startpt+i*4+2]){
				temptime+=pmdata[j];
				j++;
			}
			temptime-=markers[startpt+i*4+2];
			photons=1;
			avgperiod+=(markers[startpt+i*4]+markers[startpt+i*4+2]);
		}
		while(temptime<markers[startpt+(numpts-1)*4]&&j<pmdata.length){
			temptime+=pmdata[j];
			photons++;
			j++;
		}
		if(j>=pmdata.length){
			data[numpts-1]=photons;
		}else{
			data[numpts-1]=photons-1;
		}
		float[][] tmdata=new float[3][];
		tmdata[0]=new float[numpts/2];
		tmdata[1]=new float[numpts/2];
		tmdata[2]=new float[1];
		avgperiod/=numpts-1;
		tmdata[2][0]=(float)(pmfreq/(2.0*avgperiod));
		for(int i=0;i<numpts/2;i++){
			tmdata[0][i]=data[2*i];
			tmdata[1][i]=data[2*i+1];
		}
		return tmdata;
	}

	public int calc_wlswitch_offset(int[] pmdata,double swfreq,int pmfreq){
		float phase=calc_wlswitch_phase(pmdata,swfreq,pmfreq);
		float phasepercent=phase/360.0f;
		double offset=((double)phasepercent*(double)pmfreq)/swfreq;
		return (int)offset;
	}

	public float[] wlswitch_phase_hist(int[] pmdata,double swfreq,int pmfreq){
		double period=pmfreq/swfreq;
		double currtime=0.0;
		float[] hist=new float[100];
		for(int i=0;i<pmdata.length;i++){
			currtime+=pmdata[i];
			double pos=currtime%period;
			int pos2=(int)(99.0*pos/period);
			hist[pos2]+=1.0;
		}
		return hist;
	}

	public float[] wlswitch_phase_hist(int[] pmdata,double swfreq,int pmfreq,double ill_delay){
		double period=pmfreq/swfreq;
		double frac_delay=ill_delay*swfreq;
		double currtime=0.0;
		float[] hist=new float[100];
		for(int i=0;i<pmdata.length;i++){
			currtime+=pmdata[i];
			double pos=currtime%period;
			double frac=pos/period;
			if((frac>frac_delay&&frac<(0.5-frac_delay))||(frac>(0.5+frac_delay)&&frac<(1.0-frac_delay))){
				int pos2=(int)(99.0*frac);
				hist[pos2]+=1.0;
			}
		}
		return hist;
	}

	public float calc_wlswitch_phase(int[] pmdata,double swfreq,int pmfreq){
		double swtime2=2.0e7/swfreq;
		double currtime=0.0;
		double G=0.0;
		double S=0.0;
		for(int i=0;i<pmdata.length;i++){
			currtime+=pmdata[i];
			G+=Math.cos(2.0*Math.PI*(currtime/swtime2));
			S+=Math.sin(2.0*Math.PI*(currtime/swtime2));
			if(currtime>swtime2){
				currtime%=swtime2;
			}
		}
		float phase=(float)((Math.atan2(S,G)*180.0)/(Math.PI));
		phase-=90.0f;
		if(phase<0){
			phase+=360.0f;
		}
		if(phase>=360.0f){
			phase-=360.0f;
		}
		return phase;
	}

	public int calc_wlswitch_period(int[] pmdata,int pmfreq,int minperiod,int maxperiod){
		// here we use correlation and a grid search to find the period
		// of switching with pmfreq accuracy
		int total_time=0;
		for(int i=0;i<pmdata.length;i++){
			total_time+=pmdata[i];
		}
		int endtime=0;
		int endcounter=0;
		for(int i=(pmdata.length-1);i>=0;i--){
			endtime+=pmdata[i];
			endcounter++;
			if(endtime>maxperiod){
				break;
			}
		}
		int analysislength=pmdata.length-endcounter;
		int maxcorrsum=0;
		int maxcorrperiod=minperiod;
		for(int i=minperiod;i<=maxperiod;i++){
			int corrsum=0;
			for(int j=0;j<analysislength;j++){
				int thistime=0;
				int thiscounter=0;
				while(thistime<i){
					thistime+=pmdata[j+thiscounter];
					thiscounter++;
				}
				if(thistime==i){
					corrsum++;
				}
			}
			if(corrsum>maxcorrsum){
				maxcorrsum=corrsum;
				maxcorrperiod=i;
			}
		}
		return maxcorrperiod;
	}

	public double calc_wlswitch_freq(int[] pmdata,int pmfreq,double startfreq,double endfreq){
		int minperiod=(int)(pmfreq/endfreq);
		int maxperiod=1+(int)(pmfreq/startfreq);
		int period=calc_wlswitch_period(pmdata,pmfreq,minperiod,maxperiod);
		return (double)pmfreq/(double)period;
	}

	public float[] pm2pch(int[] pmdata,double sfreq,int pmfreq){
		float[] tmdata=pm2tm(pmdata,sfreq,pmfreq);
		return create_histogram(tmdata);
	}

	public float[][] pm2pch(int[] pmdata1,int[] pmdata2,double sfreq,int pmfreq){
		float[] tmdata1=pm2tm(pmdata1,sfreq,pmfreq);
		float[] tmdata2=pm2tm(pmdata2,sfreq,pmfreq);
		return create_2Dhistogram(tmdata1,tmdata2);
	}

	public float[] create_histogram(float[] tmdata){
		double[] hist=new double[65536];
		int max=0;
		for(int i=0;i<tmdata.length;i++){
			int intsignal=(int)tmdata[i];
			if(intsignal<=65535){
				hist[intsignal]+=1.0;
				if(intsignal>max){
					max=intsignal;
				}
			}
		}
		float[] newhist=new float[max+1];
		for(int i=0;i<=max;i++){
			newhist[i]=(float)hist[i];
		}
		return newhist;
	}

	public float[][] create_2Dhistogram(float[] tmdata1,float[] tmdata2){
		double[][] hist=new double[4096][4096];
		int xmax=0;
		int ymax=0;
		int length=tmdata1.length;
		if(tmdata2.length<length){
			length=tmdata2.length;
		}
		for(int i=0;i<length;i++){
			int intsignal1=(int)tmdata1[i];
			int intsignal2=(int)tmdata2[i];
			if(intsignal1<=4095&&intsignal2<=4095){
				hist[intsignal1][intsignal2]+=1.0;
				if(intsignal1>xmax){
					xmax=intsignal1;
				}
				if(intsignal2>ymax){
					ymax=intsignal2;
				}
			}
		}
		float[][] newhist=new float[xmax+1][ymax+1];
		for(int i=0;i<=xmax;i++){
			for(int j=0;j<=ymax;j++){
				newhist[i][j]=(float)hist[i][j];
			}
		}
		return newhist;
	}

	public float[][] create_2Dhistogram256(float[] tmdata1,float[] tmdata2){
		double[][] hist=new double[256][256];
		int xmax=0;
		int ymax=0;
		int length=tmdata1.length;
		if(tmdata2.length<length){
			length=tmdata2.length;
		}
		for(int i=0;i<length;i++){
			int intsignal1=(int)tmdata1[i];
			int intsignal2=(int)tmdata2[i];
			if(intsignal1<256&&intsignal2<256){
				hist[intsignal1][intsignal2]+=1.0;
				if(intsignal1>xmax){
					xmax=intsignal1;
				}
				if(intsignal2>ymax){
					ymax=intsignal2;
				}
			}
		}
		float[][] newhist=new float[xmax+1][ymax+1];
		for(int i=0;i<=xmax;i++){
			for(int j=0;j<=ymax;j++){
				newhist[i][j]=(float)hist[i][j];
			}
		}
		return newhist;
	}

	public int[] tm2pm(float[] tmtraj,double xinc,double pmfreq){
		// lets try to evenly space the photons over each bin
		int binclocks=(int)(xinc*pmfreq);
		int totpts=0;
		for(int i=0;i<tmtraj.length;i++)
			totpts+=(int)tmtraj[i];
		int[] pmdata=new int[totpts-1];
		// the very first photon is implicit
		int counter=0;
		for(int i=0;i<tmtraj.length;i++){
			int photons=(int)tmtraj[i];
			int spacing=(int)((double)binclocks/(double)photons);
			int thisspacing=0;
			for(int j=1;j<photons;j++){
				if(counter>=totpts)
					break;
				pmdata[counter]=spacing;
				thisspacing+=spacing;
				counter++;
			}
			// now put the first photon of the next bin
			if(counter>=totpts)
				break;
			pmdata[counter]=binclocks-thisspacing;
			counter++;
		}
		return pmdata;
	}

	public float[] tm2pmf(float[] tmtraj,double xinc,double pmfreq){
		// lets try to evenly space the photons over each bin
		int binclocks=(int)(xinc*pmfreq);
		int totpts=0;
		for(int i=0;i<tmtraj.length;i++)
			totpts+=(int)tmtraj[i];
		float[] pmdata=new float[totpts-1];
		// the very first photon is implicit
		int counter=0;
		for(int i=0;i<tmtraj.length;i++){
			int photons=(int)tmtraj[i];
			int spacing=(int)((double)binclocks/(double)photons);
			int thisspacing=0;
			for(int j=1;j<photons;j++){
				if(counter>=(totpts-1))
					break;
				pmdata[counter]=spacing;
				thisspacing+=spacing;
				counter++;
			}
			// now put the first photon of the next bin
			if(counter>=(totpts-1))
				break;
			pmdata[counter]=binclocks-thisspacing;
			counter++;
		}
		return pmdata;
	}

}
