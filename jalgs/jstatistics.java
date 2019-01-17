/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.awt.Polygon;
import java.awt.Rectangle;
import java.util.Arrays;

public class jstatistics{
	public static final String[] stats={"Avg","Sum","Max","Min","Variance","Median","Mode","StDev","StErr","RelErr","Count","Percentile","ConditionalAvg","Not0Avg","Not0StDev","Not0StErr","Not0Sum",
			"Not0Count","Not0Min","Identity","Zero","MaxPos","AvgPos","NaN"};
	public static final String[] stats2={"avg","sum","max","min","variance","median","mode","stdev","sterr","relerr","count","percentile","conditionalavg","not0avg","not0stdev","not0sterr","not0sum",
			"not0count","not0min","identity","zero","maxpos","avgpos","nan"};

	public static float getstatistic(String stat,Object data,int width,int height,Rectangle r,float[] extras){
		if(data instanceof float[]){
			float[] tempdata=(float[])data;
			if(r!=null)
				tempdata=getrect((float[])data,width,height,r);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof short[]){
			short[] tempdata=(short[])data;
			if(r!=null)
				tempdata=getrect((short[])data,width,height,r);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof byte[]){
			byte[] tempdata=(byte[])data;
			if(r!=null)
				tempdata=getrect((byte[])data,width,height,r);
			return getstat(stat,tempdata,extras);
		}else{
			int[] tempdata=(int[])data;
			if(r!=null)
				tempdata=getrect((int[])data,width,height,r);
			return getstat(stat,tempdata,extras);
		}
	}
	
	public static float[] getspectrum(String stat,Object[] data,int width,int height,Rectangle r,float[] extras){
		float[] spectrum=new float[data.length];
		for(int i=0;i<data.length;i++){
			spectrum[i]=getstatistic(stat,data[i],width,height,r,extras);
		}
		return spectrum;
	}


	public static float getstatistic(String stat,Object data,int width,int height,boolean[] mask,float[] extras){
		if(data instanceof float[]){
			float[] tempdata=(float[])data;
			if(mask!=null)
				tempdata=getmask((float[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof short[]){
			short[] tempdata=(short[])data;
			if(mask!=null)
				tempdata=getmask((short[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof byte[]){
			byte[] tempdata=(byte[])data;
			if(mask!=null)
				tempdata=getmask((byte[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}else{
			int[] tempdata=(int[])data;
			if(mask!=null)
				tempdata=getmask((int[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}
	}
	
	public static float[] getspectrum(String stat,Object[] data,int width,int height,boolean[] mask,float[] extras){
		float[] spectrum=new float[data.length];
		for(int i=0;i<data.length;i++){
			spectrum[i]=getstatistic(stat,data[i],width,height,mask,extras);
		}
		return spectrum;
	}
	
	/***************
	 * 
	 * @param stat
	 * @param data
	 * @param width
	 * @param height
	 * @param mask
	 * @param lims: int array with xmin,xmax,ymin,ymax inclusive
	 * @param extras
	 * @return
	 */
	public static float getstatistic(String stat,Object data,int width,int height,boolean[] mask,int[] lims,float[] extras){
		if(data instanceof float[]){
			float[] tempdata=(float[])data;
			if(mask!=null)
				tempdata=getmask((float[])data,width,height,mask,lims);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof short[]){
			short[] tempdata=(short[])data;
			if(mask!=null)
				tempdata=getmask((short[])data,width,height,mask,lims);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof byte[]){
			byte[] tempdata=(byte[])data;
			if(mask!=null)
				tempdata=getmask((byte[])data,width,height,mask,lims);
			return getstat(stat,tempdata,extras);
		}else{
			int[] tempdata=(int[])data;
			if(mask!=null)
				tempdata=getmask((int[])data,width,height,mask,lims);
			return getstat(stat,tempdata,extras);
		}
	}
	
	public static float[] getspectrum(String stat,Object[] data,int width,int height,boolean[] mask,int[] lims,float[] extras){
		float[] spectrum=new float[data.length];
		for(int i=0;i<data.length;i++){
			spectrum[i]=getstatistic(stat,data[i],width,height,mask,lims,extras);
		}
		return spectrum;
	}

	public static float getstatistic(String stat,Object data,int width,int height,Polygon mask,float[] extras){
		if(data instanceof float[]){
			float[] tempdata=(float[])data;
			if(mask!=null)
				tempdata=getmask((float[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof short[]){
			short[] tempdata=(short[])data;
			if(mask!=null)
				tempdata=getmask((short[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}else if(data instanceof byte[]){
			byte[] tempdata=(byte[])data;
			if(mask!=null)
				tempdata=getmask((byte[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}else{
			int[] tempdata=(int[])data;
			if(mask!=null)
				tempdata=getmask((int[])data,width,height,mask);
			return getstat(stat,tempdata,extras);
		}
	}
	
	public static float[] getspectrum(String stat,Object[] data,int width,int height,Polygon mask,float[] extras){
		float[] spectrum=new float[data.length];
		for(int i=0;i<data.length;i++){
			spectrum[i]=getstatistic(stat,data[i],width,height,mask,extras);
		}
		return spectrum;
	}

	public static float getstatistic(String stat,Object data,float[] extras){
		if(data instanceof float[]){
			float[] tempdata=(float[])data;
			return getstat(stat,tempdata,extras);
		}else if(data instanceof short[]){
			short[] tempdata=(short[])data;
			return getstat(stat,tempdata,extras);
		}else if(data instanceof byte[]){
			byte[] tempdata=(byte[])data;
			return getstat(stat,tempdata,extras);
		}else{
			int[] tempdata=(int[])data;
			return getstat(stat,tempdata,extras);
		}
	}
	
	public static float getstatistic(String stat,Object data,int start,int length,float[] extras){
		if(data instanceof float[]){
			return getstat(stat,(float[])algutils.get_subarray(data,start,length),extras);
		}else if(data instanceof short[]){
			return getstat(stat,(short[])algutils.get_subarray(data,start,length),extras);
		}else if(data instanceof byte[]){
			return getstat(stat,(byte[])algutils.get_subarray(data,start,length),extras);
		}else{
			return getstat(stat,(int[])algutils.get_subarray(data,start,length),extras);
		}
	}
	
	public static float[] getspectrum(String stat,Object[] data,float[] extras){
		float[] spectrum=new float[data.length];
		for(int i=0;i<data.length;i++){
			spectrum[i]=getstatistic(stat,data[i],extras);
		}
		return spectrum;
	}

	public static float getstat(String stat,float[] tempdata,float[] extras){
		if(tempdata==null)
			return Float.NaN;
		if(stat.equalsIgnoreCase("avg")){
			return favg(tempdata);
		}
		if(stat.equalsIgnoreCase("sum")){
			return fsum(tempdata);
		}
		if(stat.equalsIgnoreCase("max")){
			return fmax(tempdata);
		}
		if(stat.equalsIgnoreCase("min")){
			return fmin(tempdata);
		}
		if(stat.equalsIgnoreCase("variance")){
			return fvar(tempdata);
		}
		if(stat.equalsIgnoreCase("median")){
			return fmedian(tempdata);
		}
		if(stat.equalsIgnoreCase("mode")){
			return fmode(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("conditionalavg")){
			return fcondavg(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0avg")){
			return fn0avg(tempdata);
		}
		if(stat.equalsIgnoreCase("stdev")){
			return (float)Math.sqrt(fvar(tempdata));
		}
		if(stat.equalsIgnoreCase("sterr")){
			return (float)Math.sqrt(fvar(tempdata)/tempdata.length);
		}
		if(stat.equalsIgnoreCase("relerr")){
			return (float)Math.sqrt(fvar(tempdata)/tempdata.length)/favg(tempdata);
		}
		if(stat.equalsIgnoreCase("count")){
			return fcount(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("percentile")){
			return fpercentile(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0stdev")){
			return (float)Math.sqrt(fn0var(tempdata));
		}
		if(stat.equalsIgnoreCase("not0sterr")){
			return (float)Math.sqrt(fn0var(tempdata)/fn0count(tempdata));
		}
		if(stat.equalsIgnoreCase("not0count")){
			return fn0count(tempdata);
		}
		if(stat.equalsIgnoreCase("not0sum")){
			return fn0sum(tempdata);
		}
		if(stat.equalsIgnoreCase("not0min")){
			return fn0min(tempdata);
		}
		if(stat.equalsIgnoreCase("identity")){
			return 1.0f;
		}
		if(stat.equalsIgnoreCase("zero")){
			return 0.0f;
		}
		if(stat.equalsIgnoreCase("maxpos")){
			return fmaxpos(tempdata);
		}
		if(stat.equalsIgnoreCase("avgpos")){
			return favgpos(tempdata);
		}
		if(stat.equalsIgnoreCase("nan")){
			return Float.NaN;
		}
		return 0.0f;
	}

	public static float getstat(String stat,short[] tempdata,float[] extras){
		if(tempdata==null)
			return Float.NaN;
		if(stat.equalsIgnoreCase("avg")){
			return savg(tempdata);
		}
		if(stat.equalsIgnoreCase("sum")){
			return ssum(tempdata);
		}
		if(stat.equalsIgnoreCase("max")){
			return smax(tempdata);
		}
		if(stat.equalsIgnoreCase("min")){
			return smin(tempdata);
		}
		if(stat.equalsIgnoreCase("variance")){
			return svar(tempdata);
		}
		if(stat.equalsIgnoreCase("median")){
			return smedian(tempdata);
		}
		if(stat.equalsIgnoreCase("mode")){
			return smode(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("conditionalavg")){
			return scondavg(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0avg")){
			return sn0avg(tempdata);
		}
		if(stat.equalsIgnoreCase("stdev")){
			return (float)Math.sqrt(svar(tempdata));
		}
		if(stat.equalsIgnoreCase("sterr")){
			return (float)Math.sqrt(svar(tempdata)/tempdata.length);
		}
		if(stat.equalsIgnoreCase("relerr")){
			return (float)Math.sqrt(svar(tempdata)/tempdata.length)/savg(tempdata);
		}
		if(stat.equalsIgnoreCase("count")){
			return scount(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("percentile")){
			return spercentile(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0stdev")){
			return (float)Math.sqrt(sn0var(tempdata));
		}
		if(stat.equalsIgnoreCase("not0sterr")){
			return (float)Math.sqrt(sn0var(tempdata)/sn0count(tempdata));
		}
		if(stat.equalsIgnoreCase("not0count")){
			return sn0count(tempdata);
		}
		if(stat.equalsIgnoreCase("not0sum")){
			return sn0sum(tempdata);
		}
		if(stat.equalsIgnoreCase("not0min")){
			return sn0min(tempdata);
		}
		if(stat.equalsIgnoreCase("identity")){
			return 1.0f;
		}
		if(stat.equalsIgnoreCase("zero")){
			return 0.0f;
		}
		if(stat.equalsIgnoreCase("maxpos")){
			return smaxpos(tempdata);
		}
		if(stat.equalsIgnoreCase("avgpos")){
			return savgpos(tempdata);
		}
		if(stat.equalsIgnoreCase("nan")){
			return Float.NaN;
		}
		return 0.0f;
	}

	public static float getstat(String stat,int[] tempdata,float[] extras){
		if(tempdata==null)
			return Float.NaN;
		if(stat.equalsIgnoreCase("avg")){
			return iavg(tempdata);
		}
		if(stat.equalsIgnoreCase("sum")){
			return isum(tempdata);
		}
		if(stat.equalsIgnoreCase("max")){
			return imax(tempdata);
		}
		if(stat.equalsIgnoreCase("min")){
			return imin(tempdata);
		}
		if(stat.equalsIgnoreCase("variance")){
			return ivar(tempdata);
		}
		if(stat.equalsIgnoreCase("median")){
			return imedian(tempdata);
		}
		if(stat.equalsIgnoreCase("mode")){
			return imode(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("conditionalavg")){
			return icondavg(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0avg")){
			return in0avg(tempdata);
		}
		if(stat.equalsIgnoreCase("stdev")){
			return (float)Math.sqrt(ivar(tempdata));
		}
		if(stat.equalsIgnoreCase("sterr")){
			return (float)Math.sqrt(ivar(tempdata)/tempdata.length);
		}
		if(stat.equalsIgnoreCase("relerr")){
			return (float)Math.sqrt(ivar(tempdata)/tempdata.length)/iavg(tempdata);
		}
		if(stat.equalsIgnoreCase("count")){
			return icount(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("percentile")){
			return ipercentile(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0stdev")){
			return (float)Math.sqrt(in0var(tempdata));
		}
		if(stat.equalsIgnoreCase("not0sterr")){
			return (float)Math.sqrt(in0var(tempdata)/in0count(tempdata));
		}
		if(stat.equalsIgnoreCase("not0count")){
			return in0count(tempdata);
		}
		if(stat.equalsIgnoreCase("not0sum")){
			return in0sum(tempdata);
		}
		if(stat.equalsIgnoreCase("not0min")){
			return in0min(tempdata);
		}
		if(stat.equalsIgnoreCase("identity")){
			return 1.0f;
		}
		if(stat.equalsIgnoreCase("zero")){
			return 0.0f;
		}
		if(stat.equalsIgnoreCase("maxpos")){
			return imaxpos(tempdata);
		}
		if(stat.equalsIgnoreCase("avgpos")){
			return iavgpos(tempdata);
		}
		if(stat.equalsIgnoreCase("nan")){
			return Float.NaN;
		}
		return 0.0f;
	}

	public static float getstat(String stat,byte[] tempdata,float[] extras){
		if(tempdata==null)
			return Float.NaN;
		if(stat.equalsIgnoreCase("avg")){
			return bavg(tempdata);
		}
		if(stat.equalsIgnoreCase("sum")){
			return bsum(tempdata);
		}
		if(stat.equalsIgnoreCase("max")){
			return bmax(tempdata);
		}
		if(stat.equalsIgnoreCase("min")){
			return bmin(tempdata);
		}
		if(stat.equalsIgnoreCase("variance")){
			return bvar(tempdata);
		}
		if(stat.equalsIgnoreCase("median")){
			return bmedian(tempdata);
		}
		if(stat.equalsIgnoreCase("mode")){
			return bmode(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("conditionalavg")){
			return bcondavg(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0avg")){
			return bn0avg(tempdata);
		}
		if(stat.equalsIgnoreCase("stdev")){
			return (float)Math.sqrt(bvar(tempdata));
		}
		if(stat.equalsIgnoreCase("sterr")){
			return (float)Math.sqrt(bvar(tempdata)/tempdata.length);
		}
		if(stat.equalsIgnoreCase("relerr")){
			return (float)Math.sqrt(bvar(tempdata)/tempdata.length)/bavg(tempdata);
		}
		if(stat.equalsIgnoreCase("count")){
			return bcount(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("percentile")){
			return bpercentile(tempdata,extras);
		}
		if(stat.equalsIgnoreCase("not0stdev")){
			return (float)Math.sqrt(bn0var(tempdata));
		}
		if(stat.equalsIgnoreCase("not0sterr")){
			return (float)Math.sqrt(bn0var(tempdata)/bn0count(tempdata));
		}
		if(stat.equalsIgnoreCase("not0count")){
			return bn0count(tempdata);
		}
		if(stat.equalsIgnoreCase("not0sum")){
			return bn0sum(tempdata);
		}
		if(stat.equalsIgnoreCase("not0min")){
			return bn0min(tempdata);
		}
		if(stat.equalsIgnoreCase("identity")){
			return 1.0f;
		}
		if(stat.equalsIgnoreCase("zero")){
			return 0.0f;
		}
		if(stat.equalsIgnoreCase("maxpos")){
			return bmaxpos(tempdata);
		}
		if(stat.equalsIgnoreCase("avgpos")){
			return bavgpos(tempdata);
		}
		if(stat.equalsIgnoreCase("nan")){
			return Float.NaN;
		}
		return 0.0f;
	}

	public static float favg(float[] data){
		double tempavg=0.0;
		int len=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					tempavg+=(double)data[i];
					len++;
				}
			}
		}
		tempavg/=(double)len;
		return (float)tempavg;
	}

	public static float iavg(int[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			tempavg+=(double)data[i]/(double)(data.length);
		}
		return (float)tempavg;
	}

	public static float savg(short[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			tempavg+=(double)temp/(double)(data.length);
		}
		return (float)tempavg;
	}

	public static float bavg(byte[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			tempavg+=(double)temp/(double)(data.length);
		}
		return (float)tempavg;
	}

	public static float fmedian(float[] values){
		int length=values.length;
		if(length==2){
			return 0.5f*(values[0]+values[1]);
		}
		// sort the vector
		float[] tempdata=values.clone();
		Arrays.sort(tempdata);
		if((length%2.0f)==0.0f){
			return 0.5f*(tempdata[length/2]+tempdata[length/2-1]);
		}else{
			return tempdata[(int)(length/2.0f)];
		}
	}

	public static float imedian(int[] values){
		int length=values.length;
		if(length==2){
			return 0.5f*(values[0]+values[1]);
		}
		// sort the vector
		float[] tempdata=new float[length];
		for(int i=0;i<length;i++){
			tempdata[i]=(values[i]);
		}
		Arrays.sort(tempdata);
		if((length%2.0f)==0.0f){
			return 0.5f*(tempdata[length/2]+tempdata[length/2-1]);
		}else{
			return tempdata[(int)(length/2.0f)];
		}
	}

	public static float smedian(short[] values){
		int length=values.length;
		if(length==2){
			return 0.5f*(values[0]&0xffff+values[1]&0xffff);
		}
		// sort the vector
		float[] tempdata=new float[length];
		for(int i=0;i<length;i++){
			tempdata[i]=values[i]&0xffff;
		}
		Arrays.sort(tempdata);
		if((length%2.0f)==0.0f){
			return 0.5f*(tempdata[length/2]+tempdata[length/2-1]);
		}else{
			return tempdata[(int)(length/2.0f)];
		}
	}

	public static float bmedian(byte[] values){
		int length=values.length;
		if(length==2){
			return 0.5f*(values[0]&0xff+values[1]&0xff);
		}
		// sort the vector
		float[] tempdata=new float[length];
		for(int i=0;i<length;i++){
			tempdata[i]=values[i]&0xff;
		}
		Arrays.sort(tempdata);
		if((length%2.0f)==0.0f){
			return 0.5f*(tempdata[length/2]+tempdata[length/2-1]);
		}else{
			return tempdata[(int)(length/2.0f)];
		}
	}

	public static float fpercentile(float[] values,float[] percentile){
		int length=values.length;
		float[] temppercentile;
		if(percentile==null){
			temppercentile=new float[1];
			temppercentile[0]=0.95f*(length-1);
		}else{
			temppercentile=new float[percentile.length];
			for(int i=0;i<percentile.length;i++){
				temppercentile[i]=percentile[i]*0.01f*(length-1);
			}
		}
		// sort the vector
		float[] tempdata=values.clone();
		Arrays.sort(tempdata);
		float[] percentile2=new float[temppercentile.length];
		for(int i=0;i<temppercentile.length;i++){
			if(temppercentile[i]<=0.0f){
				percentile2[i]=tempdata[0];
				continue;
			}
			if(temppercentile[i]>=length-1){
				percentile2[i]=tempdata[length-1];
				continue;
			}
			int prev=(int)temppercentile[i];
			float rem=temppercentile[i]-prev;
			percentile2[i]=tempdata[prev]+rem*(tempdata[prev+1]-tempdata[prev]);
		}
		if(percentile!=null)
			for(int i=0;i<percentile.length;i++)
				percentile[i]=percentile2[i];
		return percentile2[0];
	}

	public static float ipercentile(int[] values,float[] percentile){
		int length=values.length;
		float[] temppercentile;
		if(percentile==null){
			temppercentile=new float[1];
			temppercentile[0]=0.95f*(length-1);
		}else{
			temppercentile=new float[percentile.length];
			for(int i=0;i<percentile.length;i++){
				temppercentile[i]=percentile[i]*0.01f*(length-1);
			}
		}
		// sort the vector
		float[] tempdata=new float[length];
		for(int i=0;i<length;i++)
			tempdata[i]=(values[i]);
		Arrays.sort(tempdata);
		float[] percentile2=new float[temppercentile.length];
		for(int i=0;i<temppercentile.length;i++){
			if(temppercentile[i]<=0.0f){
				percentile2[i]=tempdata[0];
			}
			if(temppercentile[i]>=length-1){
				percentile2[i]=tempdata[length-1];
			}
			int prev=(int)temppercentile[i];
			float rem=temppercentile[i]-prev;
			percentile2[i]=tempdata[prev]+rem*(tempdata[prev+1]-tempdata[prev]);
		}
		if(percentile!=null)
			for(int i=0;i<percentile.length;i++)
				percentile[i]=percentile2[i];
		return percentile2[0];
	}

	public static float spercentile(short[] values,float[] percentile){
		int length=values.length;
		float[] temppercentile;
		if(percentile==null){
			temppercentile=new float[1];
			temppercentile[0]=0.95f*(length-1);
		}else{
			temppercentile=new float[percentile.length];
			for(int i=0;i<percentile.length;i++){
				temppercentile[i]=percentile[i]*0.01f*(length-1);
			}
		}
		// sort the vector
		float[] tempdata=new float[length];
		for(int i=0;i<length;i++)
			tempdata[i]=values[i]&0xffff;
		Arrays.sort(tempdata);
		float[] percentile2=new float[temppercentile.length];
		for(int i=0;i<temppercentile.length;i++){
			if(temppercentile[i]<=0.0f){
				percentile2[i]=tempdata[0];
			}
			if(temppercentile[i]>=length-1){
				percentile2[i]=tempdata[length-1];
			}
			int prev=(int)temppercentile[i];
			float rem=temppercentile[i]-prev;
			percentile2[i]=tempdata[prev]+rem*(tempdata[prev+1]-tempdata[prev]);
		}
		if(percentile!=null)
			for(int i=0;i<percentile.length;i++)
				percentile[i]=percentile2[i];
		return percentile2[0];
	}

	public static float bpercentile(byte[] values,float[] percentile){
		int length=values.length;
		float[] temppercentile;
		if(percentile==null){
			temppercentile=new float[1];
			temppercentile[0]=0.95f*(length-1);
		}else{
			temppercentile=new float[percentile.length];
			for(int i=0;i<percentile.length;i++){
				temppercentile[i]=percentile[i]*0.01f*(length-1);
			}
		}
		// sort the vector
		float[] tempdata=new float[length];
		for(int i=0;i<length;i++)
			tempdata[i]=values[i]&0xff;
		Arrays.sort(tempdata);
		float[] percentile2=new float[temppercentile.length];
		for(int i=0;i<temppercentile.length;i++){
			if(temppercentile[i]<=0.0f){
				percentile2[i]=tempdata[0];
			}
			if(temppercentile[i]>=length-1){
				percentile2[i]=tempdata[length-1];
			}
			int prev=(int)temppercentile[i];
			float rem=temppercentile[i]-prev;
			percentile2[i]=tempdata[prev]+rem*(tempdata[prev+1]-tempdata[prev]);
		}
		if(percentile!=null)
			for(int i=0;i<percentile.length;i++)
				percentile[i]=percentile2[i];
		return percentile2[0];
	}

	public static float fmax(float[] data){
		float tempmax=0.0f;
		tempmax=data[0];
		if(Float.isNaN(tempmax)) tempmax=0.0f; //this is a pretty big assumption
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]>tempmax){
						tempmax=data[i];
					}
				}
			}
		}
		return tempmax;
	}

	public static float imax(int[] data){
		float tempmax=0.0f;
		tempmax=data[0];
		for(int i=0;i<data.length;i++){
			float temp=data[i];
			if(temp>tempmax){
				tempmax=temp;
			}
		}
		return tempmax;
	}

	public static float smax(short[] data){
		float tempmax=0.0f;
		tempmax=data[0]&0xffff;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp>tempmax){
				tempmax=temp;
			}
		}
		return tempmax;
	}

	public static float bmax(byte[] data){
		float tempmax=0.0f;
		tempmax=data[0]&0xff;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp>tempmax){
				tempmax=temp;
			}
		}
		return tempmax;
	}
	
	public static float fmaxpos(float[] data){
		float tempmax=0.0f;
		tempmax=data[0];
		int pos=0;
		if(Float.isNaN(tempmax)) tempmax=0.0f; //this is a pretty big assumption
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]>tempmax){
						tempmax=data[i];
						pos=i;
					}
				}
			}
		}
		return pos;
	}

	public static float imaxpos(int[] data){
		float tempmax=0.0f;
		tempmax=data[0];
		int pos=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i];
			if(temp>tempmax){
				tempmax=temp;
				pos=i;
			}
		}
		return pos;
	}

	public static float smaxpos(short[] data){
		float tempmax=0.0f;
		tempmax=data[0]&0xffff;
		int pos=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp>tempmax){
				tempmax=temp;
				pos=i;
			}
		}
		return pos;
	}

	public static float bmaxpos(byte[] data){
		float tempmax=0.0f;
		tempmax=data[0]&0xff;
		int pos=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp>tempmax){
				tempmax=temp;
				pos=0;
			}
		}
		return pos;
	}
	
	public static float favgpos(float[] data){
		double sum=0.0;
		double sumpos=0.0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					sumpos+=data[i]*(double)i;
					sum+=data[i];
				}
			}
		}
		return (float)(sumpos/sum);
	}
	
	public static float savgpos(short[] data){
		double sum=0.0;
		double sumpos=0.0;
		for(int i=0;i<data.length;i++){
			double temp=(double)(data[i]&0xffff);
			sumpos+=(double)temp*(double)i;
			sum+=(double)temp;
		}
		return (float)(sumpos/sum);
	}
	
	public static float bavgpos(byte[] data){
		double sum=0.0;
		double sumpos=0.0;
		for(int i=0;i<data.length;i++){
			double temp=(double)(data[i]&0xff);
			sumpos+=(double)temp*(double)i;
			sum+=(double)temp;
		}
		return (float)(sumpos/sum);
	}
	
	public static float iavgpos(int[] data){
		double sum=0.0;
		double sumpos=0.0;
		for(int i=0;i<data.length;i++){
			double temp=(double)data[i];
			sumpos+=(double)temp*(double)i;
			sum+=(double)temp;
		}
		return (float)(sumpos/sum);
	}

	public static float fmin(float[] data){
		float tempmin=0.0f;
		tempmin=data[0];
		if(Float.isNaN(tempmin)) tempmin=0.0f;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]<tempmin){
						tempmin=data[i];
					}
				}
			}
		}
		return tempmin;
	}

	public static float imin(int[] data){
		float tempmin=0.0f;
		tempmin=(data[0]);
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp<tempmin){
				tempmin=temp;
			}
		}
		return tempmin;
	}

	public static float smin(short[] data){
		float tempmin=0.0f;
		tempmin=data[0]&0xffff;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp<tempmin){
				tempmin=temp;
			}
		}
		return tempmin;
	}

	public static float bmin(byte[] data){
		float tempmin=0.0f;
		tempmin=data[0]&0xff;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp<tempmin){
				tempmin=temp;
			}
		}
		return tempmin;
	}

	public static float fsum(float[] data){
		double tempsum=0.0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					tempsum+=data[i];
				}
			}
		}
		return (float)tempsum;
	}

	public static float isum(int[] data){
		double tempsum=0.0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			tempsum+=temp;
		}
		return (float)tempsum;
	}

	public static float ssum(short[] data){
		double tempsum=0.0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			tempsum+=temp;
		}
		return (float)tempsum;
	}

	public static float bsum(byte[] data){
		double tempsum=0.0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			tempsum+=temp;
		}
		return (float)tempsum;
	}

	public static float fvar(float[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					tempavg+=data[i];
					avgsq+=data[i]*data[i];
					length++;
				}
			}
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float ivar(int[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			tempavg+=temp;
			avgsq+=temp*temp;
			length++;
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float svar(short[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			tempavg+=temp;
			avgsq+=temp*temp;
			length++;
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float bvar(byte[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			tempavg+=temp;
			avgsq+=temp*temp;
			length++;
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float fn0var(float[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]!=0.0f){
						tempavg+=data[i];
						avgsq+=data[i]*data[i];
						length++;
					}
				}
			}
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float in0var(int[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp!=0.0f){
				tempavg+=temp;
				avgsq+=temp*temp;
				length++;
			}
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float sn0var(short[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp!=0.0f){
				tempavg+=temp;
				avgsq+=temp*temp;
				length++;
			}
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float bn0var(byte[] data){
		double tempavg=0.0;
		double avgsq=0.0;
		int length=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp!=0.0f){
				tempavg+=temp;
				avgsq+=temp*temp;
				length++;
			}
		}
		tempavg/=length;
		avgsq/=length;
		avgsq-=tempavg*tempavg;
		avgsq*=((double)length/(double)(length-1));
		return (float)avgsq;
	}

	public static float fmode(float[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=fmax(data);
			histend=fmin(data);
		}
		int[] histogram=new int[histbins];
		for(int j=0;j<histbins;j++){
			histogram[j]=0;
		}
		for(int i=0;i<data.length;i++){
			int histval=(int)(((data[i]-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		int histmax=histogram[0];
		int histmaxval=0;
		for(int j=1;j<histbins;j++){
			if(histogram[j]>histmax){
				histmax=histogram[j];
				histmaxval=j;
			}
		}
		return ((histmaxval+0.5f)/histbins)*(histend-histstart)+histstart;
	}

	public static float imode(int[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=imax(data);
			histend=imin(data);
		}
		int[] histogram=new int[histbins];
		for(int j=0;j<histbins;j++){
			histogram[j]=0;
		}
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			int histval=(int)(((temp-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		int histmax=histogram[0];
		int histmaxval=0;
		for(int j=1;j<histbins;j++){
			if(histogram[j]>histmax){
				histmax=histogram[j];
				histmaxval=j;
			}
		}
		return ((histmaxval+0.5f)/histbins)*(histend-histstart)+histstart;
	}

	public static float smode(short[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=smax(data);
			histend=smin(data);
		}
		int[] histogram=new int[histbins];
		for(int j=0;j<histbins;j++){
			histogram[j]=0;
		}
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			int histval=(int)(((temp-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		int histmax=histogram[0];
		int histmaxval=0;
		for(int j=1;j<histbins;j++){
			if(histogram[j]>histmax){
				histmax=histogram[j];
				histmaxval=j;
			}
		}
		return ((histmaxval+0.5f)/histbins)*(histend-histstart)+histstart;
	}

	public static float bmode(byte[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=bmax(data);
			histend=bmin(data);
		}
		int[] histogram=new int[histbins];
		for(int j=0;j<histbins;j++){
			histogram[j]=0;
		}
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			int histval=(int)(((temp-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		int histmax=histogram[0];
		int histmaxval=0;
		for(int j=1;j<histbins;j++){
			if(histogram[j]>histmax){
				histmax=histogram[j];
				histmaxval=j;
			}
		}
		return ((histmaxval+0.5f)/histbins)*(histend-histstart)+histstart;
	}

	public static float fcondavg(float[] data,float[] extras){
		float upper,lower;
		if(extras!=null){
			upper=extras[0];
			lower=extras[1];
		}else{
			return favg(data);
		}
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]>=lower&&data[i]<=upper){
						tempavg+=data[i];
						counter++;
					}
				}
			}
		}
		return (float)(tempavg/counter);
	}

	public static float icondavg(int[] data,float[] extras){
		float upper,lower;
		if(extras!=null){
			upper=extras[0];
			lower=extras[1];
		}else{
			return iavg(data);
		}
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp>=lower&&temp<=upper){
				tempavg+=temp;
				counter++;
			}
		}
		return (float)(tempavg/counter);
	}

	public static float scondavg(short[] data,float[] extras){
		float upper,lower;
		if(extras!=null){
			upper=extras[0];
			lower=extras[1];
		}else{
			return savg(data);
		}
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp>=lower&&temp<=upper){
				tempavg+=temp;
				counter++;
			}
		}
		return (float)(tempavg/counter);
	}

	public static float bcondavg(byte[] data,float[] extras){
		float upper,lower;
		if(extras!=null){
			upper=extras[0];
			lower=extras[1];
		}else{
			return bavg(data);
		}
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp>=lower&&temp<=upper){
				tempavg+=temp;
				counter++;
			}
		}
		return (float)(tempavg/counter);
	}

	public static float fn0avg(float[] data){
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]!=0.0f){
						tempavg+=data[i];
						counter++;
					}
				}
			}
		}
		return (float)(tempavg/counter);
	}

	public static float in0avg(int[] data){
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp!=0.0f){
				tempavg+=temp;
				counter++;
			}
		}
		return (float)(tempavg/counter);
	}

	public static float sn0avg(short[] data){
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp!=0.0f){
				tempavg+=temp;
				counter++;
			}
		}
		return (float)(tempavg/counter);
	}

	public static float bn0avg(byte[] data){
		double tempavg=0.0;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp!=0.0f){
				tempavg+=temp;
				counter++;
			}
		}
		return (float)(tempavg/counter);
	}

	public static float fn0sum(float[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]!=0.0f){
						tempavg+=data[i];
					}
				}
			}
		}
		return (float)tempavg;
	}

	public static float in0sum(int[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp!=0.0f){
				tempavg+=temp;
			}
		}
		return (float)tempavg;
	}

	public static float sn0sum(short[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp!=0.0f){
				tempavg+=temp;
			}
		}
		return (float)tempavg;
	}

	public static float bn0sum(byte[] data){
		double tempavg=0.0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp!=0.0f){
				tempavg+=temp;
			}
		}
		return (float)tempavg;
	}

	public static float fcount(float[] data,float[] extras){
		if(extras==null)
			return data.length;
		int counter=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]>=extras[0]&&data[i]<=extras[1]){
						counter++;
					}
				}
			}
		}
		return counter;
	}

	public static float icount(int[] data,float[] extras){
		if(extras==null)
			return data.length;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp>=extras[0]&&temp<=extras[1]){
				counter++;
			}
		}
		return counter;
	}

	public static float scount(short[] data,float[] extras){
		if(extras==null)
			return data.length;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp>=extras[0]&&temp<=extras[1]){
				counter++;
			}
		}
		return counter;
	}

	public static float bcount(byte[] data,float[] extras){
		if(extras==null)
			return data.length;
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp>=extras[0]&&temp<=extras[1]){
				counter++;
			}
		}
		return counter;
	}

	public static float fn0count(float[] data){
		int counter=0;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]!=0.0f){
						counter++;
					}
				}
			}
		}
		return counter;
	}

	public static float in0count(int[] data){
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp!=0.0f){
				counter++;
			}
		}
		return counter;
	}

	public static float sn0count(short[] data){
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp!=0.0f){
				counter++;
			}
		}
		return counter;
	}

	public static float bn0count(byte[] data){
		int counter=0;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp!=0.0f){
				counter++;
			}
		}
		return counter;
	}

	public static float fn0min(float[] data){
		boolean first=true;
		float min=0.0f;
		for(int i=0;i<data.length;i++){
			if(!Float.isInfinite(data[i])){
				if(!Float.isNaN(data[i])){
					if(data[i]!=0.0f){
						if(first){
							min=data[i];
							first=false;
						} else {
							if(data[i]<min) min=data[i];
						}
					}
				}
			}
		}
		return min;
	}

	public static float in0min(int[] data){
		boolean first=true;
		float min=0.0f;
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			if(temp!=0.0f){
				if(first){
					min=temp;
					first=false;
				} else {
					if(temp<min) min=temp;
				}
			}
		}
		return min;
	}

	public static float sn0min(short[] data){
		boolean first=true;
		float min=0.0f;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp!=0.0f){
				if(first){
					min=temp;
					first=false;
				} else {
					if(temp<min) min=temp;
				}
			}
		}
		return min;
	}

	public static float bn0min(byte[] data){
		boolean first=true;
		float min=0.0f;
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp!=0.0f){
				if(first){
					min=temp;
					first=false;
				} else {
					if(temp<min) min=temp;
				}
			}
		}
		return min;
	}

	public static float[] getrect(float[] data,int width,int height,Rectangle r2){
		Rectangle r=r2.intersection(new Rectangle(0,0,width,height));
		float[] retarray=new float[r.width*r.height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				retarray[(i-r.y)*r.width+j-r.x]=data[i*width+j];
			}
		}
		return retarray;
	}

	public static int[] getrect(int[] data,int width,int height,Rectangle r2){
		Rectangle r=r2.intersection(new Rectangle(0,0,width,height));
		int[] retarray=new int[r.width*r.height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				retarray[(i-r.y)*r.width+j-r.x]=data[i*width+j];
			}
		}
		return retarray;
	}

	public static short[] getrect(short[] data,int width,int height,Rectangle r2){
		Rectangle r=r2.intersection(new Rectangle(0,0,width,height));
		short[] retarray=new short[r.width*r.height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				retarray[(i-r.y)*r.width+j-r.x]=data[i*width+j];
			}
		}
		return retarray;
	}

	public static byte[] getrect(byte[] data,int width,int height,Rectangle r2){
		Rectangle r=r2.intersection(new Rectangle(0,0,width,height));
		byte[] retarray=new byte[r.width*r.height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				retarray[(i-r.y)*r.width+j-r.x]=data[i*width+j];
			}
		}
		return retarray;
	}

	public static float[] getmask(float[] data,int width,int height,boolean[] mask){
		int npts=getmaskarea(mask);
		if(npts==0)
			return null;
		float[] retarray=new float[npts];
		int counter=0;
		for(int i=0;i<width*height;i++){
			if(mask[i]){
				retarray[counter]=data[i];
				counter++;
			}
		}
		return retarray;
	}
	
	public static float[] getmask(float[] data,int width,int height,boolean[] mask,int[] lims){
		int npts=getmaskarea(mask,width,height,lims);
		if(npts==0)
			return null;
		float[] retarray=new float[npts];
		int counter=0;
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
    			if(mask[j+i*width]){
    				retarray[counter]=data[j+i*width];
    				counter++;
    			}
			}
		}
		return retarray;
	}

	public static int[] getmask(int[] data,int width,int height,boolean[] mask){
		int npts=getmaskarea(mask);
		if(npts==0)
			return null;
		int[] retarray=new int[npts];
		int counter=0;
		for(int i=0;i<width*height;i++){
			if(mask[i]){
				retarray[counter]=data[i];
				counter++;
			}
		}
		return retarray;
	}
	
	public static int[] getmask(int[] data,int width,int height,boolean[] mask,int[] lims){
		int npts=getmaskarea(mask,width,height,lims);
		if(npts==0)
			return null;
		int[] retarray=new int[npts];
		int counter=0;
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
    			if(mask[j+i*width]){
    				retarray[counter]=data[j+i*width];
    				counter++;
    			}
			}
		}
		return retarray;
	}

	public static short[] getmask(short[] data,int width,int height,boolean[] mask){
		int npts=getmaskarea(mask);
		if(npts==0)
			return null;
		short[] retarray=new short[npts];
		int counter=0;
		for(int i=0;i<width*height;i++){
			if(mask[i]){
				retarray[counter]=data[i];
				counter++;
			}
		}
		return retarray;
	}
	
	/************************
	 * 
	 * @param data
	 * @param width
	 * @param height
	 * @param mask
	 * @param lims: int array with xmin,xmax,ymin,ymax inclusive
	 * @return
	 */
	public static short[] getmask(short[] data,int width,int height,boolean[] mask,int[] lims){
		int npts=getmaskarea(mask,width,height,lims);
		if(npts==0)
			return null;
		short[] retarray=new short[npts];
		int counter=0;
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
    			if(mask[j+i*width]){
    				retarray[counter]=data[j+i*width];
    				counter++;
    			}
			}
		}
		return retarray;
	}

	public static byte[] getmask(byte[] data,int width,int height,boolean[] mask){
		int npts=getmaskarea(mask);
		if(npts==0)
			return null;
		byte[] retarray=new byte[npts];
		int counter=0;
		for(int i=0;i<width*height;i++){
			if(mask[i]){
				retarray[counter]=data[i];
				counter++;
			}
		}
		return retarray;
	}
	
	public static byte[] getmask(byte[] data,int width,int height,boolean[] mask,int[] lims){
		int npts=getmaskarea(mask,width,height,lims);
		if(npts==0)
			return null;
		byte[] retarray=new byte[npts];
		int counter=0;
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
    			if(mask[j+i*width]){
    				retarray[counter]=data[j+i*width];
    				counter++;
    			}
			}
		}
		return retarray;
	}

	public static float[] getmask(float[] data,int width,int height,Polygon poly){
		Rectangle r=poly.getBounds();
		int npts=getmaskarea(poly);
		if(npts==0)
			return null;
		float[] retarray=new float[npts];
		int counter=0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					retarray[counter]=data[j+i*width];
					counter++;
				}
			}
		}
		return retarray;
	}

	public static int[] getmask(int[] data,int width,int height,Polygon poly){
		Rectangle r=poly.getBounds();
		int npts=getmaskarea(poly);
		if(npts==0)
			return null;
		int[] retarray=new int[npts];
		int counter=0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					retarray[counter]=data[j+i*width];
					counter++;
				}
			}
		}
		return retarray;
	}

	public static short[] getmask(short[] data,int width,int height,Polygon poly){
		Rectangle r=poly.getBounds();
		int npts=getmaskarea(poly);
		if(npts==0)
			return null;
		short[] retarray=new short[npts];
		int counter=0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					retarray[counter]=data[j+i*width];
					counter++;
				}
			}
		}
		return retarray;
	}

	public static byte[] getmask(byte[] data,int width,int height,Polygon poly){
		Rectangle r=poly.getBounds();
		int npts=getmaskarea(poly);
		if(npts==0)
			return null;
		byte[] retarray=new byte[npts];
		int counter=0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					retarray[counter]=data[j+i*width];
					counter++;
				}
			}
		}
		return retarray;
	}
	
	public static boolean[] poly2mask(Polygon poly,int width,int height){
		Rectangle r=poly.getBounds();
		boolean[] mask=new boolean[width*height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					mask[j+i*width]=true;
				}
			}
		}
		return mask;
	}
	
	public static float[] poly2fmask(Polygon poly,int width,int height){
		Rectangle r=poly.getBounds();
		float[] mask=new float[width*height];
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					mask[j+i*width]=1.0f;
				}
			}
		}
		return mask;
	}

	public static int getmaskarea(Polygon poly){
		Rectangle r=poly.getBounds();
		int npts=0;
		for(int i=r.y;i<(r.y+r.height);i++){
			for(int j=r.x;j<(r.x+r.width);j++){
				if(poly.contains(j,i)){
					npts++;
				}
			}
		}
		return npts;
	}

	public static int getmaskarea(boolean[] mask){
		int npts=0;
		for(int i=0;i<mask.length;i++){
			if(mask[i])
				npts++;
		}
		return npts;
	}
	
	public static int getmaskarea(boolean[] mask,int width,int height,int[] lims){
		int npts=0;
		for(int i=lims[2];i<=lims[3];i++){
			for(int j=lims[0];j<=lims[1];j++){
				if(mask[j+i*width])
					npts++;
			}
		}
		return npts;
	}

	public static int[] histogram(Object data,float[] extras){
		if(data instanceof float[])
			return fhistogram((float[])data,extras);
		if(data instanceof short[])
			return shistogram((short[])data,extras);
		if(data instanceof byte[])
			return bhistogram((byte[])data,extras);
		return null;
	}

	public static int[] fhistogram(float[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=fmin(data);
			histend=fmax(data);
		}
		int[] histogram=new int[histbins];
		for(int i=0;i<data.length;i++){
			int histval=(int)(((data[i]-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		return histogram;
	}

	public static int[] ihistogram(int[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=imin(data);
			histend=imax(data);
		}
		int[] histogram=new int[histbins];
		for(int i=0;i<data.length;i++){
			float temp=(data[i]);
			int histval=(int)(((temp-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		return histogram;
	}

	public static int[] shistogram(short[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=smin(data);
			histend=smax(data);
		}
		int[] histogram=new int[histbins];
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xffff;
			int histval=(int)(((temp-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		return histogram;
	}

	public static int[] bhistogram(byte[] data,float[] extras){
		int histbins=10;
		float histstart,histend;
		if(extras!=null){
			histbins=(int)extras[0];
			histstart=extras[1];
			histend=extras[2];
		}else{
			histstart=bmin(data);
			histend=bmax(data);
		}
		int[] histogram=new int[histbins];
		for(int i=0;i<data.length;i++){
			float temp=data[i]&0xff;
			int histval=(int)(((temp-histstart)/(histend-histstart))*histbins);
			if(histval>=0&&histval<histbins){
				histogram[histval]++;
			}
		}
		return histogram;
	}

	public static float[] histxaxis(float[] extras,int position){
		int histbins=(int)extras[0];
		float histstart=extras[1];
		float histend=extras[2];
		float[] xvals=new float[histbins];
		float binsize=(histend-histstart)/histbins;
		for(int i=0;i<histbins;i++){
			if(position==0){
				xvals[i]=histstart+i*binsize;
			}else{
				if(position==1){
					xvals[i]=histstart+i*binsize+0.5f*binsize;
				}else{
					xvals[i]=histstart+(i+1)*binsize;
				}
			}
		}
		return xvals;
	}

	public static float[] getminmax(float[] data){
		float[] minmax={data[0],data[0]};
		for(int i=1;i<data.length;i++){
			if(data[i]<minmax[0]){
				minmax[0]=data[i];
			}
			if(data[i]>minmax[1]){
				minmax[1]=data[i];
			}
		}
		return minmax;
	}

	public static float[] getminmax(int[] data){
		float[] minmax={(data[0]),(data[0])};
		for(int i=1;i<data.length;i++){
			float temp=(data[i]);
			if(temp<minmax[0]){
				minmax[0]=temp;
			}
			if(temp>minmax[1]){
				minmax[1]=temp;
			}
		}
		return minmax;
	}

	public static float[] getminmax(short[] data){
		float[] minmax={data[0]&0xffff,data[0]&0xffff};
		for(int i=1;i<data.length;i++){
			float temp=data[i]&0xffff;
			if(temp<minmax[0]){
				minmax[0]=temp;
			}
			if(temp>minmax[1]){
				minmax[1]=temp;
			}
		}
		return minmax;
	}

	public static float[] getminmax(byte[] data){
		float[] minmax={data[0]&0xff,data[0]&0xff};
		for(int i=1;i<data.length;i++){
			float temp=data[i]&0xff;
			if(temp<minmax[0]){
				minmax[0]=temp;
			}
			if(temp>minmax[1]){
				minmax[1]=temp;
			}
		}
		return minmax;
	}
}
