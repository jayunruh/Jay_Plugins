/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jguis;

// ImageJ
import ij.ImagePlus;

public class TurboRegJ_{
	// this is an adaptor class to provide access to TurboReg_ without scripting
	public static final int AFFINE=6;
	public static final int RIGID_BODY=3;
	public static final int SCALED_ROTATION=4;
	public static final int TRANSLATION=2;

	public Object tr;

	/*
	 * public TurboRegJ_(){ tr=new TurboReg_(); }
	 */

	public TurboRegJ_(Object tr){
		this.tr=tr;
	}

	public double[][] getSourcePoints(){
		// return tr.getSourcePoints();
		return (double[][])jutils.getReflectionField(tr,"sourcePoints");
	}

	public double[][] getTargetPoints(){
		// return tr.getTargetPoints();
		return (double[][])jutils.getReflectionField(tr,"targetPoints");
	}

	public void setSourcePoints(double[][] sourcePoints){
		double[][] temp=getSourcePoints();
		for(int i=0;i<sourcePoints.length;i++){
			for(int j=0;j<sourcePoints[i].length;j++){
				temp[i][j]=sourcePoints[i][j];
			}
		}
	}

	public void setTargetPoints(double[][] targetPoints){
		double[][] temp=getTargetPoints();
		for(int i=0;i<targetPoints.length;i++){
			for(int j=0;j<targetPoints[i].length;j++){
				temp[i][j]=targetPoints[i][j];
			}
		}
	}

	public ImagePlus initAlignment(final ImagePlus source,final ImagePlus target,final int transformation){
		int[] sourceCrop={0,0,source.getWidth(),source.getHeight()};
		Object[] args={source,sourceCrop,target,sourceCrop,transformation,new Boolean(false)};
		return (ImagePlus)jutils.runReflectionMethod(tr,"alignImages",args);
	}

	public ImagePlus transformImage(ImagePlus source,int width,int height,int transformation){
		Object[] args={source,width,height,transformation,new Boolean(false)};
		return (ImagePlus)jutils.runReflectionMethod(tr,"transformImage",args);
	}

}
