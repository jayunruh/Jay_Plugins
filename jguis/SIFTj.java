package jguis;

import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import mpicbg.ij.InverseTransformMapping;
import mpicbg.ij.Mapping;
import mpicbg.ij.SIFT;
import mpicbg.imagefeatures.Feature;
import mpicbg.imagefeatures.FloatArray2DSIFT;
import mpicbg.models.AbstractAffineModel2D;
import mpicbg.models.AffineModel2D;
import mpicbg.models.PointMatch;
import mpicbg.models.RigidModel2D;
import mpicbg.models.SimilarityModel2D;
import mpicbg.models.TranslationModel2D;

public class SIFTj{
	//this is my wrapper for the mpi-cbg SIFT code
	//followed the Fiji code for "SIFT_Align.java"
	public FloatArray2DSIFT.Param siftp;
	public float rod=0.92f;
	public float maxEpsilon=25.0f;
	public float minInlierRatio=0.05f;
	public static String[] modelNames={"Translation","Rigid","Similarity","Affine"};
	public int modelIndex=1;
	public boolean interpolate=true;
	private List<Feature> fs1=new ArrayList<Feature>();
	private List<Feature> fs2=new ArrayList<Feature>();

	public SIFTj(float[] params,boolean interpolate){
		//params are initialsigma, steps, minoctavesize,maxoctavesize, fdsize,fdbins,model(0-3),interpolate
		siftp=new FloatArray2DSIFT.Param();
		siftp.initialSigma=params[0];
		siftp.steps=(int)params[1];
		siftp.minOctaveSize=(int)params[2];
		siftp.maxOctaveSize=(int)params[3];
		siftp.fdSize=(int)params[4];
		siftp.fdBins=(int)params[5];
		modelIndex=(int)params[6];
		this.interpolate=interpolate;
	}
	
	public Object[] runSIFT(ImagePlus source){
		return runSIFT(source,null);
	}
	
	public Object[] runSIFT(ImagePlus source,ImagePlus simp){
		//this version aligns a secondary image with the same number of frames
		int nchans=source.getNChannels();
		int nslices=source.getNSlices();
		int nframes=source.getNFrames();
		int currchan=source.getC()-1;
		int currslice=source.getZ()-1;
		if(nframes==1){
			nframes=nslices;
			currslice=0;
			nslices=1;
		}
		int nchans2=0,nslices2=0,nframes2=0;
		ImageStack sstack=null,saligned=null;
		if(simp!=null){
			nchans2=simp.getNChannels();
			nslices2=simp.getNSlices();
			nframes2=simp.getNFrames();
    		if(nframes2==1){
    			nframes2=nslices2;
    			nslices2=1;
    		}
    		if(nframes2!=nframes){
    			IJ.error("secondary stack length doesn't match the first");
    			return null;
    		}
    		sstack=simp.getStack();
    		saligned=new ImageStack(sstack.getWidth(),sstack.getHeight());
		}
		fs1.clear();
		fs2.clear();
		ImageStack stack=source.getStack();
		ImageStack aligned=new ImageStack(stack.getWidth(),stack.getHeight());
		//ImageProcessor firstSlice=stack.getProcessor(1);
		ImageProcessor firstSlice=getProc(stack,currchan,currslice,0,nchans,nslices,nframes);
		aligned.addSlice(null,firstSlice.duplicate());
		ImageProcessor ip1;
		//ImageProcessor ip2=stack.getProcessor(1);
		ImageProcessor ip2=getProc(stack,currchan,currslice,0,nchans,nslices,nframes);
		FloatArray2DSIFT sift=new FloatArray2DSIFT(siftp);
		SIFT ijSIFT=new SIFT(sift);
		ijSIFT.extractFeatures(ip2,fs2);
		AbstractAffineModel2D model;
		if(modelIndex==0) model=new TranslationModel2D();
		else if(modelIndex==1) model=new RigidModel2D();
		else if(modelIndex==2) model=new SimilarityModel2D();
		else model=new AffineModel2D();
		Mapping mapping=new InverseTransformMapping<AbstractAffineModel2D<?>>(model);
		//String[] outtrans=new String[stack.getSize()];
		//outtrans[0]=model.createAffine().toString();
		double[][] outtrans=new double[nframes][6];
		model.createAffine().getMatrix(outtrans[0]);
		float[][] transtraj=new float[2][stack.getSize()];
		for(int i=1;i<nframes;i++){
			ip1=ip2;
			//ip2=stack.getProcessor(i+1);
			ip2=getProc(stack,currchan,currslice,i,nchans,nslices,nframes);
			fs1.clear();
			fs1.addAll(fs2);
			fs2.clear();
			ijSIFT.extractFeatures(ip2,fs2);
			Vector<PointMatch> candidates=FloatArray2DSIFT.createMatches(fs2,fs1,1.5f,null,Float.MAX_VALUE,rod);
			Vector<PointMatch> inliers=new Vector<PointMatch>();
			AbstractAffineModel2D<?> currentModel;
			if(modelIndex==0) currentModel=new TranslationModel2D();
			else if(modelIndex==1) currentModel=new RigidModel2D();
			else if(modelIndex==2) currentModel=new SimilarityModel2D();
			else currentModel=new AffineModel2D();
			boolean modelFound;
			try{
				modelFound=currentModel.filterRansac(candidates,inliers,1000,maxEpsilon,minInlierRatio);
			}catch(Exception e){
				modelFound=false;
				System.err.println(e.getMessage());
				IJ.log("model not found for slice "+i);
				IJ.log(e.getMessage());
			}
			if(modelFound){
				model.concatenate(currentModel);
				AffineTransform acmodel=currentModel.createAffine();
				//outtrans[i]=acmodel.toString();
				acmodel.getMatrix(outtrans[i]);
				//currentModel.toArray(new double[6]);
				//transtraj[0][i]=(float)acmodel.getTranslateX();
				//transtraj[1][i]=(float)acmodel.getTranslateY();
				try{
					double[] translation=model.applyInverse(new double[]{0.0,0.0});
					transtraj[0][i]=(float)translation[0];
					transtraj[1][i]=(float)translation[1];
				}catch(Exception e){}
			}
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
    				//ImageProcessor originalSlice=stack.getProcessor(i+1);
    				ImageProcessor originalSlice=getProc(stack,k,j,i,nchans,nslices,nframes);
    				originalSlice.setInterpolationMethod(ImageProcessor.BILINEAR);
    				ImageProcessor alignedSlice=originalSlice.createProcessor(stack.getWidth(),stack.getHeight());
    				if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
    				else mapping.map(originalSlice,alignedSlice);
    				aligned.addSlice(null,alignedSlice);
				}
			}
			if(simp!=null){
				for(int j=0;j<nslices2;j++){
					for(int k=0;k<nchans2;k++){
	    				//ImageProcessor originalSlice=stack.getProcessor(i+1);
	    				ImageProcessor originalSlice=getProc(sstack,k,j,i,nchans2,nslices2,nframes2);
	    				originalSlice.setInterpolationMethod(ImageProcessor.BILINEAR);
	    				ImageProcessor alignedSlice=originalSlice.createProcessor(sstack.getWidth(),sstack.getHeight());
	    				if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
	    				else mapping.map(originalSlice,alignedSlice);
	    				saligned.addSlice(null,alignedSlice);
					}
				}
			}
			IJ.showProgress(i,nframes);
		}
		ImagePlus aimp=jutils.create_hyperstack("SIFT Aligned",aligned,source);
		ImagePlus saimp=null;
		if(simp!=null){
			saimp=jutils.create_hyperstack("SIFT Aligned Secondary",saligned,source);
		}
		return new Object[]{aimp,outtrans,transtraj,saimp};
	}
	
	public Object[] runSIFT(ImagePlus source,ImagePlus simp,int nthreads){
		//this version aligns a secondary image with the same number of frames
		//this version uses multithreaded transformation
		int nchans=source.getNChannels();
		int nslices=source.getNSlices();
		int nframes=source.getNFrames();
		int currchan=source.getC()-1;
		int currslice=source.getZ()-1;
		if(nframes==1){
			nframes=nslices;
			currslice=0;
			nslices=1;
		}
		int nchans2=0,nslices2=0,nframes2=0;
		ImageStack sstack=null,saligned=null;
		if(simp!=null){
			nchans2=simp.getNChannels();
			nslices2=simp.getNSlices();
			nframes2=simp.getNFrames();
    		if(nframes2==1){
    			nframes2=nslices2;
    			nslices2=1;
    		}
    		if(nframes2!=nframes){
    			IJ.error("secondary stack length doesn't match the first");
    			return null;
    		}
    		sstack=simp.getStack();
    		saligned=new ImageStack(sstack.getWidth(),sstack.getHeight());
		}
		fs1.clear();
		fs2.clear();
		ImageStack stack=source.getStack();
		ImageStack aligned=new ImageStack(stack.getWidth(),stack.getHeight());
		//ImageProcessor firstSlice=stack.getProcessor(1);
		ImageProcessor firstSlice=getProc(stack,currchan,currslice,0,nchans,nslices,nframes);
		aligned.addSlice(null,firstSlice.duplicate());
		ImageProcessor ip1;
		//ImageProcessor ip2=stack.getProcessor(1);
		ImageProcessor ip2=getProc(stack,currchan,currslice,0,nchans,nslices,nframes);
		FloatArray2DSIFT sift=new FloatArray2DSIFT(siftp);
		SIFT ijSIFT=new SIFT(sift);
		ijSIFT.extractFeatures(ip2,fs2);
		AbstractAffineModel2D model;
		if(modelIndex==0) model=new TranslationModel2D();
		else if(modelIndex==1) model=new RigidModel2D();
		else if(modelIndex==2) model=new SimilarityModel2D();
		else model=new AffineModel2D();
		Mapping mapping=new InverseTransformMapping<AbstractAffineModel2D<?>>(model);
		//String[] outtrans=new String[stack.getSize()];
		//outtrans[0]=model.createAffine().toString();
		double[][] outtrans=new double[nframes][6];
		model.createAffine().getMatrix(outtrans[0]);
		float[][] transtraj=new float[2][stack.getSize()];
		ExecutorService executor=Executors.newFixedThreadPool(nthreads);
		for(int i=1;i<nframes;i++){
			ip1=ip2;
			//ip2=stack.getProcessor(i+1);
			ip2=getProc(stack,currchan,currslice,i,nchans,nslices,nframes);
			fs1.clear();
			fs1.addAll(fs2);
			fs2.clear();
			ijSIFT.extractFeatures(ip2,fs2);
			Vector<PointMatch> candidates=FloatArray2DSIFT.createMatches(fs2,fs1,1.5f,null,Float.MAX_VALUE,rod);
			Vector<PointMatch> inliers=new Vector<PointMatch>();
			AbstractAffineModel2D<?> currentModel;
			if(modelIndex==0) currentModel=new TranslationModel2D();
			else if(modelIndex==1) currentModel=new RigidModel2D();
			else if(modelIndex==2) currentModel=new SimilarityModel2D();
			else currentModel=new AffineModel2D();
			boolean modelFound;
			try{
				modelFound=currentModel.filterRansac(candidates,inliers,1000,maxEpsilon,minInlierRatio);
			}catch(Exception e){
				modelFound=false;
				System.err.println(e.getMessage());
				IJ.log("model not found for slice "+i);
				IJ.log(e.getMessage());
			}
			if(modelFound){
				model.concatenate(currentModel);
				AffineTransform acmodel=currentModel.createAffine();
				//outtrans[i]=acmodel.toString();
				acmodel.getMatrix(outtrans[i]);
				//transtraj[0][i]=(float)acmodel.getTranslateX();
				//transtraj[1][i]=(float)acmodel.getTranslateY();
				try{
					double[] translation=model.applyInverse(new double[]{0.0,0.0});
					transtraj[0][i]=(float)translation[0];
					transtraj[1][i]=(float)translation[1];
				}catch(Exception e){}
			}
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
    				//ImageProcessor originalSlice=stack.getProcessor(i+1);
    				ImageProcessor originalSlice=getProc(stack,k,j,i,nchans,nslices,nframes);
    				originalSlice.setInterpolationMethod(ImageProcessor.BILINEAR);
    				ImageProcessor alignedSlice=originalSlice.createProcessor(stack.getWidth(),stack.getHeight());
    				Runnable worker=new SIFTTransformRunnable((AbstractAffineModel2D)model.copy(),originalSlice,alignedSlice,interpolate);
    				executor.execute(worker);
    				//if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
    				//else mapping.map(originalSlice,alignedSlice);
    				aligned.addSlice(null,alignedSlice);
				}
			}
			if(simp!=null){
				for(int j=0;j<nslices2;j++){
					for(int k=0;k<nchans2;k++){
	    				//ImageProcessor originalSlice=stack.getProcessor(i+1);
	    				ImageProcessor originalSlice=getProc(sstack,k,j,i,nchans2,nslices2,nframes2);
	    				originalSlice.setInterpolationMethod(ImageProcessor.BILINEAR);
	    				ImageProcessor alignedSlice=originalSlice.createProcessor(sstack.getWidth(),sstack.getHeight());
	    				Runnable worker=new SIFTTransformRunnable((AbstractAffineModel2D)model.copy(),originalSlice,alignedSlice,interpolate);
	    				executor.execute(worker);
	    				//if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
	    				//else mapping.map(originalSlice,alignedSlice);
	    				saligned.addSlice(null,alignedSlice);
					}
				}
			}
			IJ.showProgress(i,nframes);
		}
		executor.shutdown();
		while(!executor.isTerminated()){}
		ImagePlus aimp=jutils.create_hyperstack("SIFT Aligned",aligned,source);
		ImagePlus saimp=null;
		if(simp!=null){
			saimp=jutils.create_hyperstack("SIFT Aligned Secondary",saligned,source);
		}
		return new Object[]{aimp,outtrans,transtraj,saimp};
	}
	
	public static int getPos(int channel,int slice,int frame,int nchans,int nslices,int nframes){
		//maps a channel slice and frame to stack position
		return 1+channel+slice*nchans+frame*nslices*nchans;
	}
	
	public static ImageProcessor getProc(ImageStack stack,int channel,int slice,int frame,int nchans,int nslices,int nframes){
		return stack.getProcessor(getPos(channel,slice,frame,nchans,nslices,nframes));
	}
	
	public static ImagePlus doTransformations(ImagePlus source,float[][][] trans,boolean interpolate,boolean global){
		//here we take a stack of affine matrix values and apply them to an image
		//this version aligns a secondary image with the same number of frames
		int nchans=source.getNChannels();
		int nslices=source.getNSlices();
		int nframes=source.getNFrames();
		int currchan=source.getC()-1;
		int currslice=source.getZ()-1;
		if(nframes==1){
			nframes=nslices;
			currslice=0;
			nslices=1;
		}
		ImageStack stack=source.getStack();
		ImageStack aligned=new ImageStack(stack.getWidth(),stack.getHeight());
		//ImageProcessor firstSlice=stack.getProcessor(1);
		ImageProcessor firstSlice=getProc(stack,currchan,currslice,0,nchans,nslices,nframes);
		aligned.addSlice(null,firstSlice.duplicate());
		ImageProcessor ip1;
		//ImageProcessor ip2=stack.getProcessor(1);
		ImageProcessor ip2=getProc(stack,currchan,currslice,0,nchans,nslices,nframes);
		//FloatArray2DSIFT sift=new FloatArray2DSIFT(siftp);
		//SIFT ijSIFT=new SIFT(sift);
		//ijSIFT.extractFeatures(ip2,fs2);
		AbstractAffineModel2D model=new AffineModel2D();
		Mapping mapping=new InverseTransformMapping<AbstractAffineModel2D<?>>(model);
		if(global){
			AffineModel2D currentModel=new AffineModel2D();
			final double[] t={(double)trans[0][0][0],(double)trans[0][1][0],(double)trans[0][0][1],(double)trans[0][1][1],(double)trans[0][0][2],(double)trans[0][1][2]};
			currentModel.set(t[0],t[1],t[2],t[3],t[4],t[5]);
			model.set(currentModel);
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
    				//ImageProcessor originalSlice=stack.getProcessor(i+1);
    				ImageProcessor originalSlice=getProc(stack,k,j,0,nchans,nslices,nframes);
    				originalSlice.setInterpolationMethod(ImageProcessor.BILINEAR);
    				ImageProcessor alignedSlice=originalSlice.createProcessor(stack.getWidth(),stack.getHeight());
    				if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
    				else mapping.map(originalSlice,alignedSlice);
    				aligned.addSlice(null,alignedSlice);
				}
			}
		}
		for(int i=1;i<nframes;i++){
			ip1=ip2;
			//ip2=stack.getProcessor(i+1);
			ip2=getProc(stack,currchan,currslice,i,nchans,nslices,nframes);
			AffineModel2D currentModel=new AffineModel2D();
			final double[] t={(double)trans[i][0][0],(double)trans[i][1][0],(double)trans[i][0][1],(double)trans[i][1][1],(double)trans[i][0][2],(double)trans[i][1][2]};
			currentModel.set(t[0],t[1],t[2],t[3],t[4],t[5]);
			if(!global) model.concatenate(currentModel);
			else model.set(currentModel);
			for(int j=0;j<nslices;j++){
				for(int k=0;k<nchans;k++){
    				//ImageProcessor originalSlice=stack.getProcessor(i+1);
    				ImageProcessor originalSlice=getProc(stack,k,j,i,nchans,nslices,nframes);
    				originalSlice.setInterpolationMethod(ImageProcessor.BILINEAR);
    				ImageProcessor alignedSlice=originalSlice.createProcessor(stack.getWidth(),stack.getHeight());
    				if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
    				else mapping.map(originalSlice,alignedSlice);
    				aligned.addSlice(null,alignedSlice);
				}
			}
			IJ.showProgress(i,nframes);
		}
		ImagePlus aimp=jutils.create_hyperstack("Transformation Aligned",aligned,source);
		return aimp;
	}
}

class SIFTTransformRunnable implements Runnable{
	AbstractAffineModel2D model;
	ImageProcessor originalSlice,alignedSlice;
	boolean interpolate;
	
	public SIFTTransformRunnable(AbstractAffineModel2D model,ImageProcessor originalSlice,ImageProcessor alignedSlice,boolean interpolate){
		this.model=model;
		this.originalSlice=originalSlice;
		this.alignedSlice=alignedSlice;
		this.interpolate=interpolate;
	}
	
	public void run(){
		Mapping mapping=new InverseTransformMapping<AbstractAffineModel2D<?>>(model);
		if(interpolate) mapping.mapInterpolated(originalSlice,alignedSlice);
		else mapping.map(originalSlice,alignedSlice);
	}
	
}
