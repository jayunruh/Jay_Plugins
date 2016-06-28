package jguis;

import jalgs.algutils;
import jalgs.jstatistics;
import jalgs.stitching;

import java.awt.Rectangle;

import mpicbg.stitching.PairWiseStitchingImgLib;
import mpicbg.stitching.PairWiseStitchingResult;
import mpicbg.stitching.StitchingParameters;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;

public class test_stitching{
	
	public static float[][] stitch2DSnakeCollection(ImageStack input,float overlap,int scandir,int ximgs,int yimgs,float ftoler){
		//scan directions are 0:upleft,1:upright,2:lowleft,3:lowright
		int xstart=0; int ystart=yimgs-1; int xscandir=1; int yscandir=1;
		if(scandir==0){ystart=0; yscandir=1;}
		if(scandir==1){xstart=ximgs-1; xscandir=-1; ystart=0; yscandir=1;}
		if(scandir==3){xstart=ximgs-1; xscandir=-1;}
		float[][] coords=stitching.getTileCoords(ximgs,yimgs,input.getWidth(),input.getHeight(),overlap,xstart,ystart,xscandir,yscandir,1.0f);
		int[][] pairs=stitching.getPairs(ximgs,yimgs,xscandir,yscandir,false);
		return stitch2DCollection(input,coords,pairs,ftoler);
	}
	
	public static float[][] stitch2DRasterCollection(ImageStack input,float overlap,int scandir,int ximgs,int yimgs,float ftoler){
		//scan directions are 0:upleft,1:upright,2:lowleft,3:lowright
		int xstart=0; int ystart=yimgs-1; int xscandir=1; int yscandir=1;
		if(scandir==0){ystart=0; yscandir=1;}
		if(scandir==1){xstart=ximgs-1; xscandir=-1; ystart=0; yscandir=1;}
		if(scandir==3){xstart=ximgs-1; xscandir=-1;}
		float[][] coords=stitching.getTileCoordsRaster(ximgs,yimgs,input.getWidth(),input.getHeight(),overlap,xstart,ystart,xscandir,yscandir,1.0f);
		int[][] pairs=stitching.getPairsRaster(ximgs,yimgs,xscandir,yscandir,false);
		return stitch2DCollection(input,coords,pairs,ftoler);
	}
	
	public static float[][] stitch2DCollection(ImageStack input, float[][] coords, int[][] pairs,float ftoler){
		int width=input.getWidth(); int height=input.getHeight();
		int slices=input.getSize();
		//offsets values are paira,pairb,origxshift,origyshift,corrxshift,corryshift,phasecorr,crosscorr,finalxshift,finalyshift,alignedx,alignedy
		float[][] offsets=new float[pairs.length][12];
		StitchingParameters params=new StitchingParameters();
		params.subpixelAccuracy=true;
		params.computeOverlap=true;
		params.checkPeaks=5;
		params.dimensionality=2;
		params.fusionMethod=0; //?
		params.addTilesAsRois=false;
		params.virtual=false;
		params.invertX=false;
		params.invertY=false;
		//params.ignoreZStage=true; //ignore z stage position
		params.channel1=0; params.channel2=0;
		params.timeSelect=0;
		//params.regThreshold=0.3;
		//get the pairwise cross-correlation values
		for(int i=0;i<pairs.length;i++){
			int pair1=pairs[i][0]; int pair2=pairs[i][1];
    		Roi roi1=getRoi(width,height,coords[0][pair2]-coords[0][pair1],coords[1][pair2]-coords[1][pair1]);
    		Roi roi2=getRoi(width,height,coords[0][pair1]-coords[0][pair2],coords[1][pair1]-coords[1][pair2]);
    		float predxshift=coords[0][pair2]-coords[0][pair1];
    		float predyshift=coords[1][pair2]-coords[1][pair1];
    		ImagePlus imp1=new ImagePlus("temp1",input.getProcessor(pair1+1));
    		ImagePlus imp2=new ImagePlus("temp2",input.getProcessor(pair2+1));
    		PairWiseStitchingResult result=PairWiseStitchingImgLib.stitchPairwise(imp1,imp2,roi1,roi2,0,0,params);
    		float[] offset=result.getOffset();
    		offsets[i][0]=(float)pair1; offsets[i][1]=(float)pair2;
    		offsets[i][2]=predxshift; offsets[i][3]=predyshift;
    		offsets[i][4]=offset[0]; offsets[i][5]=offset[1]; 
    		offsets[i][6]=result.getPhaseCorrelation(); offsets[i][7]=result.getCrossCorrelation(); 
			IJ.log("image"+pair1+ " <- " + "image"+pair2 + ": " + 
					result.getOffset()[0]+ " , "+result.getOffset()[1]+" , " + " correlation (R)=" + result.getCrossCorrelation());
    	}
		//convert these values into comparepair objects
		float[][] pairStats=new float[pairs.length][3];
		//Vector<ComparePair> pairs2=new Vector<ComparePair>();
		for(int i=0;i<pairs.length;i++){
			int pair1=pairs[i][0]; int pair2=pairs[i][1];
			//ImagePlus imp11=new ImagePlus("temp1",input.getProcessor(pair1+1));
			//ImagePlus imp21=new ImagePlus("temp2",input.getProcessor(pair2+1));
			//ImageCollectionElement ice1=new ImageCollectionElement(null,pair1);
			//ice1.setImagePlus(imp11);
			//ImageCollectionElement ice2=new ImageCollectionElement(null,pair2);
			//ice1.setImagePlus(imp21);*/
			//ImagePlusTimePoint imp1=new ImagePlusTimePoint(imp11, pair1, 1, new TranslationModel2D(), null);
			//ImagePlusTimePoint imp2=new ImagePlusTimePoint(imp21, pair2, 1, new TranslationModel2D(), null);
			//ComparePair pair=new ComparePair(imp1,imp2);
			pairStats[i][0]=offsets[i][4];
			pairStats[i][1]=offsets[i][5];
			pairStats[i][2]=offsets[i][7];
			//pair.setRelativeShift(new float[]{offsets[i][4],offsets[i][5]});
			//pair.setCrossCorrelation(offsets[i][7]);
			//pairs2.add(pair);
		}
		float[] medvals=checkCoords(pairStats,coords,pairs,width,height,ftoler);
		float[][] outcoords=preAlignJay(pairStats,0,pairs,slices,(float)params.regThreshold);
		for(int i=0;i<pairs.length;i++){
			//add the checked coordinates and the final coordinates
			//offsets[i][8]=pairs2.get(i).getRelativeShift()[0];
			//offsets[i][9]=pairs2.get(i).getRelativeShift()[1];
			offsets[i][8]=pairStats[i][0];
			offsets[i][9]=pairStats[i][1];
			offsets[i][10]=outcoords[i][0];
			offsets[i][11]=outcoords[i][1];
		}
		/*float[][] optimized=globOpt(pairs2,pairs2.get(0).getTile1(),params);
		for(int i=0;i<pairs.length;i++){
			offsets[i][8]=optimized[i][0];
			offsets[i][9]=optimized[i][1];
		}*/
		/*List<ImagePlusTimePoint> optimized=GlobalOptimization.optimize(pairs2,pairs2.get(0).getTile1(),params);
		if(optimized==null) return null;
		float[][] positions=new float[optimized.size()][2];
		int[] ids=new int[optimized.size()];
		for(int i=0;i<optimized.size();i++){
			//optimized.get(index)
			TranslationModel2D m=(TranslationModel2D)optimized.get(i).getModel();
			double[] temp=new double[2];
			m.applyInPlace(temp);
			positions[i][0]=(float)temp[0];
			positions[i][1]=(float)temp[1];
			ids[i]=optimized.get(i).getImpId();
			IJ.log(""+ids[i]+" , "+positions[i][0]+" , "+positions[i][1]);
		}
		float[][] positions2=new float[optimized.size()][2];
		for(int i=0;i<positions2.length;i++) positions2[ids[i]]=positions[i];
		for(int i=0;i<pairs.length;i++){
			int pair1=pairs[i][0]; int pair2=pairs[i][1];
			offsets[i][8]=positions[pair2][0]-positions[pair1][0];
			offsets[i][9]=positions[pair2][1]-positions[pair1][1];
		}*/
		return offsets;
	}
	
	public static Roi getRoi(int width,int height,float xshift,float yshift){
		//this returns an roi surrounding the overlapped region on the reference image
		//assume both images have the same size
		int x=0; int y=0; int w=0; int h=0;
		int ixshift=(int)Math.round(xshift);
		int iyshift=(int)Math.round(yshift);
		if(ixshift>=0 && iyshift>=0){
			x=ixshift; y=iyshift; w=width-ixshift; h=height-iyshift;
		}
		if(ixshift>=0 && iyshift<0){
			x=ixshift; y=0; w=width-ixshift; h=height+iyshift;
		}
		if(ixshift<0 && iyshift<0){
			x=0; y=0; w=width+ixshift; h=height+iyshift;
		}
		if(ixshift<0 && iyshift>=0){
			x=0; y=iyshift; w=width+ixshift; h=height-iyshift;
		}
		return new Roi(new Rectangle(x,y,w,h));
	}
	
/*	public static float[][] globOpt(Vector<ComparePair> pairs,ImagePlusTimePoint fixedImage,StitchingParameters params){
		TileConfigurationStitching tc=null;
		boolean redo;
		do{
			redo=false;
    		//assume that all pairs are good (original filtered roughly on cross correlation)
    		ArrayList<Tile<?>> tiles=new ArrayList<Tile<?>>();
    		for(int i=0;i<pairs.size();i++){
    			ComparePair pair=pairs.get(i);
    			if(pair.getCrossCorrelation()>=params.regThreshold && pair.getIsValidOverlap()){
        			Tile t1=pair.getTile1();
        			Tile t2=pair.getTile2();
        			Point p1,p2;
        			p1=new Point(new double[]{0.0,0.0});
        			p2=new Point(new double[]{-pair.getRelativeShift()[0],-pair.getRelativeShift()[1]});
        			t1.addMatch(new PointMatchStitching(p1,p2,pair.getCrossCorrelation(),pair));
        			t2.addMatch(new PointMatchStitching(p2,p1,pair.getCrossCorrelation(),pair));
        			t1.addConnectedTile(t2);
        			t2.addConnectedTile(t1);
        			if(!tiles.contains(t1)) tiles.add(t1);
        			if(!tiles.contains(t2)) tiles.add(t2);
        			pair.setIsValidOverlap(true); //redundant
    			} else {
    				pair.setIsValidOverlap(false);
    			}
    		}
    		if(tiles.size()==0){
    			IJ.log("no connected tiles found");
    			return null;
    		}
    		tc=new TileConfigurationStitching();
    		tc.addTiles(tiles);
    		//find a connected tile to fix
    		if(fixedImage.getConnectedTiles().size()>0) tc.fixTile(fixedImage);
    		else{
    			for(int i=0;i<tiles.size();++i){
    				if(tiles.get(i).getConnectedTiles().size()>0){
    					tc.fixTile(tiles.get(i));
    					break;
    				}
    			}
    		}
    		try{
    			IJ.log("prealigning");
    			tc.preAlign();
    			IJ.log("optimizing");
    			tc.optimize(10,1000,200);
    			double avgError=tc.getError();
    			double maxError=tc.getMaxError();
    			if (((avgError*params.relativeThreshold<maxError && maxError>0.95) || avgError>params.absoluteThreshold)){
    				//if error above threshold (but not too high), throw out the worst tile and try again 
    				double longestDisplacement=0;
    				PointMatch worstMatch=null;
    				for(Tile t:tc.getTiles()){
    					for(PointMatch p:(Set<PointMatch>)t.getMatches()){
    						if(p.getDistance()>longestDisplacement){
    							longestDisplacement=p.getDistance();
    							worstMatch=p;
    						}
    					}
    				}
    				ComparePair pair=((PointMatchStitching)worstMatch).getPair();
    				((PointMatchStitching)worstMatch).getPair().setIsValidOverlap(false);
					redo = true;
					IJ.log("removing worst tile and trying again");
    			}
    		} catch(Exception e){
    			IJ.error("error in optimization");
    		}
		} while(redo);
		float[][] temp=new float[pairs.size()][2];
		for(int i=0;i<pairs.size();i++){
			temp[i]=pairs.get(i).getRelativeShift();
		}
		return temp;
	}*/
	
	public static float[] checkCoords(float[][] pairStats,float[][] guesscoords,int[][] pairs,int width,int height,float ftoler){
		//this subroutine uses the median displacements for vertical and horizontal pairs
		//if values are more than ftoler*dimension different than the median, they are converted to the median
		int halfwidth=width/2;
		int halfheight=height/2;
		float[][] vertvals=new float[2][pairs.length];
		int numvert=0;
		float[][] horvals=new float[2][pairs.length];
		int numhor=0;
		float[][] diagvals=new float[2][pairs.length];
		int numdiag=0;
		int[][] pairtype=new int[pairs.length][2];
		for(int i=0;i<pairs.length;i++){
			//ComparePair pair=pairs2.get(i);
			float[] pair=pairStats[i];
			int pair1=pairs[i][0]; int pair2=pairs[i][1];
			float xdiff=guesscoords[pair2][0]-guesscoords[pair1][0];
			float ydiff=guesscoords[pair2][1]-guesscoords[pair1][1];
			float axdiff=Math.abs(xdiff);
			float aydiff=Math.abs(ydiff);
			if(axdiff>halfwidth && aydiff<halfheight){ //horizontal pairs
				pairtype[i][0]=0; pairtype[i][1]=(ydiff>0)?1:-1; //forward or backwards pair
				float mult=(float)pairtype[i][1];
				horvals[0][numhor]=mult*pair[0]; horvals[1][numhor]=mult*pair[1];
				numhor++; 
			}
			if(axdiff<halfwidth && aydiff>halfheight){ //vertical pairs
				pairtype[i][0]=1; pairtype[i][1]=(xdiff>0)?1:-1; //up or down pair
				float mult=(float)pairtype[i][1];
				vertvals[0][numvert]=mult*pair[0]; vertvals[1][numvert]=mult*pair[1]; 
				numvert++;
			}
			if(axdiff>halfwidth && aydiff>halfheight){ //diagonal pairs
				pairtype[i][0]=2; pairtype[i][1]=(xdiff>0)?1:-1; //down-right or up-left pairs
				float mult=(float)pairtype[i][1];
				diagvals[0][numdiag]=mult*pair[0]; diagvals[1][numdiag]=mult*pair[1]; 
				numdiag++; 
			}
		}
		//get the median statistics
		//assume we have enough horizontal and vertical pairs
		float vmedxdiff=jstatistics.getstatistic("median",algutils.get_subarray(vertvals[0],0,numvert),null);
		float vmedydiff=jstatistics.getstatistic("median",algutils.get_subarray(vertvals[1],0,numvert),null);
		float hmedxdiff=jstatistics.getstatistic("median",algutils.get_subarray(horvals[0],0,numhor),null);
		float hmedydiff=jstatistics.getstatistic("median",algutils.get_subarray(horvals[1],0,numhor),null);
		IJ.log("horizontal median shift");
		IJ.log(""+hmedxdiff+" , "+hmedydiff);
		IJ.log("vertical median shift");
		IJ.log(""+vmedxdiff+" , "+vmedydiff);
		float dmedxdiff=0.0f; float dmedydiff=0.0f;
		if(numdiag>1){ //we might not have diagonal pairs
			dmedxdiff=jstatistics.getstatistic("median",algutils.get_subarray(diagvals[0],0,numdiag),null);
			dmedydiff=jstatistics.getstatistic("median",algutils.get_subarray(diagvals[1],0,numdiag),null);
			IJ.log("diagonal median shift");
			IJ.log(""+dmedxdiff+" , "+dmedydiff);
		}
		float htoler=ftoler*(float)width; //pixel tolerance for difference from median
		float vtoler=ftoler*(float)height;
		int hcounter=0; int vcounter=0; int dcounter=0;
		for(int i=0;i<pairs.length;i++){
			//ComparePair pair=pairs2.get(i);
			float[] pair=pairStats[i];
			if(pairtype[i][0]==0){ //horizontal
				if(Math.abs(horvals[0][hcounter]-hmedxdiff)>htoler){
					float mult=(float)pairtype[i][1];
					pair[0]=mult*hmedxdiff;
					pair[1]=mult*hmedydiff;
					pair[2]=0.0f;
				}
				hcounter++;
			}
			if(pairtype[i][0]==1){ //vertical
				if(Math.abs(vertvals[0][vcounter]-vmedydiff)>vtoler){
					float mult=(float)pairtype[i][1];
					pair[0]=mult*vmedxdiff;
					pair[1]=mult*vmedydiff;
					pair[2]=0.0f;
				}
				vcounter++;
			}
			if(pairtype[i][0]==2){ //diagonal
				if(Math.abs(diagvals[0][dcounter]-dmedxdiff)>htoler){
					float mult=(float)pairtype[i][1];
					pair[0]=mult*dmedxdiff;
					pair[1]=mult*dmedydiff;
					pair[2]=0.0f;
				}
				dcounter++;
			}
		}
		return new float[]{hmedxdiff,hmedydiff,vmedxdiff,vmedydiff,dmedxdiff,dmedydiff};
	}
	
	public static float[][] preAlignJay(float[][] pairStats,int fixedindex,int[][] pairs,int nimages, float corrthresh){
		//start with the fixed pair and propagate backwards and forwards until all are aligned
		int naligned=1;
		int nunaligned=nimages-1;
		float[][] outcoords=new float[nimages][];
		outcoords[fixedindex]=new float[]{0.0f,0.0f};
		int prevnunaligned=nunaligned;
		while(nunaligned>0){
			for(int i=0;i<outcoords.length;i++){
				if(outcoords[i]!=null){
					//we found an aligned one, see if it has neighbors with valid phase correlation
					//if it does, add those to the aligned pool
					int[][] targets=findConnectedImages(pairs,i,nimages);
					if(targets==null) continue;
					for(int j=0;j<targets.length;j++){
						if(outcoords[targets[j][0]]==null){
    						//ComparePair pair=pairs2.get(targets[j][1]);
    						float[] pair=pairStats[targets[j][i]];
    						int targetmember=pairs[targets[j][1]][1];
    						float[] targetshift=new float[]{pair[0],pair[1]};
    						if(targetmember==i){//the query is the second member of the pair
    							targetmember=pairs[targets[j][1]][0];
    							targetshift[0]*=-1.0f;
    							targetshift[1]*=-1.0f;
    						}
    						if(pair[2]>corrthresh){
    							float[] querycoords=outcoords[i];
    							float[] newcoords={querycoords[0]+targetshift[0],querycoords[1]+targetshift[1]};
    							outcoords[targets[j][0]]=newcoords;
    							naligned++; nunaligned--;
    						}
						}
					}
				}
			}
			if(nunaligned>0 && nunaligned==prevnunaligned){
				//we have isolated images, allow them to be aligned by pairs with invalid phase correlation
				//use the first valid shift
				for(int i=0;i<outcoords.length;i++){
					if(outcoords[i]==null){
						int[][] targets=findConnectedImages(pairs,i,nimages);
						for(int j=0;j<targets.length;j++){
							if(outcoords[targets[j][0]]!=null){
								//ComparePair pair=pairs2.get(targets[j][1]);
								float[] pair=pairStats[targets[j][1]];
								int targetmember=pairs[targets[j][1]][1];
	    						float[] targetshift=new float[]{pair[0],pair[1]};
	    						targetshift[0]*=-1.0f; //this time we are aligning the query to the target
	    						targetshift[1]*=-1.0f;
	    						if(targetmember==i){//the query is the second member of the pair
	    							targetmember=pairs[targets[j][1]][0];
	    							targetshift[0]*=-1.0f;
	    							targetshift[1]*=-1.0f;
	    						}
	    						float[] targetcoords=outcoords[targetmember];
	    						float[] newcoords={targetcoords[0]+targetshift[0],targetcoords[1]+targetshift[1]};
	    						outcoords[i]=newcoords;
	    						naligned++; nunaligned--;
	    						break;
							}
						}
					}
				}
			}
		}
		return outcoords;
	}
	
	/*public void globOptIterJay(Vector<ComparePair> pairs2,float[][] coords,int fixedindex,int[][] pairs,int nimages,float ftoler){
		//here we perform one iteration of my poor-man's global optimization
		//each set of coordinates is shifted in the direction of the weighted average of the surrounding shifts
		//don't allow shifts greater than ftoler*dimension_size
		//ignore the fixed image for now, and shift everything by its values at the end
		//start with the image with highest overall weighted distance
		float maxcorr=0.0f;
		for(int i=0;i<pairs.length;i++){
			float corr=pairs2.get(i).getCrossCorrelation();
			if(corr>maxcorr) corr=maxcorr;
		}
		float[][] stats=new float[nimages][3];
		for(int i=0;i<nimages;i++){
			int[][] connected=findConnectedImages(pairs,i,nimages);
			for(int j=0;j<connected.length;j++){
				ComparePair pair=pairs2.get(connected[j][1]);
				float corr=pair.getCrossCorrelation();
				float currxshift=coords[connected[j][0]][0]-coords[connected[]]
			}
		}
	}*/
	
	public static int[][] findConnectedImages(int[][] pairs,int query,int nimages){
		int[][] targets=new int[nimages][2];  //first value is other image index, second is pair index
		int ntargets=0;
		for(int i=0;i<pairs.length;i++){
			if(pairs[i][0]==query){
				targets[ntargets][0]=pairs[i][1]; targets[ntargets][1]=i; ntargets++;
			}
			if(pairs[i][1]==query){
				targets[ntargets][0]=pairs[i][0]; targets[ntargets][1]=i; ntargets++;
			}
		}
		if(ntargets==0) return null;
		int[][] targets2=new int[ntargets][2];
		for(int i=0;i<ntargets;i++) targets2[i]=targets[i];
		return targets2;
	}
	
}
