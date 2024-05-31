package jguis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.io.File;

//import ij.IJ;
//import ij.ImagePlus;
//import jalgs.algutils;
import jalgs.interpolation;
import jalgs.jdataio;

public class pdb_volume_measurements {
	
	public static String[] atomnames={"ALA_CB", "ALA_C", "ALA_CA", "ALA_O", "ALA_N", "ARG_C", "ARG_CA", "ARG_CB", "ARG_CG", 
			"ARG_CD", "ARG_CZ", "ARG_N", "ARG_NE", "ARG_NH1", "ARG_NH2", "ARG_O", "ASN_C", "ASN_CA", 
			"ASN_CB", "ASN_CG", "ASN_N", "ASN_ND2", "ASN_O", "ASN_OD1", "ASP_C", "ASP_CA", "ASP_CB", 
			"ASP_CG", "ASP_N", "ASP_O", "ASP_OD1", "ASP_OD2", "CYS_C", "CYS_CA", "CYS_CB", "CYS_N", 
			"CYS_O", "CYS_SG", "GLN_C", "GLN_CA", "GLN_CB", "GLN_CG", "GLN_CD", "GLN_N", "GLN_NE2", 
			"GLN_O", "GLN_OE1", "GLU_C", "GLU_CA", "GLU_CB", "GLU_CG", "GLU_CD", "GLU_N", "GLU_O", 
			"GLU_OE1", "GLU_OE2", "GLY_C", "GLY_CA", "GLY_O", "GLY_N", "HIS_C", "HIS_CA", "HIS_CB", 
			"HIS_CG", "HIS_CD2", "HIS_CE1", "HIS_N", "HIS_ND1", "HIS_NE2", "HIS_O", "HYP_C", "HYP_CA", 
			"HYP_CB", "HYP_CG", "HYP_CD", "HYP_N", "HYP_O", "HYP_OD1", "ILE_C", "ILE_CA", "ILE_CB", 
			"ILE_CG1", "ILE_CG2", "ILE_CD1", "ILE_N", "ILE_O", "LEU_C", "LEU_CA", "LEU_CB", "LEU_CG", 
			"LEU_CD1", "LEU_CD2", "LEU_N", "LEU_O", "LYS_C", "LYS_CA", "LYS_CB", "LYS_CG", "LYS_CD", 
			"LYS_CE", "LYS_NZ", "LYS_N", "LYS_O", "MET_C", "MET_CA", "MET_CB", "MET_CG", "MET_CE", 
			"MET_N", "MET_O", "MET_SD", "MSE_C", "MSE_CA", "MSE_CB", "MSE_CG", "MSE_CE", "MSE_N", 
			"MSE_O", "MSE_SE", "UNK_C", "UNK_CA", "UNK_N", "UNK_O", "ACE_C", "ACE_CH3", "ACE_O", 
			"NME_N", "NME_C", "NH2_N", "PCA_C", "PCA_CA", "PCA_CB", "PCA_CG", "PCA_CD", "PCA_N", 
			"PCA_O", "PCA_OE", "PHE_C", "PHE_CA", "PHE_CB", "PHE_CG", "PHE_CD1", "PHE_CD2", "PHE_CE1", 
			"PHE_CE2", "PHE_CZ", "PHE_N", "PHE_O", "PRO_C", "PRO_CA", "PRO_CB", "PRO_CG", "PRO_CD", 
			"PRO_N", "PRO_O", "SER_C", "SER_CA", "SER_CB", "SER_N", "SER_O", "SER_OG", "THR_C", 
			"THR_CA", "THR_CB", "THR_CG2", "THR_N", "THR_O", "THR_OG1", "TRP_C", "TRP_CA", "TRP_CB", 
			"TRP_CG", "TRP_CD1", "TRP_CD2", "TRP_CE2", "TRP_CE3", "TRP_CH2", "TRP_CZ2", "TRP_CZ3", 
			"TRP_N", "TRP_NE1", "TRP_O", "TYR_C", "TYR_CA", "TYR_CB", "TYR_CG", "TYR_CD1", "TYR_CD2", 
			"TYR_CE1", "TYR_CE2", "TYR_CZ", "TYR_N", "TYR_O", "TYR_OH", "VAL_C", "VAL_CA", "VAL_CB", 
			"VAL_CG1", "VAL_CG2", "VAL_N", "VAL_O"};
	
	public static float[] fivals={0.4395f, -0.1002f, -0.1571f, -0.0233f, -0.6149f, -0.1002f, -0.1571f, 0.3212f, 0.3212f, 
			0.0116f, 0.5142f, -0.6149f, -0.1425f, -0.5995f, -0.5995f, -0.0233f, -0.1002f, -0.1571f, 
			0.0348f, -0.1002f, -0.6149f, -0.7185f, -0.0233f, -0.0233f, -0.1002f, -0.1571f, 0.0348f, 
			-0.1002f, -0.6149f, -0.0233f, -0.4087f, -0.4087f, -0.1002f, -0.1571f, 0.0116f, -0.6149f, 
			-0.0233f, 0.511f, -0.1002f, -0.1571f, 0.3212f, 0.0348f, -0.1002f, -0.6149f, -0.7185f, 
			-0.0233f, -0.0233f, -0.1002f, -0.1571f, 0.3212f, 0.0348f, -0.1002f, -0.6149f, -0.0233f, 
			-0.4087f, -0.4087f, -0.1002f, -0.2018f, -0.0233f, -0.6149f, -0.1002f, -0.1571f, 0.0348f, 
			0.2361f, 0.5185f, 0.1443f, -0.6149f, -0.266f, -0.266f, -0.0233f, -0.1002f, -0.1571f, 
			0.3212f, -0.0504f, 0.0116f, -0.5113f, -0.0233f, -0.4603f, -0.1002f, -0.1571f, -0.0015f, 
			0.4562f, 0.642f, 0.642f, -0.6149f, -0.0233f, -0.1002f, -0.1571f, 0.3212f, 0.066f, 0.642f, 
			0.642f, -0.6149f, -0.0233f, -0.1002f, -0.1571f, 0.3212f, 0.4562f, 0.4562f, 0.0116f, 
			-0.8535f, -0.6149f, -0.0233f, -0.1002f, -0.1571f, 0.3212f, 0.0116f, 0.1023f, -0.6149f, 
			-0.0233f, 0.5906f, -0.1002f, -0.1571f, 0.3212f, 0.0116f, 0.1023f, -0.6149f, -0.0233f, 
			0.6601f, -0.1002f, -0.1571f, -0.6149f, -0.0233f, -0.1002f, 0.0099f, -0.0233f, -0.6149f, 
			0.1023f, -0.7185f, -0.1002f, -0.1571f, 0.3212f, 0.0348f, -0.1002f, -0.6149f, -0.0233f, 
			-0.0233f, -0.1002f, -0.1571f, 0.3212f, 0.1492f, 0.305f, 0.305f, 0.305f, 0.305f, 0.305f, 
			-0.6149f, -0.0233f, -0.1002f, -0.1571f, 0.3212f, 0.3212f, 0.0116f, -0.5113f, -0.0233f, 
			-0.1002f, -0.1571f, 0.0116f, -0.6149f, -0.0233f, -0.4603f, -0.1002f, -0.1571f, -0.0514f, 
			0.4395f, -0.6149f, -0.0233f, -0.4603f, -0.1002f, -0.1571f, 0.0348f, 0.1492f, 0.5185f, 
			0.1492f, 0.1539f, 0.305f, 0.305f, 0.305f, 0.305f, -0.6149f, 0.0223f, -0.0233f, -0.1002f, 
			-0.1571f, 0.3212f, 0.1492f, 0.305f, 0.305f, 0.305f, 0.305f, 0.1539f, -0.6149f, -0.0233f, 
			-0.1163f, -0.1002f, -0.1571f, -0.0015f, 0.642f, 0.642f, -0.6149f, -0.0233f};
	
	public static String aanames="AVLIMFYWCGPSTNQRHKDE";
	
	public static String[] aanames2= {"ALA","VAL","LEU","ILE","MET","PHE","TYR","TRP","CYS",
			"GLY","PRO","SER","THR","ASN","GLN","ARG","HIS","LYS","ASP","GLU"};
	
	String[] pdbcolnames= {"type","atom","atype","resname","chain","residue","coords",
			"occ","temp","element"};
	
	public static HashMap<String,String> getAAHashMap(){
		HashMap<String, String> aamap = new HashMap<String, String>();
		for(int i=0;i<aanames2.length;i++) {
			aamap.put(aanames2[i],aanames.substring(i,i+1));
		}
		return aamap;
	}
	
	public static HashMap<String,Float> getfiHashMap(){
		HashMap<String, Float> fimap = new HashMap<String, Float>();
		for(int i=0;i<atomnames.length;i++) {
			fimap.put(atomnames[i],fivals[i]);
		}
		return fimap;
	}
	
	public static Object[] read_pdb(String fname){
		//returns an array of objects with pdb file columns
		//columns are 'type','atom','atype','resname','chain','residue','x','y','z','occ','temp','element'
		jdataio jdio=new jdataio();
		String[] lines=null;
		try{
			BufferedReader b=new BufferedReader(new FileReader(fname));
			lines=jdio.readstringfilelines(b);
			System.out.println("read "+lines.length+" lines from file");
			System.out.println("line1:"+lines[0]);
			b.close();
		} catch(IOException e){
			System.err.println("File Reader:"+e.getMessage());
			return null;
		}
		String[] types=new String[lines.length];
		int[] atoms=new int[lines.length];
		String[] atypes=new String[lines.length];
		String[] resnames=new String[lines.length];
		String[] chains=new String[lines.length];
		int[] residues=new int[lines.length];
		float[][] coords=new float[lines.length][];
		float[] occs=new float[lines.length];
		float[] temps=new float[lines.length];
		String[] elements=new String[lines.length];
		int na=0;
		for(int i=0;i<lines.length;i++) {
			String type=lines[i].substring(0,6).strip();
			if(type.equals("ATOM") || type.equals("HETATM")) {
				try {
					types[na]=type;
					atoms[na]=Integer.parseInt(lines[i].substring(6,11).strip());
					atypes[na]=lines[i].substring(12,16).strip();
					resnames[na]=lines[i].substring(17,20).strip();
					chains[na]=lines[i].substring(21,22).strip();
					residues[na]=Integer.parseInt(lines[i].substring(22,26).strip());
					coords[na]= new float[] {Float.parseFloat(lines[i].substring(30,38).strip()),
								Float.parseFloat(lines[i].substring(38,46).strip()),
								Float.parseFloat(lines[i].substring(46,54).strip())};
					occs[na]=Float.parseFloat(lines[i].substring(54,60).strip());
					temps[na]=Float.parseFloat(lines[i].substring(60,66).strip());
					elements[na]=lines[i].substring(76).strip();
				} catch(NumberFormatException e) {
					System.err.println("bad number in line"+i);
				}
				na++;
			}
		}
		//now get rid of all the ones we don't need
		types=(String[])getSubarray(types,na);
		atoms=(int[])getSubarray(atoms,na);
		atypes=(String[])getSubarray(atypes,na);
		resnames=(String[])getSubarray(resnames,na);
		chains=(String[])getSubarray(chains,na);
		residues=(int[])getSubarray(residues,na);
		coords=(float[][])getSubarray(coords,na);
		occs=(float[])getSubarray(occs,na);
		temps=(float[])getSubarray(temps,na);
		elements=(String[])getSubarray(elements,na);
		return new Object[]{types,atoms,atypes,resnames,chains,residues,coords,occs,temps,elements};
	}
	
	public static Object getSubarray(Object arr,int len) {
		if(arr instanceof String[]) {
			String[] arr2=new String[len];
			System.arraycopy(arr, 0, arr2, 0, len);
			return arr2;
		} else if(arr instanceof int[]) {
			int[] arr2=new int[len];
			System.arraycopy(arr, 0, arr2, 0, len);
			return arr2;
		} else if(arr instanceof float[]) {
			float[] arr2=new float[len];
			System.arraycopy(arr, 0, arr2, 0, len);
			return arr2;
		} else {
			//must be a 2D float
			float[][] arr2=new float[len][];
			for(int i=0;i<len;i++) {
				arr2[i]=((float[][])arr)[i];
			}
			return arr2;
		}
	}
	
	public static String getSequence(Object[] pdbcols) {
		int[] residues=(int[])pdbcols[5];
		System.out.println("pdb length ="+residues.length);
		String[] atype=(String[])pdbcols[2];
		String[] resname=(String[])pdbcols[3];
		String seq="";
		HashMap<String,String> aamap=getAAHashMap();
		int pos=0;
		while(pos<residues.length){
			if(atype[pos].equals("CA")){
				String tresname=resname[pos];
				if(aamap.containsKey(tresname)) {
					seq+=aamap.get(tresname);
				} else {
					seq+="X";
				}
			}
			pos++;
		}
		System.out.println("sequence:"+seq);
		return seq;
	}
	
	public static float getDist2(float[] pos,float[] coords){
		return (coords[0]-pos[0])*(coords[0]-pos[0])+(coords[1]-pos[1])*(coords[1]-pos[1])+(coords[2]-pos[2])*(coords[2]-pos[2]);
	}

	public static float fauchere(float[] fi,float[] pos,float[][] coords,float max_dist2){
	  float sum=0.0f;
	  for(int i=0;i<coords.length;i++){
	    float dist2=getDist2(pos,coords[i]);
	    if(dist2<=max_dist2){
	      sum+=100*fi[i]*Math.exp(-Math.sqrt(dist2));
	    }
	  }
	  return sum;
	}
	
	public static float[][] getExtremes(float[][] coords){
		float[] mins=coords[0].clone();
		float[] maxs=coords[0].clone();
		for(int i=1;i<coords.length;i++) {
			for(int j=0;j<3;j++) {
				mins[j]=Math.min(mins[j],coords[i][j]);
				maxs[j]=Math.max(maxs[j],coords[i][j]);
			}
		}
		return new float[][] {mins,maxs};
	}
	
	public static Object[] mlp(String pdbname,float pad,float max_dist, float resolution){
		//typical pad and max_dist values are 5-8 angstroms and resolution is 1 angstrom
		//returns the matrix and the origin coordinates, span pixels, and the pdb columns
		Object[] pdbcols=read_pdb(pdbname);
		float[][] coords=(float[][])pdbcols[6];
		String[] atypes=(String[])pdbcols[2];
		String[] resnames=(String[])pdbcols[3];
		float max_dist2=max_dist*max_dist;
		//first find the min and max coordinate values (plus pad)
		float[][] minmax=getExtremes(coords);
		for(int i=0;i<3;i++) {
			minmax[0][i]-=pad;
			minmax[1][i]+=pad;
		}
		int[] spanpix=new int[3];
		for(int i=0;i<3;i++) {
			spanpix[i]=(int)((minmax[1][i]-minmax[0][i])/resolution)+1;
		}
		float[] finums=new float[coords.length];
		HashMap<String, Float> fimap = getfiHashMap();
		for(int i=0;i<atypes.length;i++) {
			String aname=resnames[i]+"_"+atypes[i];
			if(fimap.containsKey(aname)) {
				finums[i]=fimap.get(aname);
			} else {
				finums[i]=0.0f;
			}
		}
		float[][][] matrix=new float[spanpix[2]][spanpix[1]][spanpix[0]];
		for(int i=0;i<spanpix[2];i++) {
			float zpos=i*resolution+minmax[0][2];
			for(int j=0;j<spanpix[1];j++) {
				float ypos=j*resolution+minmax[0][1];
				for(int k=0;k<spanpix[0];k++) {
					float xpos=k*resolution+minmax[0][0];
					float[] pos= {xpos,ypos,zpos};
					matrix[i][j][k]=fauchere(finums,pos,coords,max_dist2);
				}
			}
		}
		return new Object[] {matrix,minmax[0],spanpix,pdbcols};
	}
	
	public static int findMatch(String[] arr,String query) {
		for(int i=0;i<arr.length;i++) {
			if(query==arr[i]) return i;
		}
		return -1;
	}
	
	public static int findMatch(int[] arr,int query) {
		for(int i=0;i<arr.length;i++) {
			if(query==arr[i]) return i;
		}
		return -1;
	}
	
	public static float[][] flatten3D2D(float[][][] arr) {
		int width=arr[0][0].length;
		int height=arr[0].length;
		float[][] outpix=new float[arr.length][width*height];
		for(int i=0;i<arr.length;i++) {
			for(int j=0;j<height;j++) {
				for(int k=0;k<width;k++) {
					outpix[i][k+j*width]=arr[i][j][k];
				}
			}
		}
		return outpix;
	}
	
	public static byte[][] flatten3D2D(byte[][][] arr) {
		int width=arr[0][0].length;
		int height=arr[0].length;
		byte[][] outpix=new byte[arr.length][width*height];
		for(int i=0;i<arr.length;i++) {
			for(int j=0;j<height;j++) {
				for(int k=0;k<width;k++) {
					outpix[i][k+j*width]=arr[i][j][k];
				}
			}
		}
		return outpix;
	}
	
	public static float[] getCrossLinkMLPs(float[][][] mlpmat1,float resolution,float[] spanstart,int[] spanpix,
			Object[] pdbcols,int[] selres) {
		//here we get the MLP values for cross-link residues
		float[][] mlpmat=flatten3D2D(mlpmat1);
		float[][] coords=(float[][])pdbcols[6];
		String[] atypes=(String[])pdbcols[2];
		String[] resnames=(String[])pdbcols[3];
		int[] resnums=(int[])pdbcols[5];
		float[] resmlp=new float[selres.length];
		for(int i=0;i<selres.length;i++) {
			int residx=findMatch(resnums,selres[i]);
			if(residx<0) {
				resmlp[i]=Float.NaN;
				continue;
			}
			String resname=resnames[residx];
			int atmidx=residx;
			if(resname.equals("LYS")) {
				while(resnums[atmidx]==selres[i]) {
					if(atypes[atmidx].equals("NZ")) {
						float[] tcoords=coords[atmidx].clone();
						tcoords[0]=(tcoords[0]-spanstart[0])/resolution;
						tcoords[1]=(tcoords[1]-spanstart[1])/resolution;
						tcoords[2]=(tcoords[2]-spanstart[2])/resolution;
						System.out.println("respix:"+tcoords[0]+","+tcoords[1]+","+tcoords[2]);
						resmlp[i]=interpolation.interp3D(mlpmat, spanpix[0], spanpix[1], 
								tcoords[0], tcoords[1], tcoords[2]);
						break;
					} else {
						atmidx+=1;
					}
				}
			} else if(resnums[i]==1) {
				//in this case we have a N term crosslink
				while(resnums[atmidx]==selres[i]){
					if(atypes[atmidx].equals("N")) {
						float[] tcoords=coords[atmidx].clone();
						tcoords[0]=(tcoords[0]-spanstart[0])/resolution;
						tcoords[1]=(tcoords[1]-spanstart[1])/resolution;
						tcoords[2]=(tcoords[2]-spanstart[2])/resolution;
						System.out.println("respix:"+tcoords[0]+","+tcoords[1]+","+tcoords[2]);
						resmlp[i]=interpolation.interp3D(mlpmat, spanpix[0], spanpix[1], 
								tcoords[0], tcoords[1], tcoords[2]);
						break;
					} else {
						atmidx+=1;
					}
				}
			} else {
				resmlp[i]=Float.NaN;
			}
		}
		return resmlp;
	}
	
	public static float[] getCrossLinkMLPs(String pdbname,int[] selres){
		//here we get the MLP values for cross-link residues
		float resolution=1.0f;
		Object[] temp=mlp(pdbname,5.0f,5.0f,resolution);
		Object[] pdbcols=(Object[])temp[3];
		float[][][] mlpmat1=(float[][][])temp[0];
		float[] spanstart=(float[])temp[1];
		int[] spanpix=(int[])temp[2];
		return getCrossLinkMLPs(mlpmat1,resolution,spanstart,spanpix,pdbcols,selres);
	}
	
	public static Object[] getMolecularSurface(String pdbname,float resolution,float vdwrad){
		/*
		 * returns a 3D surface mask for a molecule ,the spanstart, the spanpix, and the pdbcolumns
		 * typical resolution and van der waals radii are 0.25 and 2.6 angstroms
		 */
		Object[] pdbcols=read_pdb(pdbname);
		float[][] coords=(float[][])pdbcols[6];
		String[] atypes=(String[])pdbcols[2];
		String[] resnames=(String[])pdbcols[3];
		float srad=vdwrad/resolution;
		float[][] minmax=getExtremes(coords);
		for(int i=0;i<3;i++) {
			minmax[0][i]-=vdwrad;
			minmax[1][i]+=vdwrad;
		}
		int[] spanpix=new int[3];
		for(int i=0;i<3;i++) {
			spanpix[i]=(int)((minmax[1][i]-minmax[0][i])/resolution)+1;
		}
		byte[][] mask=new byte[spanpix[2]][spanpix[1]*spanpix[0]];
		for(int i=0;i<coords.length;i++) {
			int x=(int)((coords[i][0]-minmax[0][0])/resolution);
			int y=(int)((coords[i][1]-minmax[0][1])/resolution);
			int z=(int)((coords[i][2]-minmax[0][2])/resolution);
			mask[z][x+y*spanpix[0]]=(byte)255;
		}
		float[][] result=getEDT(mask,spanpix);
		for(int i=0;i<spanpix[2];i++) {
			for(int j=0;j<mask[i].length;j++) {
				if(result[i][j]>=(-srad)) mask[i][j]=(byte)255;
			}
		}
		return new Object[] {mask,minmax[0],spanpix,pdbcols};
	}
	
	public static float[][] getEDT(byte[][] mask,int[] spanpix){
		EDT edt=new EDT(spanpix[0],spanpix[1],spanpix[2]);
		return edt.getEDT(mask);
	}
	
	public static float[] getSurfaceCrossLinkDists(byte[][] mask,float[] spanstart,int[] spanpix,
			Object[] pdbcols,int[] selres,float resolution) {
		//invert the mask to get inside surface distance
		/*byte[][] mask=algutils.clone_multidim_array(mask1);
		for(int i=0;i<mask.length;i++) {
			for(int j=0;j<mask[i].length;j++) {
				if(mask[i][j]==(byte)255) mask[i][j]=(byte)0;
				else mask[i][j]=(byte)255;
			}
		}*/
		float[][] edt=getEDT(mask,spanpix);
		float[][] coords=(float[][])pdbcols[6];
		String[] atypes=(String[])pdbcols[2];
		String[] resnames=(String[])pdbcols[3];
		int[] resnums=(int[])pdbcols[5];
		float[] resdist=new float[selres.length];
		for(int i=0;i<selres.length;i++) {
			int residx=findMatch(resnums,selres[i]);
			if(residx<0) {
				resdist[i]=Float.NaN;
				continue;
			}
			String resname=resnames[residx];
			int atmidx=residx;
			if(resname.equals("LYS")) {
				while(resnums[atmidx]==selres[i]) {
					if(atypes[atmidx].equals("NZ")) {
						float[] tcoords=coords[atmidx].clone();
						tcoords[0]=(tcoords[0]-spanstart[0])/resolution;
						tcoords[1]=(tcoords[1]-spanstart[1])/resolution;
						tcoords[2]=(tcoords[2]-spanstart[2])/resolution;
						System.out.println("respix:"+tcoords[0]+","+tcoords[1]+","+tcoords[2]);
						resdist[i]=interpolation.interp3D(edt, spanpix[0], spanpix[1], 
								tcoords[0], tcoords[1], tcoords[2])*resolution;
						break;
					} else {
						atmidx+=1;
					}
				}
			} else if(resnums[i]==1) {
				//in this case we have a N term crosslink
				while(resnums[atmidx]==selres[i]){
					if(atypes[atmidx].equals("N")) {
						float[] tcoords=coords[atmidx].clone();
						tcoords[0]=(tcoords[0]-spanstart[0])/resolution;
						tcoords[1]=(tcoords[1]-spanstart[1])/resolution;
						tcoords[2]=(tcoords[2]-spanstart[2])/resolution;
						System.out.println("respix:"+tcoords[0]+","+tcoords[1]+","+tcoords[2]);
						resdist[i]=interpolation.interp3D(edt, spanpix[0], spanpix[1], 
								tcoords[0], tcoords[1], tcoords[2])*resolution;
						break;
					} else {
						atmidx+=1;
					}
				}
			} else {
				resdist[i]=Float.NaN;
			}
		}
		return resdist;
	}
	
	public static float[] getSurfaceCrossLinkDists(String pdbname,int[] selres) {
		float resolution=0.25f;
		Object[] temp=getMolecularSurface(pdbname,resolution,2.6f);
		byte[][] mask=(byte[][])temp[0];
		float[] spanstart=(float[])temp[1];
		int[] spanpix=(int[])temp[2];
		Object[] pdbcols=(Object[])temp[3];
		return getSurfaceCrossLinkDists(mask,spanstart,spanpix,pdbcols,selres,resolution);
	}
	
	public static void runBatch(String resfile,String pdbfolder,String outfile) {
		/*
		 * here we walk though a csv file with pdb ids and resa, resb values and xl dists
		 * get lipophilicity and surface distances for all pairs and write them into the file
		 */
		List<List<String>> table=table_tools.getTableFromFile(new File(resfile), ",", false);
		List<String> collabels=table.get(0);
		table.remove(0);
		//the columns are index,UPID,resa,resb,dist,conf_a,conf_b,conf_min
		float[] idxs=table_tools.get_column_array(table, 0);
		String[] upids=table_tools.get_listtable_column(table, 1);
		float[] resas=table_tools.get_column_array(table, 2);
		float[] resbs=table_tools.get_column_array(table, 3);
		float[] dists=table_tools.get_column_array(table, 4);
		//float[] confas=table_tools.get_column_array(table, 5);
		//float[] confbs=table_tools.get_column_array(table, 6);
		//float[] confmins=table_tools.get_column_array(table, 7);
		int ncols=table.get(0).size();
		float[][] outarrs=new float[4][upids.length];
		//if the output file partially exists read it in and initialize the values
		int startidx=0;
		if((new File(outfile)).exists()) {
			List<List<String>> table2=table_tools.getTableFromFile(new File(outfile), ",", false);
			table2.remove(0);
			for(int i=0;i<4;i++) {
				float[] temp=table_tools.get_column_array(table2,ncols+i);
				startidx=temp.length;
				for(int j=0;j<temp.length;j++) {
					outarrs[i][j]=temp[j];
				}
			}
		}
		String[] newcollabels= {"a_surf_dist","b_surf_dist","a_lipophilicity","b_lipophilicity"};
		String lastpdbname="";
		Object[] mlp=new Object[] {};
		Object[] surf=new Object[] {};
		for(int i=startidx;i<upids.length;i++) {
			String pdbname="AF-"+upids[i]+"-F1-model_v4.pdb";
			String pdbpath=pdbfolder+File.separator+pdbname;
			System.out.println("analyzing:"+(int)idxs[i]+","+pdbpath);
			if(dists[i]==0.0f) {
				outarrs[0][i]=Float.NaN; outarrs[1][i]=Float.NaN; 
				outarrs[2][i]=Float.NaN; outarrs[3][i]=Float.NaN;
				System.out.println("distance is 0");
				continue;
			}
			if(!(new File(pdbpath)).exists()) {
				outarrs[0][i]=Float.NaN; outarrs[1][i]=Float.NaN; 
				outarrs[2][i]=Float.NaN; outarrs[3][i]=Float.NaN;
				System.out.println("file not found");
				continue;
			}
			int[] xlres=new int[] {(int)resas[i],(int)resbs[i]};
			if(pdbname!=lastpdbname) {
				mlp=mlp(pdbpath,5.0f,5.0f,1.0f);
			}
			float[] spanstart=(float[])mlp[1];
			int[] spanpix=(int[])mlp[2];
			float[][][] mlpmat=(float[][][])mlp[0];
			Object[] pdbcols=(Object[])mlp[3];
			float[] mlpvals=getCrossLinkMLPs(mlpmat,1.0f,spanstart,spanpix,
					pdbcols,xlres);
			if(pdbname!=lastpdbname) {
				surf=getMolecularSurface(pdbpath,0.25f,2.6f);
			}
			byte[][] surfmask=(byte[][])surf[0];
			spanstart=(float[])surf[1];
			spanpix=(int[])surf[2];
			pdbcols=(Object[])surf[3];
			float[] surfdists=getSurfaceCrossLinkDists(surfmask,spanstart,spanpix,
					pdbcols,xlres,0.25f);

			outarrs[0][i]=surfdists[0];
			outarrs[1][i]=surfdists[1];
			outarrs[2][i]=mlpvals[0];
			outarrs[3][i]=mlpvals[1];
			if(i%5==0) {
				writePartialTable(table,collabels,newcollabels,outarrs,outfile,i);
			}
			lastpdbname=pdbname;
		}
		writePartialTable(table,collabels,newcollabels,outarrs,outfile,-1);
	}
	
	public static void writePartialTable(List<List<String>> table1,List<String> collabels1,
			String[] newcollabels,float[][] outarrs,String outfile,int sublen) {
		List<List<String>> table=table_tools.copylisttable(table1);
		List<String> collabels=table_tools.copylist(collabels1);
		for(int i=0;i<outarrs.length;i++) {
			table_tools.add_listtable_column(table, outarrs[i], -1);
		}
		for(int i=0;i<newcollabels.length;i++) collabels.add(newcollabels[i]);
		if(sublen<0) {
			table_tools.writeTableToFile(outfile,table_tools.list2stringarray(collabels),
					table,1);
		} else {
			table_tools.writeTableToFile(outfile,table_tools.list2stringarray(collabels),
					table.subList(0, sublen),1);
		}
	}
	
	public static void main(String[] args) {
		/*
		 * args should be a csv with idx,UPID,resa,resb, and dist values for crosslinks
		 * then a path to the pdb folder
		 * then a path for output
		 */
		runBatch(args[0],args[1],args[2]);
		//this was the previous testing paradigm
		/*boolean writemaps=false;
		String pdbname=args[0];
		String parent=(new File(pdbname)).getParent();
		int cl1=Integer.parseInt(args[1]);
		int cl2=Integer.parseInt(args[2]);
		int[] xlres=new int[] {cl1,cl2};
		System.out.println("path="+pdbname);
		Object[] mlp=mlp(pdbname,5.0f,5.0f,1.0f);
		float[] spanstart=(float[])mlp[1];
		System.out.println("mlp start coords = "+spanstart[0]+","+spanstart[1]+","+spanstart[2]);
		int[] spanpix=(int[])mlp[2];
		System.out.println("mlp span pix = "+spanpix[0]+","+spanpix[1]+","+spanpix[2]);
		float[][][] mlpmat=(float[][][])mlp[0];
		if(writemaps) {
			float[][] mlpmat2=flatten3D2D(mlpmat);
			ImagePlus imp=new ImagePlus("mlp_map",jutils.array2stack(mlpmat2,spanpix[0],spanpix[1]));
			IJ.saveAsTiff(imp,parent+File.separator+"mlp.tif");
		}
		Object[] pdbcols=(Object[])mlp[3];ls 0
		getSequence(pdbcols);
		float[] mlpvals=getCrossLinkMLPs(mlpmat,1.0f,spanstart,spanpix,
				pdbcols,xlres);
		Object[] surf=getMolecularSurface(pdbname,0.25f,2.6f);
		byte[][] surfmask=(byte[][])surf[0];

		spanstart=(float[])surf[1];
		System.out.println("surface start coords = "+spanstart[0]+","+spanstart[1]+","+spanstart[2]);
		spanpix=(int[])surf[2];
		System.out.println("surface span pix = "+spanpix[0]+","+spanpix[1]+","+spanpix[2]);
		if(writemaps) {
			ImagePlus imp=new ImagePlus("surface_mask",jutils.array2stack(surfmask,spanpix[0],spanpix[1]));
			IJ.saveAsTiff(imp,parent+File.separator+"surface.tif");
		}
		pdbcols=(Object[])surf[3];
		float[] surfdists=getSurfaceCrossLinkDists(surfmask,spanstart,spanpix,
				pdbcols,xlres,0.25f);
		System.out.println("mlp0,"+mlpvals[0]+",mlp1,"+mlpvals[1]+",sdist0,"+surfdists[0]+",sdist1,"+surfdists[1]);*/
	}

}