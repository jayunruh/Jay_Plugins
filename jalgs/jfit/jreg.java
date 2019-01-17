package jalgs.jfit;

import jalgs.matrixsolve;

public class jreg{
	//uses the kabsch algorithm

	public Object[] rigid_body_fiducials(float[][] f1,float[][] f2){
		// start by calculating the centroids
		float[] centroid1=centroidND(f1);
		float[] centroid2=centroidND(f2);
		int ndims=f1.length;
		if(ndims>3)
			return null;
		int nfids=f1[0].length;
		float[][] disp1=new float[ndims][nfids];
		float[][] disp2=new float[ndims][nfids];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<nfids;j++){
				disp1[i][j]=f1[i][j]-centroid1[i];
				disp2[i][j]=f2[i][j]-centroid2[i];
			}
		}
		float[][] covar=new float[ndims][ndims];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<ndims;j++){
				for(int k=0;k<nfids;k++){
					covar[i][j]+=disp1[i][k]*disp2[j][k];
				}
			}
		}
		// now perform the svd
		float[] w=new float[ndims];
		float[][] v=new float[ndims][ndims];
		float[][] u=jsvd.svdcmp(covar,w,v);
		// v=mattrans(v);
		// get the rotation matrix
		float[][] diag=new float[ndims][ndims];
		for(int i=0;i<ndims;i++)
			diag[i][i]=1.0f;
		diag[ndims-1][ndims-1]=matdet(matmult(v,u));
		float[][] r=matmult(v,diag);
		r=matmult(r,mattrans(u));
		float[] rcentroid1=new float[ndims];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<ndims;j++){
				rcentroid1[i]+=r[i][j]*centroid1[j];
			}
		}
		float[] t=new float[ndims];
		for(int i=0;i<ndims;i++)
			t[i]=centroid2[i]-rcentroid1[i];
		//r is the rotation matrix (ndims x ndims)
		//t is the translation vector (ndims length)
		//w is from svd
		//centroids are simply the centroids of the two data sets
		return new Object[]{r,t,w,centroid1,centroid2};
	}

	public Object[] scaled_rotation_fiducials(float[][] f1,float[][] f2){
		//the f1 and f2 matrices are ndims x npts
		// start by calculating the centroids
		float[] centroid1=centroidND(f1);
		float[] centroid2=centroidND(f2);
		int ndims=f1.length;
		if(ndims>3)
			return null;
		int nfids=f1[0].length;
		float[][] disp1=new float[ndims][nfids];
		float[][] disp2=new float[ndims][nfids];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<nfids;j++){
				disp1[i][j]=f1[i][j]-centroid1[i];
				disp2[i][j]=f2[i][j]-centroid2[i];
			}
		}
		float[][] covar=new float[ndims][ndims];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<ndims;j++){
				for(int k=0;k<nfids;k++){
					covar[i][j]+=disp1[i][k]*disp2[j][k];
				}
			}
		}
		// now perform the svd
		float[] w=new float[ndims];
		float[][] v=new float[ndims][ndims];
		float[][] u=jsvd.svdcmp(covar,w,v);
		// v=mattrans(v);
		// get the rotation matrix
		float[][] diag=new float[ndims][ndims];
		for(int i=0;i<ndims;i++)
			diag[i][i]=1.0f;
		diag[ndims-1][ndims-1]=matdet(matmult(v,u));
		float[][] r=matmult(v,diag);
		r=matmult(r,mattrans(u));
		float[] rcentroid1=new float[ndims];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<ndims;j++){
				rcentroid1[i]+=r[i][j]*centroid1[j];
			}
		}
		float[] t=new float[ndims];
		for(int i=0;i<ndims;i++)
			t[i]=centroid2[i]-rcentroid1[i];
		//estimate the scaling from the ratio of the average displacements from the centroid
		float ad1=0.0f;
		float ad2=0.0f;
		for(int i=0;i<nfids;i++){
			for(int j=0;j<ndims;j++){
				ad1+=disp1[j][i]*disp1[j][i];
				ad2+=disp2[j][i]*disp2[j][i];
			}
		}
		float s=(float)Math.sqrt(ad2/ad1);
		//r is the rotation matrix (ndims x ndims)
		//t is the translation vector (ndims length)
		//s is the scaling factor
		return new Object[]{r,t,new Float(s),centroid1,centroid2};
	}
	
	/**************************
	 * here we use least squares to find the best affine transformation
	 * @param f1: reference point set
	 * @param f2: point set to transform
	 * @return
	 */
	public Object[] affine_transformation_fiducials(float[][] f1,float[][] f2){
		// start by calculating the centroids
		float[] centroid1=centroidND(f1);
		float[] centroid2=centroidND(f2);
		int ndims=f1.length;
		if(ndims>3)
			return null;
		int nfids=f1[0].length;
		float[][] disp1=new float[ndims][nfids];
		float[][] disp2=new float[ndims][nfids];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<nfids;j++){
				disp1[i][j]=f1[i][j]-centroid1[i];
				disp2[i][j]=f2[i][j]-centroid2[i];
			}
		}
		double sumx2=0.0; double sumxy=0.0; double sumy2=0.0; double sumz2=0.0; double sumxz=0.0; double sumyz=0.0;
		double sumx=0.0; double sumy=0.0; double sumz=0.0;
		double osumx=0.0; double osumy=0.0; double osumz=0.0;
		double osumx2=0.0; double osumy2=0.0; double osumz2=0.0;
		double osumxy=0.0; double osumxz=0.0; double osumyz=0.0;
		double osumyx=0.0; double osumzx=0.0; double osumzy=0.0;
		for(int i=0;i<nfids;i++){
			sumx+=disp2[0][i]; sumy+=disp2[1][i]; sumz+=disp2[2][i];
			sumx2+=disp2[0][i]*disp2[0][i]; sumy2+=disp2[1][i]*disp2[1][i]; sumz2+=disp2[2][i]*disp2[2][i];
			sumxy+=disp2[0][i]*disp2[1][i]; sumxz+=disp2[0][i]*disp2[2][i]; sumyz+=disp2[1][i]*disp2[2][i];
			osumx+=disp1[0][i]; osumy+=disp1[1][i]; osumz+=disp1[2][i];
			osumx2+=disp1[0][i]*disp2[0][i]; osumy2+=disp1[1][i]*disp2[1][i]; osumz2+=disp1[2][i]*disp2[2][i]; 
			osumxy+=disp1[0][i]*disp2[1][i]; osumxz+=disp1[0][i]*disp2[2][i]; osumyz+=disp1[1][i]*disp2[2][i]; 
			osumyx+=disp2[0][i]*disp1[1][i]; osumzx+=disp2[0][i]*disp1[2][i]; osumzy+=disp2[1][i]*disp1[2][i]; 
		}
		//our general approach is to list the unknowns of the affine transform linearly: a11, a12, a13, a14, a21, a22, a23, a24, a31...		
		//then multiply (from the left) by a 12 x 12 matrix with lots of zeros and solve with matrix inversion
		//to simplify, every third line is collected together to make the matrix block diagonal
		//then the a vectors are a1i,... ; a2i...; and a3i,..
		//the right hand side vector is collected every third just like the left hand side matrix
		//in this way, we do three 4x4 solutions rather than one 12x12 solution
		double[][] mat={{sumx2,sumxy,sumxz,sumx},{sumxy,sumy2,sumyz,sumy},{sumxz,sumyz,sumz2,sumz},{sumx,sumy,sumz,1.0}};
		double[] vec1={osumx2,osumxy,osumxz,osumx};
		double[] vec2={osumyx,osumy2,osumyz,osumy};
		double[] vec3={osumzx,osumzy,osumz2,osumz};
		matrixsolve ms=new matrixsolve();
		double[] a1i=new double[4];
		ms.gjsolve(mat,vec1,a1i,4);
		double[] a2i=new double[4];
		ms.gjsolve(mat,vec2,a2i,4);
		double[] a3i=new double[4];
		ms.gjsolve(mat,vec3,a3i,4);
		double[][] result={a1i,a2i,a3i,{0,0,0,1}};
		return new Object[]{result,centroid1,centroid2};
	}

	public float centroid1D(float[] data){
		float sum=data[0];
		for(int i=1;i<data.length;i++)
			sum+=data[i];
		return sum/data.length;
	}

	public float[] centroidND(float[][] data){
		float[] centroid=new float[data.length];
		for(int i=0;i<data.length;i++)
			centroid[i]=centroid1D(data[i]);
		return centroid;
	}

	public float[][] matmult(float[][] a,float[][] b){
		float[][] res=new float[a.length][b[0].length];
		for(int i=0;i<a.length;i++){
			for(int j=0;j<b.length;j++){
				for(int k=0;k<a[0].length;k++){
					res[i][j]+=a[i][k]*b[k][j];
				}
			}
		}
		return res;
	}

	public float[][] mattrans(float[][] a){
		float[][] res=new float[a[0].length][a.length];
		for(int i=0;i<a.length;i++){
			for(int j=0;j<a[0].length;j++){
				res[j][i]=a[i][j];
			}
		}
		return res;
	}

	public float matdet(float[][] mat){
		if(mat.length==2)
			return matdet2(mat);
		if(mat.length==3)
			return matdet3(mat);
		return 0.0f;
	}

	public float matdet2(float[][] mat){
		return mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
	}

	public float matdet3(float[][] mat){
		float temp=mat[0][0]*mat[1][1]*mat[2][2];
		temp+=mat[0][1]*mat[1][2]*mat[2][0];
		temp+=mat[0][2]*mat[1][0]*mat[2][1];
		temp-=mat[0][2]*mat[1][1]*mat[2][0];
		temp-=mat[0][1]*mat[1][0]*mat[2][2];
		temp-=mat[0][0]*mat[1][2]*mat[2][1];
		return temp;
	}

}
