package jalgs.jfit;

public class jreg{

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
		return new Object[]{r,t,w,centroid1,centroid2};
	}

	public Object[] scaled_rotation_fiducials(float[][] f1,float[][] f2){
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
		float[][] rotdisp1=matmult(r,disp1);
		float s=0.0f;
		float den=0.0f;
		for(int i=0;i<nfids;i++){
			for(int k=0;k<ndims;k++){
				s+=rotdisp1[k][i]*disp2[k][i];
				den+=disp1[k][i]*disp1[k][i];
			}
		}
		s/=den;
		float[] rcentroid1=new float[ndims];
		for(int i=0;i<ndims;i++){
			for(int j=0;j<ndims;j++){
				rcentroid1[i]+=s*r[i][j]*centroid1[j];
			}
		}
		float[] t=new float[ndims];
		for(int i=0;i<ndims;i++)
			t[i]=centroid2[i]-rcentroid1[i];
		return new Object[]{r,t,new Float(s)};
	}

	public float centroid1D(float[] data){
		float sum=data[0];
		for(int i=1;i<data.length;i++)
			sum+=data[i];
		return sum/(float)data.length;
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
