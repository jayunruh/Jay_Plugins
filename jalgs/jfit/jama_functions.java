package jalgs.jfit;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class jama_functions{
	//code from http://the-lost-beauty.blogspot.com/2009/04/moore-penrose-pseudoinverse-in-jama.html
	
	public static double[][] pinv(double[][] input){
		Matrix x=new Matrix(input);
		Matrix inv=pinv(x);
		return inv.getArray();
	}
	
	public static double[][] pinv_jsvd(double[][] input){
		Matrix x=new Matrix(input);
		Matrix inv=pinv_jsvd(x);
		return inv.getArray();
	}
	
	 public static double MACHEPS = 2E-16;

	 /**
	  * Updates MACHEPS for the executing machine.
	  */
	 public static void updateMacheps() {
	  MACHEPS = 1;
	  do
	   MACHEPS /= 2;
	  while (1 + MACHEPS / 2 != 1);
	 }

	 /**
	  * Computes the Moore-Penrose pseudoinverse using the SVD method.
	  * 
	  * Modified version of the original implementation by Kim van der Linde.
	  */
	 public static Matrix pinv(Matrix x) {
	  int rows = x.getRowDimension();
	  int cols = x.getColumnDimension();
	  if (rows < cols) {
	   Matrix result = pinv(x.transpose());
	   if (result != null)
	    result = result.transpose();
	   return result;
	  }
	  SingularValueDecomposition svdX = new SingularValueDecomposition(x);
	  if (svdX.rank() < 1)
	   return null;
	  double[] singularValues = svdX.getSingularValues();
	  double tol = Math.max(rows, cols) * singularValues[0] * MACHEPS;
	  double[] singularValueReciprocals = new double[singularValues.length];
	  for (int i = 0; i < singularValues.length; i++)
	   if (Math.abs(singularValues[i]) >= tol)
	    singularValueReciprocals[i] =  1.0 / singularValues[i];
	  double[][] u = svdX.getU().getArray();
	  double[][] v = svdX.getV().getArray();
	  int min = Math.min(cols, u[0].length);
	  double[][] inverse = new double[cols][rows];
	  for (int i = 0; i < cols; i++)
	   for (int j = 0; j < u.length; j++)
	    for (int k = 0; k < min; k++)
	     inverse[i][j] += v[i][k] * singularValueReciprocals[k] * u[j][k];
	  //return new Matrix(inverse);
	  return new Matrix(inverse);
	 }
	 
	 public static Matrix pinv_jsvd(Matrix x){
		  int rows = x.getRowDimension(); //this is m, m>n
		  int cols = x.getColumnDimension(); //this is n
		  if (rows < cols) {
		   Matrix result = pinv_jsvd(x.transpose());
		   if (result != null)
		    result = result.transpose();
		   return result;
		  }
		  //these are the diagonals of the m x n sigma matrix
		  double[] singularValues=new double[cols];
		  double[][] v=new double[cols][cols]; //this is n by n
		  double[][] u=jsvd.svdcmp(x.getArray(),singularValues,v);
		  //SingularValueDecomposition svdX = new SingularValueDecomposition(x);
		  //if (svdX.rank() < 1)
		   //return null;
		  //double[] singularValues = svdX.getSingularValues();
		  double tol = Math.max(rows, cols) * singularValues[0] * MACHEPS;
		  double[] singularValueReciprocals = new double[singularValues.length];
		  for (int i = 0; i < singularValues.length; i++)
		   if (Math.abs(singularValues[i]) >= tol)
		    singularValueReciprocals[i] =  1.0 / singularValues[i];
		  //double[][] u = svdX.getU().getArray();
		  //double[][] v = svdX.getV().getArray();
		  int min = Math.min(cols, u[0].length);
		  double[][] inverse = new double[cols][rows];
		  for (int i = 0; i < cols; i++)
		   for (int j = 0; j < u.length; j++)
		    for (int k = 0; k < min; k++)
		     inverse[i][j] += v[i][k] * singularValueReciprocals[k] * u[j][k];
		  //return new Matrix(inverse);
		  return new Matrix(inverse);
	 }

}
