/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jfit;

public class jsvd{
	// this code is adapted from numerical recipes

	public static float[][] svdcmp(float[][] a1,float w[],float[][] v){
		// returns the u matrix
		//here A = UWV* (W is the diagonal)
		//alternatively can have M = USV*
		int m=a1.length;
		int n=a1[0].length;
		float[][] a=copyarr(a1);
		int flag,i,its,j,jj,k,l,nm;
		nm=0;
		l=0;
		float anorm,c,f,g,h,s,scale,x,y,z;
		float[] rv1=new float[n];
		g=scale=anorm=0.0f;
		for(i=1;i<=n;i++){
			l=i+1;
			rv1[i-1]=scale*g;
			g=0.0f;
			s=0.0f;
			scale=0.0f;
			if(i<=m){
				for(k=i;k<=m;k++)
					scale+=Math.abs(a[k-1][i-1]);
				if(scale!=0.0f){
					for(k=i;k<=m;k++){
						a[k-1][i-1]/=scale;
						s+=a[k-1][i-1]*a[k-1][i-1];
					}
					f=a[i-1][i-1];
					g=-jsvd.sign((float)Math.sqrt(s),f);
					h=f*g-s;
					a[i-1][i-1]=f-g;
					for(j=l;j<=n;j++){
						s=0.0f;
						for(k=i;k<=m;k++)
							s+=a[k-1][i-1]*a[k-1][j-1];
						f=s/h;
						for(k=i;k<=m;k++)
							a[k-1][j-1]+=f*a[k-1][i-1];
					}
					for(k=i;k<=m;k++)
						a[k-1][i-1]*=scale;
				}
			}
			w[i-1]=scale*g;
			g=0.0f;
			s=0.0f;
			scale=0.0f;
			if(i<=m&&i!=n){
				for(k=l;k<=n;k++)
					scale+=Math.abs(a[i-1][k-1]);
				if(scale!=0.0f){
					for(k=l;k<=n;k++){
						a[i-1][k-1]/=scale;
						s+=a[i-1][k-1]*a[i-1][k-1];
					}
					f=a[i-1][l-1];
					g=-sign((float)Math.sqrt(s),f);
					h=f*g-s;
					a[i-1][l-1]=f-g;
					for(k=l;k<=n;k++)
						rv1[k-1]=a[i-1][k-1]/h;
					for(j=l;j<=m;j++){
						s=0.0f;
						for(k=l;k<=n;k++)
							s+=a[j-1][k-1]*a[i-1][k-1];
						for(k=l;k<=n;k++)
							a[j-1][k-1]+=s*rv1[k-1];
					}
					for(k=l;k<=n;k++)
						a[i-1][k-1]*=scale;
				}
			}
			anorm=Math.max(anorm,(Math.abs(w[i-1])+Math.abs(rv1[i-1])));
		}
		for(i=n;i>=1;i--){
			if(i<n){
				if(g!=0.0f){
					for(j=l;j<=n;j++)
						v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
					for(j=l;j<=n;j++){
						s=0.0f;
						for(k=l;k<=n;k++)
							s+=a[i-1][k-1]*v[k-1][j-1];
						for(k=l;k<=n;k++)
							v[k-1][j-1]+=s*v[k-1][i-1];
					}
				}
				for(j=l;j<=n;j++){
					v[i-1][j-1]=0.0f;
					v[j-1][i-1]=0.0f;
				}
			}
			v[i-1][i-1]=1.0f;
			g=rv1[i-1];
			l=i;
		}
		for(i=Math.min(m,n);i>=1;i--){
			l=i+1;
			g=w[i-1];
			for(j=l;j<=n;j++)
				a[i-1][j-1]=0.0f;
			if(g!=0.0f){
				g=1.0f/g;
				for(j=l;j<=n;j++){
					s=0.0f;
					for(k=l;k<=m;k++)
						s+=a[k-1][i-1]*a[k-1][j-1];
					f=(s/a[i-1][i-1])*g;
					for(k=i;k<=m;k++)
						a[k-1][j-1]+=f*a[k-1][i-1];
				}
				for(j=i;j<=m;j++)
					a[j-1][i-1]*=g;
			}else{
				for(j=i;j<=m;j++)
					a[j-1][i-1]=0.0f;
			}
			++a[i-1][i-1];
		}
		for(k=n;k>=1;k--){
			for(its=1;its<=30;its++){
				flag=1;
				for(l=k;l>=1;l--){
					nm=l-1;
					if(Math.abs(rv1[l-1])+anorm==anorm){
						flag=0;
						break;
					}
					if(Math.abs(w[nm-1])+anorm==anorm)
						break;
				}
				if(flag!=0){
					c=0.0f;
					s=1.0f;
					for(i=l;i<=k;i++){
						f=s*rv1[i-1];
						rv1[i-1]=c*rv1[i-1];
						if(Math.abs(f)+anorm==anorm)
							break;
						g=w[i-1];
						h=pythag(f,g);
						w[i-1]=h;
						h=1.0f/h;
						c=g*h;
						s=-f*h;
						for(j=1;j<=m;j++){
							y=a[j-1][nm-1];
							z=a[j-1][i-1];
							a[j-1][nm-1]=y*c+z*s;
							a[j-1][i-1]=z*c-y*s;
						}
					}
				}
				z=w[k-1];
				if(l==k){
					if(z<0.0f){
						w[k-1]=-z;
						for(j=1;j<=n;j++)
							v[j-1][k-1]=-v[j-1][k-1];
					}
					break;
				}
				if(its==30)
					return null;
				x=w[l-1];
				nm=k-1;
				y=w[nm-1];
				g=rv1[nm-1];
				h=rv1[k-1];
				f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0f*h*y);
				g=pythag(f,1.0f);
				f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
				c=1.0f;
				s=1.0f;
				for(j=l;j<=nm;j++){
					i=j+1;
					g=rv1[i-1];
					y=w[i-1];
					h=s*g;
					g=c*g;
					z=pythag(f,h);
					rv1[j-1]=z;
					c=f/z;
					s=h/z;
					f=x*c+g*s;
					g=g*c-x*s;
					h=y*s;
					y*=c;
					for(jj=1;jj<=n;jj++){
						x=v[jj-1][j-1];
						z=v[jj-1][i-1];
						v[jj-1][j-1]=x*c+z*s;
						v[jj-1][i-1]=z*c-x*s;
					}
					z=pythag(f,h);
					w[j-1]=z;
					if(z!=0.0f){
						z=1.0f/z;
						c=f*z;
						s=h*z;
					}
					f=c*g+s*y;
					x=c*y-s*g;
					for(jj=1;jj<=m;jj++){
						y=a[jj-1][j-1];
						z=a[jj-1][i-1];
						a[jj-1][j-1]=y*c+z*s;
						a[jj-1][i-1]=z*c-y*s;
					}
				}
				rv1[l-1]=0.0f;
				rv1[k-1]=f;
				w[k-1]=x;
			}
		}
		return a;
	}
	
	public static double[][] svdcmp(double[][] a1,double w[],double[][] v){
		// returns the u matrix
		//here A = UWV* (W is the diagonal)
		//alternatively can have M = USV*
		int m=a1.length;
		int n=a1[0].length;
		double[][] a=copyarr(a1);
		int flag,i,its,j,jj,k,l,nm;
		nm=0;
		l=0;
		double anorm,c,f,g,h,s,scale,x,y,z;
		double[] rv1=new double[n];
		g=scale=anorm=0.0;
		for(i=1;i<=n;i++){
			l=i+1;
			rv1[i-1]=scale*g;
			g=0.0;
			s=0.0;
			scale=0.0;
			if(i<=m){
				for(k=i;k<=m;k++)
					scale+=Math.abs(a[k-1][i-1]);
				if(scale!=0.0){
					for(k=i;k<=m;k++){
						a[k-1][i-1]/=scale;
						s+=a[k-1][i-1]*a[k-1][i-1];
					}
					f=a[i-1][i-1];
					g=-jsvd.sign(Math.sqrt(s),f);
					h=f*g-s;
					a[i-1][i-1]=f-g;
					for(j=l;j<=n;j++){
						s=0.0;
						for(k=i;k<=m;k++)
							s+=a[k-1][i-1]*a[k-1][j-1];
						f=s/h;
						for(k=i;k<=m;k++)
							a[k-1][j-1]+=f*a[k-1][i-1];
					}
					for(k=i;k<=m;k++)
						a[k-1][i-1]*=scale;
				}
			}
			w[i-1]=scale*g;
			g=0.0;
			s=0.0;
			scale=0.0;
			if(i<=m&&i!=n){
				for(k=l;k<=n;k++)
					scale+=Math.abs(a[i-1][k-1]);
				if(scale!=0.0){
					for(k=l;k<=n;k++){
						a[i-1][k-1]/=scale;
						s+=a[i-1][k-1]*a[i-1][k-1];
					}
					f=a[i-1][l-1];
					g=-sign(Math.sqrt(s),f);
					h=f*g-s;
					a[i-1][l-1]=f-g;
					for(k=l;k<=n;k++)
						rv1[k-1]=a[i-1][k-1]/h;
					for(j=l;j<=m;j++){
						s=0.0;
						for(k=l;k<=n;k++)
							s+=a[j-1][k-1]*a[i-1][k-1];
						for(k=l;k<=n;k++)
							a[j-1][k-1]+=s*rv1[k-1];
					}
					for(k=l;k<=n;k++)
						a[i-1][k-1]*=scale;
				}
			}
			anorm=Math.max(anorm,(Math.abs(w[i-1])+Math.abs(rv1[i-1])));
		}
		for(i=n;i>=1;i--){
			if(i<n){
				if(g!=0.0){
					for(j=l;j<=n;j++)
						v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
					for(j=l;j<=n;j++){
						s=0.0;
						for(k=l;k<=n;k++)
							s+=a[i-1][k-1]*v[k-1][j-1];
						for(k=l;k<=n;k++)
							v[k-1][j-1]+=s*v[k-1][i-1];
					}
				}
				for(j=l;j<=n;j++){
					v[i-1][j-1]=0.0;
					v[j-1][i-1]=0.0;
				}
			}
			v[i-1][i-1]=1.0;
			g=rv1[i-1];
			l=i;
		}
		for(i=Math.min(m,n);i>=1;i--){
			l=i+1;
			g=w[i-1];
			for(j=l;j<=n;j++)
				a[i-1][j-1]=0.0;
			if(g!=0.0){
				g=1.0/g;
				for(j=l;j<=n;j++){
					s=0.0;
					for(k=l;k<=m;k++)
						s+=a[k-1][i-1]*a[k-1][j-1];
					f=(s/a[i-1][i-1])*g;
					for(k=i;k<=m;k++)
						a[k-1][j-1]+=f*a[k-1][i-1];
				}
				for(j=i;j<=m;j++)
					a[j-1][i-1]*=g;
			}else{
				for(j=i;j<=m;j++)
					a[j-1][i-1]=0.0;
			}
			++a[i-1][i-1];
		}
		for(k=n;k>=1;k--){
			for(its=1;its<=30;its++){
				flag=1;
				for(l=k;l>=1;l--){
					nm=l-1;
					if(Math.abs(rv1[l-1])+anorm==anorm){
						flag=0;
						break;
					}
					if(Math.abs(w[nm-1])+anorm==anorm)
						break;
				}
				if(flag!=0){
					c=0.0;
					s=1.0;
					for(i=l;i<=k;i++){
						f=s*rv1[i-1];
						rv1[i-1]=c*rv1[i-1];
						if(Math.abs(f)+anorm==anorm)
							break;
						g=w[i-1];
						h=pythag(f,g);
						w[i-1]=h;
						h=1.0/h;
						c=g*h;
						s=-f*h;
						for(j=1;j<=m;j++){
							y=a[j-1][nm-1];
							z=a[j-1][i-1];
							a[j-1][nm-1]=y*c+z*s;
							a[j-1][i-1]=z*c-y*s;
						}
					}
				}
				z=w[k-1];
				if(l==k){
					if(z<0.0){
						w[k-1]=-z;
						for(j=1;j<=n;j++)
							v[j-1][k-1]=-v[j-1][k-1];
					}
					break;
				}
				if(its==30)
					return null;
				x=w[l-1];
				nm=k-1;
				y=w[nm-1];
				g=rv1[nm-1];
				h=rv1[k-1];
				f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
				g=pythag(f,1.0);
				f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
				c=1.0;
				s=1.0;
				for(j=l;j<=nm;j++){
					i=j+1;
					g=rv1[i-1];
					y=w[i-1];
					h=s*g;
					g=c*g;
					z=pythag(f,h);
					rv1[j-1]=z;
					c=f/z;
					s=h/z;
					f=x*c+g*s;
					g=g*c-x*s;
					h=y*s;
					y*=c;
					for(jj=1;jj<=n;jj++){
						x=v[jj-1][j-1];
						z=v[jj-1][i-1];
						v[jj-1][j-1]=x*c+z*s;
						v[jj-1][i-1]=z*c-x*s;
					}
					z=pythag(f,h);
					w[j-1]=z;
					if(z!=0.0){
						z=1.0/z;
						c=f*z;
						s=h*z;
					}
					f=c*g+s*y;
					x=c*y-s*g;
					for(jj=1;jj<=m;jj++){
						y=a[jj-1][j-1];
						z=a[jj-1][i-1];
						a[jj-1][j-1]=y*c+z*s;
						a[jj-1][i-1]=z*c-y*s;
					}
				}
				rv1[l-1]=0.0;
				rv1[k-1]=f;
				w[k-1]=x;
			}
		}
		return a;
	}

	public static float pythag(float a,float b){
		float absa=Math.abs(a);
		float absb=Math.abs(b);
		if(absa>absb)
			return absa*(float)Math.sqrt(1.0+(absb/absa)*(absb/absa));
		else
			return (absb==0.0f)?0.0f:absb*(float)Math.sqrt(1.0+(absa/absb)*(absa/absb));
	}
	
	public static double pythag(double a,double b){
		double absa=Math.abs(a);
		double absb=Math.abs(b);
		if(absa>absb)
			return absa*Math.sqrt(1.0+(absb/absa)*(absb/absa));
		else
			return (absb==0.0)?0.0:absb*Math.sqrt(1.0+(absa/absb)*(absa/absb));
	}

	public static float sign(float a,float b){
		return Math.signum(b)*Math.abs(a);
	}
	
	public static double sign(double a,double b){
		return Math.signum(b)*Math.abs(a);
	}

	public static float[][] copyarr(float[][] input){
		float[][] temp=new float[input.length][];
		for(int i=0;i<input.length;i++){
			temp[i]=input[i].clone();
		}
		return temp;
	}
	
	public static double[][] copyarr(double[][] input){
		double[][] temp=new double[input.length][];
		for(int i=0;i<input.length;i++){
			temp[i]=input[i].clone();
		}
		return temp;
	}

}
