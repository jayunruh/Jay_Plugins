/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

public class point3D implements Cloneable{
	public int x,y,z;
	public int rx,ry,rz;

	public point3D(int x1,int y1,int z1){
		x=x1;
		y=y1;
		z=z1;
		rx=x;
		ry=y;
		rz=z;
	}

	public void moveto(int x1,int y1,int z1){
		x=x1;
		y=y1;
		z=z1;
		rx=x;
		ry=y;
		rz=z;
	}

	public void reset(){
		rx=x;
		ry=y;
		rz=z;
	}

	public void translate(int transx,int transy,int transz){
		x+=transx;
		y+=transy;
		z+=transz;
		rx+=transx;
		ry+=transy;
		rz+=transz;
	}

	public void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		double radx=Math.toRadians(degx1);
		double rady=Math.toRadians(degy1);
		double radz=Math.toRadians(degz1);
		setrotation(radx,rady,radz,centerx,centery,centerz);
	}

	public void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		setrotation(radx1,rady1,radz1,centerx,centery,centerz);
	}

	public void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		setrotation(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
	}
	
	public void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz){
		addrotation(cosx1,cosy1,cosz1,sinx1,siny1,sinz1,centerx,centery,centerz);
	}

	public void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz){
		double radx=Math.toRadians(degx1);
		double rady=Math.toRadians(degy1);
		double radz=Math.toRadians(degz1);
		addrotation(radx,rady,radz,centerx,centery,centerz);
	}

	public void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz){
		addrotation(radx1,rady1,radz1,centerx,centery,centerz);
	}

	public void transform_perspective(double horizon_dist,int centerx,int centery,int centerz){
		if(horizon_dist>0.0){
			double tempx=rx-centerx;
			double tempy=ry-centery;
			double tempz=rz-centerz;
			double temphordist=(tempz+horizon_dist)/horizon_dist;
			if(temphordist<=0){
				tempx=0;
				tempy=0;
			}else{
				tempx*=temphordist;
				tempy*=temphordist;
			}
			rx=centerx+(int)tempx;
			ry=centery+(int)tempy;
		}
	}

	public void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz){
		// here the horizon is in the foreground
		if(horizon_dist>0.0){
			double tempx=rx-centerx;
			double tempy=ry-centery;
			double tempz=centerz-rz;
			double temphordist=(tempz+horizon_dist)/horizon_dist;
			if(temphordist<=0){
				tempx=0;
				tempy=0;
			}else{
				tempx*=temphordist;
				tempy*=temphordist;
			}
			rx=centerx+(int)tempx;
			ry=centery+(int)tempy;
		}
	}

	public void setrotation(double dx,double dy,double dz,int centerx,int centery,int centerz){
		// rotate about the x, y, and z axes in order
		double tempx=x-centerx;
		double tempy=y-centery;
		double tempz=z-centerz;
		if(dz!=0.0){
			double sinval=Math.sin(-dz);
			double cosval=Math.cos(-dz);
			double tempx1=tempx*cosval-tempy*sinval;
			double tempy1=tempx*sinval+tempy*cosval;
			tempx=tempx1;
			tempy=tempy1;
		}
		if(dy!=0.0){
			double sinval=Math.sin(dy);
			double cosval=Math.cos(dy);
			double tempx1=tempx*cosval+tempz*sinval;
			double tempz1=-tempx*sinval+tempz*cosval;
			tempx=tempx1;
			tempz=tempz1;
		}
		if(dx!=0.0){
			double sinval=Math.sin(dx);
			double cosval=Math.cos(dx);
			double tempy1=tempy*cosval-tempz*sinval;
			double tempz1=tempy*sinval+tempz*cosval;
			tempy=tempy1;
			tempz=tempz1;
		}
		rx=centerx+(int)tempx;
		ry=centery+(int)tempy;
		rz=centerz+(int)tempz;
	}

	public void setrotation(double cosdx,double cosdy,double cosdz,double sindx,double sindy,double sindz,int centerx,int centery,int centerz){
		// rotate about the x, y, and z axes in order
		double tempx=x-centerx;
		double tempy=y-centery;
		double tempz=z-centerz;
		if(sindz!=0.0){
			double sinval=sindz;
			double cosval=cosdz;
			double tempx1=tempx*cosval-tempy*sinval;
			double tempy1=tempx*sinval+tempy*cosval;
			tempx=tempx1;
			tempy=tempy1;
		}
		if(sindy!=0.0){
			double sinval=sindy;
			double cosval=cosdy;
			double tempx1=tempx*cosval+tempz*sinval;
			double tempz1=-tempx*sinval+tempz*cosval;
			tempx=tempx1;
			tempz=tempz1;
		}
		if(sindx!=0.0){
			double sinval=sindx;
			double cosval=cosdx;
			double tempy1=tempy*cosval-tempz*sinval;
			double tempz1=tempy*sinval+tempz*cosval;
			tempy=tempy1;
			tempz=tempz1;
		}
		rx=centerx+(int)tempx;
		ry=centery+(int)tempy;
		rz=centerz+(int)tempz;
	}

	public void addrotation(double dx,double dy,double dz,int centerx,int centery,int centerz){
		// rotate about the x, y, and z axes in order
		double tempx=rx-centerx;
		double tempy=ry-centery;
		double tempz=rz-centerz;
		if(dz!=0.0){
			double sinval=Math.sin(-dz);
			double cosval=Math.cos(-dz);
			double tempx1=tempx*cosval-tempy*sinval;
			double tempy1=tempx*sinval+tempy*cosval;
			tempx=tempx1;
			tempy=tempy1;
		}
		if(dy!=0.0){
			double sinval=Math.sin(dy);
			double cosval=Math.cos(dy);
			double tempx1=tempx*cosval+tempz*sinval;
			double tempz1=-tempx*sinval+tempz*cosval;
			tempx=tempx1;
			tempz=tempz1;
		}
		if(dx!=0.0){
			double sinval=Math.sin(dx);
			double cosval=Math.cos(dx);
			double tempy1=tempy*cosval-tempz*sinval;
			double tempz1=tempy*sinval+tempz*cosval;
			tempy=tempy1;
			tempz=tempz1;
		}
		rx=centerx+(int)tempx;
		ry=centery+(int)tempy;
		rz=centerz+(int)tempz;
	}
	
	public void addrotation(double cosdx,double cosdy,double cosdz,double sindx,double sindy,double sindz,int centerx,int centery,int centerz){
		// rotate about the x, y, and z axes in order
		double tempx=rx-centerx;
		double tempy=ry-centery;
		double tempz=rz-centerz;
		if(sindz!=0.0){
			double sinval=sindz;
			double cosval=cosdz;
			double tempx1=tempx*cosval-tempy*sinval;
			double tempy1=tempx*sinval+tempy*cosval;
			tempx=tempx1;
			tempy=tempy1;
		}
		if(sindy!=0.0){
			double sinval=sindy;
			double cosval=cosdy;
			double tempx1=tempx*cosval+tempz*sinval;
			double tempz1=-tempx*sinval+tempz*cosval;
			tempx=tempx1;
			tempz=tempz1;
		}
		if(sindx!=0.0){
			double sinval=sindx;
			double cosval=cosdx;
			double tempy1=tempy*cosval-tempz*sinval;
			double tempz1=tempy*sinval+tempz*cosval;
			tempy=tempy1;
			tempz=tempz1;
		}
		rx=centerx+(int)tempx;
		ry=centery+(int)tempy;
		rz=centerz+(int)tempz;
	}
	
	public void add_rot_about_vector(double theta,double ux,double uy,double uz,int centerx,int centery,int centerz){
		add_rot_about_vector(Math.cos(theta),Math.sin(theta),ux,uy,uz,centerx,centery,centerz);
	}
	
	public void add_rot_about_vector(double cosval,double sinval,double ux,double uy,double uz,int centerx,int centery,int centerz){
		double tempx=rx-centerx;
		double tempy=ry-centery;
		double tempz=rz-centerz;
		double tempx1=tempx*(cosval+ux*ux*(1.0-cosval))+tempy*(ux*uy*(1.0-cosval)-uz*sinval)+tempz*(ux*uy*(1.0-cosval)+uy*sinval);
		double tempy1=tempx*(ux*uy*(1.0-cosval)+uz*sinval)+tempy*(uy*uy*(1.0-cosval)+cosval)+tempz*(uz*uy*(1.0-cosval)-ux*sinval);
		double tempz1=tempx*(ux*uz*(1.0-cosval)-uy*sinval)+tempy*(uz*uy*(1.0-cosval)+ux*sinval)+tempz*(uz*uz*(1.0-cosval)+cosval);
		rx=centerx+(int)tempx1;
		ry=centery+(int)tempy1;
		rz=centerz+(int)tempz1;
	}
	
	public void set_rot_about_vector(double theta,double ux,double uy,double uz,int centerx,int centery,int centerz){
		set_rot_about_vector(Math.cos(theta),Math.sin(theta),ux,uy,uz,centerx,centery,centerz);
	}
	
	public void set_rot_about_vector(double cosval,double sinval,double ux,double uy,double uz,int centerx,int centery,int centerz){
		double tempx=x-centerx;
		double tempy=y-centery;
		double tempz=z-centerz;
		double tempx1=tempx*(cosval+ux*ux*(1.0-cosval))+tempy*(ux*uy*(1.0-cosval)-uz*sinval)+tempz*(ux*uy*(1.0-cosval)+uy*sinval);
		double tempy1=tempx*(ux*uy*(1.0-cosval)+uz*sinval)+tempy*(uy*uy*(1.0-cosval)+cosval)+tempz*(uz*uy*(1.0-cosval)-ux*sinval);
		double tempz1=tempx*(ux*uz*(1.0-cosval)-uy*sinval)+tempy*(uz*uy*(1.0-cosval)+ux*sinval)+tempz*(uz*uz*(1.0-cosval)+cosval);
		rx=centerx+(int)tempx1;
		ry=centery+(int)tempy1;
		rz=centerz+(int)tempz1;
	}
	
	public void set_obs_rot(double dx,double dy,double dz,double ux,double uy,double uz,int centerx,int centery,int centerz){
		//this rotates the object from 0,0,0 so that it is observed from vector, u
		//rotate about the cross product between the oberver and the z axis
		double[] cross={uy,-ux,0.0};
		double norm=Math.sqrt(cross[0]*cross[0]+cross[1]*cross[1]);
		cross[0]/=norm; cross[1]/=norm;
		double sintheta=norm;
		double costheta=uz;
		//double theta=Math.acos(uz);
		set_rot_about_vector(costheta,-sintheta,cross[0],cross[1],cross[2],centerx,centery,centerz);
	}
	
	public void transform(double[][] transmat,int centerx,int centery,int centerz){
		double tempx=x-centerx;
		double tempy=y-centery;
		double tempz=z-centerz;
		double tempx2=transmat[0][0]*tempx+transmat[0][1]*tempy+transmat[0][2]*tempz;
		double tempy2=transmat[1][0]*tempx+transmat[1][1]*tempy+transmat[1][2]*tempz;
		double tempz2=transmat[2][0]*tempx+transmat[2][1]*tempy+transmat[2][2]*tempz;
		rx=(int)(tempx2+centerx);
		ry=(int)(tempy2+centery);
		rz=(int)(tempz2+centerz);
	}
	
	public void addtransform(double[][] transmat,int centerx,int centery,int centerz){
		double tempx=rx-centerx;
		double tempy=ry-centery;
		double tempz=rz-centerz;
		double tempx2=transmat[0][0]*tempx+transmat[0][1]*tempy+transmat[0][2]*tempz;
		double tempy2=transmat[1][0]*tempx+transmat[1][1]*tempy+transmat[1][2]*tempz;
		double tempz2=transmat[2][0]*tempx+transmat[2][1]*tempy+transmat[2][2]*tempz;
		rx=(int)(tempx2+centerx);
		ry=(int)(tempy2+centery);
		rz=(int)(tempz2+centerz);
	}
	
	public point3D clone(){
		point3D temp=new point3D(x,y,z);
		temp.rx=rx;
		temp.ry=ry;
		temp.rz=rz;
		return temp;
	}

}
