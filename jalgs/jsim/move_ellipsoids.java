/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class move_ellipsoids{
	//here we implement diffusion for a rod-like object
	//diffusion is separated into radial and axial directions depending on the geometry of the rod
	//also have to handle rotational diffusion around the axial axes
	rngs random;
	public float fboxpixels,Dax,Drad,Drot,a,b,Dax2,Drad2,Drot2,a2,b2;
	public sim_setup set;

	/*public move_ellipsoids(sim_setup set1,float[] Dvals){
		set=set1;
		random=set1.random;
		fboxpixels=set1.fboxpixels;
		this.Dax=Dvals[0]; this.Drad=Dvals[1]; this.Drot=Dvals[2];
	}*/
	
	public move_ellipsoids(sim_setup set1,float a,float b,float visc){
		//note that a and b are in microns and viscosity is in centipoises
		set=set1;
		random=set1.random;
		fboxpixels=set1.fboxpixels;
		float[] Dvals=get_perrin_D(a,b,visc);
		this.a=a; this.b=b;
		Dax=Dvals[0]; Drad=Dvals[1]; Drot=Dvals[2];
		Dax2=Dax*(set.frametime2/(set.pixelsize*set.pixelsize));
		Drad2=Drad*(set.frametime2/(set.pixelsize*set.pixelsize));
		Drot2=Drot*set.frametime2;
		a2=a/set.pixelsize;
		b2=b/set.pixelsize;
	}
	
	public sim_particle[] create_particles(){
		//first calculate the center
		float xc=(float)random.unidev(fboxpixels,0.0f);
		float yc=(float)random.unidev(fboxpixels,0.0f);
		float zc=(float)random.unidev(fboxpixels,0.0f);
		//now get a random vector
		double[] randvec=set.random.random_3D_vector();
		float halfdist=a2;
		float[] coords1={xc+halfdist*(float)randvec[0],yc+halfdist*(float)randvec[1],zc+halfdist*(float)randvec[2]};
		float[] coords2={xc-halfdist*(float)randvec[0],yc-halfdist*(float)randvec[1],zc-halfdist*(float)randvec[2]};
		handle_boundary(coords1,coords2);
		sim_particle[] sps=new sim_particle[2];
		sps[0]=new sim_particle(coords1,new int[1],false,0);
		sps[1]=new sim_particle(coords2,new int[1],false,0);
		return sps;
	}
	
	public void do_move_particle(sim_particle sp1,sim_particle sp2){
		float[] coords1=sp1.coords;
		float[] coords2=sp2.coords;
		float[] axvec=get_ax_vec(coords1,coords2);
		//now translate along the axial direction
		float xdev=(float)random.gasdev(0.0,Math.sqrt(2.0*Dax2));
		translate(coords1,axvec,xdev);
		translate(coords2,axvec,xdev);
		//now calculate the radial direction (perpendicular to the axial direction)
		//use the cross product with the z axis (garanteed to be in the xy plane)
		float[] rad1=cross_prod(axvec,new float[]{0.0f,0.0f,1.0f});
		float ydev=(float)random.gasdev(0.0,Math.sqrt(2.0*Drad2));
		translate(coords1,rad1,ydev);
		translate(coords2,rad1,ydev);
		//now use the cross product between the previous two vectors to get the final one
		float[] rad2=cross_prod(axvec,rad1);
		float zdev=(float)random.gasdev(0.0,Math.sqrt(2.0*Drad2));
		translate(coords1,rad2,zdev);
		translate(coords2,rad2,zdev);
		//now rotate about the two radial vectors
		float thetadev=(float)random.gasdev(0.0f,Math.sqrt(2.0*Drot2));
		float costhetadev=(float)Math.cos(thetadev);
		float sinthetadev=(float)Math.sin(thetadev);
		//calculate the center point
		float[] center=get_center(coords1,coords2);
		rot_about_vector(coords1,center,rad1,costhetadev,sinthetadev);
		rot_about_vector(coords2,center,rad1,costhetadev,sinthetadev);
		rot_about_vector(rad2,new float[]{0.0f,0.0f,0.0f},rad1,costhetadev,sinthetadev);
		float thetadev2=(float)random.gasdev(0.0f,Math.sqrt(2.0*Drot2));
		float costhetadev2=(float)Math.cos(thetadev2);
		float sinthetadev2=(float)Math.sin(thetadev2);
		rot_about_vector(coords1,center,rad1,costhetadev2,sinthetadev2);
		rot_about_vector(coords2,center,rad1,costhetadev2,sinthetadev2);
		//now handle periodic boundary conditions
		handle_boundary(coords1,coords2);
	}
	
	public void handle_boundary(float[] coords1,float[] coords2){
		//handle periodic boundary conditions for the center, not the inviduals
		//that way they don't get separated
		float[] center=get_center(coords1,coords2);
		if(center[0]<0.0f){
			coords1[0]+=fboxpixels; coords2[0]+=fboxpixels;
		}
		if(center[0]>fboxpixels){
			coords1[0]-=fboxpixels; coords2[0]-=fboxpixels;
		}
		if(center[1]<0.0f){
			coords1[1]+=fboxpixels; coords2[1]+=fboxpixels;
		}
		if(center[1]>fboxpixels){
			coords1[1]-=fboxpixels; coords2[1]-=fboxpixels;
		}
		if(center[2]<0.0f){
			coords1[2]+=fboxpixels; coords2[2]+=fboxpixels;
		}
		if(center[2]>fboxpixels){
			coords1[2]-=fboxpixels; coords2[2]-=fboxpixels;
		}
	}
	
	public static float[] get_center(float[] coords1,float[] coords2){
		return new float[]{0.5f*(coords1[0]+coords2[0]),0.5f*(coords1[1]+coords2[1]),0.5f*(coords1[2]+coords2[2])};
	}
	
	public static float get_3D_dist(float[] coords1,float[] coords2){
		return (float)Math.sqrt((coords2[0]-coords1[0])*(coords2[0]-coords1[0])+(coords2[1]-coords1[1])*(coords2[1]-coords1[1])+(coords2[2]-coords1[2])*(coords2[2]-coords1[2]));
	}
	
	public static float[] get_perrin_D(float a1,float b1,float visc1){
		//a and b are in microns
		//visc is in centipoises (mPa*s)
		float a=a1/1000000.0f; //convert to meters
		float b=b1/1000000.0f; //convert to meters
		float visc=visc1/1000.0f; //convert to Pa*s
		float k=(float)(1.38e-23);
		float T=298;
		float p=a/b; //a is the long axis
		float vol=(4.0f/3.0f)*a*b*b*(float)Math.PI;
		float xi=(float)Math.sqrt(Math.abs(p*p-1.0f))/p;
		float S=(float)Math.log((1.0f+xi)/(1.0f-xi))/xi;
		float fsphere=6.0f*visc*(float)Math.PI*(float)Math.pow(3.0f*vol/(4.0f*(float)Math.PI),1.0f/3.0f);
		float fp=2.0f*(float)Math.pow(p,2.0f/3.0f)/S;
		float ftot=fp*fsphere;
		float temp=1.0f/(p*p);
		float feq=(4.0f/3.0f)*(temp-p*p)/(2.0f-S*(2.0f-temp));
		float frsphere=6.0f*visc*vol;
		float Dreq=k*T/(feq*frsphere);
		float const1=p*p-1.0f;
		float const2=(2.0f*p*p-1.0f)/(float)Math.pow(const1,1.5);
		const2*=(float)Math.log((p+Math.sqrt(const1))/(p-Math.sqrt(const1)));
		const2-=2.0f*p/const1;
		float Rtax=(8.0f*b/3.0f)*(1.0f/const2);
		float Dtax=k*T/(6.0f*(float)Math.PI*visc*Rtax);
		const2=(2.0f*p*p-3.0f)/(float)Math.pow(const1,1.5);
		const2*=(float)Math.log(p+Math.sqrt(const1));
		const2+=p/const1;
		float Rteq=(8.0f*b/3.0f)*(1.0f/const2);
		float Dteq=k*T/(6.0f*(float)Math.PI*visc*Rteq);
		//now convert to microns squared per second
		Dteq*=(float)1.0E12;
		Dtax*=(float)1.0E12;
		return new float[]{Dteq,Dtax,Dreq};
	}
	
	public void translate(float[] coords,float[] direction,float dist){
		for(int i=0;i<coords.length;i++){
			coords[i]+=dist*direction[i];
		}
	}
	
	public void rot_about_vector(float[] coords,float[] center,float[] axis,float cosval,float sinval){
		//double cosval=Math.cos(theta);
		//double sinval=Math.sin(theta);
		float rx=coords[0]; float ry=coords[1]; float rz=coords[2];
		float ux=axis[0]; float uy=axis[1]; float uz=axis[2];
		double tempx=rx-center[0];
		double tempy=ry-center[1];
		double tempz=rz-center[2];
		double tempx1=tempx*(cosval+ux*ux*(1.0-cosval))+tempy*(ux*uy*(1.0-cosval)-uz*sinval)+tempz*(ux*uy*(1.0-cosval)+uy*sinval);
		double tempy1=tempx*(ux*uy*(1.0-cosval)+uz*sinval)+tempy*(uy*uy*(1.0-cosval)+cosval)+tempz*(uz*uy*(1.0-cosval)-ux*sinval);
		double tempz1=tempx*(ux*uz*(1.0-cosval)-uy*sinval)+tempy*(uz*uy*(1.0-cosval)+ux*sinval)+tempz*(uz*uz*(1.0-cosval)+cosval);
		rx=center[0]+(float)tempx1;
		ry=center[1]+(float)tempy1;
		rz=center[2]+(float)tempz1;
		coords[0]=rx; coords[1]=ry; coords[2]=rz;
	}
		
	public float[] get_ax_vec(float[] coords1,float[] coords2){
		float[] coords3=new float[coords1.length];
		for(int i=0;i<coords1.length;i++) coords3[i]=coords2[i]-coords1[i];
		get_norm_vec(coords3);
		return coords3;
	}
	
	public void get_norm_vec(float[] coords){
		float length=0.0f;
		for(int i=0;i<coords.length;i++){
			length+=coords[i]*coords[i];
		}
		length=(float)Math.sqrt(length);
		for(int i=0;i<coords.length;i++){
			coords[i]/=length;
		}
	}
	
	public float[] cross_prod(float[] u,float[] v){
		float[] temp={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
		get_norm_vec(temp);
		return temp;
	}

}
