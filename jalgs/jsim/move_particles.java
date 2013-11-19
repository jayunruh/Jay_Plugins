/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jsim;

public class move_particles implements Cloneable{
	rngs random;
	public int confineindex,boundaryindex;
	public float D2,fboxpixels,yflowrate2,gridsizepixels,jumpprob,dfraction;
	public boolean restrict;
	public sim_setup set;

	public move_particles(sim_species ss1,sim_setup set1){
		set=set1;
		random=set1.random;
		confineindex=set1.confineindex;
		boundaryindex=set1.boundaryindex;
		D2=ss1.D2;
		fboxpixels=set1.fboxpixels;
		yflowrate2=set1.yflowrate2;
		gridsizepixels=set1.gridsizepixels;
		jumpprob=set1.jumpprob;
		dfraction=set1.dfraction;
		restrict=set1.restrict;
	}

	public Object Clone() throws CloneNotSupportedException{
		return super.clone();
	}

	public void do_move_particle(sim_particle sp){
		// now move the particle according to its identity
		// confine index of 1 is xy, 2 is xz, 3 is x line
		float tempd=D2;
		if(restrict&&sp.stuck){
			tempd*=dfraction;
		}
		float xdev=yflowrate2+(float)random.gasdev(0.0,Math.sqrt(2.0*tempd));
		sp.coords[0]+=xdev;
		float ydev=0.0f;
		if(confineindex<2){
			ydev=(float)random.gasdev(0.0,Math.sqrt(2.0*tempd));
			sp.coords[1]+=ydev;
		}
		float zdev=0.0f;
		if(confineindex!=1&&confineindex!=3){
			zdev=(float)random.gasdev(0.0,Math.sqrt(2.0*tempd));
			sp.coords[2]+=zdev;
		}

		// now implement the periodic boundary conditions
		if(boundaryindex==0){
			// here we have periodic boundaries
			if(sp.coords[0]>fboxpixels){
				sp.coords[0]-=fboxpixels;
			}
			if(sp.coords[0]<0.0f){
				sp.coords[0]+=fboxpixels;
			}
			if(confineindex!=2){
				if(sp.coords[1]>fboxpixels){
					sp.coords[1]-=fboxpixels;
				}
				if(sp.coords[1]<0.0f){
					sp.coords[1]+=fboxpixels;
				}
			}
			if(confineindex!=1){
				if(sp.coords[2]>fboxpixels){
					sp.coords[2]-=fboxpixels;
				}
				if(sp.coords[2]<0.0f){
					sp.coords[2]+=fboxpixels;
				}
			}
		}else{
			if(boundaryindex==1){
				// here we have reflective boundaries
				if(sp.coords[0]>fboxpixels){
					sp.coords[0]=2.0f*fboxpixels-sp.coords[0];
				}
				if(sp.coords[0]<0.0f){
					sp.coords[0]=-sp.coords[0];
				}
				if(confineindex!=2){
					if(sp.coords[1]>fboxpixels){
						sp.coords[1]=2.0f*fboxpixels-sp.coords[1];
					}
					if(sp.coords[1]<0.0f){
						sp.coords[1]=-sp.coords[1];
					}
				}
				if(confineindex!=1){
					if(sp.coords[2]>=fboxpixels){
						sp.coords[2]=2.0f*fboxpixels-sp.coords[2];
					}
					if(sp.coords[2]<0.0f){
						sp.coords[2]=-sp.coords[2];
					}
				}
			}else{
				// here we create a new particle
				float[] tempcoords={sp.coords[0],sp.coords[1],sp.coords[2]};
				float[] dev={xdev,ydev,zdev};
				boolean temp=false;
				boolean[] outs=new boolean[3];
				boolean[] ins=new boolean[3];
				if(sp.coords[0]>fboxpixels){
					outs[0]=true;
					temp=true;
				}
				if(sp.coords[0]<0.0f){
					temp=true;
					ins[0]=true;
				}
				if(confineindex!=2){
					if(sp.coords[1]>fboxpixels){
						outs[1]=true;
						temp=true;
					}
					if(sp.coords[1]<0.0f){
						temp=true;
						ins[1]=true;
					}
				}
				if(confineindex!=1){
					if(sp.coords[2]>=fboxpixels){
						outs[2]=true;
						temp=true;
					}
					if(sp.coords[2]<0.0f){
						temp=true;
						ins[2]=true;
					}
				}
				if(temp){
					sp.coords=generate_new_particle(tempcoords,dev,outs,ins);
				}
			}
		}
		// finally implement the spatial restrictions
		if(restrict){
			handle_stuck(sp,xdev,ydev);
		}
	}

	private float[] generate_new_particle(float[] coords,float[] dev,boolean[] outs,boolean[] ins){
		float[] fcoords={coords[0]-dev[0],coords[1]-dev[1],coords[2]-dev[2]};
		// find out how far the particle was from each edge (before it went out)
		float distance=0.0f;
		if(outs[0]){
			distance=fboxpixels-fcoords[0];
		}else{
			if(ins[1]){
				distance=fcoords[1];
			}else{
				if(ins[0]){
					distance=fcoords[0];
				}else{
					if(outs[1]){
						distance=fboxpixels-fcoords[1];
					}else{
						if(ins[2]){
							distance=fcoords[2];
						}else{
							distance=fboxpixels-fcoords[2];
						}
					}
				}
			}
		}
		// now generate a new particle with that distance from the edge and move
		// it one step
		// repeat until it stays inside
		do{
			int temp=0;
			if(confineindex==0){
				temp=(int)random.unidev(5.99999,0.0);
			}
			if(confineindex==1){
				temp=(int)random.unidev(3.99999,0.0);
			}
			if(confineindex==2){
				temp=(int)random.unidev(5.99999,2.0);
				if(temp<4){
					temp-=2;
				}
			}
			if(temp==0){
				fcoords[0]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
				if(confineindex!=2){
					fcoords[1]=distance;
				}
				if(confineindex!=1){
					fcoords[2]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
				}
			}else{
				if(temp==1){
					fcoords[0]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
					if(confineindex!=2){
						fcoords[1]=fboxpixels-distance;
					}
					if(confineindex!=1){
						fcoords[2]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
					}
				}else{
					if(temp==2){
						fcoords[1]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
						fcoords[0]=distance;
						if(confineindex!=1){
							fcoords[2]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
						}
					}else{
						if(temp==3){
							fcoords[1]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
							fcoords[0]=fboxpixels-distance;
							if(confineindex!=1){
								fcoords[2]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
							}
						}else{
							if(temp==4){
								fcoords[0]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
								if(confineindex!=2){
									fcoords[1]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
								}
								fcoords[2]=distance;
							}else{
								fcoords[0]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
								if(confineindex!=2){
									fcoords[1]=(float)random.unidev((double)(fboxpixels-distance),(double)distance);
								}
								fcoords[2]=fboxpixels-distance;
							}
						}
					}
				}
			}
			fcoords[0]+=yflowrate2+(float)random.gasdev(0.0,Math.sqrt(2.0*D2));
			if(confineindex!=2){
				fcoords[1]+=(float)random.gasdev(0.0,Math.sqrt(2.0*D2));
			}
			if(confineindex!=1){
				fcoords[2]+=(float)random.gasdev(0.0,Math.sqrt(2.0*D2));
			}
		}while((fcoords[0]<fboxpixels&&fcoords[0]>0.0f)&&((fcoords[1]<fboxpixels&&fcoords[1]>=0.0f)&&(fcoords[2]<fboxpixels&&fcoords[2]>=0.0f)));
		return fcoords;
	}

	private void handle_stuck(sim_particle sp,float xdev,float ydev){
		boolean newstuck=sp.stuck;
		boolean oldstuck=sp.stuck;
		float fx=sp.coords[0]-xdev;
		float fy=sp.coords[1]-ydev;
		int gridunitx=(int)(sp.coords[0]/gridsizepixels);
		int fgridunitx=(int)(fx/gridsizepixels);
		int gridunity=(int)(sp.coords[1]/gridsizepixels);
		int fgridunity=(int)(fy/gridsizepixels);
		if(confineindex<2){
			if((gridunitx!=fgridunitx&&gridunity==fgridunity)||(gridunitx==fgridunitx&&gridunity!=fgridunity)){
				// we crossed a boundary--check if that was okay--if not,
				// reflect
				newstuck=!oldstuck;
				if(jumpprob==0.0f||(float)random.unidev(1.0,0.0)>jumpprob){
					if(fgridunitx!=gridunitx){
						if(gridunitx>fgridunitx){
							sp.coords[0]=2.0f*(gridsizepixels*(float)gridunitx)-sp.coords[0];
						}else{
							sp.coords[0]=2.0f*(gridsizepixels*(float)fgridunitx)-sp.coords[0];
						}
					}else{
						if(gridunity>fgridunity){
							sp.coords[1]=2.0f*(gridsizepixels*(float)gridunity)-sp.coords[1];
						}else{
							sp.coords[1]=2.0f*(gridsizepixels*(float)fgridunity)-sp.coords[1];
						}
					}
					newstuck=oldstuck;
				}
			}
		}else{
			if(gridunitx!=fgridunitx){
				// we crossed a boundary--check if that was okay--if not,
				// reflect
				newstuck=!oldstuck;
				if(jumpprob==0.0f||(float)random.unidev(1.0,0.0)>jumpprob){
					if(gridunitx>fgridunitx){
						sp.coords[0]=2.0f*(gridsizepixels*(float)gridunitx)-sp.coords[0];
					}else{
						sp.coords[0]=2.0f*(gridsizepixels*(float)fgridunitx)-sp.coords[0];
					}
					newstuck=oldstuck;
				}
			}
		}
		sp.stuck=newstuck;
	}

}
