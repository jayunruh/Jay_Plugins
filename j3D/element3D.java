/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package j3D;

import java.awt.*;

public abstract class element3D implements Cloneable{

	public Color color;

	public boolean thick;

	public abstract int getzpos();

	public abstract void translate(int transx,int transy,int transz);

	public abstract void moveto(int ptx,int pty,int ptz);

	public abstract void rotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz);

	public abstract void rotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz);

	public abstract void rotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz);
	
	public abstract void rotatecossinabout(double cosval1,double sinval1,double ux1,double uy1,double uz1,int centerx1,int centery1,int centerz1);

	public abstract void addrotatedeg(double degx1,double degy1,double degz1,int centerx,int centery,int centerz);

	public abstract void addrotaterad(double radx1,double rady1,double radz1,int centerx,int centery,int centerz);
	
	public abstract void addrotatecossin(double cosx1,double cosy1,double cosz1,double sinx1,double siny1,double sinz1,int centerx,int centery,int centerz);

	public abstract void transform_perspective(double horizon_dist,int centerx,int centery,int centerz);

	public abstract void transform_negative_perspective(double horizon_dist,int centerx,int centery,int centerz);
	
	public abstract void transform(double[][] transmat,int centerx,int centery,int centerz);

	public abstract void drawelement(Graphics g);
	
	public abstract element3D clone();

}
