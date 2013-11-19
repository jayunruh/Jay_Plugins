/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs.jseg;

import java.util.List;

public interface track_interface{

	public int getNFrames();

	public List<float[]> getNextFrameParams();

	public void put_assignments(int[] assignments);

	public void show_progress(int progress,int nframes);

}
