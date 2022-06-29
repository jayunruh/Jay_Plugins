/*******************************************************************************
 * Copyright (c) 2014 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Frame;
import ij.plugin.*;
import ij.text.*;
import java.util.*;
import jguis.*;

public class close_all_tables_jru_v1 implements PlugIn {

	public void run(String arg) {
		Object[] windowlist=jutils.getTableWindowList(false);
		Frame[] frames=(Frame[])windowlist[0];
		for(int i=0;i<frames.length;i++){
			((TextWindow)frames[i]).close();
		}
	}

}
