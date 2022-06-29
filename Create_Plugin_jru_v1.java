/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import java.util.Calendar;
//import fiji.scripting.*;

public class Create_Plugin_jru_v1 implements PlugIn {

	public void run(String arg) {
		String[] plugintypes={"Blank","Image","Plot","Table","Macro","Plot3D"};
		GenericDialog gd=new GenericDialog("Options");
		gd.addChoice("Plugin Type",plugintypes,plugintypes[1]);
		gd.showDialog(); if(gd.wasCanceled()){return;}
		int pt=gd.getNextChoiceIndex();
		if(pt==4){
			Editor ed = (Editor)IJ.runPlugIn("ij.plugin.frame.Editor", "");
			if (ed==null) return;
			ed.create("My_Macro.ijm", "");
			return;
		}
		String importtext="";
		String codetext="\n";
		if(pt==1){
			codetext="\t\tImagePlus imp=WindowManager.getCurrentImage();\n";
			codetext+="\t\tint width=imp.getWidth(); int height=imp.getHeight();\n";
		}
		if(pt==2){
			codetext="\t\tImageWindow iw=WindowManager.getCurrentWindow();\n";
			codetext+="\t\tfloat[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,\"getXValues\");\n";
			codetext+="\t\tfloat[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,\"getYValues\");\n";
			codetext+="\t\tint[] npts=(int[])jutils.runPW4VoidMethod(iw,\"getNpts\");\n";
			importtext="import jguis.*;\n";
		}
		if(pt==3){
			codetext="\t\tTextWindow[] tw=jutils.selectTables(false,1,new String[]{\"Table\"});\n";
			codetext+="\t\tif(tw==null || tw.length<1) return;\n";
			codetext+="\t\tTextPanel tp=tw[0].getTextPanel();\n";
			codetext+="\t\tString[] col_labels=table_tools.getcollabels(tp);\n";
			codetext+="\t\tList<List<String>> listtable=table_tools.table2listtable(tp);\n";
			importtext="import ij.text.*;\nimport java.util.*;\nimport jguis.*;\n";
		}
		if(pt==5){
			codetext="\t\tImageWindow iw=WindowManager.getCurrentWindow();\n";
			codetext+="\t\tfloat[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,\"getXValues\");\n";
			codetext+="\t\tfloat[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,\"getYValues\");\n";
			codetext+="\t\tfloat[][][] zvals=(float[][][])jutils.runPW4VoidMethod(iw,\"getZValues\");\n";
			codetext+="\t\tint[][] npts=(int[][])jutils.runPW4VoidMethod(iw,\"getNpts\");\n";
			codetext+="\t\t//delete the following lines if not a trajectory\n";
			codetext+="\t\tfloat[][] tzvals=zvals[0];\n";
			codetext+="\t\tint[] tnpts=npts[0];\n";
			importtext="import jguis.*;\n";
		}
		String year=""+Calendar.getInstance().get(Calendar.YEAR);
		String pluginName = "My_Plugin.java";
		String className = "My_Plugin";
		String text = "";
		text+="/*******************************************************************************\n";
		text+=" * Copyright (c) "+year+" Jay Unruh, Stowers Institute for Medical Research.\n";
		text+=" * All rights reserved. This program and the accompanying materials\n";
		text+=" * are made available under the terms of the GNU Public License v2.0\n";
		text+=" * which accompanies this distribution, and is available at\n";
		text+=" * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html\n";
		text+=" ******************************************************************************/\n";
		text += "import ij.*;\n";
		text += "import ij.process.*;\n";
		text += "import ij.gui.*;\n";
		if(pt==3) text += "import java.awt.Frame;\n";
		else text += "import java.awt.*;\n";
		text += "import ij.plugin.*;\n";
		text+=importtext;
		text += "\n";
		text += "public class "+className+" implements PlugIn {\n";
		text += "\n";
		text += "\tpublic void run(String arg) {\n";
		text+=codetext;
		text += "\t}\n";
		text += "\n";
		text += "}\n";
		/*try{
			Class<?> temp=Class.forName("fiji.scripting.TextEditor");
			TextEditor te=Script_Editor.getInstance();
			if(te==null || !te.isVisible()){
				te=new TextEditor(pluginName,text);
			} else {
				te.createNewDocument();
				EditorPane ep=te.getEditorPane();
				ep.setText(text);
				ep.setFileName(pluginName);
				//ep.setLanguageByExtension(".java");
			}
		}catch (Exception e){*/
			Editor ed = (Editor)IJ.runPlugIn("ij.plugin.frame.Editor", "");
			if (ed==null) return;
			ed.create(pluginName, text);
		//}
	}

}
