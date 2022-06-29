package jguis;

import ij.IJ;

public class ijpluginrunner{
	ClassLoader ijcl;
	Object plugin;
	
	public ijpluginrunner(String classname){
		ijcl=IJ.getClassLoader();
		plugin=loadClass(classname);
	}
	
	public Object getPlugin(){
		return plugin;
	}
	
	public Object loadClass(String classname){
		try{
			return ijcl.loadClass(classname).newInstance();
		} catch(Exception e){
			return null;
		}
	}
	
	public Object runPluginMethod(String method,Object[] args){
		//run the method using reflection
		return jutils.runReflectionMethod(plugin,method,args);
	}
	
	public static Object getClass(String classname){
		try{
			return IJ.getClassLoader().loadClass(classname).newInstance();
		} catch(Exception e){
			return null;
		}
	}

}
