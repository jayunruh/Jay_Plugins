package jguis;

import ij.IJ;
import j3D.element3D;
import j3D.line3D;
import j3D.renderer;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.jocl.*;

import jalgs.*;

public class opencl_maxproj{
	boolean valid=false;
	public int slices,projtype,newsize,width,height,ch,xrem,yrem,zrem;
	public float zratio,thresh;
	public line3D[] lines;
	public element3D[] elements;
	public renderer r;
	public cl_mem[] memObjects;
	public cl_command_queue commandQueue;
	public cl_kernel kernel;
	public cl_context context;
	public cl_program program;
	
	public opencl_maxproj(Object[] image,int width,int height,float zratio,int newsize,int chs,int ch,float thresh,int projtype,line3D[] lines,element3D[] elements,renderer r){
		//String javaLibPath = System.getProperty("java.library.path");
  	  //System.out.println(javaLibPath);
		this.projtype=projtype;
		this.newsize=newsize;
		this.width=width;
		this.height=height;
		this.zratio=zratio;
		this.ch=ch;
		this.thresh=thresh;
        slices=image.length/chs;
        
        this.lines=lines;
        this.elements=elements;
        this.r=r;
		
        String programSource="";
        jdataio jdio=new jdataio();
        BufferedReader b;
        try{
        	b=new BufferedReader(new FileReader("maxproj2.cl"));
        	programSource=jdio.readstringfile(b);
        	//IJ.log(programSource);
        	b.close();
        }catch(IOException e){
        	e.printStackTrace();
        }

        if(programSource==null) return;
        long numBytes[] = new long[1];

        //float[] maxproj=new float[newsize*newsize];
        //float[] depth=new float[newsize*newsize];

        // Obtain the platform IDs and initialize the context properties
        IJ.log("Obtaining platform...");
        cl_platform_id platforms[] = new cl_platform_id[1];
        CL.clGetPlatformIDs(platforms.length, platforms, null);
        cl_context_properties contextProperties = new cl_context_properties();
        contextProperties.addProperty(CL.CL_CONTEXT_PLATFORM, platforms[0]);
        
        // Create an OpenCL context on a GPU device
        context = CL.clCreateContextFromType(contextProperties, CL.CL_DEVICE_TYPE_GPU, null, null, null);
        if (context == null){
            // If no context for a GPU device could be created,
            // try to create one for a CPU device.
            context = CL.clCreateContextFromType(contextProperties, CL.CL_DEVICE_TYPE_CPU, null, null, null);
            if (context == null){
          	  IJ.log("Unable to create a context");
                return;
            }
        }

        // Enable exceptions and subsequently omit error checks in this sample
        CL.setExceptionsEnabled(true);
        
        // Get the list of GPU devices associated with the context
        CL.clGetContextInfo(context, CL.CL_CONTEXT_DEVICES, 0, null, numBytes); 
        
        // Obtain the cl_device_id for the first device
        int numDevices = (int) numBytes[0] / Sizeof.cl_device_id;
        cl_device_id devices[] = new cl_device_id[numDevices];
        /*for(int i=0;i<devices.length;i++){
        	IJ.log(devices[i].toString());
        }*/
        CL.clGetContextInfo(context, CL.CL_CONTEXT_DEVICES, numBytes[0],Pointer.to(devices), null);

        // Create a command-queue
        commandQueue = CL.clCreateCommandQueue(context, devices[0], 0, null);

        // Allocate the memory objects for the input- and output data
        memObjects = new cl_mem[5];
        //memObjects[0] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,Sizeof.cl_float * image2.length, inputarr, null);
        memObjects[0]=CL.clCreateBuffer(context,CL.CL_MEM_READ_ONLY,Sizeof.cl_float*width*height*image.length,null,null);
        memObjects[2]=CL.clCreateBuffer(context,CL.CL_MEM_READ_WRITE,Sizeof.cl_float*newsize*newsize*6,null,null);
        memObjects[1] = CL.clCreateBuffer(context, CL.CL_MEM_READ_WRITE,Sizeof.cl_float * newsize*newsize, null, null);
        memObjects[3] = CL.clCreateBuffer(context, CL.CL_MEM_READ_WRITE,Sizeof.cl_float * newsize*newsize, null, null);
        float[] limits={0.0f,(float)(width-1),0.0f,(float)(height-1),0.0f,(float)(slices-1)};
        memObjects[4] = CL.clCreateBuffer(context, CL.CL_MEM_READ_WRITE,Sizeof.cl_float * 6, null, null);
        
        // Create the program from the source code
        program = CL.clCreateProgramWithSource(context,1, new String[]{ programSource }, null, null);
        
        // Build the program
        CL.clBuildProgram(program, 0, null, null, null, null);
        
        // Create the kernel
        kernel = CL.clCreateKernel(program, "maxproj2", null);
        
        // Set the arguments for the kernel
        //loop over the input array dimensions
        for(int i=0;i<image.length;i++){
      	  CL.clEnqueueWriteBuffer(commandQueue,memObjects[0],CL.CL_TRUE,i*width*height*Sizeof.cl_float,width*height*Sizeof.cl_float,Pointer.to(algutils.convert_arr_float2(image[i])),0,null,null);
        }
        //new enqueue the line data
        //make sure the line data is not null before running the kernel
        if(lines!=null) CL.clEnqueueWriteBuffer(commandQueue,memObjects[2],CL.CL_TRUE,0,newsize*newsize*6*Sizeof.cl_float,Pointer.to(lines2array(lines)),0,null,null);
        //and the limits data
        CL.clEnqueueWriteBuffer(commandQueue,memObjects[4],CL.CL_TRUE,0,6*Sizeof.cl_float,Pointer.to(limits),0,null,null);
        //and the limit data
        CL.clSetKernelArg(kernel, 0,Sizeof.cl_mem, Pointer.to(memObjects[0]));
        CL.clSetKernelArg(kernel, 1,Sizeof.cl_mem, Pointer.to(memObjects[2]));
        CL.clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(memObjects[1]));
        CL.clSetKernelArg(kernel, 3, Sizeof.cl_mem, Pointer.to(memObjects[3]));
        CL.clSetKernelArg(kernel, 4, Sizeof.cl_mem, Pointer.to(memObjects[4]));
        CL.clSetKernelArg(kernel, 5, Sizeof.cl_int, Pointer.to(new int[]{width}));
        CL.clSetKernelArg(kernel, 6, Sizeof.cl_int, Pointer.to(new int[]{height}));
        CL.clSetKernelArg(kernel, 7, Sizeof.cl_int, Pointer.to(new int[]{slices}));
        CL.clSetKernelArg(kernel, 8, Sizeof.cl_int, Pointer.to(new int[]{chs}));
        CL.clSetKernelArg(kernel, 9, Sizeof.cl_int, Pointer.to(new int[]{ch}));
        CL.clSetKernelArg(kernel, 10, Sizeof.cl_int, Pointer.to(new int[]{projtype}));
        CL.clSetKernelArg(kernel, 11, Sizeof.cl_float, Pointer.to(new float[]{zratio}));
        CL.clSetKernelArg(kernel, 12, Sizeof.cl_float, Pointer.to(new float[]{thresh}));
        CL.clSetKernelArg(kernel, 13, Sizeof.cl_int, Pointer.to(new int[]{newsize}));
        CL.clSetKernelArg(kernel, 14, Sizeof.cl_int, Pointer.to(new int[]{0}));
	}
	
	public opencl_maxproj(Object[] image,int width,int height,float zratio,int newsize,int chs,int ch,float thresh,int projtype){
		//String javaLibPath = System.getProperty("java.library.path");
  	  //System.out.println(javaLibPath);
		this(image,width,height,zratio,newsize,chs,ch,thresh,projtype,null,null,null);
        lines=new line3D[newsize*newsize];
        elements=new element3D[newsize*newsize];
        xrem=(int)(0.5f*(newsize-width));
        yrem=(int)(0.5f*(newsize-height));
        float zsize=zratio*(float)slices;
        zrem=(int)(0.5f*(newsize-(int)zsize));
        for(int i=0;i<newsize;i++){
      	  for(int j=0;j<newsize;j++){
      		  lines[j+i*newsize]=new line3D(j-xrem,i-yrem,-zrem,j-xrem,i-yrem, newsize-1-zrem,null);
      		  elements[j+i*newsize]=lines[j+i*newsize];
      	  }
        }
        r=new renderer(newsize,newsize);
		r.centerz=newsize/2-zrem;
		r.centerx=newsize/2-xrem;
		r.centery=newsize/2-yrem;
		r.setelementarray(elements);
		setLines(lines);
	}
	
	public void setRotation(float rx,float ry,float rz){
		//units are degrees
		r.setrotation(rx,ry,rz);
		CL.clEnqueueWriteBuffer(commandQueue,memObjects[2],CL.CL_TRUE,0,newsize*newsize*6*Sizeof.cl_float,Pointer.to(lines2array(lines)),0,null,null);
	}
	
	public void setLines(line3D[] lines2){
		setLines(lines2,0.0f,0.0f,0.0f);
	}
	
	public void setLines(line3D[] lines2,float xshift,float yshift,float zshift){
		CL.clEnqueueWriteBuffer(commandQueue,memObjects[2],CL.CL_TRUE,0,newsize*newsize*6*Sizeof.cl_float,Pointer.to(lines2array(lines2,xshift,yshift,zshift)),0,null,null);
	}
	
	public void setLimits(float[] limits){
		//limits are xmin,xmax,ymin,ymax,zmin,zmax
		float[] lim2=limits.clone();
		if(lim2[0]<0.0f) lim2[0]=0.0f;
		if(lim2[2]<0.0f) lim2[2]=0.0f;
		if(lim2[4]<0.0f) lim2[4]=0.0f;
		if(lim2[1]>=(float)(width-1)) lim2[1]=(float)(width-1)-0.001f;
		if(lim2[3]>=(float)(height-1)) lim2[3]=(float)(height-1)-0.001f;
		if(lim2[5]>=(float)(slices-1)) lim2[5]=(float)(slices-1)-0.001f;
		CL.clEnqueueWriteBuffer(commandQueue,memObjects[4],CL.CL_TRUE,0,6*Sizeof.cl_float,Pointer.to(limits),0,null,null);
	}
	
	public void setProjType(int projtype){
		this.projtype=projtype;
		CL.clSetKernelArg(kernel, 10, Sizeof.cl_int, Pointer.to(new int[]{projtype}));
	}
	
	public void setThresh(float thresh){
		this.thresh=thresh;
		CL.clSetKernelArg(kernel, 12, Sizeof.cl_float, Pointer.to(new float[]{thresh}));
	}
	
	public void setChannel(int ch){
		this.ch=ch;
        CL.clSetKernelArg(kernel, 9, Sizeof.cl_int, Pointer.to(new int[]{ch}));
	}
	
	public float[] createProjection(){
		 long global_work_size[] = new long[]{newsize*newsize};
		 float[] maxproj=new float[newsize*newsize];
		 CL.clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, global_work_size, null, 0, null, null);
		 if(projtype==3){
        	//use the kernel with the calcsurf value set to 1 to apply the reflectivity
            CL.clSetKernelArg(kernel, 14, Sizeof.cl_int, Pointer.to(new int[]{1}));
            CL.clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, global_work_size, null, 0, null, null);
            CL.clSetKernelArg(kernel, 14, Sizeof.cl_int, Pointer.to(new int[]{0}));
            CL.clEnqueueReadBuffer(commandQueue, memObjects[1], CL.CL_TRUE, 0,newsize*newsize*Sizeof.cl_float, Pointer.to(maxproj), 0, null, null);
		 } else {
    	  	CL.clEnqueueReadBuffer(commandQueue, memObjects[1], CL.CL_TRUE, 0,newsize*newsize*Sizeof.cl_float, Pointer.to(maxproj), 0, null, null);
		 }
		 return maxproj;
	}
	
	public float[] getDepth(){
		float[] depth=new float[newsize*newsize];
		CL.clEnqueueReadBuffer(commandQueue, memObjects[3], CL.CL_TRUE, 0,newsize*newsize*Sizeof.cl_float, Pointer.to(depth), 0, null, null);
		return depth;
	}
	
	public void dispose(){
        // Release kernel, program, and memory objects
        CL.clReleaseMemObject(memObjects[0]);
        CL.clReleaseMemObject(memObjects[1]);
        CL.clReleaseMemObject(memObjects[2]);
        CL.clReleaseMemObject(memObjects[3]);
        CL.clReleaseMemObject(memObjects[4]);
        CL.clReleaseKernel(kernel);
        CL.clReleaseProgram(program);
        CL.clReleaseCommandQueue(commandQueue);
        CL.clReleaseContext(context);
	}
	
    public static float[] lines2array(line3D[] lines){
  	  float[] arr=new float[lines.length*6];
  	  for(int i=0;i<lines.length;i++){
  		  int off=i*6;
  		  arr[off]=lines[i].pt1.rx;
  		  arr[off+1]=lines[i].pt1.ry;
  		  arr[off+2]=lines[i].pt1.rz;
  		  arr[off+3]=lines[i].pt2.rx;
  		  arr[off+4]=lines[i].pt2.ry;
  		  arr[off+5]=lines[i].pt2.rz;
  	  }
  	  return arr;
    }
    
    public static float[] lines2array(line3D[] lines,float xshift,float yshift,float zshift){
    	  float[] arr=new float[lines.length*6];
    	  for(int i=0;i<lines.length;i++){
    		  int off=i*6;
    		  arr[off]=lines[i].pt1.rx-xshift;
    		  arr[off+1]=lines[i].pt1.ry-yshift;
    		  arr[off+2]=lines[i].pt1.rz-zshift;
    		  arr[off+3]=lines[i].pt2.rx-xshift;
    		  arr[off+4]=lines[i].pt2.ry-yshift;
    		  arr[off+5]=lines[i].pt2.rz-zshift;
    	  }
    	  return arr;
      }
}
