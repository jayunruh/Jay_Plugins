package jguis;

import java.util.*;

import Glacier2.CannotCreateSessionException;
import Glacier2.PermissionDeniedException;
import ij.IJ;
import omero.*;
import omero.api.*;
import omero.model.*;
import omero.sys.ParametersI;
import pojos.DatasetData;
import pojos.ImageData;

public class Ome_file_reader{
	private client client;
	private client unsecureClient;
	private ServiceFactoryPrx entry;
	
	public static void main(String[] args){
		new Ome_file_reader();
	}
	
	public Ome_file_reader(){
    	client=new client("columbus01.sgc.loc",4063);
		//unsecureClient=client.createClient(false);
		try{
			//entry=client.createSession("HDAC5","Parama");
			entry=client.createSession("Polyploidy","tamara");
			System.out.println("session created");
			IContainerPrx proxy=entry.getContainerService();
			ParametersI param = new ParametersI();
			long userId = entry.getAdminService().getEventContext().userId;
			System.out.println("user id = "+userId);
			param.exp(omero.rtypes.rlong(userId));

			//indicate to load the images
			param.leaves();
			//List<IObject> results = proxy.loadContainerHierarchy(Dataset.class.getName(), new ArrayList<Long>(), param);
			//List<IObject> results = proxy.loadContainerHierarchy(Screen.class.getName(), new ArrayList<Long>(), param);
			//List<IObject> results = proxy.loadContainerHierarchy(Project.class.getName(), new ArrayList<Long>(), param);
			//List<Image> results=proxy.getImages(Image.class.getName(),Arrays.asList(imageId),new ParametersI());
			GregorianCalendar gc = new GregorianCalendar();
	    	gc = new GregorianCalendar(gc.get(Calendar.YEAR), 
					gc.get(Calendar.MONTH), 
					gc.get(Calendar.DAY_OF_MONTH), 0, 0, 0);
			long startTime=gc.getTime().getTime();
			System.out.println("start time = "+startTime);
			ParametersI po=new ParametersI();
			po.leaves();
			//start time is milliseconds since Jan 1 1970
			po.startTime(omero.rtypes.rtime(startTime-1L));
			List<Image> results=proxy.getImagesByOptions(po);
			//You can directly interact with the IObject or the Pojos object.
			//Follow interaction with the Pojos.
			//Iterator<IObject> i = results.iterator();
			System.out.println("results length = "+results.size());
			for(int i=0;i<results.size();i++){
				Image temp=results.get(i);
				System.out.println(temp.getName().getValue());
			}
			/*DatasetData dataset;
			Set<ImageData> images;
			Iterator<ImageData> j;
			ImageData image;
			int counter=0;
			while (i.hasNext()) {
			  dataset = new DatasetData((Dataset) i.next());
			  images = dataset.getImages();
			  System.out.println("dataset "+counter+" length = "+images.size());
			  j = images.iterator();
			  int counter2=0;
			  while (j.hasNext()) {
				 
			    image = j.next();
			    System.out.println("image "+counter2+" "+image.getId());
			    counter2++;
			  }
			  counter++;
			}*/
		}catch(CannotCreateSessionException e){
			System.out.println(e.toString());
		}catch(PermissionDeniedException e){
			System.out.println(e.toString());
		}catch(ServerError e){
			System.out.println(e.toString());
		}
		dispose();
	}
	
	public void dispose(){
		client.closeSession();
	}
	
}
