package jguis;

import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.type.NativeType;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.NumericType;

public class fijiutils{
	//these are utilities that reference fiji and related libraries to keep them isolated from my jutils code

	/*************
	 * wrap an imageplus in an imglib2 img (RAI)
	 * @param imp
	 * @return
	 */
	public static <T extends NumericType<T> & NativeType<T>> Img<T> wrapImpImg(ImagePlus imp){
	//public static <T extends Type<T>> Img<T> wrapImpImg(ImagePlus imp){
		Img<T> wrapped=ImagePlusAdapter.wrap(imp);
		return wrapped;
	}
	
	//public static <T extends NumericType<T> & NativeType<T>> boolean copyRAI(final Img<T> source, final IterableInterval<T> dest){
	public static <T extends Type<T>> boolean copyRAI(final Img<T> source, final IterableInterval<T> dest){
		//imglib2 doesn't raster through images so have to create a localizing cursor and a random access point
		//Cursor< T > cursorSource = source.localizingCursor();
		RandomAccess< T > randomAccess=source.randomAccess();
	    //RandomAccess< T > randomAccess = dest.randomAccess();
	    Cursor<T> cursorDest=dest.localizingCursor();
	    //while(cursorSource.hasNext()){
	    	//cursorSource.fwd();
	    	//randomAccess.setPosition(cursorSource);
	    	//randomAccess.get().set(cursorSource.get());
	    	//cursorDest.
	    //}
	    while(cursorDest.hasNext()){
	    	cursorDest.fwd();
	    	randomAccess.setPosition(cursorDest);
	    	cursorDest.get().set(randomAccess.get());
	    }
	    return true;
	}

}
