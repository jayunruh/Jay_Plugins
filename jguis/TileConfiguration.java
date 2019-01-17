/**
 * License: GPL
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package jguis;

import ij.IJ;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;
import java.util.concurrent.ExecutionException;

import mpicbg.models.ErrorStatistic;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;
import mpicbg.models.PointMatch;
import mpicbg.models.Tile;
import mpicbg.models.TileUtil;
import mpicbg.util.RealSum;


/**
 * A configuration of tiles.
 *
 * Add all tiles that build a common interconnectivity graph to one
 * configuration, fix at least one of the tiles and optimize the configuration.
 *
 * @author Stephan Saalfeld <saalfelds@janelia.hhmi.org>
 */
public class TileConfiguration
{

	final static protected DecimalFormat decimalFormat = new DecimalFormat();
	final static protected DecimalFormatSymbols decimalFormatSymbols = new DecimalFormatSymbols();

	final protected HashSet< Tile< ? > > tiles = new HashSet< Tile< ? > >();
	final public HashSet< Tile< ? > > getTiles(){ return tiles; }

	final protected HashSet< Tile< ? > > fixedTiles = new HashSet< Tile< ? > >();
	final public HashSet< Tile< ? > > getFixedTiles(){ return fixedTiles; }

	protected double minError = Double.MAX_VALUE;
	final public double getMinError() {	return minError; }

	protected double maxError = 0.0;
	final public double getMaxError() { return maxError; }

	protected double error = Double.MAX_VALUE;
	final public double getError() { return error; }

	public TileConfiguration()
	{
		decimalFormatSymbols.setGroupingSeparator( ',' );
		decimalFormatSymbols.setDecimalSeparator( '.' );
		decimalFormat.setDecimalFormatSymbols( decimalFormatSymbols );
		decimalFormat.setMaximumFractionDigits( 3 );
		decimalFormat.setMinimumFractionDigits( 3 );
	}

	protected void println( final String s ){ System.out.println( s ); }

	/**
	 * Cleanup.
	 */
	public void clear()
	{
		tiles.clear();
		fixedTiles.clear();

		minError = Double.MAX_VALUE;
		maxError = 0.0;
		error = Double.MAX_VALUE;
	}

	/**
	 * Add a single {@link Tile}.
	 *
	 * @param t
	 */
	public void addTile( final Tile< ? > t ){ tiles.add( t ); }

	/**
	 * Add a {@link Collection} of {@link Tile Tiles}.
	 *
	 * @param t
	 */
	public void addTiles( final Collection< ? extends Tile< ? > > t ){ tiles.addAll( t ); }

	/**
	 * Add all {@link Tile Tiles} of another {@link TileConfiguration}.
	 *
	 * @param t
	 */
	public void addTiles( final TileConfiguration t ){ tiles.addAll( t.tiles ); }

	/**
	 * Fix a single {@link Tile}.
	 *
	 * @param t
	 */
	public void fixTile( final Tile< ? > t ){ fixedTiles.add( t ); }

	/**
	 * Apply the model of each {@link Tile} to all its
	 * {@link PointMatch PointMatches}.
	 */
	protected void apply()
	{
		for ( final Tile< ? > t : tiles )
			t.apply();
	}

	/**
	 * Estimate min/max/average displacement of all
	 * {@link PointMatch PointMatches} in all {@link Tile Tiles}.
	 */
	protected void updateErrors()
	{
		double cd = 0.0;
		minError = Double.MAX_VALUE;
		maxError = 0.0;
		for ( final Tile< ? > t : tiles )
		{
			t.updateCost();
			final double d = t.getDistance();
			if ( d < minError ) minError = d;
			if ( d > maxError ) maxError = d;
			cd += d;
		}
		cd /= tiles.size();
		error = cd;
	}

	/**
	 * Update all {@link PointMatch Correspondences} in all {@link Tile Tiles}
	 * and estimate the average displacement.
	 */
	protected void update()
	{
		double cd = 0.0;
		minError = Double.MAX_VALUE;
		maxError = 0.0;
		for ( final Tile< ? > t : tiles )
		{
			t.update();
			final double d = t.getDistance();
			if ( d < minError ) minError = d;
			if ( d > maxError ) maxError = d;
			cd += d;
		}
		cd /= tiles.size();
		error = cd;

	}

	/**
	 * Minimize the displacement of all {@link PointMatch Correspondence pairs}
	 * of all {@link Tile Tiles}
	 *
	 * @param maxAllowedError do not accept convergence if error is > max_error
	 * @param maxIterations stop after that many iterations even if there was
	 *   no minimum found
	 * @param maxPlateauwidth convergence is reached if the average absolute
	 *   slope in an interval of this size and half this size is smaller than
	 *   0.0001 (in double accuracy).  This is assumed to prevent the algorithm
	 *   from stopping at plateaus smaller than this value.
	 */
	public void optimizeSilently(
			final ErrorStatistic observer,
			final double maxAllowedError,
			final int maxIterations,
			final int maxPlateauwidth ) throws NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		int i = 0;
		boolean proceed = i < maxIterations;
		/* initialize the configuration with the current model of each tile */
		apply();
		while ( proceed ){
			for ( final Tile< ? > tile : tiles ){
				if ( fixedTiles.contains( tile ) ) continue;
				tile.fitModel();
				tile.apply();
			}
			updateErrors();
			observer.add( error );
			//IJ.log(""+error);
			if ( i > maxPlateauwidth ){
				proceed = error > maxAllowedError;
				int d = maxPlateauwidth;
				while ( !proceed && d >= 1 ){
					try{
						proceed |= Math.abs( observer.getWideSlope( d ) ) > 0.0001;
					}
					catch ( final Exception e ) { e.printStackTrace(); }
					d /= 2;
				}
			}
			IJ.log(""+i+" , "+observer.mean+" , "+observer.max);
			proceed &= ++i < maxIterations;
		}
	}


	/**
	 * Computes a pre-alignemnt of all non-fixed {@link Tile}s by propagating the pairwise
	 * models. This does not give a correct registration but a very good starting point
	 * for the global optimization. This is necessary for models where the global optimization
	 * is not guaranteed to converge like the {@link HomographyModel2D}, {@link RigidModel3D}, ...
	 *
	 * @return - a list of {@link Tile}s that could not be pre-aligned
	 * @throws NotEnoughDataPointsException
	 * @throws {@link IllDefinedDataPointsException}
	 */
	public List< Tile< ? > > preAlign() throws NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		// first get order all tiles by
		// a) unaligned
		// b) aligned - which initially only contains the fixed ones
		final ArrayList< Tile< ? > > unAlignedTiles = new ArrayList< Tile< ? > >();
		final ArrayList< Tile< ? > > alignedTiles = new ArrayList< Tile< ? > >();

		// if no tile is fixed, take another */
		if ( getFixedTiles().size() == 0 )
		{
			final Iterator< Tile< ? > > it = this.getTiles().iterator();
			alignedTiles.add( it.next() );
			while ( it.hasNext() )
				unAlignedTiles.add( it.next() );
		}
		else
		{
			for ( final Tile< ? > tile : this.getTiles() )
			{
				if ( this.getFixedTiles().contains( tile ) )
					alignedTiles.add( tile );
				else
					unAlignedTiles.add( tile );
			}
		}

		// we go through each fixed/aligned tile and try to find a pre-alignment
		// for all other unaligned tiles
		for ( final ListIterator< Tile< ?> > referenceIterator = alignedTiles.listIterator(); referenceIterator.hasNext(); )
		{
			// once all tiles are aligned we can quit this loop
			if ( unAlignedTiles.size() == 0 )
				break;

			// get the next reference tile (either a fixed or an already aligned one
			final Tile< ? > referenceTile = referenceIterator.next();

			// transform all reference points into the reference coordinate system
			// so that we get the direct model even if we are not anymore at the
			// level of the fixed tile
			referenceTile.apply();

			// now we go through the unaligned tiles to see if we can align it to the current reference tile one
			for ( final ListIterator< Tile< ?> > targetIterator = unAlignedTiles.listIterator(); targetIterator.hasNext(); )
			{
				// get the tile that we want to preregister
				final Tile< ? > targetTile = targetIterator.next();

				// target tile is connected to reference tile
				if ( referenceTile.getConnectedTiles().contains( targetTile ) )
				{
					// extract all PointMatches between reference and target tile and fit a model only on these
					final ArrayList< PointMatch > pm = getConnectingPointMatches( targetTile, referenceTile );

					// are there enough matches?
					if ( pm.size() > targetTile.getModel().getMinNumMatches() )
					{
						// fit the model of the targetTile to the subset of matches
						// mapping its local coordinates target.p.l into the world
						// coordinates reference.p.w
						// this will give us an approximation for the global optimization
						targetTile.getModel().fit( pm );

						// now that we managed to fit the model we remove the
						// Tile from unaligned tiles and add it to aligned tiles
						targetIterator.remove();

						// now add the aligned target tile to the end of the reference list
						int countFwd = 0;

						while ( referenceIterator.hasNext() )
						{
							referenceIterator.next();
							++countFwd;
						}
						referenceIterator.add( targetTile );

						// move back to the current position
						// (+1 because it add just behind the current position)
						for ( int j = 0; j < countFwd + 1; ++j )
							referenceIterator.previous();
					}
				}

			}
		}

		return unAlignedTiles;
	}

	/**
	 * Returns an {@link ArrayList} of {@link PointMatch} that connect the targetTile and the referenceTile. The order of the
	 * {@link PointMatch} is PointMatch.p1 = target, PointMatch.p2 = reference. A {@link Model}.fit() will then solve the fit
	 * so that target.p1.l is mapped to reference.p2.w.
	 *
	 * @param targetTile - the {@link Tile} for which a {@link Model} can fit
	 * @param referenceTile - the {@link Tile} to which target will map
	 *
	 * @return - an {@link ArrayList} of all {@link PointMatch} that target and reference share
	 */
	public ArrayList<PointMatch> getConnectingPointMatches( final Tile<?> targetTile, final Tile<?> referenceTile )
	{
		final Set< PointMatch > referenceMatches = referenceTile.getMatches();
		final ArrayList< Point > referencePoints = new ArrayList<Point>( referenceMatches.size() );

		// add all points from the reference tile so that we can search for them
		for ( final PointMatch pm : referenceMatches )
			referencePoints.add( pm.getP1() );

		// the result arraylist containing only the pointmatches from the target file
		final ArrayList< PointMatch > connectedPointMatches = new ArrayList<PointMatch>();

		// look for all PointMatches where targetTile.PointMatch.Point2 == referenceTile.PointMatch.Point1
		// i.e. a PointMatch of the target tile that links a Point which is part of reference tile
		for ( final PointMatch pm : targetTile.getMatches() )
			if ( referencePoints.contains( pm.getP2() ) )
				connectedPointMatches.add( pm );

		return connectedPointMatches;
	}

}
