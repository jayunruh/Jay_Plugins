Jay_Plugins
===========

ImageJ plugins from Jay

This code is licensed under the GPL.  Please read the license.txt file for more details.  I recommend using the Fiji update site or my website (http://research.stowers.org/imagejplugins/) for download and update.

This code is organized roughly according to function.  The top level files are the actual ImageJ plugins.  ImageJ only allows 100 plugins per jar file, so Jay_Plugins2 accomodates the overflow.

The jalgs package includes most of the algorithms and tries to avoid gui functions and even dependencies.  There are a few awt dependencies (for the Polygon class) as well as a few files which depend on the jama library.  Within that package, algutils contains lots of static utility methods, jdataio contains lots of io methods, and jstatistics is my master statistics class.  jfft contains my fft codes along with associated things like convolution and image filters.  jfit contains linear and non-linear least squares codes along with utilities for global fitting and error analysis by support plane or monte carlo analysis.  The jseg package has codes for image segmentation.  The most used code is the findblobs3 class which creates indexed objects images and performs typical measurements and morphological operations on them.  The jsim package has my diffusion and fcs simulation codes.

The j3D package contains codes for transformations on 3D objects.  The renderer class aggregates the objects for collective transformation.

The jguis package contains lots of graphical interface codes as well as file readers and codes that interface with ImageJ.  This contains all of my plotting tools, segmentation interfaces, fitting interfaces, etc.  The jutils class contains a large number of static utilities for ImageJ operations.

The jhcsfcs package contains little used high content screening methods for FCS data.