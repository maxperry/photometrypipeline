#Photometry processing library (misc function that can be reused)

import astropy.io.fits as pf
import fnmatch
import os
from numpy import sqrt
import numpy as np



    # make a circle for identifying sources in images
    
"""
NAME:
	cicle
PURPOSE:
	Create x and y array that represent circle coordinates given the center and radius of circle
INPUTS:
	xcenter, ycenter - Center coordinates of circle
	radius			 - Radius in same units as center coordinates
OUTPUTS:
	Returns circle coordinates
EXAMPLE:
	cir = circle(0,0,12)
"""
def circle( xcenter, ycenter, radius ):
    points = np.linspace( 0., 2.*np.pi, 100 )
    x = xcenter + radius * np.cos(points)
    y = ycenter + radius * np.sin(points)
    return np.transpose([x,y])

"""
NAME:
	nearest
PURPOSE:
	Find coordinates (Cartesian) that are within the specified distance of a reference coordinate
INPUTS:
	x, y      - reference coordinate
	xarr,yarr - array of coordinates to compare to reference coordinate
	maxdist   - maximum distance that coordinates can be from reference coordinates
OUTPUTS:
	Returns mask of values that are within maximum distance from reference
EXAMPLE:
	match = nearest(ra[0]*cos(dec[0]*pi/180.),dec[0],refra*cos(refdec*pi/180.),refdec, 1./3600.)
"""
def nearest(x, y, xarr, yarr, maxdist):
	
	dist = sqrt( (x-xarr)**2 + (y-yarr)**2)
	good = (dist < maxdist)
	
	return good

"""
NAME:
	choosefiles
PURPOSE:
	Returns files in directory "loc" that contain "selection" text
INPUTS:
	selection - text required to be in filename
OPTIONAL KEYWORD INPUTS:
	loc - directory where to look, default is current working directory
OUTPUT:
	Matches to search
EXAMPLE:
	matches = choosefiles('coadd*.fits', loc='.')

Written by John Capone (jicapone@astro.umd.edu)
"""
#
def choosefiles( selection, loc='.' ):
	matches = []
	for files in os.listdir(loc):
		if fnmatch.fnmatch( files, selection ):
			matches.append(files)
	return matches

"""
NAME:
	weightedge
PURPOSE:
	To find the innermost corner of a weighted image for a good crop
	Looks for point where values are roughly constant between rows or columns
	Can change this with scale keyword	
INPUTS:
	array - array to search
	itarray - array of indices to search through	
OPTIONAL KEYWORD INPUTS:
	column - set True if iterating through x-axis
	row    - set True if iterating through y-axis
	scale  - set value to look for scaling to previous row or column
OUTPUT:
	Returns the column or row of "nonzero" component on weighted file
EXAMPLE:
	leftedge = weightedge(data, range(xaxis), column=1, scale=0.99)
	
Written by Vicki Toy (vtoy@astro.umd.edu)
"""
def weightedge(array, itarray, scale=1, column=None, row=None):
	oldsum = 0
	for i in itarray:
		if column is not None:
			newsum = sum(array[:,i])
		elif row is not None:
			newsum = sum(array[i,:])
			
		if scale*newsum >= oldsum:
			oldsum = newsum
		else:
			return i-1		

'''
Pared down version of hextract.pro from the IDL Astro Library
Translated to python by Vicki Toy (vtoy@astro.umd.edu)

NAME:
	hextractlite
PURPOSE:
	Extracts subimage from an array and updates astrometry in FITS file.  Saves to newfile name
INPUTS:
	newfile     - file name to save altered fits file to
	data        - array to be altered
	fitsheader  - header to be altered
	x1,x2,y1,y2 - respectively: first and last x pixel, first and last y pixel to be
				  extracted from data array.  Need to be integer scalars
OUTPUTS:
	Saves altered data to newfile
EXAMPLE:
	hextractlite(newfile, data, fitsheader, 100,500,20,1900)
'''

def hextractlite(newfile, data, fitsheader, x1, x2, y1, y2):

	#Truncates the coordinates
	x1 = int(x1)
	y1 = int(y1)
	x2 = int(x2)
	y2 = int(y2)

	fitsheader.add_history('HEXTRACT: Original image size was '+str(fitsheader['NAXIS2'])+' by '+str(fitsheader['NAXIS1']))
	fitsheader.add_history('Extracted Image: ['+ str(y1)+':'+str(y2+1)+','+str(x1)+':'+ str(x2+1)+']')

	#Naxis altered to reflect the new subarray
	fitsheader.update('naxis1', x2-x1+1)
	fitsheader.update('naxis2', y2-y1+1)
	
	#Crpix shifted by new origin
	oldcrpix1 = fitsheader['crpix1']
	oldcrpix2 = fitsheader['crpix2']
	fitsheader.update('crpix1', oldcrpix1-x1)
	fitsheader.update('crpix2', oldcrpix2-y1)

	#Extract subarray
	#Python has opposite coordinate convention
	newdata = data[y1:y2+1, x1:x2+1]	
	
	#Saves new fits file and header to specified fits file name
	pf.writeto(newfile, newdata, fitsheader, clobber=True)

"""
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Lance Simms, Stanford University 2009
#NAME:
#   Djs_Iterstat.py
# 
# PURPOSE:
#   Compute the mean, median and/or sigma of data with iterative sigma clipping.
#
# CALLING SEQUENCE:
#   djs_iterstat, image, [sigrej=, maxiter=, mean=, median=, sigma=, mask=]
#
# INPUTS:
#   image:      Input data
#
# OPTIONAL INPUTS:
#   SigRej:     Sigma for rejection# default to 3.0
#   MaxIter:    Maximum number of sigma rejection iterations# default to 10
#   Min:        The minimum data value to keep; default to min of array
#   Max:        The maximum data value to keep; default to max of array
#   RejVal:     A value that was used as a flag for bad data; ie. -1
#   BinData:   1 to bin the data to find a mode; 0 to not 
# OUTPUTS:
#
# OPTIONAL OUTPUTS:
#   Mean:       Computed mean
#   Median:     Computed median
#   Sigma:      Computed sigma
#   Mode:       Computed Mode
#   Mask:       Mask set to 1 for good points, and 0 for rejected points
#
# PROCEDURES CALLED:
#
# COMMENTS:
#   This routine is based upon Mark Dickinson's IRAF (!) script ITERSTAT.
#   It iteratively rejects outliers as determined by SIGREJ.  It stops
#   when one of the following conditions is met:
#   (1) The maximum number of iterations, as set by MAXITER, is reached.
#   (2) No new pixels are rejected, as compared to the previous iteration.
#   (3) At least 2 pixels remain from which to compute statistics.  If not,
#       then the returned values are based upon the previous iteration.
# 
# REVISION HISTORY:
#   IDL 16-Jun-1999  Written by David Schlegel, Princeton
#   IDL 11-Sep-2000  Speeded up by Hogg and Eisenstein
#   IDL 18-Sep-2000  Note change in MASK values to =1 for good (unrejected) points.
#   PYTHON 12-Jul-2008 Adapted to Scipy by Lance Simms, Stanford
#
#EXAMPLE:
#from Djs_IterStat import *
#FMean, FSig, FMedian, FMask = Djs_Iterstat(10+cos(frange(10000)/1000),RejVal=-1, SigRej=1.3)
#mplot.plot(FMask)
#
############################################################
"""
from numpy import *

def djs_iterstat(InputArr, SigRej=3.0, MaxIter=10, Mask=0,\
                 Max='', Min='', RejVal='', BinData=0):
 
  NGood    = InputArr.size  
  ArrShape = InputArr.shape
  if NGood == 0: 
    print 'No data points given'
    return 0, 0, 0, 0, 0
  if NGood == 1:
    print 'Only one data point; cannot compute stats'
    return 0, 0, 0, 0, 0

  #Determine Max and Min
  if Max == '':
    Max = InputArr.max()
  if Min == '':
    Min = InputArr.min()
 
  if unique(InputArr).size == 1:
    return 0, 0, 0, 0, 0
 
  Mask  = zeros(ArrShape, dtype=byte)+1
  #Reject those above Max and those below Min
  Mask[InputArr > Max] = 0
  Mask[InputArr < Min] = 0
  if RejVal != '' :  Mask[InputArr == RejVal]=0
  FMean = sum(1.*InputArr*Mask) / NGood
  FSig  = sqrt(sum((1.*InputArr-FMean)**2*Mask) / (NGood-1))

  NLast = -1
  Iter  =  0
  NGood = sum(Mask)
  if NGood < 2:
    return -1, -1, -1, -1, -1

  while (Iter < MaxIter) and (NLast != NGood) and (NGood >= 2) :

    LoVal = FMean - SigRej*FSig
    HiVal = FMean + SigRej*FSig
    NLast = NGood

    Mask[InputArr < LoVal] = 0
    Mask[InputArr > HiVal] = 0
    NGood = sum(Mask)

    if NGood >= 2:
      FMean = sum(1.*InputArr*Mask) / NGood
      FSig  = sqrt(sum((1.*InputArr-FMean)**2*Mask) / (NGood-1))
      SaveMask = Mask.copy()
    else:
      SaveMask = Mask.copy()

    Iter = Iter+1
  if sum(SaveMask) > 2:
    FMedian = median(InputArr[SaveMask == 1])
    if BinData == 1:
      HRange  = InputArr[SaveMask==1].max()-InputArr[SaveMask==1].min()
      bins_In = arange(HRange)+InputArr[SaveMask==1].min()
      Bins, N = histOutline.histOutline(InputArr[SaveMask == 1], binsIn = bins_In)
      FMode   = Bins[(where(N == N.max()))[0]].mean()
    else: 
      FMode = 0
  else:
    FMedian = FMean
    FMode   = FMean
    

  return FMean, FSig, FMedian, FMode, SaveMask 

