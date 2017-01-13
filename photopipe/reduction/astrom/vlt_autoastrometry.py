#!/usr/bin/env python
#  

# NAME:
#	vlt_autoastrometry.py
#
# PURPOSE:
#	A fast astrometric solver that uses a combination of pair-distance matching and asterism 
#   matching to solve a field without knowing the position angle.  The full list of steps is:
#    	1 - Extract all stars from the image.  Attempts to filter out cosmic rays,
#         	galaxies, bad columns, and so on.
#     	2 - Query a catalog to get a list of known star positions.
#     	3 - Trim the catalogs if necessary and try to optimize the search strategy
#     	    based on the relative catalog densities and areas.
#     	4 - Calculate the distance from every star to its neighbors, both in the
#      		image and the catalog.  Stars too far apart are ignored.
#         	If the distance between a star and a neighbor matches the distance between
#         	some catalog star and one of its neighbors, calculate that PA too and store
#         	as a potential match.
#     	5 - Look for stars that seem to have a lot of distance matches with catalog
#         	stars.  If the PA's also match, store the list as an asterism match.
#     	6 - Prune out asterisms whose matches seem least useful.
#     	7 - Calculate analytically the PA shift and offset by averaging asterism centers.
#     	8 - Update the header and write to disk.
#
# INPUT FILES:
# 	Sextractor input files:
#   	You can modify the following sextractor inputs if you choose.  This should
#     	rarely be necessary, if ever.
#       	sex.config - Overall input file
#        	sex.conv   - Convolution matrix for detection
#     These files are written to disk using internal defaults if not present.
#
# OUTPUT FILES:
# 	Explanation of output files:
#    	DS9 region files -
#      		cat.wcs.reg        - WCS positions of all objects in the catalog.
#      		det.im.reg         - XY positions of all non-flagged detected objects.
#      		matchlines.im.reg  - spoke lines for matched asterisms from XY positions.
#      		matchlines.wcs.reg - spoke lines for matched asterisms from WCS positions.
#    	Text output -
#      		cat.txt            - objects in the catalog as a text file (RA, dec, mag)
#      		det.init.txt       - objects in the image as a text file ('RA', 'dec', mag)
#      		det.final.txt      - objects in the image as a text file (RA, dec, mag)
#      		match.list         - list of matched stars: (X Y RA dec)
#    	Image output -
#      		a[image.fits]      - output image with astrometry (if not specified with -o)
#      		temp.fits          - image with provisional WCS from your -px and -pa inputs
#
# CATALOG HELP:
#	Leave the catalog field blank will use SDSS if available and USNO otherwise.
#   The catalog query uses wcstools (tdc-www.harvard.edu/wcstools).  However, you
#   can also use your own catalog file on disk if you prefer using -c [filename]
#   The default format is a text file with the first three columns indicating
#   ra, dec, magnitude.  However, you can change the order of the columns by
#   adding, e.g.
#		#:1,2,6
#   to the first line.  In this case, this would indicate that the RA is in the
#   1st column, dec in the 2nd, and magnitude in the 6th.   The mag column can be
#   omitted completely, although if the catalog is not the same depth as the
#   image this may compromise the search results.
#
# EXAMPLES:
# 	For an image with provisional WCS header but possibly incorrect PA:
#   	autoastrometry image.fits
#   For an image with provisional WCS, definitely correct PA:
#       autoastrometry image.fits -upa 0.5
#   For an image with no WCS (or bad WCS) and a pixel scale of 0.35"/pixel:
#       autoastrometry image.fits -px 0.35
#   For an image with no WCS, pixel scale of 0.35", and PA of 121 degrees:
#       autoastrometry image.fits -px 0.35 -pa 121 -upa 0.5
#   For an image with no WCS, pixel scale of 0.35, and positive parity:
#       autoastrometry image.fits -px 0.35 -inv
#    Use a catalog on disk instead of wcstools query:    
#       autoastrometry image.fits -c catalog.txt
#    Widen the catalog search to 12x12 arcminutes if the pointing is bad:
#       autoastrometry image.fits -b 720
#    Specify seeing (7 pixels) to better exclude cosmic rays and galaxies:
#       autoastrometry image.fits -s 7
#    Combine lots of options for a maximally directed solution:
#       autoastrometry image.fits -px 0.35 -pa 121 -upa 0.5 -inv -b 600 -s 7
#    (Substitute 'autoastrometry' with 'python autoastrometry.py' if not aliased.)
#
# INSTALLATION:
#	Save this file anywhere on disk, and call it from the command line
#		"python vlt_autoastrometry.py"
#   Required python packages:  numpy, pyfits, optparse, shutil, urllib, and optionally ephem
#		Additionally need homemade packages (astrometrydist, astrometrystats, astrometrysources)
#   You must also have sextractor installed: if the path is nonstandard, edit the global 
#	variable below to specify. 
#   For help, type "python vlt_autoastrometry.py -h"
#
# TROUBLESHOOTING:
# 	Supplying the correct pixel scale (within 1%) and correct parity is critical if the 
#   image does not already contain this information in the FITS header. If you have 
#   difficulty solving a field correctly, double-check these values. If still having
#   trouble, try opening temp.fits and an archival image of the field (from DSS, etc.)
#   and loading the .reg files in DS9.  The problem might be in the telescope pointing/header
#   info (in this case, increase the boxsize) or good matching stars may be thrown away or
#   confused by artifacts (in this case, specify a seeing value).  If the PA is known, 
#   restricting it can also help (try -upa 0.5); by default all orientations are searched.
#   If still having issues, e-mail dperley@astro.berkeley.edu for help.
#
#   AUTHOR: Daniel Perley (dperley@astro.caltech.edu)
#   last significant modifications 2012-04-23
#	 
#	modified by Vicki Toy (vtoy@astro.umd.edu) 11/19/2013
#
# 	4/23/12: program can actually be overwhelmed by too many good matches (too high maxrad).
# 	need to fix this.
# 	some possible future improvements:
# 		verification to assess low-confidence solutions
# 		full automatic retry mode (parity detection, etc.)
# 		dealing with unknown pixel scale
# 		run wcstools for distortion parameters
# 		merge catalog check with catalog search to save a query
# 		improve the CR rejection further... maybe think about recognizing elliptical "seeing"?


sexpath = ''  # if "sex" works in any directory, leave blank
defaulttolerance = 0.01  # these defaults should generally not be altered.
defaultpatolerance = 1.4   
defaultminfwhm = 1.5
defaultmaxfwhm = 25

fastmatch = 1
showmatches = 0

import numpy
from optparse import OptionParser
import shutil
import pyfits
import os, sys
import urllib
import astrometrydist
import astrometrystats
import astrometrysources
import warnings

try:
   import ephem
except:
   pass

####################################
"""
Saves object file with RA, DEC, and mag
"""
def writetextfile(filename, objlist):
    out = open(filename,'w')
    for ob in objlist:
      out.write("%11.7f %11.7f %5.2f\n" % (ob.ra, ob.dec, ob.mag))
    out.close()

####################################
"""
Creates region file for DS9 to box objects from list
"""
def writeregionfile(filename, objlist, color="green",sys=''):
    if sys == '': sys = 'wcs'
    out = open(filename,'w')
    i = -1
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    if sys == 'wcs': 
      out.write('fk5\n')
      for ob in objlist:
        i += 1
        out.write("point(%.7f,%.7f) # point=boxcircle text={%i}\n" % (ob.ra, ob.dec, i))
    if sys == 'img': 
      out.write('image\n')
      for ob in objlist:
        i += 1
        out.write("point(%.3f,%.3f) # point=boxcircle text={%i}\n" % (ob.x, ob.y, i))
    out.close()

##########################
"""
NAME:
	autoastrometry
PURPOSE:
	Runs sextractor on image and compares with a catalog to find the rotation and shift.
	Uses this to correct WCS coordinates.
INPUTS:
	filename - fits file that want to correct
	pixelscale - pixel scale
	pa - position angle
	inv - inversion, used to tell the parity of the WCS coordinates
	uncpa - uncertainty in position angle
	userra - user supplied right ascension
	userdec - user supplied declination
	minfwhm - minimum full width half max used to eliminate non-star sources
	maxfwhm - maximum full width half max used to eliminate non-star sources
	maxellip - maximum ellipticity used to eliminate non-star sources
	boxsize - size of desired field to look in catalog
	maxrad - maximum radius to look for distances between stars when performing astrometry
	tolerance - used to find the expected number of false multiples
	catalog - catalog to look in
	nosolve - specifies that don't want to solve the WCS correction, just want catalog to be output
	overwrite - allows program to overwrite files
	saturation - saturation level for sextractor
	quiet - doesn't print as much information
OUTPUTS:
	Saves file with prefix 'a' when astrometry properly calculated
	
	Returns the number of matches, sky offset position angle, the stddev of the postion angle,
	the shift in both RA and DEC in arcseconds as well as the total stddev of the shift in arcseconds
EXAMPLE:
	(nmatch, skyoffpa, stdevpa, raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec) = autoastrometry(filename,
		pixelscale=-1,pa=-999,inv=0,uncpa=-1,userra=-999, userdec=-999, minfwhm=1.5,maxfwhm=20,maxellip=0.5,
		boxsize=-1,maxrad=-1,tolerance=0.010,catalog='',nosolve=0,overwrite=False, outfile='', saturation=-1, 
		quiet=False)	
NOTE:
    ASSUMES RA is in direction of y and dec in direction of x
"""

def autoastrometry(filename,pixelscale=-1,pa=-999,inv=0,uncpa=-1,userra=-999, userdec=-999, minfwhm=1.5,maxfwhm=20,maxellip=0.5,boxsize=-1,maxrad=-1,tolerance=0.010,catalog='',nosolve=0,overwrite=False, outfile='', saturation=-1, quiet=False):
  
    # Get some basic info from the header  
    try: 
       fits = pyfits.open(filename)
       fits.verify('silentfix')
    except:
       print 'Error opening', filename
       if os.path.isfile(filename)==False: print 'File does not exist.'
       return -1
    h = fits[0].header  #ideally check for primary extension, or even iterate
    sfilename = filename
    
    #Position angle set to 0 if pixel scale set
    if pixelscale > 0 and pa == -999: pa = 0

	#If pixel scale set and position angle set (default is -999), 
	#calculate CD*, CR* header keywords and write as temp.fits
    if pixelscale > 0 and pa > -360:
 
       parad = pa * numpy.pi / 180.
       pxscaledeg = pixelscale / 3600.
       
       if inv > 0: 
          parity = -1
       else:
          parity = 1
          
       if  360. > userra >= 0.:
          ra = userra
       else:
          try:
             ra = astrometrystats.rasex2deg(h['RA'])
          except:
             ra = astrometrystats.rasex2deg(h['CRVAL1'])
             
       if 360. > userdec >= 0.:
          dec = userdec
       else:
          try:
             dec = astrometrystats.decsex2deg(h['DEC'])
          except:
             dec = astrometrystats.decsex2deg(h['CRVAL2'])
       
       epoch = float(h.get('EPOCH', 2000))
       equinox = float(h.get('EQUINOX', epoch)) #If RA and DEC are not J2000 then convert

       if abs(equinox-2000) > 0.5:
           print 'Converting equinox from', equinox, 'to J2000'
           try:
              j2000 = ephem.Equatorial(ephem.Equatorial(str(ra/15), str(dec), epoch=str(equinox)),epoch=ephem.J2000)
              [ra, dec] = [astrometrystats.rasex2deg(j2000.ra), astrometrystats.decsex2deg(j2000.dec)]
           except:
              print 'PyEphem is not installed but is required to precess this image.'
              return -1
       h.set("CD1_1",  pxscaledeg * numpy.cos(parad)*parity)
       h.set("CD1_2",  pxscaledeg * numpy.sin(parad))
       h.set("CD2_1", -pxscaledeg * numpy.sin(parad)*parity)
       h.set("CD2_2",  pxscaledeg * numpy.cos(parad))
       h.set("CRPIX1", h['NAXIS1']/2)
       h.set("CRPIX2", h['NAXIS2']/2)
       h.set("CRVAL1", ra)
       h.set("CRVAL2", dec)
       h.set("CTYPE1","RA---TAN")
       h.set("CTYPE2","DEC--TAN")
       h.set("EQUINOX",2000.0)
       
       if os.path.isfile('temp.fits'): os.remove('temp.fits')
       fits[0].header = h
       fits.writeto('temp.fits',output_verify='silentfix') #,clobber=True
       fits.close()
       fits = pyfits.open('temp.fits')
       h = fits[0].header
       sfilename = 'temp.fits'

    #Read the WCS values from file header (even if we put it there in the first place)
    try:
        # no longer drawing RA and DEC from here.
        nxpix = h['NAXIS1']
        nypix = h['NAXIS2']
    except:
        print 'Cannot find necessary WCS header keyword NAXIS*'
        sys.exit(1)
    try: 
        cra  =  float(h['CRVAL1'])
        cdec = float(h['CRVAL2'])
        crpix1 = float(h['CRPIX1'])  
        crpix2 = float(h['CRPIX2'])
        cd11 = float(h['CD1_1'])
        cd22 = float(h['CD2_2'])
        cd12 = float(h['CD1_2']) # deg / pix
        cd21 = float(h['CD2_1'])
    except:
        print 'Cannot find necessary WCS header keyword CRVAL*, CRPIX*, or CD*_*'
        print 'Must specify pixel scale (-px VAL) or provide provisional basic WCS info via CD matrix.'
        sys.exit(1)
        
	#This section deals with manipulating WCS coordinates, for thorough description see:
	#iraf.noao.edu/iraf/ftp/misc/fitswcs_draft.ps

	#Determine parity from sign of determinant of CD matrix   
    if cd11 * cd22 < 0 or cd12 * cd21 > 0:
       parity = -1
    else:
       parity = 1

    #Calculates CDELTA1 (xscale) and CDELTA2 (yscale) which is how much RA or DEC 
    #changes when you move along a column or row
    xscale = numpy.sqrt(cd11**2 + cd21**2)
    yscale = numpy.sqrt(cd12**2 + cd22**2)
    
    #Calculates CROTA2 based from transformations between CD matrix and CDELT values
    initpa = -parity * numpy.arctan2(cd21 * yscale, cd22 * xscale) * 180 / numpy.pi
    
    #Find field width based on largest dimension and calculates the area of the field
    #as well as the pixel scale in arcseconds
    xscale = abs(xscale)
    yscale = abs(yscale)
    fieldwidth = max(xscale * nxpix, yscale * nypix) * 3600. /2.0   
    area_sqdeg = xscale * nxpix * yscale * nypix
    
    area_sqmin = area_sqdeg * 3600. 
    area_sqsec = area_sqmin * 3600. 
    pixscale = numpy.sqrt(xscale*yscale) * 3600.

	#Finds center pixel in each dimension
    centerx = nxpix/2
    centery = nypix/2
    
    #Calculate how the center pixel relates the to the header value's crpix1/2.  
    #Theoretically, centerx and crpix1 should be the same.  But most of the time, they are not. 
    #This is due to a number of reasons, ranging from the algorithm used to calculate the 
    #astrometry to simply cropping the image at some point in the reduction.  
    #So, the crpix header value may be hundreds of pixels off from the actual center of the actual field.
    centerdx = centerx - crpix1
    centerdy = centery - crpix2
    
    #Calculate the RA and DEC at center of field to correct initial guess
    centerra  = cra  - centerdx*xscale*numpy.cos(initpa*numpy.pi/180.) + centerdy*yscale*numpy.sin(initpa*numpy.pi/180.)
    centerdec = cdec + parity*centerdx*xscale*numpy.sin(-initpa*numpy.pi/180.) + centerdy*yscale*numpy.cos(initpa*numpy.pi/180.)
    # this has only been checked for a PA of zero.
    
    if quiet == False:
       print 'cra=%10.6f, centerdx=%10.6f, xscale=%10.6f, centerdy=%10.6f, yscale=%10.6f' % (cra,centerdx,xscale,centerdy,yscale)
       print 'Initial WCS info:'
       print '   pixel scale:     x=%.4f"/pix,   y=%.4f"/pix' % (xscale*3600, yscale*3600)
       print '   position angle: PA=%.2f' % initpa
       if parity ==  1: print '   normal parity'
       if parity == -1: print '   inverse parity'
       print '   center:        RA=%10.6f, dec=%9.6f' % (centerra, centerdec)
       print '   field width: %10.6f' % (fieldwidth)
    
    #Run sextract (runs sextractor) to produce image star catalog    
    goodsexlist = astrometrysources.sextract(sfilename, nxpix, nypix, 3, 12, minfwhm=minfwhm, maxfwhm=maxfwhm, maxellip=maxellip, saturation=saturation, sexpath=sexpath, quiet=quiet)
	
	#If there are less than 4 good objects, ends program and writes images to txt and region files
    ngood = len(goodsexlist)
    if ngood < 4:
       print 'Only', ngood, 'good stars were found in the image.  The image is too small or shallow, the detection'
       print 'threshold is set too high, or stars and cosmic rays are being confused.'      
       #Saves text file that contains RA, DEC, and mag of sextractor list
       writetextfile('det.init.txt', goodsexlist)
       writeregionfile('det.im.reg', goodsexlist, 'red', 'img')
       return -1
       
	#Finds source number density
    density = len(goodsexlist) / area_sqmin
    
    if quiet == False:
    	print 'Source density of %f4 /arcmin^2' % density
    
    #If set to only solve for catalog and not astrometry, save good list
    if nosolve == 1: 
       if catalog == '': catalog = 'det.ref.txt'
       #Saves text file that contains RA, DEC, and mag of sextractor list
       writetextfile(catalog, goodsexlist)
       return

    #If no catalog specified, find catalog with > 15 entries for center RA and DEC encircled by radius of 180 arcseconds
    #If no catalog found after that, end program
    if catalog == '':
    
    	#trycats = ['tmpsc']
        trycats = ['sdss','tmpsc', 'ub2', 'tmc']
        for trycat in trycats:
        
            if trycat == 'sdss':
            	testqueryurl = "http://cas.sdss.org/dr7/en/tools/search/x_radial.asp?ra="+str(centerra)+"&dec="+str(centerdec)+"&radius="+str(3)+"&entries=top&topnum=50000&format=csv"
            	commentlen = 4
            else:
            	testqueryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + trycat +  "&ra=" + str(centerra) + "&dec=" + str(centerdec) + "&system=J2000&rad=" + str(-180)
            
            check = urllib.urlopen(testqueryurl)
            checklines = check.readlines()
            check.close()

            # Find comment line for those pulled by Harvard using comment identifier
            if trycat != 'sdss':
                indices = [i for i, s in enumerate(checklines) if '--' in s]
                commentlen = indices[0] + 1
            
            if len(checklines) > commentlen:
                catalog = trycat
                if quiet == False:
                	print testqueryurl
                	print 'Using catalog', catalog
                break
        if (catalog == ''):
            print 'No catalog is available.  Check your internet connection.'
            return -1
	
    #Load in reference star catalog with boxsize
    if (boxsize == -1):
        boxsize = fieldwidth/2.
    
    decwidth  = xscale * nxpix*3600./2.
    rawidth = yscale * nypix*3600./2.

    catlist = astrometrysources.getcatalog(catalog, centerra, centerdec, boxsize,rawidth,decwidth)
    
    ncat = len(catlist)
    catdensity = ncat / (2*boxsize/60.)**2
    
    if quiet == False:
    	print ncat, 'good catalog objects.'
    	print 'Source density of %f4 /arcmin^2' % catdensity
    
    #Throws up warning if very few catalog objects, stops program if no catalog objects found
    if 0 < ncat < 5:
       print 'Only ', ncat, ' catalog objects in the search zone.  Increase the magnitude threshold or box size.'

    if ncat == 0 :
       print
       print 'No objects found in catalog.'
       print 'The web query failed, all stars were excluded by the FHWM clip, or the image'
       print 'is too small.  Check input parameters or your internet connection.'
       return -1   
  
    #If this image is actually shallower than reference catalog, trim the reference catalog down
    if ncat > 16 and catdensity > 3 * density:
        print 'Image is shallow.  Trimming reference catalog...'
        while catdensity > 3 * density:
            catlist = catlist[0:len(catlist)*4/5]
            ncat = len(catlist)
            catdensity = ncat / (2*boxsize/60.)**2

    #If the image is way deeper than USNO, trim the image catalog down
    #if ngood > 8 and density > 4 * catdensity:
    #    print 'Image is deep.  Trimming image catalog...'
    #    while density > 4 * catdensity and ngood > 8:
    #        goodsexlist = goodsexlist[0:len(goodsexlist)*4/5]
    #        ngood = len(goodsexlist)
    #        density = ngood / area_sqmin

    #If too many objects, do some more trimming
    if ngood*ncat > 120*120*4:
        print 'Image and/or catalog still too deep.  Trimming...'
        while ngood*ncat > 120*120*4:
            if density > catdensity: 
                goodsexlist = goodsexlist[0:len(goodsexlist)*4/5]
                ngood = len(goodsexlist)
                density = ngood / area_sqmin
            else:
                catlist = catlist[0:len(catlist)*4/5]
                ncat = len(catlist)
                catdensity = ncat / (2*boxsize/60.)**2   

    #Remove fainter object in close pairs for both lists
    goodsexlist = astrometrydist.tooclose(goodsexlist, minsep=3, quiet=quiet)
    catlist = astrometrydist.tooclose(catlist, minsep=3, quiet=quiet)

    #Saves text file that contains RA, DEC, and mag of sextractor list
    writetextfile('det.init.txt', goodsexlist)
    writeregionfile('det.im.reg', goodsexlist, 'red', 'img')
    writetextfile('cat.txt', catlist)
    writeregionfile('cat.wcs.reg', catlist, 'green', 'wcs')   

    ##### The catalogs have now been completed. Now start getting into the actual astrometry. #####
    
    #Maximum radius (in arcseconds) calculated by looking at the radius of 15 object in the sparsest dataset
    #Must at least 60 arcseconds or 75% of the field width (whichever is smaller)
    minrad = 5.0
    if (maxrad == -1):  
    	#60 arcsec/arcmin *sqrt(15 stars/(stars/arcsec^2)) / 2 <-- density is for a box, so half the length of the box is the radius                                          
        maxrad = 30.0*(15.0 /min(density,catdensity))**0.5  
        maxrad = max(maxrad, 60.0)
        
        if maxrad == 60.0:                           
             minrad = 10.0   # in theory could scale this up further to reduce #comparisons i.e. instead of 15 stars, x stars
        maxrad = min(maxrad, fieldwidth*3./4)       
        
    #Finds the number of objects expected within circular area (chooses smaller of the entire field or of area 
    #in between min and max radius). NOTE: density is per arcmin^2, while the radii are in arcsec, hence the conversion factor. 
    circlearea     = (numpy.pi*(maxrad/60.)**2 - numpy.pi*(minrad/60)**2) #in arcmin^2
    circdensity    = density * min([area_sqmin, circlearea])
    circcatdensity = catdensity * circlearea	#Finds number of catalog objects expected within circular area
    catperimage    = catdensity * area_sqmin	#Finds number of catalog objects expected within field
    
    if quiet == False:
    	print 'After trimming: '
    	print '   ', len(goodsexlist), 'detected objects (%.2f/arcmin^2, %.1f/searchzone)' % (density, circdensity)
    	print '   ', len(catlist),     'catalog objects (%.2f/arcmin^2, %.1f/searchzone)' % (catdensity, circcatdensity)
	
	#Sets position angle tolerance and calculates the expected number of false multiples
    patolerance = defaultpatolerance
    expectfalsepairs = ngood * ncat * circdensity**1 * circcatdensity**1 * tolerance**1 * (patolerance/360.)**0
    expectfalsetrios = ngood * ncat * circdensity**2 * circcatdensity**2 * tolerance**2 * (patolerance/360.)**1
    expectfalsequads = ngood * ncat * circdensity**3 * circcatdensity**3 * tolerance**3 * (patolerance/360.)**2
    expectfalsequint = ngood * ncat * circdensity**4 * circcatdensity**4 * tolerance**4 * (patolerance/360.)**3

	#Guess that 30% of the sextractor sources overlap with stars, finds estimate of how many real matches we expect
    overlap1 = 0.3 * min(1,catdensity/density) # fraction of stars in image that are also in catalog - a guess
    truematchesperstar = (circdensity * overlap1) # but how many matches >3 and >4?  some annoying binomial thing
    
    #Default required match is 3 (triangle), but can require more if we expect a log of false triples or less if not many sources in either image
    reqmatch = 3
    if expectfalsetrios > 30 and truematchesperstar >= 4: reqmatch = 4   
       #should check that this will actually work for the catalog, too.
    if catperimage <= 6 or ngood <= 6: reqmatch = 2 
    if catperimage <= 3 or ngood <= 3: reqmatch = 1
        #for an extremely small or shallow image
        
    #Calculates the matched stars between sextractor and catalog
    if quiet == False:
    	print 'Pair comparison search radius: %.2f"'%maxrad
    	print 'Using reqmatch =', reqmatch
    	
    (primarymatchs, primarymatchc, mpa) = astrometrydist.distmatch(goodsexlist, catlist, maxrad, minrad, reqmatch, patolerance, uncpa, showmatches=showmatches, fastmatch=fastmatch, quiet=quiet)
  
    #Quits program if no matches or too few matches found (gives different error readouts)
    nmatch = len(primarymatchs)
    if nmatch == 0:
        print ' No valid matches found!'
        if quiet == False:
           print ' Possible issues:'
           print '  - The specified pixel scale (or PA or parity) is incorrect.  Double-check the input value.'
           print '  - The field is outside the catalog search region.  Check header RA/DEC or increase search radius.'
           print '  - The routine is flooded by bad sources.  Specify or check the input seeing.'
           print '  - The routine is flagging many real stars.  Check the input seeing.'
           print ' You can display a list of detected/catalog sources using det.im.reg and cat.wcs.reg.'
        return -1
    if nmatch <= 2:
        print 'Warning: only', nmatch, 'match(es).  Astrometry may be unreliable.'
        if quiet == False:
           print '   Check the pixel scale and parity and consider re-running.'
        #return -1
        warning = 1

    #We now have the PA and a list of stars that are almost certain matches.
    offpa = astrometrystats.median(mpa)  #get average PA from the excellent values
    stdevpa = astrometrystats.stdev(mpa)
    
    skyoffpa = -parity*offpa # This appears to be necessary for the printed value to agree with our normal definition.
    
    if quiet == False:
    	print 'PA offset:'
    	print '  dPA = %.3f  (unc. %.3f)' % (skyoffpa, stdevpa)

    # Rotate the image to the new, correct PA
    #  NOTE: when CRPIX don't match CRVAL this shifts the center and screws things up.  
    #  I don't understand why they don't always match.  [[I think this was an equinox issue.
    #  should be solved now, but be alert for further problems.]]
    
    #Greisen et al.:
    #WCS_i = SUM[j] (CD_ij)(p_j - CRPIX_j)      i.e.
    # RA - CRVAL1 = CD1_1 (x - CRPIX1) + CD1_2 (y - CRPIX2)
    #dec - CRVAL2 = CD2_1 (x - CRPIX1) + CD2_2 (y - CRPIX2)   [times a projection scale...]

    #Rotate CD matrix to account for additional rotation calculated from catalog star matching
    rot = offpa * numpy.pi/180
      #...the image itself
    h.set("CD1_1", numpy.cos(rot)*cd11 - numpy.sin(rot)*cd21 )
    h.set("CD1_2", numpy.cos(rot)*cd12 - numpy.sin(rot)*cd22 )  # a parity issue may be involved here?
    h.set("CD2_1", numpy.sin(rot)*cd11 + numpy.cos(rot)*cd21 )
    h.set("CD2_2", numpy.sin(rot)*cd12 + numpy.cos(rot)*cd22 )
      #...the coordinates (so we don't have to resex) 
    for i in range(len(goodsexlist)):  #do all of them, though this is not necessary
        goodsexlist[i].rotate(offpa,cra,cdec)

    #Saves text file that contains RA, DEC, and mag of sextractor list
    writetextfile('det.wcs.txt', goodsexlist)
    
    #Calculate shift in RA and DEC for each object compared to its catalog match
    imraoffset = []
    imdecoffset = []
    for i in range(len(primarymatchs)):
        imraoffset.append(goodsexlist[primarymatchs[i]].ra - catlist[primarymatchc[i]].ra)
        imdecoffset.append(goodsexlist[primarymatchs[i]].dec - catlist[primarymatchc[i]].dec)
    
    #Find the median and standard deviation of offsets    
    raoffset = -astrometrystats.median(imraoffset)
    decoffset = -astrometrystats.median(imdecoffset)
    rastd = astrometrystats.stdev(imraoffset)*numpy.cos(cdec*numpy.pi/180)  # all of these are in degrees
    decstd = astrometrystats.stdev(imdecoffset)
    stdoffset = numpy.sqrt(rastd**2 + decstd**2)
    
    #Change from degrees to arcseconds
    raoffsetarcsec = raoffset*3600*numpy.cos(cdec*numpy.pi/180)
    decoffsetarcsec = decoffset*3600
    totoffsetarcsec = (raoffsetarcsec**2 + decoffset**2)**0.5
    stdoffsetarcsec = stdoffset*3600
    
    if quiet == False:
    	print 'Spatial offset:'
    	print '  dra = %.2f",  ddec = %.2f"  (unc. %.3f")' % (raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec)
    
    #If standard deviation of total offset is larger than 10 arcseconds end program
    warning = 0
    if (stdoffset*3600 > 10.0):
        print 'WARNING: poor solution - some matches may be bad.  Check pixel scale?'
        return -1
        warning = 1
    
    #Shift center pixel to match catalog values
    h.set("CRVAL1", cra + raoffset)
    h.set("CRVAL2", cdec + decoffset)
    
    #Add keywords to header to detail changes made to WCS coordinates
    try:
       oldcat = h['ASTR_CAT']
       h.set("OLD_CAT",oldcat, "Earlier reference catalog")
    except:
       pass
    h.set("ASTR_CAT", catalog, "Reference catalog for vlt_autoastrometry")
    h.set("ASTR_UNC", stdoffsetarcsec, "Astrometric scatter vs. catalog (arcsec)")
    h.set("ASTR_SPA", stdevpa, "Measured uncertainty in PA (degrees)")
    h.set("ASTR_DPA", skyoffpa, "Change in PA (degrees)")
    h.set("ASTR_OFF", totoffsetarcsec, "Change in center position (arcsec)")
    h.set("ASTR_NUM", len(primarymatchs), "Number of matches")

    #Write out a match list to allow doing a formal fit with WCStools.
    outmatch = open('match.list','w')
    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        ci = primarymatchc[i]
        outmatch.write("%s %s  %s %s\n" % (goodsexlist[si].x, goodsexlist[si].y, catlist[ci].ra, catlist[ci].dec))
    outmatch.close()                                                     
    
    # Could repeat with scale adjustment
    # Could then go back to full good catalog and match all sources
    
    #Create new file with new header
    if overwrite: outfile = filename
    if outfile == '': 
        slashpos = filename.rfind('/')
        dir = filename[0:slashpos+1]
        fil = filename[slashpos+1:]
        outfile = dir+'a'+fil # alternate behavior would always output to current directory
    try:
        os.remove(outfile)
    except:
        pass
    fits[0].header = h
    fits.writeto(outfile,output_verify='silentfix') #,clobber=True
    
    if quiet == False:
    	print 'Written to '+outfile

    fits.close()
    
    #Return relevant offsets
    return (nmatch, skyoffpa, stdevpa, raoffsetarcsec, decoffsetarcsec, stdoffsetarcsec)

######################################################################
def usage():
    (xdir,xname) = os.path.split(sys.argv[0])
    print "Usage:  %s filename(s) [-x pixelscale -p PA -inv -b boxsize -s seeing -upa PAunc -l saturation]" % xname
    print "     or %s -h for instructions and more options." % xname

####################################
"""
Parses out keywords from argv input and runs autoastrometry.  If multiple input files and
set to solve, then prints out return values of autoastrometry
"""

def main():
    
    files=[]

    if (len(sys.argv)==1):
        usage()
        sys.exit(1)

    #defaults
    nosolve = 0
    overwrite = 0
    outfile = ''
    warnings.filterwarnings('ignore')
    
    #Sets up option flag parser with help files.  If unsure how to run run program with -h flag
    
    parser = OptionParser()
    
    parser.add_option("-l",	"--sat", help="saturation level (do not use stars exceeding)", 
    					default=-1, action="store", type="float", dest="sat")

    parser.add_option("-x", "--pix", help="pixel scale in arcsec/pixel (within 1%)", 
    					default=-1, action="store", type="float", dest="pix")    

    parser.add_option("-b", "--box", help="half-width of box for reference catalog query (arcsec)", 
    					default=-1, action="store", type="float", dest="box") 

    parser.add_option("-s", "--see", help="approximate seeing in pixels for CR/star/galaxy ID'ing", 
    					default=-1, action="store", type="float", dest="see") 
    					
    parser.add_option("-u", "--upa", help="Uncertainty of the position angle (degrees)", 
    					default=-1, action="store", type="float", dest="upa") 
    					    					
    parser.add_option("-c", "--catalog", help="Catalog to use (tmpsc, ub2, tmc, sdss, or file: see --catalog)", 
    					default='', action="store", dest="cat")  					    					

    parser.add_option("-t", "--toler", help="Amount of slack allowed in match determination", 
    					default=defaulttolerance, action="store", type="float", dest="tol") 
    					    					
    parser.add_option("-e", help="Maximum ellipticity", 
    					default=-1, action="store", type="float", dest="mel") 
   					    					    					    					
    parser.add_option("-p", "--pa", help="The position angle in degrees.  Not usually needed.", 
    					default=-999, action="store", type="float", dest="pa") 

    parser.add_option("-i", "--inv", help="Reverse(=positive) parity.", 
    					default=False, action="store_true", dest="inv") 
    					
    parser.add_option("-q", "--quiet", help="Quiet", 
    					default=False, action="store_true", dest="quiet")     					

    parser.add_option("-m", help="Maximum distance to look for star pairs (in radius).", 
    					default=-1, action="store", type="float", dest="maxrad") 
    					
    parser.add_option("-r", "--ra", help="Right ascension", 
    					default='', action="store", dest="userra")     					

    parser.add_option("-d", "--dec", help="Declination", 
    					default='', action="store", dest="userdec")   
    					
    parser.add_option("-o", "--output", help="Specify output file (no argument overwrites input file)", 
    					default='', action="store", dest="outfile")       					

    parser.add_option("-n", "--nosolve", help="Do not attempt to solve astrometry; just write catalog.", 
    					default='', action="store", dest="nosolvecat")
    					
    					
    (options, args) = parser.parse_args()
    
    #If nosolve set, then set nosolve=1
    if options.nosolvecat != '': 
    	options.cat = options.nosolvecat
    	nosolve = 1
    if options.inv:
    	inv = 1
    else:
    	inv = 0
    
    #Turns userra and userdec into degrees or sets number default (parser stored as strings because it is in sexagesimal)
    if options.userra == '':
    	userra = -999
    else:
    	userra = astrometrystats.rasex2deg(options.userra)
    	
    if options.userdec == '':
    	userdec = -999
    else:
    	userdec = astrometrystats.decsex2deg(options.userdec)
    
    #Finds input files and stores in variable files.  Assumes input files are argv that come before flagged keywords
    #If no files present prints warning and quits
    isys=0
    for value in sys.argv[1:]:
    	if value[0] == '-':
    		break
    	isys += 1

    files = sys.argv[1:isys+1]
    
    if len(files) == 0:
    	print 'No files selected!'
    	return
    
    if (options.see == -1):
    	minfwhm = defaultminfwhm
    	maxfwhm = defaultmaxfwhm
    else:
    	minfwhm = 0.7 * options.see
    	maxfwhm = 2.0 * options.see
    
    #Find directory where this python code is located
    call = sys.argv[0]
    place = call.find('vlt_autoastrometry.py')
    propath = call[:place]
    

    #Copies parameter file and configuration file from default files located in program folder   
    shutil.copyfile(propath+'tempsource.param', 'tempsource.param')
    
    if not os.path.exists('sex.config'): 
    	shutil.copyfile(propath+'sex.config', 'sex.config')
    
    if not os.path.exists('sex.conv'):    
    	shutil.copyfile(propath+'sex.conv', 'sex.conv')    
    	
    if not os.path.exists('default.nnw'):    
    	shutil.copyfile(propath+'default.nnw', 'default.nnw') 

    nimage = len(files)
    failures = []
    questionable = []
    multiinfo = []    

    for file in files:
    	if options.quiet == False:
        	print 'Processing', file
        	print 'userra and dec'
        	print userra, userdec
        
        fitinfo = autoastrometry(file,pixelscale=options.pix,pa=options.pa,inv=inv,
        						uncpa=options.upa,minfwhm=minfwhm,maxfwhm=maxfwhm,
        						maxellip=options.mel,boxsize=options.box, maxrad=options.maxrad, 
        						userra=userra, userdec=userdec, tolerance=options.tol, 
        						catalog=options.cat, nosolve=nosolve, overwrite=overwrite, 
        						outfile=options.outfile, saturation=options.sat, quiet=options.quiet)
		
        if nosolve: continue
        
        #Save output of autoastrometry to multiinfo
        if type(fitinfo)==int: 
           fitinfo = (0,0,0,0,0,0)
        
        multiinfo.append(fitinfo)

		#If number of matches (fitsinfo[0]) from autoastrometry is 0, list in failures
		#If std is greater than 2 then list in questionable
        if (fitinfo[0] == 0):   #number of matches
            failures.append(file)
        if (fitinfo[5] > 2):    #stdev of offset
            questionable.append(file)
    
    #Returns failed and questionable images.  Also prints all values out from multiinfo       
    if nimage > 1 and nosolve==0:

        if len(failures) == 0 and len(questionable) == 0 and options.quiet == False:
            print 'Successfully processed all images!'
        else:
            print 'Finished processing all images.'
        
        if len(questionable) > 0:
            print 'The following images solved but have questionable astrometry: '
            print '    ',
            for f in questionable: print f
            
        if len(failures) > 0:
            print 'The following images failed to solve: '
            print '    ',
            for f in failures: print f

        print "%25s " %'Filename', 
        print "%6s %8s (%6s)  %7s %7s (%6s)" % ('#match', 'dPA ', 'stdev', 'dRA', 'dDec', 'stdev')
        for i in range(len(files)):
            info = multiinfo[i]
            print "%25s " % files[i], 
            if info[0] > 0:
               print "%6d %8.3f (%6.3f)  %7.3f %7.3f (%6.3f)" % info
            else:
               print "failed to solve"
               
	#Removes temp.param file               
    try:
       os.remove('tempsource.param')
    except:
       print 'Could not remove tempsource.param for some reason'
    
"""    
        if (arg.find("-examp") == 0):
            examplehelp()
            sys.exit(1)
        if (arg.find("-troub") == 0):
            troublehelp()
            sys.exit(1)
        if (arg.find("-catal") == 0):
            cataloghelp()
            sys.exit(1)
        if (arg.find("-output") == 0):
            outputhelp()
            sys.exit(1)
        if (arg.find("-input") == 0):
            inputhelp()
            sys.exit(1)
        if (arg.find("-algor") == 0):
            inputhelp()
            sys.exit(1) 
"""

######################################################################
# Running as executable
if __name__=='__main__':
    main()

######################################################################

