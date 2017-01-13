import sys
import numpy
import os
import astrometrystats
import urllib
import pyfits as pf
from astropy import wcs

########################################
class Obj:
    ra = 0.0
    dec = 0.0
    mag = 0.0
    
    ra_rad = 0.0
    dec_rad = 0.0
    
    def __init__(self, inra, indec, inmag):
        self.ra = inra
        self.dec = indec
        self.ra_rad = inra * numpy.pi/180
        self.dec_rad = indec * numpy.pi/180
        self.mag = inmag
    
    def rotate(self, dpa_deg, ra0, dec0):
        dpa_rad = dpa_deg * numpy.pi/180
        sindpa = numpy.sin(dpa_rad)
        cosdpa = numpy.cos(dpa_rad)
        rascale = numpy.cos(dec0*numpy.pi/180)
        
        #this is only valid for small fields away from the pole.
        x = (self.ra  - ra0 ) * rascale
        y = (self.dec - dec0)
        
        xrot = cosdpa * x - sindpa * y
        yrot = sindpa * x + cosdpa * y
        
        self.ra   = (xrot / rascale) + ra0
        self.dec  =  yrot + dec0 
        self.ra_rad  = self.ra  * numpy.pi/180
        self.dec_rad =  self.dec * numpy.pi/180

#########################################
class SexObj(Obj):
    x = 0.
    y = 0.
    mag = 0.0
    magerr = 0.0
    ellip = 0.0
    fwhm = 0.0
    flag = 0
    
    def __init__(self, row):
    
        self.x = row[0]
        self.y = row[1]
        self.ra = row[2]
        self.dec = row[3]
        self.mag = row[4]
        self.magerr = row[5]
        self.ellip = row[6]
        self.fwhm = row[7]
        self.flag = row[8]

        self.ra_rad = self.ra * numpy.pi/180
        self.dec_rad = self.dec * numpy.pi/180

####################################
"""
NAME:
	sextract
PURPOSE:
	Runs sextractor with sex.config configuration file and saturation level.  Removes any borders, corners,
	values above max ellipticity (galaxies) and beyond fwhm constraints (to try to remove galaxies and cosmic rays)
	Returns list with sextractor outputs.
INPUTS:
	sexfilename - name of fits file to run sextract on
	nxpix - number of pixels in the x direction
	nypix - number of pixels in the y direction
	border - sets number of pixels in border to remove (default=3)
	corner - sets number of pixels in corner to remove (default=12)
	minfwhm - minimum acceptable full width half max used to remove cosmic rays (default=1.5)
	maxfwhm - maximum acceptable full width half max used to remove galaxies (default=25)
	maxellip - maximum ellipticity, used to remove galaxies (default=0.5)
	saturation - saturation level that sextractor will use (default=-1)
	sexpath - location of where sextractor can be run from (default='')
OUTPUTS:
	Returns good sextractor list of objects in Obj class ordered by magnitude
EXAMPLE:
	goodsexlist = sextract(slist, 1024, 1024, border=3, corner=12, minfwhm=1.5, maxfwhm=25, maxellip=0.5, saturation=55000, sexpath='')
"""

def sextract(sexfilename, nxpix, nypix, border=3, corner=12, minfwhm=1.5, maxfwhm=25, maxellip=0.5, saturation=-1, sexpath='', quiet=False):

    if maxellip == -1: maxellip = 0.5
    
    if saturation > 0: 
       sexsaturation = saturation
    else:
       sexsaturation = 9e5

	#Runs sextractor and reads in sextractor output file (temp.cat) 
    try:
       # Sextract the image !
       if quiet == False:
       	print sexpath+"sex " + sexfilename + " -c sex.config -SATUR_LEVEL "+str(sexsaturation)
       os.system(sexpath+"sex " + sexfilename + " -c sex.config -SATUR_LEVEL "+str(sexsaturation))
    except:
       print ' Error: Problem running sextractor'
       print ' Check that program is installed and runs at command line using ' + sexpath+'sex'
       sys.exit(1)

	# Read in the sextractor catalog
    try:
    	cat = numpy.loadtxt('temp.cat', dtype='float', comments='#')
    except:
    	print 'Cannot load sextractor output file!'
    	sys.exit(1)   	
    
    if len(cat) == 0:
    	print 'Sextractor catalog is empty: try a different catalog?'	
    	sys.exit(1)

    minx = border
    miny = border
    maxx = nxpix - border    # This should be generalized
    maxy = nypix - border
    
    #Separate columns from sextractor output
    x 		= cat[:,0]
    y 		= cat[:,1]
    ra 		= cat[:,2]
    dec 	= cat[:,3]
    mag 	= cat[:,4]
    magerr 	= cat[:,5]
    ellip 	= cat[:,6]
    fwhm 	= cat[:,7]
    flag 	= cat[:,8]
    a_imag  = cat[:,9]
    b_imag  = cat[:,10]

    #Initial filtering creates mask that will remove borders, corners,
    #values above max ellipticity (galaxies), and within fwhm constraints
    mask = (ellip <= maxellip) & (fwhm >= minfwhm) & (fwhm <= maxfwhm) \
    		& (x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy) \
    		& (x+y >= corner) & (x+nypix-y >= corner) \
		& (nxpix-x >= corner) & (nxpix-x+nypix-y >= corner) \
		& (a_imag > 1) & (b_imag > 1)
 	
    #Removes flagged values if saturation level set
    if saturation > 0:
    	mask = mask & (flag == 0)   
    
    #VLT Removed code from autoastrometry.py line 387-432 that rules out false detections
    #If too many false positives need to rewrite section in more coherent way
    #OC: Might also be good to screen for false detections created by bad columns/rows
    
    fwhmlist = fwhm[mask]
	
    #Calculates the 20% value of the sextractor masked FWHM and the mode 
    #to distinguish stars from galaxies (removed by ellipticity) and cosmic rays
    if len(fwhmlist) > 5:
    	sfwhmlist = sorted(fwhmlist)
    	fwhm20 	  = sfwhmlist[len(fwhmlist)/5]	#percentile values
    	fwhmmode  = astrometrystats.most(sfwhmlist, vmax=0, vmin=0)
    else:
    	fwhmmode  = minfwhm
    	fwhm20    = minfwhm
    
	#Creates new FWHM cutoff to distinguish stars from cosmic rays
	#This method done from trial and error from original program to calculate typical cutoffs between stars and cosmic
	#OC: if CR's are bigger and more common than stars, this is dangerous...
    refinedminfwhm = numpy.median([0.75*fwhmmode,0.9*fwhm20,minfwhm]) 
    if quiet == False:
    	print 'Refined min FWHM:', refinedminfwhm, 'pix'
    
    #Masks out fwhm with more stringent fwhm condition and creates newmask
    #that combines the new masked out fwhm values
    refwhmmask = fwhm > refinedminfwhm    
    newmask = refwhmmask & mask
    
    #Create array with bad values removed and then sort by magnitude (column 4)
    goodsext = cat[newmask]
    sortedgoodsext = goodsext[goodsext[:,4].argsort()]
    
    if quiet == False:
    	print len(fwhmlist), 'objects detected in image ('+ str(len(fwhmlist)-len(sortedgoodsext)) +' discarded)'
    
    #Save RA, DEC, and mag into class and save class objects into list
    goodsexlist = []
    for value in sortedgoodsext:
		#RA = column 2, DEC = column 3, mag = column 4 with columns starting at 0 
    	indObj = SexObj(value)
    	goodsexlist.append(indObj)
    
    return goodsexlist
    
####################################
"""
NAME:
	getcatalog
PURPOSE:
	Gets catalog given the center RA, DEC, and box size.  If given a standard catalog, grabs from harvard page, otherwise
	will attempt to open given file.  Saves RA, DEC, mag (in that column order) to an output array for values within magnitute
	and proper motion limits.
INPUTS:
	catalog - catalog name, usually a standard catalog
	ra - right ascension of center of field
	dec - declination of center of field
	boxsize - size of region to look for sources, maximum is fieldwidth
	rawidth - size of half of ra region covered by image
	decwidth - size of half of dec region covered by image
	minmag - minimum magnitude will accept
	maxmag - maximum magnitude will accept
	maxpm - maximum proper motion
OUTPUTS:
	Returns catalog list of objects in Obj class ordered by magnitude
EXAMPLE:
	catlist = getcatalog(catalog, ra, dec, boxsize, rawidth, decwidth, minmag=8.0, maxmag=-1, maxpm=60.)
"""
def getcatalog(catalog, ra, dec, boxsize, rawidth,decwidth,minmag=8.0, maxmag=-1, maxpm=60.):
    # Get catalog from USNO

	#If maximum magnitude is not set, assign max magnitude based on catalog max or set default
    if maxmag == -1:
        maxmag = 999 #default (custom catalog)
        if catalog == 'tmpsc': 	maxmag = 20.0
        if catalog == 'ub2': 	maxmag = 21.0
        if catalog == 'sdss': 	maxmag = 22.0
        if catalog == 'tmc': 	maxmag = 20.0
  
  	#If catalog is one of the standard 4, then set values for indices for catalogs
  	#Search catalog with given box size (if no box size selected, will be fieldwidth)
  	#Read in file to variable catlines
    if (catalog == 'tmpsc' or catalog =='ub2' or catalog=='sdss' or catalog=='tmc'):
        usercat 	= 0
        racolumn 	= 1
        deccolumn 	= 2
        magcolumn 	= 6
        if (catalog == 'tmpsc' or catalog == 'tmc'): magcolumn=3
        pmracolumn 	= 10
        pmdeccolumn = 11   
        
        #SDSS direct query has different format than other scat queries
        if catalog == 'sdss':

        	#queryurl = "http://cas.sdss.org/dr7/en/tools/search/x_radial.asp?ra="+str(ra)+"&dec="+str(dec)+"&radius="+str(boxsize/60.0)+"&entries=top&topnum=6400&format=csv"
        	#queryurl = "http://cas.sdss.org/dr7/en/tools/search/x_rect.asp?min_ra="+str(ra-boxsize/3600.)+"&max_ra="+str(ra+boxsize/3600.)+"&min_dec="+str(dec-2.*boxsize/3600.)+"&max_dec="+str(dec+2.*boxsize/3600.)+"&entries=top&topnum=6400&format=csv"
        	queryurl = "http://cas.sdss.org/dr7/en/tools/search/x_rect.asp?min_ra="+str(ra-rawidth/3600.)+"&max_ra="+str(ra+rawidth/3600.)+"&min_dec="+str(dec-decwidth/3600.)+"&max_dec="+str(dec+decwidth/3600.)+"&entries=top&topnum=6400&format=csv"
   	
        	racolumn    = 7
        	deccolumn   = 8
        	magcolumn   = 12
        	pmracolumn  = 99
        	pmdeccolumn = 99
        	print queryurl
        else:
        	queryurl = "http://tdc-www.harvard.edu/cgi-bin/scat?catalog=" + catalog +  "&ra=" + str(ra) + "&dec=" + str(dec) + "&system=J2000&dra=" + str(rawidth)+ '&ddec=' +str(decwidth) + "&sort=mag&epoch=2000.00000&nstar=6400"
        cat = urllib.urlopen(queryurl)
        catlines = cat.readlines()
        cat.close()
        
        if len(catlines) > 6400-20:
           print 'WARNING: Reached maximum catalog query size.'
           print '         Gaps may be present in the catalog, leading to a poor solution or no solution.'
           print '         Decrease the search radius.'
    else:
        usercat = 1
        try:
           cat = open(catalog,'r')
           print 'Reading user catalog ', catalog
        except:
           print 'Failed to open user catalog ', catalog
           print 'File not found or invalid online catalog.  Specify tmpsc, ub2, sdss, or tmc.'
           return []
        racolumn = 0
        deccolumn = 1   # defaults
        magcolumn = -1    #  (to override, specify in first line using format #:0,1,2)  
        catlines = cat.readlines()
        cat.close()
    
    #Finds which lines are not comments (anything before a line starting with '-----' only for regular catalogs (user catalog has no comments)
    #Then gathers ra, dec, mag (takes first given magnitude) from specified columns set for each catalog
    #If magnitude not given skip entries, also skip entries that have too faint of too bright of magnitudes
    #as well as values not within proper motion parameters
    if usercat == 1:
    	comment = False
    else:
    	comment = True
    catlist = []
    
    #SDSS direct query has different format than other scat queries
    if catalog == 'sdss':
    	comment = False
    	catlines = catlines[1:]
    
    for line in catlines:
    
    	if not comment:
    		if catalog == 'sdss':
    			cline = line.split(',')
    		else:
    			cline = line.split()

    		narg = len(cline)
    		
    		if line[0:2] == '#:':
    			inlinearg = line[2:].split(',')
    			racolumn = int(inlinearg[0])-1
    			deccolumn = int(inlinearg[1])-1
    			if len(inlinearg) > 2: magcolumn = int(inlinearg[2])-1
    			continue
    		
    		if cline[racolumn].find(':') == -1:
    			ra = float(cline[racolumn])
    		else:
    			ra = astrometrystats.rasex2deg(cline[racolumn])
    			
    		if cline[deccolumn].find(':') == -1:
    			dec = float(cline[deccolumn])
    		else:
    			dec = astrometrystats.decsex2deg(cline[deccolumn])
    		
    		if magcolumn >= 0 and narg > magcolumn:
    			try:
    				mag = float(cline[magcolumn])
    			except:
    				mag = float(cline[magcolumn][0:-2])
    		else:
				mag = maxmag #Lets all magnitudes pass for user set catalog
    		
        	if usercat == 0 and narg > pmracolumn and narg > pmdeccolumn:
        		pmra = float(cline[pmracolumn])
        		pmdec = float(cline[pmdeccolumn])
        	else:
        		pmra = pmdec = 0
        		
        	if mag > maxmag: continue #don't believe anything this faint
        	if mag < minmag: continue #ignore anything this bright
        	if abs(pmra) > maxpm or abs(pmdec) > maxpm: continue		

        	iobj = Obj(ra, dec, mag) #process the line into an object
        	catlist.append(iobj)			
			#Add to catalog array in vertical stack
        				
    	if line.find('---') != -1:
    		comment = False

    #Sort by magnitude
    catlist.sort(astrometrystats.magcomp)

    return catlist