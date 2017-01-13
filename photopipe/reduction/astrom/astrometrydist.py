import numpy
import astrometrystats

"""
NAME:
	distance
PURPOSE:
	Calculates great circle distance between two points
INPUTS:
	obj1, obj2 - Must be objects created from the Obj class specified in astrometrysources.py
OUTPUTS:
	Returns the great circle distance between two points in arcseconds.
EXAMPLE:
	separation = distance(obj1, obj2)
"""

def distance(obj1, obj2):

	#Finds difference in declination and right ascension (both in radians) between the objects    
    ddec = obj2.dec_rad - obj1.dec_rad
    dra  = obj2.ra_rad - obj1.ra_rad
    
    #Calculates the great circle distance using the Haversine formula which is best suited for small angles (i.e. small separation)
    dist_rad = 2 * numpy.arcsin(numpy.sqrt( (numpy.sin(ddec/2.))**2 + numpy.cos(obj1.dec_rad) * numpy.cos(obj2.dec_rad) * (numpy.sin(dra/2.))**2))

	#Calculates and returns distance in arcseconds
    dist_deg = dist_rad * 180. / numpy.pi
    dist_sec = dist_deg * 3600.
    return dist_sec

########################################
"""
NAME:
	tooclose
PURPOSE:
	Finds objects that are within a specified minimum separation (set in arcseconds) and
	removes the dimmer of the two objects.  Returns the same list with objects removed
INPUTS:
	glist - list of objects created from the Obj class specified in astrometrysources.py
	minsep - minimum separation allowed, will remove objects closer than this value.  Set in arcsecs and default set to 3
OUTPUTS:
	Returns same list of objects with objects within minimum separation removed.
EXAMPLE:
	goodlist = tooclose(goodlist, minsep=5) 
"""
def tooclose(glist, minsep = 3, quiet = False):
    deleteset = set()
    
    #Finds great circle distances between objects in the list and removes the dimmer object
    #if less than minimum separation
    for i in range(len(glist)):
        for j in range(i+1, len(glist)):
            if i == j: continue
            dist = distance(glist[i], glist[j])
            if dist < minsep:
                if glist[i].mag > glist[j].mag: 
                   deleteset.add(i)
                else:
                   deleteset.add(j)

    deletelist = list(deleteset)
    for d in sorted(deletelist, reverse=True):
        del glist[d]
        if quiet == False:
        	print 'deleted index '+ str(d)
        
    return glist

########################################
"""
NAME:
	quickdistance
PURPOSE:
	Calculates a quick cartesian distance in arcseconds
INPUTS:
	obj1, obj2 - Must be objects created from the Obj class specified in astrometrysources.py
	cosdec - the cos(declination) that you want to calculate the distance at.  This is usually a mean or median of a set
OUTPUTS:
	Returns the cartesian distance in arcseconds
EXAMPLE:
	qdist = quickdistance(obj1,obj2, cos(dec_med))
"""
#Non-great-circle distance is much faster
def quickdistance(obj1, obj2, cosdec):
    ddec = obj2.dec - obj1.dec
    dra  = obj2.ra  - obj1.ra
    if dra > 180: dra = 360 - dra
    
    return 3600 * numpy.sqrt(ddec**2 + (cosdec*dra)**2)

########################################
"""
NAME:
	posangle
PURPOSE:
	Calculates the spherical position angle between two objects in degrees
INPUTS:
	obj1, obj2 - Must be objects created from the Obj class specified in astrometrysources.py
OUTPUTS:
	Returns position angle in degrees
EXAMPLE:
	pa = posangle(obj1, obj2)
"""
def posangle(obj1, obj2):
    
    dra  = obj2.ra_rad - obj1.ra_rad
    
    #Spherical position angle calculation in radians
    pa_rad = numpy.arctan2(numpy.cos(obj1.dec_rad)*numpy.tan(obj2.dec_rad)-numpy.sin(obj1.dec_rad)*numpy.cos(dra), numpy.sin(dra));
    pa_deg = pa_rad * 180./numpy.pi;
    pa_deg = 90. - pa_deg  				 #defined as degrees east of north
    while pa_deg > 200: pa_deg -= 360.   # make single-valued
    while pa_deg < -160: pa_deg += 360.  # note there is a crossing point at PA=200, images at this exact PA
    
    return pa_deg                        # will have the number of matches cut by half at each comparison level

########################################
"""
NAME:
	imdistance
PURPOSE:
	Calculates the cartesian difference between two objects in the image
INPUTS:
	obj1, obj2 - Must be objects created from the Obj class specified in astrometrysources.py
OUTPUTS:
	Returns difference in pixels
EXAMPLE:
	imdis = imdistance(obj1, obj2)
"""
#Pixel distance
def imdistance(obj1, obj2):
    return ((obj1.x - obj2.x)**2 + (obj1.y - obj2.y)**2)**0.5
########################################
"""
NAME: 
	calcdist
PURPOSE:
	Calculates the cartesian distance between objects within a list and excludes values not within specified radius.  
	Saves the distance and object index in lists
INPUTS:
	slist - list of objects created from the Obj class specified in astrometrysources.py
	maxrad - maximum radius (in degrees) of circle to search for other objects
	minrad - minimum radius (in degrees) of circle to search for other objects
	rascale - the cos(declination) that you want to calculate the distance at.  This is usually a mean or median of a set
OUTPUTS:
	Returns list of distances and matchids
EXAMPLE:
	(distances, matchids) = calcdist(sourcelist, 180, 10, cos(dec_med)
"""
def calcdist(slist, maxrad, minrad, rascale):
    dists = []
    matchids = []
    for i in range(len(slist)):
        d = []
        dj = []
        for j in range(len(slist)):
            if i == j: continue
            if abs(slist[i].dec - slist[j].dec) > maxrad: continue
            if rascale*abs(slist[i].ra - slist[j].ra) > maxrad: continue
            dist = quickdistance(slist[i], slist[j], rascale)
            if dist > minrad and dist < maxrad :
               d.append(dist)
               dj.append(j)
        dists.append(d)
        matchids.append(dj)
    
    return dists, matchids

####################################
"""
NAME:
	distmatch
PURPOSE:
	Matches sextractor sources with catalog sources by measuring distances to stars within 
	the same catalog (that fall within specified radius) and taking the ratio of the distances 
	between stars for both sextractor sources and catalog sources.  Only looks at stars 
	that meet required match number (default is 3). This is like a modified triangle match.  
	Then calculates position angle between matches.  Tosses out matches if beyond position 
	angle tolerance, if there are too many matches.  Returns the sextractor source that 
	matches the catalog source and the position angle between the two.
INPUTS:
	sexlist, catlist - lists of objects created from the Obj class specified in astrometrysources.py
	maxrad - maximum radius (in degrees) of circle to search for other objects
	minrad - minimum radius (in degrees) of circle to search for other objects
	reqmatch - required number of matches for each source to be accepted (how many lines need to match catalog)
	patolerance - position angle tolerance (in degrees)
	uncpa - user set uncertainty of position angle (in degrees)
	showmatches - keyword that displays match information in the end if set to 1 (type of verbose)
	fastmatch   - keyword that exits the for loop once 16 great matches are found (saves computing time for dense fields)
OUTPUTS:
	Returns lists of sextractor sources indices that match catalog source indices and the position angle between the two
EXAMPLE:
	(primarymatchs, primarymatchc, mpa) = distmatch(slist, clist, maxrad=130, minrad=20, 
	reqmatch=3, patolerance=1.2, uncpa=-1, showmatches=0, fastmatch=1)
"""
def distmatch(sexlist, catlist, maxrad=180, minrad=10, reqmatch=3, patolerance=1.2,uncpa=-1, showmatches=0, fastmatch=1, quiet=False):
    
    if reqmatch < 2:
       print 'Warning: reqmatch >=3 suggested'
    if patolerance <= 0: 
       print 'PA tolerance cannot be negative!!!'
       patolerance = abs(patolerance)
    if uncpa < 0: uncpa = 720

    declist = []
    for s in sexlist:
       declist.append(s.dec_rad)
    avdec_rad = astrometrystats.median(declist)       		# faster distance computation
    rascale = numpy.cos(avdec_rad)          # will mess up meridian crossings, however

 	#Calculates distances between objects in same list for both image catalog and reference catalog  
    (sexdists, sexmatchids) = calcdist(sexlist, maxrad, minrad, rascale)
    (catdists, catmatchids) = calcdist(catlist, maxrad, minrad, rascale)
    
    # Now look for matches in the reference catalog to distances in the image catalog.    
    countgreatmatches = 0

    smatch = []
    cmatch = []
    mpa = []
    offset = []
    offpa = []
    nmatch = []   
    primarymatchs = []
    primarymatchc = []
    
    #For each sextractor source A
    for si in range(len(sexdists)):
        sexdistarr = sexdists[si]
        sexidarr = sexmatchids[si]
        if len(sexdistarr) < 2: continue
        
        #For each catalog source B
        for ci in range(len(catdists)):
            catdistarr = catdists[ci]
            catidarr = catmatchids[ci]
            if len(catdistarr) < 2: continue
            match = 0
            smatchin = []
            cmatchin = []
            
            #For each sextractor source A, use the distances to other sextractor sources
            #that are within the radius range
            for sj in range(len(sexdistarr)):
                sexdist = sexdistarr[sj]
                newmatch = 1
                
                #For each catalog source B, use the distances to other catalog sources
                #that are within the radius range and compare the ratio to 
                #(sextractor source A - sextractor sources) to (catalog source B - catalog sources)
                #If within 5% then save the match (only count 1 match per catalog row) 
                #VLT loosened tolerance from 3% to 5%
                for cj in range(len(catdistarr)):
                    catdist = catdistarr[cj]
                    if abs((sexdist/catdist)-1.0) < 0.05:

                        match += newmatch
                        newmatch = 0 #further matches before the next sj loop indicate degeneracies
                        smatchin.append(sexmatchids[si][sj])
                        cmatchin.append(catmatchids[ci][cj])                 

			#If number of matches is above or equal to required number of matches for sextractor source A
			#(i.e. sextractor source A's distance to other sources makes it a likely candidate to be a catalog source)
			#Find difference in position angle between matched angle between sextractor source A and matched source with
			#catalog source B and matched catalog source.
			#Remove values that are larger than user specified value or if above position angle tolerance.
            if match >= reqmatch: 
                
                dpa = []
                
                # Here, dpa[n] is the mean rotation of the PA from the primary star of this match
                #  to the stars in its match RELATIVE TO those same angles for those same stars
                #  in the catalog.  Therefore it is a robust measurement of the rotation.
                for i in range(len(smatchin)):
                    ddpa = posangle(sexlist[si],sexlist[smatchin[i]]) - posangle(catlist[ci],catlist[cmatchin[i]])
                    while ddpa > 200: ddpa  -= 360.
                    while ddpa < -160: ddpa += 360.
                    dpa.append(ddpa)

                #If user was confident the initial PA was right, remove bad PA's right away
                for i in range(len(smatchin)-1,-1,-1):
                    if abs(dpa[i]) > uncpa: 
                        del smatchin[i]
                        del cmatchin[i]
                        del dpa[i]
                    
                if len(smatchin) < 2: continue
                
                #Finds mode of difference position angle but allowing values to be +/- patolerance    
                dpamode = astrometrystats.most(dpa, vmin=patolerance*3, vmax=patolerance*3)
                
                #Remove deviant matches by PA
                for i in range(len(smatchin)-1,-1,-1):
                    if abs(dpa[i] - dpamode) > patolerance:
                        del smatchin[i]
                        del cmatchin[i]
                        del dpa[i]
                				
                if len(smatchin) < 2: continue
                
                #Finds the number of degenerate matches
                ndegeneracies = len(smatchin)-len(astrometrystats.unique(smatchin)) + len(cmatchin)-len(astrometrystats.unique(cmatchin))
                    # this isn't quite accurate (overestimates if degeneracies are mixed up)

				#Save values
                mpa.append(dpamode)
                primarymatchs.append(si)
                primarymatchc.append(ci)
                smatch.append(smatchin)
                cmatch.append(cmatchin)
                nmatch.append(len(smatchin)-ndegeneracies)
                
                #If the number of unique sextractor sources is greater than 6, count as a great match
                if (len(smatchin)-ndegeneracies > 6): countgreatmatches += 1
        
        #Breaks loop if over 16 great matches are found and fastmatch
        if countgreatmatches > 16 and fastmatch == 1: break #save processing time
    
    #If no matches found end program and return empty lists
    nmatches = len(smatch)
    if (nmatches == 0):
        print 'Found no potential matches of any sort (including pairs).'
        print 'The algorithm is probably not finding enough real stars to solve the field.  Check seeing.'
        return [], [], []   
    
    #Get rid of matches that don't pass the reqmatch cut (2nd cut after removing bad position angles)
    for i in range(len(primarymatchs)-1,-1,-1):
        if nmatch[i] < reqmatch:
            del mpa[i]
            del primarymatchs[i]
            del primarymatchc[i]
            del smatch[i]
            del cmatch[i]
            del nmatch[i]

	#If no remaining matches then exit program
    if len(smatch) < 1:
        print 'Found no matching clusters of reqmatch =', reqmatch
        return [], [], []

    #If we still have lots of matches, get rid of those with the minimum number of submatches
    #(that is, increase reqmatch by 1)
    minmatch = min(nmatch)
    countnotmin = 0
    
    #Finds how many matches have more than the minimum require matches
    for n in nmatch:
       if n > minmatch: countnotmin += 1
    
    #If the number of matches is above 16 and there are more than 3 sources with more than the 
    #required number of matches, then delete any source with the bare minimum number of matches
    if len(nmatch) > 16 and countnotmin > 3:
        print 'Too many matches: increasing reqmatch to', reqmatch+1
        for i in range(len(primarymatchs)-1,-1,-1):
            if nmatch[i] == minmatch:
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
          
    nmatches = len(smatch) # recalculate with the new reqmatch and with prunes supposedly removed
    if quiet == False:
    	print 'Found',nmatches,'candidate matches.'

    # Kill the bad matches
    rejects = 0
    
    # Use only matches with a consistent PA (finds mode counting those within 3 patolerance)
    offpa = astrometrystats.most(mpa, vmin=3*patolerance, vmax=3*patolerance)
    
    #Removes values with rotational positional offsets that are above tolerance, then removes
    #values that are above 2 sigma of those values
    if len(smatch) > 2:    

        #Coarse iteration for anything away from the mode
        for i in range(len(primarymatchs)-1,-1,-1):
            if abs(mpa[i] - offpa) > patolerance:
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                rejects += 1
           
        medpa = astrometrystats.median(mpa)
        stdevpa = astrometrystats.stdev(mpa)
        refinedtolerance = (2.0 * stdevpa) #VLT Changed from arbitrary value of 2.2 from original script

        #Fine iteration to flag outliers now that we know most are reliable
        for i in range(len(primarymatchs)-1,-1,-1):
            if abs(mpa[i] - offpa) > refinedtolerance:
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                rejects += 1  #these aren't necessarily bad, just making more manageable.

    # New verification step: calculate distances and PAs between central stars of matches
    ndistflags = [0]*len(primarymatchs)
    for v in range(2):  #two iterations

        if len(primarymatchs) == 0: break

        #Find distances between central stars of matches and compare sextractor source distances to catalog sources
        for i in range(len(primarymatchs)):
            for j in range(len(primarymatchs)):
                if i == j: continue
                si = primarymatchs[i]
                ci = primarymatchc[i]
                sj = primarymatchs[j]
                cj = primarymatchc[j]
    
                sexdistij = distance(sexlist[si], sexlist[sj])
                catdistij = distance(catlist[ci], catlist[cj])
                
                try:
                   if abs((sexdistij/catdistij)-1.0) > 0.05: #VLT loosened tolerance from 3% to 5%
                      ndistflags[i] += 1
                except:  # (occasionally will get divide by zero)
                   pass

        #Delete bad clusters that were flagged for every match
        ntestmatches = len(primarymatchs)
        for i in range(ntestmatches-1,-1,-1):
            if ndistflags[i] == ntestmatches-1:   #if every comparison is bad, this is a bad match
                del mpa[i]
                del primarymatchs[i]
                del primarymatchc[i]
                del smatch[i]
                del cmatch[i]
                del nmatch[i]
                rejects += 1
                
    nmatches = len(primarymatchs)

    if quiet == False:
    	print 'Rejected', rejects, 'bad matches.'
    	print 'Found', nmatches, 'good matches.'

	#If no remaining matches, return empty lists
    if nmatches == 0:
       return [], [], []

    #Returns pixel scale (great circle distance in catalog [ra, dec]/cartesian distance in sextractor source [x,y])
    pixscalelist = []
    if len(primarymatchs) >= 2:
        for i in range(len(primarymatchs)-1):
            for j in range(i+1,len(primarymatchs)):
                si = primarymatchs[i]
                ci = primarymatchc[i]
                sj = primarymatchs[j]
                cj = primarymatchc[j]
                
                try:
                    pixscalelist.append(distance(catlist[ci],catlist[cj])/imdistance(sexlist[si],sexlist[sj]))
                except:
                    pass
                    
        pixelscale = astrometrystats.median(pixscalelist)
        pixelscalestd = astrometrystats.stdev(pixscalelist)
        
        if quiet == False:
        	if len(primarymatchs) >= 3:
        		print 'Refined pixel scale measurement: %.4f"/pix (+/- %.4f)' % (pixelscale, pixelscalestd)
        	else:
        		print 'Refined pixel scale measurement: %.4f"/pix' % pixelscale
           
	#If showmatches keyword set then print which objects match
    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        ci = primarymatchc[i]
        if quiet == False:
        	print  '%3i' % si, 'matches', '%3i' % ci, ' (dPA =%7.3f)' % mpa[i],
        
        #Keyword set in main program
        if showmatches:
           print
           if len(smatch[i]) < 16:
              print '  ', si, '-->', smatch[i], 
              if len(smatch[i]) >= 7: print
              print '  ', ci, '-->', cmatch[i]
           else:
              print '  ', si, '-->', smatch[i][0:10], '+', len(smatch[i])-10, 'more'
              print '  ', ci, '-->', cmatch[i][0:10], '+'#, len(cmatch[i])-10, ' more'
           if i+1 >= 10 and len(primarymatchs)-10 > 0: 
              print (len(primarymatchs)-10), 'additional matches not shown.'
              break
        else:
           if quiet == False:
              print ':', str(len(smatch[i])).strip(), 'rays'
      
	#Create region files for DS9 with the sextractor sources (matchlines.im.reg) and catalog (matchlines.wcs.reg)
    out = open('matchlines.im.reg','w')
    i = -1
    color='red'
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    out.write('image\n')
    for i in range(len(primarymatchs)):
        si = primarymatchs[i]
        for j in range(len(smatch[i])):
           sj = smatch[i][j] 
           out.write("line(%.3f,%.3f,%.3f,%.3f) # line=0 0\n" % (sexlist[si].x, sexlist[si].y, sexlist[sj].x, sexlist[sj].y))
    out.close()

    out = open('matchlines.wcs.reg','w')
    i = -1
    color='green'
    out.write('# Region file format: DS9 version 4.0\nglobal color='+color+' font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    out.write('fk5\n')
    for i in range(len(primarymatchs)):
        ci = primarymatchc[i]
        for j in range(len(smatch[i])):
           cj = cmatch[i][j]
           out.write("line(%.5f,%.5f,%.5f,%.5f) # line=0 0\n" % (catlist[ci].ra, catlist[ci].dec, catlist[cj].ra, catlist[cj].dec))
    out.close()

    #future project: if not enough, go to the secondary offsets   
    
    #Returns sextractor sources and catalog sources that appear to have matches along with the mode of the position angle
    return (primarymatchs, primarymatchc, mpa)