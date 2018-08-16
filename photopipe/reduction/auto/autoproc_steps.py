import glob
import os
import pyfits as pf
import numpy as np
import autoproc_depend as apd
from astropy import wcs
import re
import datetime
from astropy.time import Time
import sys

inpipevar = {'autoastrocommand':'autoastrometry', 'getsedcommand':'get_SEDs', 
			'sexcommand':'sex' , 'swarpcommand':'swarp' , 'rmifiles':0,  
			'prefix':'', 'datadir':'' , 'imworkingdir':'' , 'overwrite':0 , 'verbose':1, 
			'flatfail':'' , 'fullastrofail':'' ,
			'pipeautopath':'' , 'refdatapath':'', 'defaultspath':'' }

def autopipedefaults(pipevar=inpipevar):

    """
    NAME:
        autopipedefaults
    PURPOSE:
        Sets commonly used variables for pipeautoproc to use throughout each step
        Uses pipeautoproc.par to set variables, otherwise set to default values
        Saves in a dictionary
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                    but can be set to default) 
    EXAMPLE:
        autopipedefaults(pipevar=inpipevar)
    """

    print 'Setting pipeline parameters (DEFAULTS)'
    
    path = os.path.dirname(__file__)
        
    pipevar['pipeautopath'] = path
    sfile = path+'/pipeautoproc.par'
    
    if os.path.isfile(sfile):
        f = open(sfile,'r')
        for line in f.readlines():
            line = ''.join(line.split())
            colpos = line.find(':')
            
            keyword = line[0:colpos]
            value = line[colpos+1:]
            pipevar[keyword] = value
        f.close()

    if pipevar['refdatapath'] == '': 
        pipevar['refdatapath'] = pipevar['pipeautopath']+'/refdata'

    if pipevar['defaultspath'] == '': 
        pipevar['defaultspath'] = pipevar['pipeautopath']+'/defaults'

    if pipevar['imworkingdir'] != '' and not(os.path.exists(pipevar['imworkingdir'])): 
        print 'Creating imaging working directory: ',  pipevar['imworkingdir']
        os.makedirs(pipevar['imworkingdir'])
        
def autopipeprepare(pipevar=inpipevar):

    """
    NAME:
        autopipeprepare
    PURPOSE:
        Runs pipeprepare on every valid file and saves files with prefix 'p'.  Changes 
        header with more manageable keywords and does bias/dark subtraction if bias/dark 
        master exists (compares header keywords in files and bias/dark master)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeprepare(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.pipeprepare()
    """
    
    print 'PREPARE'
    
    # Looks for existing files in given data directory using prefix
    files = glob.glob(pipevar['datadir'] + pipevar['prefix'] + '*.fits')
    pfiles = glob.glob(pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    if pipevar['verbose']: print 'Found', len(files), 'files'
    
    # Finds any master bias files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    biasfiles = glob.glob(pipevar['imworkingdir'] + 'bias*')
    biascamera = []
    if len(biasfiles) > 0:
        for bfile in biasfiles:
            head = pf.getheader(bfile)  
            camera = int(head['CAMERA'])
            biascamera += [camera]
            
    # Finds any master dark files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    darkfiles = glob.glob(pipevar['imworkingdir'] + 'dark*')
    darkcamera = []
    if len(darkfiles) > 0:
        for dfile in darkfiles:
            head = pf.getheader(dfile)  
            camera = int(head['CAMERA'])
            darkcamera += [camera]     
        
    # For each file (that doesn't have an existing p file or can be overwritten), 
    # run pipeprepare on it with output file being saved into the imworkingdir, 
    # will run bias subtraction if bias master available (checks based on how bias 
    # file and data file are named
    for file in files:

        fileroot = os.path.basename(file)
        outnameim = pipevar['imworkingdir'] + 'p' + fileroot
        
        head = pf.getheader(file)
        camera = int(head['CAMERA'])
                
        try:
            bcamloc  = biascamera.index(camera)    
            biasfile = biasfiles[bcamloc]      
        except:
            biasfile = None
        
        try:
            dcamloc  = darkcamera.index(camera)    
            darkfile = darkfiles[dcamloc]      
        except:
            darkfile = None

        if (outnameim not in pfiles) or (pipevar['overwrite'] != 0):
            apd.pipeprepare(file, outname=outnameim, biasfile=biasfile, darkfile=darkfile,
                             verbose = pipevar['verbose'])
        else:
            print 'Skipping prepare. File already exists'
            
def autopipeimflatten(pipevar=inpipevar):
    
    """
    NAME:
        autopipeflatten
    PURPOSE:
        Flatten data using flat with matching filter name
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeflatten(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.flatpipeproc()
    """
    
    print 'FLATTEN'
    
    # Finds prepared files and checks to see if there are any existing flattened files
    # Find flats in imworkingdir with name flat somewhere in a fits file name
    files  = glob.glob(pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
    ffiles = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    flats  = glob.glob(pipevar['imworkingdir'] + '*flat*.fits')

    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    # If there are flats, then grab the filter from each of them, 
    # otherwise end program
    flatfilts = []
    if len(flats) > 0:
        for flat in flats:
            head = pf.getheader(flat)
            filter = head['FILTER']
            flatfilts += [filter]
    else:
        print 'No flats found for any filter!'
        return
    
    # Create outfile name and check to see if outfile already exists.  If it doesn't or
    # overwrite enabled then take filter from file and find where the flat filter matches
    # If no flats match filter, store in pipevar.flatfail, otherwise run flatproc on file
    for file in files:
        print file
        fileroot = os.path.basename(file)
        outnameim = pipevar['imworkingdir'] + 'f' + fileroot
        
        if (outnameim not in ffiles) or (pipevar['overwrite'] != 0):
            head = pf.getheader(file)
            filter = head['FILTER']

            try:
                flatfileno = flatfilts.index(filter)
            except:
                print 'Flat field not found for '+ file +' (filter='+filter+')'
                pipevar['flatfail'] += ' ' + file
                continue
            
            flatfile = flats[flatfileno]
        
            if pipevar['verbose']: print 'Flattening', file, 'using', flatfile
            
            apd.flatpipeproc(file, flatfile, flatminval=0.3)
        
        else:
            print 'Skipping flatten. File already exists'
            
    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        
def autopipemakesky(pipevar=inpipevar):
    """
    NAME:
        autopipemakesky
    PURPOSE:
        Combine sky flats based on filter type (sigma clipping for sources)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipemakesky(pipevar=inpipevar)
    DEPENDENCIES:
        astroproc_depend.skypipecombine, astroproc_depend.medclip, Sextractor
    FUTURE IMPROVEMENTS:
        skypipecombine slow, find better algorithm
    """   
    
    print 'MAKE SKY'
    
    # Copies necessary parameter file for sextractor if not in current working directory
    if not os.path.isfile('source.param'): 
        os.system('cp ' + pipevar['defaultspath'] + '/source.param .')
    if not os.path.isfile('sex_source.config'): 
        os.system('cp ' + pipevar['defaultspath'] + '/sex_source.config .')
    if not os.path.isfile('sex.conv'): 
        os.system('cp ' + pipevar['defaultspath'] + '/sex.conv .')
    if not os.path.isfile('defaulf.nnw'): 
        os.system('cp ' + pipevar['defaultspath'] + '/default.nnw .')
        
    # Finds files with given prefix
    files  = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    filters = []
    for file in files:
        head = pf.getheader(file)
        filter = head['FILTER']
        filters += [filter]
    filters = np.array(filters)
    
    # Unique list of filters
    filterlist = set(filters)
    
    # For each unique filter, combine sky files using skycombine if more than 2 files
    # Otherwise return list of unprocessed files
    for filt in filterlist:
        skyflats = np.where(filters == filt)
        outflatname = pipevar['imworkingdir'] + 'sky-' + filt + '.fits'
        
        if len(skyflats[0]) >= 2:
        
            if os.path.isfile(outflatname) and pipevar['overwrite'] == 0:
                print 'Skipping makesky for '+filt+'. File already exists'
                continue

            files = np.array(files)
            if pipevar['verbose']:
                print filt, '-band sky flats.'
                print files[skyflats]
            
                apd.skypipecombine(files[skyflats], outflatname, file,
                    pipevar, removeobjects=True, type='sky')
        else:
            print 'Unable to produce a flat field for this setting: ' + filt
            print 'Will not be able to further process ' + str(len(skyflats)) + \
                  ' image(s) without a flat from another source:'

            for i in np.arange(len(skyflats[0])):
                print '    ' + str(files[skyflats[i]])           

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')

def autopipeskysubmed(pipevar=inpipevar):
    """
    NAME:
        autopipeskysubmed
    PURPOSE:
        Subtract median, does NOT use master sky
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeskysub(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.skypipeproc
    """
    
    print 'SKY-SUBTRACT MEDIAN ONLY'
    
    # Find data that needs to be sky subtracted
    files  = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    sfiles = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    # For each file if output files don't exist or override set check if we have master 
    # skyflat for filter, sky subtract if it exists using skypipeproc
    for file in files:
        print file
        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 's' + fileroot    
    
        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print 'Skipping sky subtraction for '+file+'. File already exists'
            continue
            
        head = pf.getheader(file)
        data = pf.getdata(file)
        
        data -= np.median(data)

        filter = head['FILTER'] 
        
        if pipevar['verbose'] > 0:
            print 'Sky subtracting median only'
            
        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by skypipeprocmed ' + date)   
             
        pf.writeto(outfile, data, head, clobber=pipevar['overwrite'])

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')

        
def autopipeskysub(pipevar=inpipevar):
    """
    NAME:
        autopipeskysub
    PURPOSE:
        Subtracts master sky flat from data and subtracts median.
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeskysub(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.skypipeproc
    """
    
    print 'SKY-SUBTRACT'
    
    # Find data that needs to be sky subtracted
    files  = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    sfiles = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    skys = glob.glob(pipevar['imworkingdir'] + '*sky-*.fits')
    
    if len(skys) == 0:
        print 'No master sky files found, cannot sky subtract'
        return
        
    # Find the associated filter of each master skyflat
    skyfilts = []
    for sky in skys:
        head = pf.getheader(sky)
        filter = head['FILTER']
        skyfilts += [filter]
    
    # For each file if output files don't exist or override set check if we have master 
    # skyflat for filter, sky subtract if it exists using skypipeproc
    for file in files:

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 's' + fileroot    
    
        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print 'Skipping sky subtraction for '+file+'. File already exists'
            continue
            
        head = pf.getheader(file)
        filter = head['FILTER'] 

        # Find corresponding master skyflat
        try:
            skyloc  = skyfilts.index(filter)    
            skyfile = skys[skyloc]      
        except:
            print 'Sky field not found for ', file
            pipevar['flatfail'] += ' ' + file
            continue
        
        if pipevar['verbose'] > 0:
            print 'Sky subtracting', file, 'using', skyfile
        
        apd.skypipeproc(file, skyfile, outfile)


    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')

def autopipecrcleanim(pipevar=inpipevar):
    
    """
    NAME:
        autopipecrcleanim
    PURPOSE:
        Removes cosmic rays
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipecrcleanim(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.cosmiczap
    FUTURE IMPROVEMENTS:
        Slow, alter cosmics.py?
        Get readnoise from header
    """
    
    print 'CRCLEAN'

    # Find data that needs to be cosmic ray zapped
    files  = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
    zfiles = glob.glob(pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
 
    # For each file check that objects meet count limits and exposure time
   	# (i.e. short exposure time with lot of counts will be ignored), also targets that are
   	# calibration files will be ignored.
   	# Run cosmiczap on the files and have output files be 'z'+file plus weight files
    
    for file in files:
    
        head = pf.getheader(file)
                
        try:
            target  = head['TARGNAME']
        except:
            print 'Requires header keywords: TARGNAME. Check file.'
            continue
        
        if 'flat' in target.lower(): continue
        if 'twilight' in target.lower(): continue

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 'z' + fileroot 
          
        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print 'Skipping crzap for '+file+'. File already exists'
            continue 
                
        if pipevar['verbose'] > 0:
            print 'Cleaning cosmic rays from', file      
        
        # Runs cosmics.py
        apd.cosmiczap(file, outfile, sigclip=6.0, maxiter=1, verbose=pipevar['verbose'])
        
        
    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')

def autopipeastrometry(pipevar=inpipevar):
    """
    NAME:
        autopipepipeastrometry
    PURPOSE:
        Calculate astrometry of image files to fix WCS coordinates (shift and rotation) 
        in header. Using fast astrometry solver (vlt_autoastrometry.py) that using 
        pair-distance matching and asterism matching.  Returns file with corrected WCS 
        coordinates saved as 'a'+fitsfile. Run Scamp for additional astrometry 
        corrections, twice, once for basic individual LOOSE correction, second correct all
        together.  Uses distortion of 3 as default, but uses 7 if distortion parameters 
        high (i.e. RATIR H2RG)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeastrometry(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.astrometry, vlt_autoastrometry.py, scamp, sextractor
    FUTURE IMPROVEMENTS:
        Better distinction between first and second scamp run
    """
    
    print 'ASTROMETRY'
    
    files = glob.glob(pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')
    
    # If no files, look for those that were not cosmic ray zapped
    if len(files) == 0:
        files = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    # Calculate relative astrometric solution
    for file in files:
        
        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 'a' + fileroot
        
        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print 'Skipping astrometry for '+file+'. File already exists' 
            continue
        
        head = pf.getheader(file)
        
        targ = head['TARGNAME']
        sat  = head['SATURATE']
        
        if 'flat' in targ: continue
        
        cmd = 'python ' + pipevar['autoastrocommand'] + ' ' + file + ' -l ' + str(sat)
        
        # Run direct astrometry
        if pipevar['verbose'] > 0:
            os.system(cmd)
            print cmd
        else:
            os.system(cmd + ' -q')
            
        if not os.path.isfile(outfile):
            pipevar['fullastrofail'] += ' ' + file

    if not os.path.isfile('astrom.param'): 
        os.system('cp ' + pipevar['defaultspath'] + '/astrom.param .')
    if not os.path.isfile('astrom.conv'): 
        os.system('cp ' + pipevar['defaultspath'] + '/astrom.conv .') 
    if not os.path.isfile('default.sex'): 
        os.system('cp ' + pipevar['defaultspath'] + '/default.sex .') 
    if not os.path.isfile('default.missfits'): 
        os.system('cp ' + pipevar['defaultspath'] + '/default.missfits .') 
    if not os.path.isfile('scamp.conf'): 
        os.system('cp ' + pipevar['defaultspath'] + '/scamp.conf .')          
    
    # Calculate astrometry again using Scamp. First identify objects using sextractor, 
    # then Scamp will solve by comparing reference catalog (currently set by default to 
    # SDSS) to sources found by sextractor. Adds WCS corrections and second astrometry 
    # parameters to header
    afiles = glob.glob(pipevar['imworkingdir'] + 'azsfp' + pipevar['prefix'] + '*.fits')
    
    # If no files, look for those that were not cosmic ray zapped
    if len(afiles) == 0:
        afiles = glob.glob(pipevar['imworkingdir'] + 'asfp' + pipevar['prefix'] + '*.fits')

    if len(afiles) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    afiletarg = []
    afilefilt = []
    
    for afile in afiles:
        head = pf.getheader(afile)
        afiletarg += [head['TARGNAME']]
        afilefilt += [head['FILTER']]
    
    atargets = set(afiletarg)
    afilters = set(afilefilt)
    
    afiletarg = np.array(afiletarg)
    afilefilt = np.array(afilefilt) 
    afiles    = np.array(afiles) 
    
    for atarg in atargets:
        for afilt in afilters:
        
            thisatarget = np.where(np.logical_and(atarg == afiletarg, afilt == afilefilt))
            atfimages = afiles[thisatarget]
            
            try:
                head = pf.getheader(atfimages[0])
            except:
                continue
                
            # If scamp has already been run, skip
            try:
                test = head['ASTIRMS1']
                print 'Skipping scamp astrometry for: ',atarg,afilt,' Files already exist' 
                continue 
            except:                    
                # Run sextractor to find sources, then use those catalogs to run scamp
                # with loose fitting constraints
                apd.astrometry(atfimages, scamprun=1, pipevar=pipevar)
            
                # Do same thing again but with more stringent scamp parameters
                apd.astrometry(atfimages, scamprun=2, pipevar=pipevar)
                   
                
    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits, zsfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')

def autopipestack(pipevar=inpipevar, customcat=None, customcatfilt=[]):
    """
    NAME:
        autopipepipestack
    PURPOSE:
    	Does zeropoint correction on each individual frame using sextractor and get_SEDs. 
        Creates flux scale (newflxsc) from how close to median of zeropoint values.  Uses
        flux scale to stack images in Swarp (has moved bad zeropoint values and bad newflxsc
        values to marked folders - badzptfit/ and badflxsc/) and calculates absolute zeropoint 
        correction of coadd. Saves zeropoint plot as zpt_(FILTER).ps
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipestack(pipevar=inpipevar)
    DEPENDENCIES:
        SWarp, get_SEDs, calc_zpt, findsexobj (sextractor)
    """
  
    print 'STACK'
    
    os.system('export CDSCLIENT=http') #Fix for problem with timeout with CDSCLIENT
    
    qtcmd = 'True'; quiet = 1
    if pipevar['verbose'] > 0: quiet = 0; qtcmd = 'False'
    
    # If swarp configuration file ('default.swarp') does not exist, move swarp 
    # output default configuration file
    if not os.path.isfile('default.swarp'): 
        os.system(pipevar['swarpcommand']+' -d > default.swarp')
    
    # Find files that have had astrometry performed on them, stop program if don't exist
    files = glob.glob(pipevar['imworkingdir'] + 'a*sfp' + pipevar['prefix'] + '*.fits')
    print pipevar['imworkingdir'] + 'a*sfp' + pipevar['prefix'] + '*.fits'
    print files
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return

    filetargs = []; fileexpos = []; filefilts = []; fileairmv = [] 
    filesatvs = []; filearms1 = []; filearms2 = []; filetime  = []
    
    # Grab information in the headers of astrometry corrected file and save to array
    for i,file in enumerate(files):
        head = pf.getheader(file)
        obstime = Time(head['DATE-OBS'], format='isot', scale='utc')
        
        
        # Strip target name of whitespace
        filetargs += [re.sub(r'\s+', '', head['TARGNAME'])]; 
        fileexpos += [head['EXPTIME']]; filefilts += [head['FILTER']]
        fileairmv += [head['AIRMASS']]; filesatvs += [head['SATURATE']]
        filearms1 += [head['ASTRRMS1']]; filearms2 += [head['ASTRRMS2']]
        filetime  += [obstime.jd]
    
    files     = np.array(files); filetargs = np.array(filetargs)
    fileexpos = np.array(fileexpos); filefilts = np.array(filefilts)
    filesatvs = np.array(filesatvs); fileairmv = np.array(fileairmv)
    filearms1 = np.array(filearms1); filearms2 = np.array(filearms2)
    filetime  = np.array(filetime)
    targets = set(filetargs)
    
    # Dictionary of corresponding columns for catalog file
    catdict = {'u': 2, 'g': 3, 'r': 4, 'i': 5, 'z': 6, 'y': 7, 
               'B': 8, 'V': 9, 'R':10, 'I':11, 'J':12, 'H':13, 'K':14,
               'ue':15,'ge':16,'re':17,'ie':18,'ze':19,'ye':20,
               'Be':21,'Ve':22,'Re':23,'Ie':24,'Je':25,'He':26,'Ke':27,'mode':28}
    
    # Finds files with same target and the filters associated with this target
    for targ in targets:

        thistarget = np.where(filetargs == targ)
        if len(thistarget) == 0: continue
        
        thistargetfilts = set(filefilts[thistarget])
        
        # Find files that have the same target and same filter and store information 
        # on the exposure times and airmass. Only use good Scamp astrometric fit files
        for filter in thistargetfilts:
            stacki = (filetargs == targ) & (filefilts == filter) &\
                     (filearms1 < 2.0e-4) & (filearms1 > 5.0e-6) &\
                     (filearms2 < 2.0e-4) & (filearms2 > 5.0e-6)
                                            
            if sum(stacki) == 0: continue

            stacklist = files[stacki]
            stackexps = fileexpos[stacki]
            stackairm = fileairmv[stacki]
            stacktime = filetime[stacki]
            
            medexp = np.median(stackexps)
            medair = np.median(stackairm)
            minair = min(stackairm)
            maxair = max(stackairm)
            totexp = sum(stackexps)
            nstack = len(stacklist)
            firsttime = Time(stacktime[0], format='jd', scale='utc').isot
            lasttime  = Time(stacktime[-1], format='jd', scale='utc').isot
            medtime   = Time(np.median(stacktime), format='jd', scale='utc').isot
                        
            zpts = []
            
            # Find stars for each individual frame and try to find matches with coadded 
            # frame with each frame optimized with PSF size
            for sfile in stacklist:
                head = pf.getheader(sfile)
                ipixscl = head['PIXSCALE']
                
                apd.findsexobj(sfile, 3.0, pipevar,pix=ipixscl,aperture=20.0, quiet=quiet)
                starfile = sfile + '.stars'
                
                svars = np.loadtxt(starfile, unpack=True)
                xim  = svars[1,:]
                yim  = svars[2,:]
                mag  = svars[3,:]
                mage = svars[4,:]
                flag = svars[5,:]
                elon = svars[8,:]
                
                # astropy does not like SWarp PV keywords or unicode, temporarily delete
                for key in head.keys():
                    if any(x in key for x in ['PV1_', 'PV2_', 'COMMENT', 'HISTORY']):
                        try:
                            del head[key]
                        except KeyError as error:
                            print error
                        
                w = wcs.WCS(head)
                wrd = w.all_pix2world(np.transpose([xim, yim]), 0)
                imfile  = sfile + '.im'
                catfile = sfile + '.cat'
                
                # Save stars from image
                np.savetxt(imfile, np.transpose([wrd[:,0],wrd[:,1],mag]))
                
                # Filter name correction:
                if filter == 'Z' or filter == 'Y': filter = filter.lower()
                
                if 'SDSS' in filter:
                    filter = filter[-1].lower()

                nocustomcat = False
                # If custom catalog provided, match the same objects as the *.im file
                # and create refmag (reference magnitudes) that have the same matched
                # indices                    
                if customcat != None and filter in customcatfilt:

                    print 'USING CUSTOM CATALOG'
                    in_data = np.loadtxt(imfile)
                    input_coords = in_data[:, :2]
                    cat_data = np.loadtxt(customcat, skiprows=1)
                    cat_coords = cat_data[:,:2]
                    
                    cat_matches, tmp = apd.identify_matches(input_coords, cat_coords)
                    
                    refmag = np.zeros(len(mag)) + 99
                    mode = np.zeros(len(mag)) + -1
                    for i, i_ind in enumerate(cat_matches):
                        if i_ind > 0:
                            #print input_coords[i], cat_coords[i_ind]
                            refmag[i] = cat_data[i_ind,catdict[filter]]
                            mode[i] = 4
                    
                    # If no matching indices, run with the standard catalogs
                    if sum(refmag<90.0) == 0: nocustomcat = True
                else:
                    nocustomcat = True
                    
                # If custom catalog not provided, catalog doesn't include filter, or 
                # no objects from catalog found in image then
                # use get_SEDs.py to make catalog using 2MASS + (SDSS or APASS or USNOB1)
                if nocustomcat:
                    # Create catalog star file 
                    # (python get_SEDs.py imfile filter catfile USNOB_THRESH alloptstars)
                    sedcmd = 'python ' + pipevar['getsedcommand'] + ' ' + imfile + ' ' +\
                         filter + ' ' + catfile + " 15 True "+ qtcmd
                
                    if pipevar['verbose'] > 0: print sedcmd
                    os.system(sedcmd)
                
                    if not os.path.isfile(catfile):
                        zpts += [float('NaN')]
                        continue
                
                    # Read in catalog file
                    cvars = np.loadtxt(catfile, unpack=True)
                    refmag = cvars[catdict[filter],:]
                    mode   = cvars[catdict['mode'],:]
                    

                
                # Find catalog filter values and only cutoff values of actual detections
                goodind = (mode != -1) & (refmag < 90.0) & (flag < 8) & (elon <= 1.5)
                
                refmag = refmag[goodind]
                obsmag = mag[goodind]
                obserr = mage[goodind]
                obswts = np.zeros(len(obserr))
                obskpm = np.zeros(len(obsmag))
                
                # Store magnitudes and weights (with minimum magnitude error of 0.01)
                for i in np.arange(len(obsmag)):
                    if obserr[i] < 0.1:
                        obskpm[i] = obsmag[i]
                        obswts[i] = 1.0/(max(obserr[i], 0.01)**2)

                if len(refmag) > 0 and len(obskpm) > 0 and len(obswts) > 0:
                    zpt, scats, rmss = apd.calc_zpt(np.array([refmag]), np.array([obskpm]), 
                                                np.array([obswts]), sigma=3.0)
                                
                    # Reload because we had to remove distortion parameters before
                    head = pf.getheader(sfile)
                    data = pf.getdata(sfile)
                    head['ABSZPT']   = (zpt[0] + 25.0, 'Relative zeropoint from calc_zpt')
                    head['ABSZPTSC'] = (scats[0], 'Robust scatter of relative zeropoint')
                    head['ABSZPRMS'] = (rmss[0], 'RMS of relative zeropoint')

                    pf.update(sfile, data, head)
                    zpts += zpt     
                else:
                    zpts = [np.inf]        
                
            # Move files with bad zeropoint calculations to folder 'badzptfit' 
            # and do not use those frames
            zpts = np.array(zpts)
            goodframes = np.isfinite(zpts)
            badframes  = ~np.isfinite(zpts)
            
            if len(zpts[badframes]) != 0:
                if not os.path.exists(pipevar['imworkingdir']+'/badzptfit'):
                    os.makedirs(pipevar['imworkingdir']+'/badzptfit')
                for file in stacklist[badframes]:
                    os.system('mv ' + file + ' ' +  pipevar['imworkingdir']+'/badzptfit/')
                zpts = zpts[goodframes]
                newstack = stacklist[goodframes]
            else:
                newstack = stacklist
            
            badnewflxsc = []
            # Add relative zeropoint values to headers and calculate flux scale. 
            # Remove unphysical fluxscale files
            medzp = np.median(zpts)
            for i,file in enumerate(newstack):
                head = pf.getheader(file)
                head['NEWFLXSC'] = (1.0/(10.0**( (zpts[i] - medzp)/2.5 )), 
                                    'Flux scaling based on median zp') 
                
                if 1.0/(10.0**( (zpts[i] - medzp)/2.5 )) < 0.1:
                    badnewflxsc += [file]
            
                data = pf.getdata(file)
                pf.update(file, data, head)
            
            removedframes = []
            # Removes files that have bad newflxsc values and removes from stack list
            if len(badnewflxsc) > 0:
                if not os.path.exists(pipevar['imworkingdir']+'/badflxsc'):
                    os.makedirs(pipevar['imworkingdir']+'/badflxsc')
                
                for ibad in badnewflxsc:    
                    os.system('mv ' + ibad +' '+ pipevar['imworkingdir']+'badflxsc/')
                
                removedframes += badnewflxsc
            
                # Remove files that have bad newflxsc values from list of stack
                bad = set(badnewflxsc)
                newstack = [x for x in newstack if x not in bad]                       

            newtextslist = ' '.join(newstack)   
            
            stackcmd = pipevar['swarpcommand']
            
            # Keywords to carry through
            stackcmd += ' -COPY_KEYWORDS OBJECT,TARGNAME,FILTER,' +\
                        'INSTRUME,PIXSCALE,WAVELENG,DATE-OBS,AIRMASS,FLATFLD,FLATTYPE '
                             	
            # Create output variables that will be used by SWarp
            outfl = pipevar['imworkingdir'] + 'coadd' + targ + '_'+ re.sub(r'[^\w]', '', medtime)+'_'+ filter + '.fits'
            outwt = pipevar['imworkingdir'] + 'coadd' + targ + '_'+ re.sub(r'[^\w]', '', medtime)+'_'+ filter + '.weight.fits'
            
            if pipevar['verbose'] > 0:
                stackcmd += ' -VERBOSE_TYPE NORMAL '
            else:
                stackcmd = stackcmd + ' -VERBOSE_TYPE QUIET '
            
            # Coadd with flux scale
            stackcmd = stackcmd + ' -SUBTRACT_BACK N -WRITE_XML N -IMAGEOUT_NAME ' +\
                       outfl + ' -WEIGHTOUT_NAME ' + outwt +\
                       ' -FSCALE_KEYWORD NEWFLXSC ' + newtextslist
                       
            if pipevar['verbose'] > 0:
                print stackcmd
            
            os.system(stackcmd)
            head   = pf.getheader(outfl)
            pixscl = head['PIXSCALE']
            
            try:
                apd.findsexobj(outfl, 10.0, pipevar, pix=pixscl, aperture=20.0,
                           wtimage=outwt, quiet=quiet)
            except:
                sys.exit('Problem opening coadd fits file, may need to coadd in smaller bin size')
                           
            head   = pf.getheader(outfl)
            cpsfdi = 1.34 * float(head['SEEPIX'])
            
            # Run sextractor again on new coadd file
            apd.findsexobj(outfl, 3.0, pipevar, pix=pixscl, aperture=cpsfdi, 
                           wtimage=outwt, quiet=quiet)
            
            head = pf.getheader(outfl)
            
            coaddvars = np.loadtxt(outfl+'.stars', unpack=True)
            xim  = coaddvars[1,:]
            yim  = coaddvars[2,:]
            mag  = coaddvars[3,:]
            mage = coaddvars[4,:]
            flag = coaddvars[5,:]
            elon = coaddvars[8,:]
                
            # astropy does not like SWarp PV keywords or unicode, temporarily delete
            for key in head.keys():
                try:
                    if any(x in key for x in ['PV1_', 'PV2_', 'COMMENT', 'HISTORY']):
                        del head[key]
                except KeyError as error:
                    print error
                        
            w = wcs.WCS(head)
            wrd = w.all_pix2world(np.transpose([xim, yim]), 0)
            imfile  = outfl + '.im'
            catfile = outfl + '.cat'

            # Save stars from image
            np.savetxt(imfile, np.transpose([wrd[:,0],wrd[:,1],mag]))

            # Filter name correction:
            if filter == 'Z' or filter == 'Y': filter = filter.lower()

            nocustomcat = False
            # If custom catalog provided, match the same objects as the *.im file
            # and create refmag (reference magnitudes) that have the same matched
            # indices                    
            if customcat != None and filter in customcatfilt:

                print 'USING CUSTOM CATALOG'
                in_data = np.loadtxt(imfile)
                input_coords = in_data[:, :2]
                cat_data = np.loadtxt(customcat, skiprows=1)
                cat_coords = cat_data[:,:2]
                    
                cat_matches, tmp = apd.identify_matches(input_coords, cat_coords)
                    
                refmag = np.zeros(len(mag)) + 99
                mode = np.zeros(len(mag)) + -1
                for i, i_ind in enumerate(cat_matches):
                    if i_ind > 0:
                        #print input_coords[i], cat_coords[i_ind]
                        refmag[i] = cat_data[i_ind,catdict[filter]]
                        mode[i] = 4
                    
                # If no matching indices, run with the standard catalogs
                if sum(refmag<90.0) == 0: nocustomcat = True
            else:
                nocustomcat = True
                    
            # If custom catalog not provided, catalog doesn't include filter, or 
            # no objects from catalog found in image then
            # use get_SEDs.py to make catalog using 2MASS + (SDSS or APASS or USNOB1)
            if nocustomcat:
                # Create catalog star file 
                # (python get_SEDs.py imfile filter catfile USNOB_THRESH alloptstars)
                sedcmd = 'python ' + pipevar['getsedcommand'] + ' ' + imfile + ' ' +\
                    filter + ' ' + catfile + " 15 True "+ qtcmd
                
                if pipevar['verbose'] > 0: print sedcmd
                os.system(sedcmd)
                
                if not os.path.isfile(catfile):
                    zpts += [float('NaN')]
                    continue
                
                # Read in catalog file
                cvars = np.loadtxt(catfile, unpack=True)
                refmag = cvars[catdict[filter],:]
                mode   = cvars[catdict['mode'],:]

            # Find catalog filter values and only cutoff values of actual detections
            goodind = (mode != -1) & (refmag < 90.0) & (flag < 8) & (elon <= 1.3)
            
            refmag = refmag[goodind]
            obsmag = mag[goodind]
            obserr = mage[goodind]
            obswts = np.zeros(len(obserr))
            obskpm = np.zeros(len(obsmag))

            # Store magnitudes and weights (with minimum magnitude error of 0.01)
            for i in np.arange(len(obsmag)):
                if obserr[i] < 0.1:
                    obskpm[i] = obsmag[i]
                    obswts[i] = 1.0/(max(obserr[i], 0.01)**2)
            
            czpts, cscats, crmss = apd.calc_zpt(np.array([refmag]), np.array([obskpm]), 
                                    np.array([obswts]), sigma=1.0,
                                    plotter=pipevar['imworkingdir']+'zpt_'+filter+'.ps')
            
            chead = pf.getheader(outfl)
            
            # Add zeropoint keywords to header
            chead['SPIX']     = (cpsfdi, 'Final aperture size')
            chead['ABSZPT']   = (czpts[0]+25.0, 'Absolute zeropoint from calc_zpt')
            chead['ABSZPTSC'] = (cscats[0], 'Robust scatter of absolute zeropoint')
            chead['ABSZPRMS'] = (crmss[0], 'RMS of absolute zeropoint')
            
            # Add summary of stack information to header
            chead['DATE1']     = (firsttime, 'First frame time')
            chead['DATEN']     = (lasttime, 'Last frame time')
            chead['DATE']     = (medtime, 'Median frame time')
            chead['NSTACK']   = nstack
            chead['AIRMASS']  = (medair, 'Median exposure airmass')
            chead['AIRMIN']   = (minair, 'Minimum exposure airmass')
            chead['AIRMAX']   = (maxair, 'Maximum exposure airmass')
            chead['EXPTIME']  = (medexp, 'Effective rescaled exposure time')
            chead['TOTALEXP'] = (totexp, 'Total summed integration time')
            chead['MAXEXP']   = (max(stackexps), 'Length of longest exposure')
            chead['MINEXP']   = (min(stackexps), 'Length of shortest exposure')
            
            for i, file in enumerate(newstack):
                chead['STACK'+str(i)] = file    
      		
      		cdata = pf.getdata(outfl)
      		pf.update(outfl, cdata, chead)
      		
      	    if len(removedframes) > 0:
      	        print 'Removed frames with bad zeropoint fits: ' 
      	        print removedframes

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits, zsfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')  		                            
        os.system('rm -f ' + pipevar['imworkingdir'] + 'a*fp' + pipevar['prefix'] + '*.im')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'a*fp' + pipevar['prefix'] + '*.stars')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'a*fp' + pipevar['prefix'] + '*.cat')
        