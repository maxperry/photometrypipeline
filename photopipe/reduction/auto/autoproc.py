import autoproc_steps as ap
import os

def autoproc(datadir=None, imdir=None, start=None, stop=None, only=None, step=None,
    nocrclean=False, nomastersky=False, redo=False, quiet=False, rmifiles=False,
    customcat=None, customcatfilt=[]):

    """
    NAME:
       autoproc
    PURPOSE:
        Fully-automated reduction of imaging.
    CALLING SEQUENCE:
        ratautoproc, [settings]
    INPUTS (all optional):
        datadir       - Location of raw data (current directory if unspecified)
        imdir	      - Location of processed data ('imredux' if unspecified)
        start	      - Start with this step, skipping previous ones
        stop          - End with this step, skipping subsequent ones
        only          - Do only this step
        step          -  (completely identical to only, takes precedence)
        redo          - Repeat step(s), overwriting any existing files
        nocrclean     - Do not zap cosmic rays
        nomastersky   - Do not create master sky, only subtract median of sky
        quiet	      - (mainly) silent output unless errors
        rmifiles      - Removes intermediate files
        customcat     - Custom catalog (txt file) to determine instrumental zeropoint 
                        corrections
                        Must be in same format as what get_SEDs.py produces 
                        (i.e. ra(deg)	dec(deg)	u	g	r	i	z	y	B	V	R	
                        I	J	H	K	u_err	g_err	r_err	i_err	y_err	B_err	
                        V_err	R_err	I_err	J_err	H_err	K_err	Mode)
                        First line will be skipped (use for headers)
                        everything but JHK are expected to be in AB magnitudes, JHK should
                        be in same units as 2MASS (Vega mags)
        customcatfilt - Filters relevant to custom catalog file (all other filters will
                        use get_SEDs.py to calculate catalog from 2MASS + (SDSS or APASS
                        or USNOB1) in that order
                      
    ADDITIONAL OPTIONS:
        If any of the following files are found in the directory where autoproc
        is run, it will change the default behavior.
       
        pipeautoproc.par - Contains various defaults (mainly directory paths)
    COMMENTS:
        This code is meant to be fully automated, producing science-quality
        output with a single command line with zero intervention.  Reduces
        RATIR and LMI data.
        
        The production of images is quite high-level and includes photometric 
        calibration of the field (although the accuracy of this has not been
        robustly tested).
        
        The program works in a series of steps following standard CCD reduction
        techniques, automatically recognizing calibration files, matching
        calibrations together and science frames with appropriate calibrations,
        etc.  Users looking for more control over the reductions can run each
        step individually with the "step" command.

        If the reduction is interrupted, the program will resume without
        incident (although any failed steps may be repeated); each task is 
        run independently of all others.
   
        The program tries very hard to correct for observer mistakes. But it's not perfect.
        If problems occur, generally the easiest fix is to delete any 
        offending files and rerun. More significant problems generally require direct 
        modification of the code. 

        Filenames for the input raw images are expected to be in 
        2*.fits format. They can be either in the working directory or 
        in a different directory specified by datadir.

        The code runs in a series of steps in the following order:

        1. "Prepare" the data by converting the multi-header extension FITS files
        to a standard single frame, and adding extra information to the header. 
        Output: p*.fits (written to ./imredux/ by default.)

        2. "Flatten" data.  Divide each image by the flatfield.  A more
        refined cropping is also done at this stage, depending on the placement
        of the image area during the run (which is variable.) Assumes master flat
        (ex. flat_H.fits in imredux/ folder)
        
        3. "Makesky" makes a master sky (ex. sky-H.fits in imredux/ folder)
        
        4. "Skysubtract" subtracts out master sky 
        Output: sp*.fits (written to ./imredux/ by default.)
        
        5. Removes cosmic rays, using the independent routines cosmics.py that
        uses Laplacian cosmic ray indentification. See that program for more information.  
        This can be a time-consuming process.
        Output: zfp*.fits

        6. Solve astrometry of the field against the best available online catalog
        (SDSS/2MASS/APASS/USNO-B1.0), using the independent vlt_autoastrometry.py code.
        (Also requires sextractor to be installed.) Uses two passes of Scamp for a 
        secondary correction.  Scamp accounts for distortion.  Uses higher distortion
        parameters if already supplied distortion keywords.
        Output: azfp*.fits

        7. Stack exposures.  A weighted, masked median is performed for each field.
        Requires swarp to be installed and properly linked.
        Output: coadd[object].[filter].fits
 
        EXAMPLES:
            1.  In a directory containing a full night worth of RATIR data, enter
                    autoproc()
                This will automatically execute all of the steps above.

            2.  Run only the "prepare" step, on data stored in a separate directory:
                    autoproc(step='prepare', datadir='raw/')

        TROUBLESHOOTING:
            If the pipeline crashes during operation:
            
                * Did an external call (autoastrometry, swarp, sex) fail?
                    Check that you have specified paths correctly and have the required 
                    software installed for astrometry/coadding
                    
                * Did it encounter some other error while processing a file?  
                    Check which file the program was working on when it crashed. If the
                    file is not essential, try deleting it and re-running the pipeline
                    starting with the current step (or delete ALL of the file's precursors
                    with the same file number and rerun the pipeline.) 
            
            If processing completed, but the results are problematic:

                * Did it return without actually processing any data? 
                    Make sure that it is in the current working directory or that you 
                    have correctly pointed to the directory containing raw data with the 
                    "datadir" keyword. If you are re-doing a step, it will not overwrite 
                    existing files by default. Set the redo keyword to overwrite files 
                    from a previously-attempted step (be sure to set "start" or "step" 
                    unless you want to restart the pipeline from the beginning.)

                * Were some files skipped?
                    If the pipeline encountered a non-fatal problem processing an
                    individual image (such as an inability to flatfield) then it will 
                    not process that file any further. For most cases a summary of the 
                    problems will be printed out at the end of processing. If a file is 
                    not being processed and you do not see it in the final summary,
                    you can simply rerun the pipeline (without deleting any files and 
                    without setting the redo flag) and it will try to repeat any failed 
                    steps of this nature. 
    
    - Pipeline adapted from Dan Perley LRIS pipeline 
      (http://www.dark-cosmology.dk/~dperley/code/code.html)
      
    - Modifications made by Vicki Toy and John Capone
    """
    
    # Load default parameters and interpret user arguments.
    pipevar = {'autoastrocommand':'autoastrometry', 'getsedcommand':'get_SEDs', 
            'sexcommand':'sex' , 'swarpcommand':'swarp' , 'rmifiles':0, 
            'prefix':'', 'datadir':'' , 'imworkingdir':'' , 'overwrite':0 , 'verbose':1,
            'flatfail':'' , 'fullastrofail':'' ,
            'pipeautopath':'' , 'refdatapath':'', 'defaultspath':'' }

    if imdir    != None: pipevar['imworkingdir'] = imdir
    
    ap.autopipedefaults(pipevar=pipevar)
    
    if redo     == True: pipevar['overwrite'] = 1
    if quiet    == True: pipevar['verbose'] = 0
    if datadir  != None: pipevar['datadir'] = datadir
    if rmifiles  == True: pipevar['rmifiles'] = 1
    
    # Step options
    steps = ['prepare', 'flatten', 'makesky', 'skysub', 'crclean', 'astrometry', 'stack']
    
    # If start is specified, truncate steps to start at specified step.
    # If invalid step end program with error
    if start != None:
        try:
            w = steps.index(start)
            steps = steps[w:]
        except:
            print "Invalid starting step: ", start
            print "Must be one of: ", steps
            return
            
    # If stop is specified, truncate steps to end at specified step.
    # If invalid step end program with error
    if stop != None:
        try:
            w = steps.index(stop)
            steps = steps[:w+1]
        except:
            print "Invalid stopping step: ", stop
            print "Must be one of: ", steps
            return

    # If step specified set only to specified step
    if step != None: only = step
    
    # If only specified (including step), set steps to run as specified step.
    # If invalid step end program with error
    if only != None:
        try:
            w = steps.index(only)
            steps = steps[w]
            if start != None: print 'Note that start is also set.'
            if stop != None:  print 'Note that stop is also set.'
        except:
            print "Invalid step: ", only
            print "Must be one of: ", steps
            return    
        
    # Check if autoastrometry, sextractor, swarp are installed and functioning before running steps
    if 'astrometry'in steps or 'stack' in steps:
        
        # Check sextractor
        if os.path.isfile('temp.txt'): os.system('rm -f temp.txt')
        os.system(pipevar['sexcommand'] + ' -d > temp.txt')
        
        if os.stat('temp.txt').st_size == 0:
            print "Error: Sextractor is not installed or not configured."
            print "       Cannot run image alignment steps. Configure or stop='crclean'"
        
        # Check autoastrometry
        if os.path.isfile('temp.txt'): os.system('rm -f temp.txt')
        os.system('python ' + pipevar['autoastrocommand'] + ' > temp.txt')

        if os.stat('temp.txt').st_size == 0:
            print "Error: Autoastrometry is not installed or not configured."
            print "       Cannot run image alignment steps. Configure or stop='crclean'" 
    
    if  'stack' in steps:

        # Check swarp
        if os.path.isfile('temp.txt'): os.system('rm -f temp.txt')
        os.system(pipevar['swarpcommand'] + ' -d > temp.txt')
        
        if os.stat('temp.txt').st_size == 0:
            print "Error: Swarp is not installed or not configured."
            print "       Cannot run image coadds. Configure or stop='astrometry'"
    
    if isinstance(steps, str): steps = [steps]       
    # Runs each processing step specified in the correct order (crclean is optional)      
    for step in steps:
        
        if step == 'prepare': ap.autopipeprepare(pipevar=pipevar)
        if step == 'flatten': ap.autopipeimflatten(pipevar=pipevar)
        if step == 'makesky' and nomastersky == False: ap.autopipemakesky(pipevar=pipevar)
        if step == 'skysub' and nomastersky == False:  ap.autopipeskysub(pipevar=pipevar)
        if step == 'skysub' and nomastersky == True:  ap.autopipeskysubmed(pipevar=pipevar)    
        if step == 'crclean' and nocrclean == False: ap.autopipecrcleanim(pipevar=pipevar)
        if step == 'astrometry': ap.autopipeastrometry(pipevar=pipevar),
        if step == 'stack'     : ap.autopipestack(pipevar=pipevar, customcat=customcat, customcatfilt=customcatfilt)

    # Prints the files that were not flat fielded due to problems with file
    if pipevar['flatfail'] != '':
        print 'Unable to flat-field the following images:'
        ffail = pipevar['flatfail'].split()
        for f in ffail: print f
    
    # Prints the files that were not astrometry corrected due to problems with the file 
    if pipevar['fullastrofail'] != '':
        print 'All astrometry failed for the following images (not stacked):'
        afail = pipevar['fullastrofail'].split()
        for f in afail: print f
    
    print 'Processing complete.'
    
    # Remove any files that were created during the reduction process
    if os.path.isfile('temp*.*'): os.system('rm -f temp*.*')
    if os.path.isfile('det.*'): os.system('rm -f det.*')
    if os.path.isfile('cat.*'): os.system('rm -f cat.*')    