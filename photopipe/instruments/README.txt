To add a new instrument:

- Only edit "specific_instruments.py"

- Create instrument with inherited class "instrument"
    - initialize instrument with instrument name and
      number of cameras
    - Define the functions that were defined in the
      instrument class as abstract methods:
    
        has_cam_bias
        has_cam_dark
        change_header_keywords
        slice
        is_cam_split
        get_cam_sat
        get_cam_gain
        get_filter
        possible_filters
        get_centered_filter
        original_file_format
        change_file_names
     
       Information of what each of these functions
       should do (input/outputs) are commented in
       the instrument class
    - In the change_header_keywords must have or create
      the following keywords for future pipeline use:
      
        CAMERA   (int)
        SIMPLE
        BITPIX
        NAXIS
        NAXIS1 
        NAXIS2
        EXPTIME  (in seconds)
        SATURATE (ADU/DN)
        GAIN     (e-/ADU)
        FILTER   (case-sensitive, currently only supports ugrizyBVRIJHK.
                  in stacking phase allows filters with "SDSS-*" to be changed
                  to lower case filter name *. To expand number of filters need
                  to include more SED fitting in code/photometry/dependencies/get_SEDs.py)
        PIXSCALE
        TARGNAME
        AIRMASS
        CD*_*    (WCS CD matrix)
        CRPIX1   (WCS info)
        CRPIX2   (WCS info)
        CTYPE1   (WCS info)
        CTYPE2   (WCS info)
        (PV*_* if there are distortion parameters)
        
      The following keywords are optional:
        
        INSTRUME
        LATITUDE
        LONGITUD
        BINNING
        BINX
        BINY
        WAVELENG
        TARGNAME
        PIXSCALE
        WAVELENG
        UTC
        OBJECT

- Add instrument with instrument name and class name to
  instrument_dict