#!/bin/bash

set -e
echo 'Installing python dependencies'

# http://stackoverflow.com/a/27776822
case "$(uname -s)" in

   Darwin)
	 # upgrade pip
	 python -m pip install --upgrade pip

	 echo '- installing numpy scipy matplotlib'
	 pip install --user numpy scipy matplotlib
	 echo '- installing pyfits'
	 pip install --user pyfits
	 echo '- installing astropy'
	 pip install --user --no-deps astropy
     ;;

   Linux)
	 # upgrade pip
	 python -m pip install --upgrade pip
	 
	 # numpy, scipy
	 echo '- installing gfortran libatlas-base-dev'
	 sudo apt-get install gfortran libatlas-base-dev
	 echo '- installing numpy scipy'
	 pip install --user numpy scipy

	 # matplotlib
	 echo '- installing matplotlib'
	 sudo apt-get install libpng12-dev libfreetype6-dev libxft-dev
	 pip install --user matplotlib

	 # astropy
	 echo '- installing astropy'
	 pip install --user --no-deps astropy

	 # pyfits
	 echo '- installing pyfits'
	 sudo apt-get install python-tk
	 pip install --user pyfits
     ;;

   CYGWIN*|MINGW32*|MSYS*)
     echo 'ERROR Windows is not supported.'
     exit 1
     ;;

   # Add here more strings to compare
   # See correspondence table at the bottom of this answer

   *)
     echo 'other OS'
     echo 'ERROR Only Linux and Mac OS X are supported.'
     exit 1
     ;;
esac

echo '- installing jsonschema'
pip install --user jsonschema


