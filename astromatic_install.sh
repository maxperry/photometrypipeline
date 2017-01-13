#!/bin/bash

set -e
echo 'Installing AstrOmatic software'

# http://stackoverflow.com/a/27776822
case "$(uname -s)" in

	Darwin)
		sudo port selfupdate
		echo '- installing atlas (will take a while...)'
		sudo port install atlas +nofortran 
		echo '- installing sextractor'
		sudo port install sextractor 
		echo '- installing swarp'
		sudo port install swarp
		echo '- installing scamp'
		sudo port install scamp
		echo '- installing missfits'
		sudo port install missfits
		echo '- installing cdsclient'
		sudo port install cdsclient
		;;

	Linux)
		ASTRDIR="$(mktemp -d -t 'astromatic')"
		
		# sextractor
		echo '- installing sextractor'
		sudo apt-get install libatlas3gf-base libatlas-headers
		sudo apt-get install sextractor

		if ! command -v wget >/dev/null 2>&1; then
			sudo apt-get install wget
		fi

		# swarp
		echo '- installing swarp'
		cd $ASTRDIR && \
			wget http://www.astromatic.net/download/swarp/swarp-2.38.0.tar.gz && \
			tar -zxvf swarp-2.38.0.tar.gz && \
			cd swarp-2.38.0 && \
			./configure && \
			make && \
			sudo make install

		# dependencies for cdsclient, clapack and scamp
		echo '- installing dependencies for cdsclient, clapack and scamp'		
		sudo apt-get install libfftw3-dev liblapack-dev libplplot-dev libshp-dev

		# cdsclient
		echo '- installing cdsclient'
		cd $ASTRDIR && \
			wget http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz && \
			tar -zxvf cdsclient.tar.gz && \
			cd cdsclient-3.83/ && \
			./configure && \
			make && \
			sudo make install

		# clapack
		echo '- installing clapack'
		cd $ASTRDIR && \
			wget http://www.netlib.org/clapack/clapack.tgz && \
			tar -zxvf clapack.tgz && \
			cd CLAPACK-3.2.1 && \
			mv make.inc.example make.inc && \
			sed -i 's/blas$(PLAT).a/libcblaswr.a -lcblas -latlas/' make.inc && \
			make cblaswrap && \
			sudo make lapacklib && \
			mv lapack_LINUX.a liblapack.a && \
			sudo cp libcblaswr.a /usr/local/lib/ && \
			sudo cp liblapack.a /usr/local/lib/


		# scamp
		echo '- installing scamp'
		cd $ASTRDIR && \
			wget http://www.astromatic.net/download/scamp/scamp-2.0.4.tar.gz && \
			tar -zxvf scamp-2.0.4.tar.gz && \
			cd scamp-2.0.4 && \
			./configure --with-atlas-incdir=/usr/include/atlas/ && \	
			sudo make && \
			sudo make install

		# missfits
		echo '- installing missfits'
		cd $ASTRDIR && \
			wget http://www.astromatic.net/download/missfits/missfits-2.8.0.tar.gz && \
			tar -zxvf missfits-2.8.0.tar.gz && \
			cd missfits-2.8.0 && \
			./configure  && \
			sudo make && \
			sudo make install && \

		# clean temp
		rm -Rf $ASTRDIR
		;;


	CYGWIN*|MINGW32*|MSYS*)
		echo ERROR Windows is not supported.
		exit 1
		;;

   # Add here more strings to compare
   # See correspondence table at the bottom of this answer

	*)
		echo 'other OS'
		echo ERROR Only Linux and Mac OS X are supported.
		exit 1
		;;
esac

