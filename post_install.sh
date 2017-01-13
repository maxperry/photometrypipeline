#!/bin/bash

set -e
echo 'Setting up environment'

##########################
# check for dependencies #
##########################

# SExtractor
command -v sex >/dev/null 2>&1 || { echo >&2 "SExtractor is required but it's not installed.  Aborting."; exit 1; }

# Swarp
command -v swarp >/dev/null 2>&1 || { echo >&2 "Swarp is required but it's not installed.  Aborting."; exit 1; }

# Scamp
command -v scamp >/dev/null 2>&1 || { echo >&2 "Scamp is required but it's not installed.  Aborting."; exit 1; }

# MissFITS
command -v missfits >/dev/null 2>&1 || { echo >&2 "MissFITS is required but it's not installed.  Aborting."; exit 1; }

# cdsclient
command -v findcat >/dev/null 2>&1 || { echo >&2 "cdsclient is required but it's not installed.  Aborting."; exit 1; }

##########################
# set up  environment #
##########################

PACKAGE_NAME="photopipe"

# set project root directory
# http://stackoverflow.com/a/246128
# export PIPE_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # allows startup.sh to be called from any directory.  will fail if last component of path to startup.sh is a symlink.
PIPE_ROOT=$(python -c "import site; print(site.getsitepackages()[0])")"/"$PACKAGE_NAME

AUTOPROC_CONFIG_PATH=$PIPE_ROOT"/reduction/auto/pipeautoproc.par"
SWARP_PATH="$(command -v swarp)"
SEX_PATH="$(command -v sex)"
SCAMP_PATH="$(command -v scamp)"

touch $AUTOPROC_CONFIG_PATH

# http://unix.stackexchange.com/a/77278
cat <<EOT >> $AUTOPROC_CONFIG_PATH
autoastrocommand	:	$PIPE_ROOT/reduction/astrom/vlt_autoastrometry.py
getsedcommand		:	$PIPE_ROOT/photometry/dependencies/get_SEDs.py
swarpcommand		:	$SWARP_PATH
sexcommand			:	$SEX_PATH
scampcommand		:	$SCAMP_PATH
prefix				:	2
datadir				:	
EOT

# add pipeline directories to python path
PYTHONPATH=$PIPE_ROOT/instruments:$PIPE_ROOT/photometry:$PIPE_ROOT/photometry/dependencies:$PIPE_ROOT/reduction:$PIPE_ROOT/reduction/dependencies:$PIPE_ROOT/reduction/auto:$PIPE_ROOT/reduction/astrom:$PYTHONPATH

# save permanently
echo '# photopipe environment' >> /etc/profile
echo 'export PIPE_ROOT='$PIPE_ROOT >> /etc/profile
echo 'export PYTHONPATH='$PYTHONPATH >> /etc/profile
