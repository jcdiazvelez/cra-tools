#!/bin/sh

# from bash or tcsh, call this script as:
# eval `/net/local/software/setup.sh`

. /net/local/software/os_arch.sh

version=1
while getopts 'v:' OPTION
do
	case $OPTION in
		v)
			version=$OPTARG
			;;
		*)
			echo "Uknown option '$OPTION'!"
			exit 1
			;;
	esac
done

if [ -z "$SROOT" ];
	then echo "NOTICE: /net/local/software is now deprecated. Please change to using the /cvmfs/icecube.wisc.edu/setup.sh script instead of /net/local/software/setup.sh. This script invocation is likely located in your .bash_profile." 1>&2
fi

SROOT=/net/local/software/$OS_ARCH
I3_PORTS=$SROOT/i3ports/v3
I3_SITE_CMAKE_DIR=$SROOT/i3ports/site_cmake

# Switcheroo!
case $version in
	2)
		if [[ -d /net/local/software/v2/$OS_ARCH ]]; then
			SROOT=/net/local/software/v2/$OS_ARCH;
		else
			echo 'echo '/net/local/software/setup.sh: No v2 software yet for $OS_ARCH, reverting to v1.' ;';
		fi
		;;
	*)
		;;
esac

PATH=$SROOT/bin:$I3_PORTS/bin:$PATH

PKG_CONFIG_PATH=$SROOT/lib/pkgconfig:$PKG_CONFIG_PATH
LD_LIBRARY_PATH=$SROOT/lib:$I3_PORTS/lib:$LD_LIBRARY_PATH
PYTHONPATH=$SROOT/lib/python2.6/site-packages:$I3_PORTS/lib/python2.6/site-packages:$PYTHONPATH
MANPATH=$SROOT/man:$SROOT/share/man:$MANPATH

# Port version specific bits
LD_LIBRARY_PATH=$I3_PORTS/lib/Minuit2-5.24.00:$I3_PORTS/lib/boost-1.38.0:$I3_PORTS/lib/log4cplus-1.0.2:$LD_LIBRARY_PATH
if [ -d $I3_PORTS/qt-4.6.0 ]; then
	LD_LIBRARY_PATH=$I3_PORTS/qt-4.6.0/lib:$LD_LIBRARY_PATH
	PATH=$I3_PORTS/qt-4.6.0/bin:$PATH
elif [ -d $I3_PORTS/qt-4.4.3 ]; then
	LD_LIBRARY_PATH=$I3_PORTS/qt-4.4.3/lib:$LD_LIBRARY_PATH
	PATH=$I3_PORTS/qt-4.4.3/bin:$PATH
fi

# ROOT specific bits
if [ -d $I3_PORTS/root-v5.30.06 ]; then
	: ${ROOTVER="5.30.06"}
elif [ -d $I3_PORTS/root-v5.30.05 ]; then
	: ${ROOTVER="5.30.05"}
elif [ -d $I3_PORTS/root-v5.28.00 ]; then
	: ${ROOTVER="5.28.00"}
elif [ -d $I3_PORTS/root-v5.24.00b ]; then
	: ${ROOTVER="5.24.00b"}
fi
: ${ROOTSYS="$I3_PORTS/root-v$ROOTVER"}
PATH=$ROOTSYS/bin:$PATH
LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

# MPI, if installed
if [ -d /usr/lib64/openmpi/bin ]; then
	PATH=/usr/lib64/openmpi/bin:$PATH
fi

# GotoBLAS
GOTO_NUM_THREADS=1

# Java is the future. Enterprise. Continuous Improvement. TPS Reports. Profit.
case $OS_ARCH in
	RHEL_6.0_amd64)
		JAVA_HOME=/usr/lib/jvm/java       ;;
	RHEL_5.0_ia32)
		JAVA_HOME=/usr/java/jdk1.5.0_12   ;;
	RHEL_5.0_amd64)
		JAVA_HOME=/usr/java/default       ;;
	RHEL_4.0_ia32)
		JAVA_HOME=/usr/java/j2sdk1.4.2_14 ;;
	RHEL_4.0_amd64)
		JAVA_HOME=/usr/java/j2sdk1.4.2    ;;
esac


for name in SROOT I3_PORTS I3_SITE_CMAKE_DIR PATH MANPATH PKG_CONFIG_PATH LD_LIBRARY_PATH PYTHONPATH ROOTSYS OS_ARCH GCC_VERSION JAVA_HOME GOTO_NUM_THREADS
do
  eval VALUE=\$$name
  case ${SHELL##*/} in 
	bash)
		echo 'export '$name=$VALUE' ;' ;;
	zsh)
		echo 'export '$name=$VALUE' ;' ;;
	tcsh)
		echo 'setenv '$name' '$VALUE' ;' ;;
	csh)
		echo 'setenv '$name' '$VALUE' ;' ;;
	*)
		echo 'echo I do not know how to deal with shell ${SHELL} ;' ;;
  esac
done

