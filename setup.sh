#!/bin/bash -e

case `uname` in
	Linux) ECHO="echo -e" ;;
	*) ECHO="echo" ;;
esac

usage() {
	$ECHO "setup.sh [options]"
	$ECHO
	$ECHO "Options:"
	$ECHO "-m            \tinstall MadGraph5_aMC@NLO v2.5.4"
	$ECHO "-p            \tinstall Pythia 8.226 (w/in MadGraph, includes LHAPDF6, HepMC, MG5aMC_PY8_interface)"
	$ECHO "-d            \tinstall Delphes 3.4.1"
	$ECHO "-n            \tnumber of cores for Madgraph (default=1)"
	$ECHO "-h            \tdisplay this message and exit"
	exit $1
}

CURDIR=`pwd`
INSTALL_MADGRAPH=""
INSTALL_PYTHIA=""
INSTALL_DELPHES=""
NUMCORES_MADGRAPH=1
# check arguments
while getopts "mpdhn:" opt; do
	case "$opt" in
	m) INSTALL_MADGRAPH=yes
	;;
	p) INSTALL_PYTHIA=yes
	;;
	d) INSTALL_DELPHES=yes
	;;
	n) NUMCORES_MADGRAPH=$OPTARG
	;;
	h) usage 0
	;;
	esac
done

if [ -z "$INSTALL_MADGRAPH" ] && [ -z "$INSTALL_DELPHES" ]; then
	usage 1
fi

if [ -n "$INSTALL_PYTHIA" ] && [ -z "$INSTALL_MADGRAPH" ]; then
	$ECHO "Pythia installation requested without MadGraph - not supported"
	exit 1
fi

# setup compiler and environment
LCGDIR=/cvmfs/sft.cern.ch/lcg/views/LCG_87/x86_64-slc6-gcc49-opt
source $LCGDIR/setup.sh

if [ -n "$INSTALL_DELPHES" ]; then
	cd $CURDIR
	# get delphes source and compile
	wget -q http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.1.tar.gz
	tar -zxf Delphes-3.4.1.tar.gz
	cd Delphes-3.4.1
	make

	# cleanup
	cd $CURDIR
	rm Delphes-3.4.1.tar.gz
fi

if [ -n "$INSTALL_MADGRAPH" ]; then
	cd $CURDIR
	$ECHO "Download MG5_aMC_v2.5.4.tar.gz"
	wget -q http://madgraph.physics.illinois.edu/Downloads/MG5_aMC_v2.5.4.tar.gz
	if [ -n "$INSTALL_PYTHIA" ]; then
		$ECHO "Download pythia8226.tgz"
#		wget -q http://home.thep.lu.se/~torbjorn/pythia8/pythia8226.tgz
		xrdcp -s root://cmseos.fnal.gov//store/user/pedrok/pythia8226a.tgz .
#		$ECHO "Download hepmc2.06.09.tgz"
#		wget -q http://hepmc.web.cern.ch/hepmc/releases/hepmc2.06.09.tgz
	fi
	tar -xzf MG5_aMC_v2.5.4.tar.gz
	cd MG5_aMC_v2_5_4

	# configure
	echo "set auto_update 0" > mgconfigscript
	echo "set automatic_html_opening False" >> mgconfigscript
	echo "set output_dependencies internal" >> mgconfigscript
	if [ ${NUMCORES_MADGRAPH} -eq 1 ]; then
		echo "set run_mode 0" >> mgconfigscript
	else
		echo "set run_mode 2" >> mgconfigscript
		echo "set nb_core" ${NUMCORES_MADGRAPH} >> mgconfigscript
	fi

	# install Pythia8 w/in MadGraph
	if [ -n "$INSTALL_PYTHIA" ]; then
		BOOSTDIR=$LCGDIR/include/boost-1_62
		ZLIBDIR=$LCGDIR/lib
#		echo "install hepmc --hepmc_tarball=$CURDIR/hepmc2.06.09.tgz" >> mgconfigscript
		echo "install lhapdf6 --with_boost=$BOOSTDIR --with_zlib=$ZLIBDIR" >> mgconfigscript
		echo "install pythia8 --with_zlib=$ZLIBDIR --pythia8_tarball=$CURDIR/pythia8226a.tgz" >> mgconfigscript
	fi

	# for consistency
	if [ -n "$INSTALL_DELPHES" ]; then
		echo "set delphes_path $CURDIR/Delphes-3.4.1" >> mgconfigscript
	fi

	# run configuration
	echo "save options" >> mgconfigscript
	./bin/mg5_aMC mgconfigscript

	# cleanup
	cd $CURDIR
	rm MG5_aMC_v2.5.4.tar.gz
#	rm hepmc2.06.09.tgz
	rm pythia8226a.tgz
fi

