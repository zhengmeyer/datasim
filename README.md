[![DOI](https://zenodo.org/badge/39930086.svg)](https://zenodo.org/badge/latestdoi/39930086)


# Requirements:
*DiFX-2.4 or above
*GSL
*GSLCBLAS
*automake autoconf libtool

# Description:

Limitations:
1. 2-bits sampling
2. Only generate single-thread multi-channel VDIF
3. Default SEFD for 3 antennas
4. VDIF packet of all antennas need to repsent the same time length

# Installation:

$ libtoolize --force
$ aclocal -I m4
$ autoheader
$ automake --force-missing --add-missing
$ autoconf
$ ./configure --prefix=$HOME
$ make
$ make install

# Usage:
1. source DiFX setup.bash file
2. source HOPS hops.bash file
3. 'datasim --help' shows a list of possible options 
4. generate v2d and VEX file: './genv2dvex.sh exp_dir input'
5. in the directory of the generated v2d and VEX file, 
   create simulated .vdif data:
   'datasim -f 100 -s 1000,2000,3000 -d 12345 datasim_7004.input'
   where -f is the source flux density, -s is a list of antenna SEFDs,
   -d is the random number seed
