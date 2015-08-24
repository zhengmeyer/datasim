#!/bin/bash
# usage ./postproc.sh 7002 M87
[ $# -eq 2 ] || { echo "Usage: ./postproc.sh expNum source"
echo "(e.g. ./postproc.sh 7002 M87)" ; exit 1 ; }

# remove directory .difx and expNum if exist
[ -d datasim_$1.difx ] && { rm -rf datasim_$1.difx; echo delete datasim_$1.difx directory ...; }
[ -d $1 ] && { rm -rf $1; echo delete directory $1; }

cp $HOME/data/rundifx/default.machines datasim_$1.machines
cp $HOME/data/rundifx/default.threads datasim_$1.threads
mpirun --machinefile datasim_$1.machines -np 5 mpifxcorr datasim_$1.input
#startdifx -n -f datasim_$1.input
[ -d datasim_$1.difx ] || { echo No .difx dir; exit 1; }
ff=`type -p fourfit`
[ -z "$ff" ] && {
    echo no fourfit to run
    exit 1
} || {
    echo Postprocessing to $1/075-0558 with
    echo $ff
    difx2mark4 -v -e $1 datasim_$1
    root=`echo $1/075-0558/$2.* | sed "s/.*$2.//"`
    [ -s "$1/075-0558/$2.$root" ] || { echo No Root; exit 1; }
    echo '' ; echo Root is $root ; echo ''
    [ -z "x-m2 -c" ] || {
    # beg hacking
    ff_opts_x="-m2 -c cf-$1"

    cat > cf-$1 <<..EOF
    optimize_closure true
    * other things...
..EOF
    echo \
    $ff $ff_opts_x $1/075-0558/$2.$root 
    $ff $ff_opts_x $1/075-0558/$2.$root
    # end hacking
  } 
#[ -z "`type psmerge `" -o -z "`type ps2pdf`" ] && {
#  echo no psmerge or ps2pdf
#  exit 2
#    } || {
#  echo Making a summary pdf file
  pushd $1/075-0558
  for d in ??.?.*.$root
        do fplot -d $d.ps $d ; ls -l $d.ps ; done
#  psmerge -o$2.ps *.$root.ps
#  ps2pdf $2.ps
#  popd
#  ls -l $1/075-0558/$2.pdf | cut -c24-
#    }
}
exit 0
# eof

