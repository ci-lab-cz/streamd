#!/bin/bash
#  args: wdir
OMP_NUM_THREADS=2
cd $wdir
# MD
>&2 echo 'Script running:***************************** MD simulation *********************************'
>&2 echo 'Run simulation:'
gmx mdrun -s $tpr -v -deffnm $deffnm_next -cpi $cpt -noappend -nsteps $new_mdsteps || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx trjcat -f $xtc $deffnm_next\.part*.xtc -o $deffnm_next\.xtc -settime << INPUT
0
c
INPUT