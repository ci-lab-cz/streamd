#!/bin/bash
#  args: wdir
cd $wdir
unset OMP_NUM_THREADS
# MD
>&2 echo 'Script running:***************************** Continue MD simulation *********************************'
>&2 echo 'Run simulation:'

gmx convert-tpr -s $tpr -until $new_mdtime_ps -o $deffnm_next\.tpr
gmx mdrun -s $deffnm_next\.tpr -v -deffnm $deffnm_next -cpi $cpt -noappend -nt $ncpu -nb $compute_device $device_param $gpu_args || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx trjcat -f $xtc $deffnm_next\.part*.xtc -o $deffnm_next\.xtc -settime -tu fs << INPUT
0
c
INPUT
rm $deffnm_next\.part*.*