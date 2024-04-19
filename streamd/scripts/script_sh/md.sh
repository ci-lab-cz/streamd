#!/bin/bash
#  args: wdir
cd $wdir
# MD
echo 'Script running:***************************** MD simulation *********************************'
echo 'Run simulation:'
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o $deffnm.tpr -maxwarn 1 || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx mdrun -deffnm $deffnm -s $deffnm.tpr || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }