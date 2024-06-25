#!/bin/bash
#  args: wdir
cd $wdir
# MD
echo 'Script running:***************************** MD simulation *********************************'
echo 'Run simulation:'
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o $deffnm.tpr -maxwarn 1 || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx mdrun -deffnm $deffnm -s $deffnm.tpr -nt $ncpu -nb $compute_device -update $compute_device -pme $compute_device -bonded $compute_device -pmefft $compute_device || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }