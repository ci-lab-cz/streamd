#!/bin/bash
#  args: wdir
OMP_NUM_THREADS=2
cd $wdir
# MD
echo 'Script running:***************************** MD simulation *********************************'
echo 'Run simulation:'
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_out.tpr -maxwarn 1 || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx mdrun -deffnm md_out -s md_out.tpr || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }