#!/bin/bash
#  args: wdir
OMP_NUM_THREADS=2
cd $wdir
# MD
>&2 echo 'Script running:***************************** MD simulation *********************************'
>&2 echo 'Run simulation:'
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_out.tpr -maxwarn 1 || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx mdrun -deffnm md_out -s md_out.tpr || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }