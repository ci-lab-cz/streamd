#!/bin/bash
#  args: wdir, deffnm, device_param, gpu_args
cd $wdir
unset OMP_NUM_THREADS
# MD
echo 'Script running:***************************** MD simulation *********************************'
echo 'Run simulation:'
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o $deffnm.tpr -maxwarn 1 || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx mdrun -deffnm $deffnm -s $deffnm.tpr -nt $ncpu $device_param $gpu_args || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }