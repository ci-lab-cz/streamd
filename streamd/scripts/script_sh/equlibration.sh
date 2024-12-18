#!/bin/bash
#  args: wdir
cd $wdir
unset OMP_NUM_THREADS
#Energy minimization
if [ ! -f em.gro ]; then
>&2 echo 'Script running:***************************** Energy minimization *********************************'
gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 2
gmx mdrun -v -deffnm em -s em.tpr -nt $ncpu -nb $compute_device $gpu_args || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

gmx energy -f em.edr -o $wdir_out_analysis/potential_$system_name.xvg <<< "Potential"
fi

# NVT
if [ ! -f nvt.gro ]; then
>&2 echo 'Script running:***************************** NVT *********************************'
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -s nvt.tpr -nt $ncpu -nb $compute_device $device_param $gpu_args || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

gmx energy -f nvt.edr -o $wdir_out_analysis/temperature_$system_name.xvg  <<< "Temperature"
fi

# NPT
if [ ! -f npt.gro ]; then
>&2 echo 'Script running:***************************** NPT *********************************'
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr  -maxwarn 1
gmx mdrun -deffnm npt -s npt.tpr -nt $ncpu -nb $compute_device $device_param $gpu_args || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

gmx energy -f npt.edr -o $wdir_out_analysis/pressure_$system_name.xvg <<< "Pressure"
gmx energy -f npt.edr -o $wdir_out_analysis/density_$system_name.xvg <<< "Density"
fi