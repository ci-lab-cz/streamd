#!/bin/bash
#  args: wdirz
OMP_NUM_THREADS=2
cd $wdir
#Energy minimization
if [ ! -f em.gro ]; then
>&2 echo 'Script running:***************************** Energy minimization *********************************'
gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -n index.ndx -o em.tpr -maxwarn 1
gmx mdrun -nt 32 -v -deffnm em -s em.tpr  || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

gmx energy -f em.edr -o potential.xvg <<< "Potential"
fi

# NVT
if [ ! -f nvt.gro ]; then
>&2 echo 'Script running:***************************** NVT *********************************'
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt -s nvt.tpr  || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

gmx energy -f nvt.edr -o temperature.xvg  <<< "Temperature"
fi

# NPT
if [ ! -f npt.gro ]; then
>&2 echo 'Script running:***************************** NPT *********************************'
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr  -maxwarn 1
gmx mdrun -deffnm npt -s npt.tpr  || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

gmx energy -f npt.edr -o pressure.xvg <<< "Pressure"
gmx energy -f npt.edr -o density.xvg <<< "Density"
fi