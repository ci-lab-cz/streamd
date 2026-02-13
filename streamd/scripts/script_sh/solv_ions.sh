#!/bin/bash
# wdir
#Solvate

if [[ -z "$box_type" ]]; then
    echo "Missing required variable: box_type"
    exit 1
fi

if [[ -z "$box_padding_nm" ]]; then
    echo "Missing required variable: box_padding_nm"
    exit 1
fi

>&2 echo 'Script running:***************************** Solvation step *********************************'
cd $wdir
gmx editconf -f complex.gro -o newbox.gro -c -d "$box_padding_nm" -bt "$box_type" || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
#gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.2 #Warning about bad box - wrong number of atoms (https://gromacs.org-gmx-users.maillist.sys.kth.narkive.com/q4NXMAoY/bad-box-error) Lena version
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

#Add ions
>&2 echo 'Script running:***************************** Ions *********************************'
#use grompp to assemble *tpr file () using any *mdp file
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr   || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
#-maxwarn 10
#FATAL error SDMSO type not found -> try renaming to SDmso

printf '%s\n' "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
