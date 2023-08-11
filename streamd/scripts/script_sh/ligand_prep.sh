#!/bin/bash
#  args: input_dirname lfile script_path resid charge molid
cd $input_dirname
parmchk2 -i $molid.mol2 -f mol2 -o $molid.frcmod || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
tleap -f tleap.in  || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
python $script_path"/"pmed_amb2gmx.py -p $molid.prmtop -x $molid.inpcrd -o $molid || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }

cp $molid.top $molid.itp
sed -i '/system/,+2 d' $molid.itp
sed -i '/molecules/,+2 d' $molid.itp
sed -i '/defaults/,+2 d' $molid.itp

# restraints
gmx genrestr -f $molid.gro -o posre_$molid.itp -fc 1000 1000 1000 <<< 2 || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }