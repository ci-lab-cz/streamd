#!/bin/bash
#  args: input_dirname script_path molid ligand_forcefield
cd "$input_dirname" || { echo "Cannot enter ligand working directory: $input_dirname" >&2; exit 1; }
: "${ligand_forcefield:=gaff}"
#case "$ligand_forcefield" in
#    gaff|gaff2)
#        ;;
#    *)
#        echo "Unsupported ligand force field: $ligand_forcefield" >&2
#        exit 1
#        ;;
#esac
parmchk2 -i "$molid.mol2" -f mol2 -o "$molid.frcmod" -s "$ligand_forcefield" || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
tleap -f tleap.in  || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
python "$script_path/pmed_amb2gmx.py" -p "$molid.prmtop" -x "$molid.inpcrd" -o "$molid" || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }

cp "$molid.top" "$molid.itp"
sed -i '/system/,+2 d' "$molid.itp"
sed -i '/molecules/,+2 d' "$molid.itp"
sed -i '/defaults/,+2 d' "$molid.itp"

# restraints
gmx make_ndx -f "$molid.gro" -o index.ndx << INPUT
2 & ! a H*
q
INPUT
printf '%s\n' "3" | gmx genrestr -f "$molid.gro" -o "posre_$molid.itp" -n index.ndx -fc 1000 1000 1000 || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
