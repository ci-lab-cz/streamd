#!/bin/bash
#  args: input_dirname lfile script_path resid charge molid ligand_forcefield >> run_md_sh.log 2>&1
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
antechamber -i "$lfile" -fi mdl -o "$molid.mol2" -fo mol2 -at "$ligand_forcefield" -c bcc -pf y -s 2 -nc "$charge" -rn "$resid" -dr "$dr" || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
