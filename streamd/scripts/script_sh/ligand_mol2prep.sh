#!/bin/bash
#  args: input_dirname lfile script_path resid charge molid >> run_md_sh.log 2>&1
cd $input_dirname
antechamber -i $lfile -fi mdl -o $molid.mol2 -fo mol2 -c bcc -pf y -s 2 -nc $charge -rn $resid || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }