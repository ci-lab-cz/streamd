#!/bin/bash
#  args: wdir index_protein_ligand dtstep wdir_out_analysis
cd $wdir

echo 'Script running:***************************** Analysis of MD simulation *********************************'

printf '%s\n' "System" | gmx trjconv -s $tpr -f $xtc -pbc nojump -o $deffnm\_noj_noPBC.xtc || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
#printf '%s\n' "Protein  System" | gmx trjconv -s $tpr -f $deffnm.xtc -o $deffnm\_noPBC.xtc -pbc mol -center
printf '%s\n%s\n' "$index_group" "System" | gmx trjconv -s $tpr -f $deffnm\_noj_noPBC.xtc -o md_centermolsnoPBC.xtc -pbc mol -center -n index.ndx || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
# use it for PBSA https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/issues/33
printf '%s\n%s\n' "$index_group" "System" | gmx trjconv -s $tpr -f md_centermolsnoPBC.xtc -fit rot+trans -o md_fit.xtc -n index.ndx || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

printf '%s\n%s\n' "$index_group" "non-Water" | gmx trjconv -s $tpr -f md_centermolsnoPBC.xtc -fit rot+trans -o md_fit_nowater.xtc -n index.ndx || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
printf '%s\n' "non-Water" | gmx convert-tpr -s $tpr -o  md_out_nowater.tpr

printf '%s\n' "System" | gmx trjconv -s $tpr -f md_fit.xtc -dt $dtstep -o md_short_forcheck.xtc || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }

printf '%s\n' "Protein" | gmx gyrate -s $tpr -f md_fit.xtc -n index.ndx -o $wdir_out_analysis/gyrate_$system_name.xvg || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}"; }
printf '%s\n' "Protein" | gmx rmsf -s $tpr -f md_fit.xtc -n index.ndx -o $wdir_out_analysis/rmsf_$system_name.xvg -oq $wdir_out_analysis/rmsf_$system_name.pdb -res || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}"; }

printf '%s\n' "System" | gmx trjconv -s $tpr -f md_fit.xtc -o frame.pdb -b 10 -e 11  -n index.ndx || { echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}"; }

rm md_centermolsnoPBC.xtc
rm $deffnm\_noj_noPBC.xtc
