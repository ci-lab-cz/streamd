#!/bin/bash
#  args: wdir molid index_protein_ligand  index_ligand_noH tu
cd $wdir
>&2 echo 'Script running:***************************** Ligand Analysis of MD simulation *********************************'

gmx rms -s $tpr -f md_fit.xtc -o rmsd_$molid\.xvg -n index.ndx  -tu $tu <<< "Backbone  $index_ligand_noH" || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}"; }
gmx rms -s em.tpr -f md_fit.xtc -o rmsd_$molid\_xtal.xvg -n index.ndx -tu $tu <<< "Backbone  $index_ligand_noH" || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}"; }


