#!/bin/bash
#  args: wdir molid index_protein_ligand  index_ligand_noH tu
cd $wdir

gmx rms -s $tpr -f md_fit.xtc -o rmsd_$molid\.xvg -n index.ndx  -tu $tu <<< "Backbone  $index_ligand_noH" || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }
gmx rms -s em.tpr -f md_fit.xtc -o rmsd_$molid\_xtal.xvg -n index.ndx -tu $tu <<< "Backbone  $index_ligand_noH" || { >&2 echo "Failed to run command  at line ${LINENO} of ${BASH_SOURCE}" && exit 1; }


