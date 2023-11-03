#!/bin/bash
# ligand prep
echo 'Script running:***************************** 1. Start Gaussian preparation *********************************'
cd $input_dirname
$activate_gaussian
echo 'Script running:***************************** 2. Gaussian preparation *********************************'
python $script_path/prepare_Gaussian_input.py -i $lfile --opt_param 1.com --charges_param 2.com --freq_param 3.com --output $molid || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
echo 'Script running:***************************** 3. Conformer Gaussian optimization *********************************'
$gaussian_version < $molid\_opt_gaussian.com > $molid\_opt_gaussian.log || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
cp ligand.chk ligand_opt.chk
echo 'Script running:***************************** 4. Charges Gaussian calculation *********************************'
$gaussian_version < $molid\_charg_gaussian.com > $molid\_charg_gaussian.log || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
echo 'Script running:***************************** 5. RESP Charges calculation *********************************'
# j5 - ignore damaged bonds during optimization
antechamber -fi gout -fo mol2 -i $molid\_charg_gaussian.log -o _$molid\_tofix.mol2 -c resp -j 5 -nc $charge -pf y -s 2 -dr no -rn $resid || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
echo 'Script running:***************************** 6. Coords and Bonds fix *********************************'
# return original coords and fix damaged bonds
python $script_path/mol2_fix_coordsbonds.py --mol $lfile --mol2 _$molid\_tofix.mol2 -o _$molid\_fixed.mol2 || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }
# better to check $molid.mol2 by visual inspection
# correct mol2 from fixed previous mol2 and j 1 - use previous bonds from fixed mol2
antechamber -fi mol2 -fo mol2 -i _$molid\_fixed.mol2 -o $molid.mol2 -nc $charge -pf y -s 2 -dr no -rn $resid -j 1 || { echo "Failed to run command  at line ${LINENO} in ${BASH_SOURCE}" && exit 1; }

#if [ $mode == 'forcefield' ]
#then
#>&2 echo 'Script running:***************************** 7. Freq Gaussian calculation *********************************'
#$gaussian_version < $molid\_freq_gaussian.com > $molid\_freq_gaussian.log
#formchk $molid\_freq_gaussian.chk > $molid\_freq.fchk
#>&2 echo 'Script running:***************************** 8. CartHess2FC.py *********************************'
## get optimized structure
#newzmat -ichk -opdb ligand_opt.chk opt.pdb
## check if there are not zero values
#CartHess2FC.py -p opt.pdb -f $molid\_freq.fchk  > $molid\.forcefield