#!/bin/bash
#input_dirname lfile script_path name
cd $input_dirname
LNAME="$(awk -F'/' '{print $NF}' <<<$lfile | cut -d'.' -f1)"
if [ $(awk -F'/' '{print $NF}' <<<$lfile |  awk -F. '{print $NF}') != gro ]
 then
  #prepare tleap - add ff from relevant conda env path
  cp $script_path"/tleap.in" .
  sed -i "1s/.*/source "$(echo $CONDA_PREFIX | sed 's/\//\\\//g')"\/dat\/leap\/cmd\/leaprc\.gaff/" tleap.in
  sed -i "s/ligand/$name/g" tleap.in
  if [ $(awk -F'/' '{print $NF}' <<<$lfile |  awk -F. '{print $NF}') != mol2 ]
  then
   charge=$( python $script_path"/"getcharge.py -i $lfile)
   antechamber -i $lfile -fi mdl -o $name.mol2 -fo mol2 -c bcc -pf y -s 2 -nc $charge -rn $name
  else
   cp $lfile $name.mol2
   >&2 echo 'Script running: Using prepared ligand *.mol2 file'
  fi
  parmchk2 -i $name.mol2 -f mol2 -o $name.frcmod
  tleap -f tleap.in
  python $script_path"/"pmed_amb2gmx.py -p $name.prmtop -x $name.inpcrd -o $name
 else
#put into dir
  cp $input_dirname\/$LNAME.gro $name.gro
  cp $input_dirname\/$LNAME.top $name.top
  >&2 echo 'Script running: Using prepared ligand *.gro file'
fi

cp $name.top $name.itp
sed -i '/system/,+2 d' $name.itp
sed -i '/molecules/,+2 d' $name.itp
sed -i '/defaults/,+2 d' $name.itp

# restraints
gmx genrestr -f $name.gro -o posre_$name.itp -fc 1000 1000 1000 <<< 2

