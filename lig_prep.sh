#!/bin/bash
#  args: input_dirname lfile script_path name
cd $input_dirname
LNAME="$(awk -F'/' '{print $NF}' <<<$lfile | cut -d'.' -f1)"
if [ $(awk -F'/' '{print $NF}' <<<$lfile |  awk -F. '{print $NF}') != gro ]
 then
  #prepare tleap - add ff from relevant conda env path
  cp $script_path"/tleap.in" .
  sed -i "1s/.*/source "$(echo $CONDA_PREFIX | sed 's/\//\\\//g')"\/dat\/leap\/cmd\/leaprc\.gaff/" tleap.in
  sed -i "s/ligand/$LNAME/g" tleap.in
  if [ $(awk -F'/' '{print $NF}' <<<$lfile |  awk -F. '{print $NF}') != mol2 ]
  then
   charge=$( python $script_path"/"getcharge.py -i $lfile)
   antechamber -i $lfile -fi mdl -o $LNAME.mol2 -fo mol2 -c bcc -pf y -s 2 -nc $charge -rn $name || { >&2 echo "Failed to run command  at line ${LINENO} in bash script" && exit 1; }
  else
   cp $lfile $LNAME.mol2
   >&2 echo 'Script running: Using prepared ligand *.mol2 file'
  fi
  parmchk2 -i $LNAME.mol2 -f mol2 -o $LNAME.frcmod || { >&2 echo "Failed to run command  at line ${LINENO} in bash script" && exit 1; }
  tleap -f tleap.in || { >&2 echo "Failed to run command  at line ${LINENO} in bash script" && exit 1; }
  python $script_path"/"pmed_amb2gmx.py -p $LNAME.prmtop -x $LNAME.inpcrd -o $LNAME || { >&2 echo "Failed to run command  at line ${LINENO} in bash script" && exit 1; }
 else
#put into dir
  cp $input_dirname\/$LNAME.gro $LNAME.gro
  cp $input_dirname\/$LNAME.top $LNAME.top
  >&2 echo 'Script running: Using prepared ligand *.gro file'
fi

cp $LNAME.top $LNAME.itp
sed -i '/system/,+2 d' $LNAME.itp
sed -i '/molecules/,+2 d' $LNAME.itp
sed -i '/defaults/,+2 d' $LNAME.itp

# restraints
gmx genrestr -f $LNAME.gro -o posre_$LNAME.itp -fc 1000 1000 1000 <<< 2 || { >&2 echo "Failed to run command  at line ${LINENO} in bash script" && exit 1; }

