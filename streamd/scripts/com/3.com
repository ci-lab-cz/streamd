$RunGauss
%Chk=ligand.chk
%Mem=120GB
%NProcShared=1
#N B3LYP/6-31G* Freq Geom=AllCheckpoint Guess=Read
Integral=(Grid=UltraFine) SCF=XQC IOp(7/33=1)
