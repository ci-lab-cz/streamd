# Data Preparation

Before running MD, ensure the system is correctly prepared to avoid failures during topology generation or simulation.

## Hydrogen Handling and Protonation
StreaMD re-adds hydrogens by default using `gmx pdb2gmx -water tip3p -ignh`. This ignores original hydrogens and recreates them based on residue names.

> CAUTION: Rename residues to match correct protonation states (CYS->CYX, HIS->HID/HIE/HIP) before running. If you prefer to keep existing hydrogens, pass `--noignh`, but you must resolve any atom name recognition issues yourself or provide prepared `protein.gro`, `topol.top`, and `posre.itp`.

## Target Preparation
- Fill missing residues and loops.
- Optional: remove explicit water and cofactors from crystal structures unless they must be retained.
- Remove co-crystallized ligands.
- Add hydrogens according to protonation states; rename histidines to HID/HIE/HIP as appropriate.
- Co-crystallized ligands can be downloaded from the PDB as `.sdf` and prepared separately.

- Chimera commands for histidine renaming:
  ```
  setattr r type HID :HIS@HD1,DD1,TD1,HND
  setattr r type HIP :HID@HE2,DE2,TE2
  setattr r type HIE :HIS@HE2
  ```
- Optional tooling: Chimera (Modeller loop/refinement, Dock Prep) or pdb4amber.

Prepared examples and scripts: https://github.com/ci-lab-cz/docking-files

## Ligand Preparation (Optional)
- Ensure ligand/cofactor coordinates are aligned to the protein if running multiple ligands.
- Ligands/cofactors must be provided as `.mol` or `.sdf` files with correct 3D coordinates, hydrogens, protonation state, stereochemistry, and desired tautomer; StreaMD does not fix these automatically.

## Docking Procedure (Optional)
If ligand poses are needed, perform docking before MD:
- Use EasyDock for automated docking: https://github.com/ci-lab-cz/easydock
- EasyDock can handle protonation and stereoisomers and produces docking poses in `.sdf`, matching the required input format.


Ready to run? Proceed to {doc}`running_md` for command examples and continuation/extension workflows.
