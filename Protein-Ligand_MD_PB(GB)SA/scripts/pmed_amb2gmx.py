import parmed as pmd
import argparse

def arg_parser():
	parser = argparse.ArgumentParser(description="Convert AMBER topology to GROMACS using ParmED")
	parser.add_argument("--prmtop", "-p", type=str,
						help="*prmtop file")
	parser.add_argument("--inpcrd", "-x", type=str,
						help="*inpcrd file")
	parser.add_argument("--output", "-o", type=str,
						help="Output name")
	return parser


parser = arg_parser()
args = parser.parse_args()

# convert AMBER topology to GROMACS, CHARMM formats
amber = pmd.load_file(args.prmtop, args.inpcrd)
# Save a GROMACS topology and GRO file
amber.save(args.output+'.top')
amber.save(args.output+'.gro')
