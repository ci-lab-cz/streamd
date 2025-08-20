import pandas as pd
from streamd.run_gbsa import parse_gmxMMPBSA_results

def test_parse_gmxMMPBSA_results(tmp_path):
    content = """GENERALIZED BORN:
Complex Energy Terms
Frame #,BOND,ANGLE
1,1.0,2.0

Receptor Energy Terms
Frame #,BOND,ANGLE
1,3.0,4.0

POISSON BOLTZMANN:
Ligand Energy Terms
Frame #,BOND,ANGLE
1,5.0,6.0
"""
    f = tmp_path / "FINAL_RESULTS_MMPBSA_test.csv"
    f.write_text(content)
    df = parse_gmxMMPBSA_results(f)
    assert set(df.columns) == {"Frame", "BOND", "ANGLE", "Region", "Method"}
    assert len(df) == 3
    gb_complex = df[(df["Region"] == "Complex") & (df["Method"] == "GB")]["BOND"].iloc[0]
    assert gb_complex == 1.0
    pb_ligand = df[(df["Region"] == "Ligand") & (df["Method"] == "PB")]["ANGLE"].iloc[0]
    assert pb_ligand == 6.0
