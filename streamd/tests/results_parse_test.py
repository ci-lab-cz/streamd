import pandas as pd
from streamd.run_gbsa import parse_gmxMMPBSA_results, parse_gmxMMPBSA_output

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


def test_parse_gmxMMPBSA_output(tmp_path):
    content = """ENTROPY RESULTS (INTERACTION ENTROPY):
Energy Method          Entropy σ(Int. Energy)    Average           SD        SEM
-------------------------------------------------------------------------------
GB                          IE          3.28       0.63         0.99       0.44

Energy Method          Entropy σ(Int. Energy)    Average           SD        SEM
-------------------------------------------------------------------------------
PB                          IE          3.28       0.63         0.99       0.44

GENERALIZED BORN:
Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
ΔTOTAL                  -51.34          0.80       2.82         0.36       1.26
Using Interaction Entropy Approximation:
ΔG binding =    -50.71 +/-    2.99

POISSON BOLTZMANN:
Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
ΔTOTAL                  -46.69          1.08       2.69         0.48       1.20
Using Interaction Entropy Approximation:
ΔG binding =    -46.06 +/-    2.87
"""
    f = tmp_path / "FINAL_RESULTS_MMPBSA_test.dat"
    f.write_text(content)
    res = parse_gmxMMPBSA_output(f)
    assert res["GBSA"]["ΔTOTAL_Average"] == "-51.34"
    assert res["PBSA"]["ΔTOTAL_Average"] == "-46.69"
