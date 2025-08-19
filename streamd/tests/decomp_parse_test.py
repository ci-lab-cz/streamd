import pandas as pd
from streamd.run_gbsa import parse_gmxMMPBSA_decomp

def test_parse_gmxMMPBSA_decomp(tmp_path):
    content = """Generalized Born Decomposition Energies
Complex:
Total Decomposition Contribution (TDC)
Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
1,R:A:ALA:1,1,2,3,4,5,6

Sidechain Decomposition Contribution (SDC)
Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
1,R:A:ALA:1,0.1,0.2,0.3,0.4,0.5,0.6

Receptor:
Backbone Decomposition Contribution (BDC)
Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
1,R:A:ALA:1,0.7,0.8,0.9,1.0,1.1,1.2

Poisson Boltzmann Decomposition Energies
Complex:
Total Decomposition Contribution (TDC)
Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
1,R:A:ALA:1,10,20,30,40,50,60
"""
    f = tmp_path / "FINAL_DECOMP_MMPBSA_test.csv"
    f.write_text(content)
    df = parse_gmxMMPBSA_decomp(f)
    assert set(df.columns) == {
        "Frame",
        "Residue",
        "Internal",
        "van der Waals",
        "Electrostatic",
        "Polar Solvation",
        "Non-Polar Solv.",
        "TOTAL",
        "Region",
        "Contribution",
        "Method",
    }
    assert len(df) == 4
    complex_tdc = df[(df["Region"] == "Complex") & (df["Contribution"] == "TDC") & (df["Method"] == "GB")].iloc[0]
    assert complex_tdc["Internal"] == 1.0
    receptor_bdc = df[(df["Region"] == "Receptor") & (df["Contribution"] == "BDC")].iloc[0]
    assert receptor_bdc["TOTAL"] == 1.2
    pb_row = df[(df["Region"] == "Complex") & (df["Method"] == "PB")].iloc[0]
    assert pb_row["Electrostatic"] == 30
