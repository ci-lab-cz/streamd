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

Ligand Energy Terms
Frame #,BOND,ANGLE
1,5.0,6.0

Delta Energy Terms
Frame #,BOND,ANGLE
1,7.0,8.0
"""
    f = tmp_path / "FINAL_RESULTS_MMPBSA_test.csv"
    f.write_text(content)
    df = parse_gmxMMPBSA_results(f)
    assert set(df.columns) == {"Frame", "BOND", "ANGLE", "Region"}
    assert len(df) == 4
    assert df[df["Region"] == "Complex"]["BOND"].iloc[0] == 1.0
    assert df[df["Region"] == "Delta"]["ANGLE"].iloc[0] == 8.0
