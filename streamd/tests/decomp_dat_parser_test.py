from streamd.run_gbsa import parse_gmxMMPBSA_decomp_dat

def test_parse_gmxMMPBSA_decomp_dat(tmp_path):
    content = (
        "Energy Decomposition Analysis (All units kcal/mol): Generalized Born model\n\n"
        "Complex:\n"
        "Total Energy Decomposition:\n"
        "Residue,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,\n"
        ",Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean\n"
        "R:A:LEU:83,1,0.1,0.01,2,0.2,0.02,3,0.3,0.03,4,0.4,0.04,5,0.5,0.05,6,0.6,0.06\n\n"
        "Sidechain Energy Decomposition:\n"
        "Residue,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,\n"
        ",Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean\n"
        "R:A:LEU:83,1.1,0.11,0.011,2.1,0.21,0.021,3.1,0.31,0.031,4.1,0.41,0.041,5.1,0.51,0.051,6.1,0.61,0.061\n\n"
        "Energy Decomposition Analysis (All units kcal/mol): Poisson Boltzmann model\n"
        "Complex:\n"
        "Total Energy Decomposition:\n"
        "Residue,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,\n"
        ",Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean\n"
        "R:A:LEU:83,10,1,0.1,20,2,0.2,30,3,0.3,40,4,0.4,50,5,0.5,60,6,0.6\n"
    )
    f = tmp_path / "FINAL_DECOMP_MMPBSA_test.dat"
    f.write_text(content)
    df = parse_gmxMMPBSA_decomp_dat(f)
    assert not df.empty
    assert "Internal Avg." in df.columns
    assert "TOTAL Std. Err. of Mean" in df.columns
    gb_total = df[(df["Region"] == "Complex") & (df["Contribution"] == "Total") & (df["Method"] == "GB")].iloc[0]
    assert gb_total["van der Waals Avg."] == 2.0
    sidechain = df[(df["Contribution"] == "Sidechain")].iloc[0]
    assert sidechain["Electrostatic Std. Dev."] == 0.31
    pb_total = df[(df["Method"] == "PB")].iloc[0]
    assert pb_total["Polar Solvation Avg."] == 40.0
