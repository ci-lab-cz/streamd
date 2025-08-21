from streamd.run_gbsa import parse_gmxMMPBSA_decomp_dat, start

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



def test_start_aggregates_decomp_dat(tmp_path):
    gb_pb_content = (
        "Energy Decomposition Analysis (All units kcal/mol): Generalized Born model\n\n"
        "Complex:\n"
        "Total Energy Decomposition:\n"
        "Residue,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,\n"
        ",Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean\n"
        "R:A:LEU:83,1,0.1,0.01,2,0.2,0.02,3,0.3,0.03,4,0.4,0.04,5,0.5,0.05,6,0.6,0.06\n\n"
        "Energy Decomposition Analysis (All units kcal/mol): Poisson Boltzmann model\n"
        "Complex:\n"
        "Total Energy Decomposition:\n"
        "Residue,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,\n"
        ",Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean\n"
        "R:A:LEU:83,10,1,0.1,20,2,0.2,30,3,0.3,40,4,0.4,50,5,0.5,60,6,0.6\n"
    )
    # results_dat = tmp_path / "FINAL_RESULTS_MMPBSA_test.dat"
    # results_dat.write_text("dummy")
    decomp_dat = tmp_path / "FINAL_DECOMP_MMPBSA_test.dat"
    decomp_dat.write_text(gb_pb_content)
    # decomp_csv = tmp_path / "FINAL_DECOMP_MMPBSA_test.csv"
    # decomp_csv.write_text(
    #     "Generalized Born Decomposition Energies\nComplex:\nTotal Decomposition Contribution (TDC)\n"
    #     "Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL\n"
    #     "1,R:A:ALA:1,1,2,3,4,5,6\n"
    # )
    mmpbsa_in = tmp_path / "mmpbsa.in"
    mmpbsa_in.write_text("&general\n/\n&decomp\n/\n")
    start(
        wdir_to_run=None,
        tpr=None,
        xtc=None,
        topol=None,
        index=None,
        out_wdir=tmp_path,
        mmpbsa_file=str(mmpbsa_in),
        ncpu=1,
        ligand_resid="UNL",
        append_protein_selection=None,
        hostfile=None,
        unique_id="test",
        bash_log="bash.log",
        gmxmmpbsa_decomp_dat_files=[str(decomp_dat)],
        clean_previous=False,
        debug=False,
    )
    assert (tmp_path / "GBSA_decomp_avg_test.csv").is_file()
    assert (tmp_path / "PBSA_decomp_avg_test.csv").is_file()