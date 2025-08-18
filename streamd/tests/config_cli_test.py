"""Tests for CLI and YAML config precedence across entry points."""

import sys
import types
import yaml


def test_run_md_cli_overrides_config(monkeypatch, tmp_path):
    """Run run_md.main ensuring CLI args override YAML and extras are ignored."""

    dummy_module = types.SimpleNamespace
    # Stub internal modules to avoid heavy dependencies
    monkeypatch.setitem(sys.modules, "streamd.analysis.md_system_analysis", dummy_module(run_md_analysis=lambda *a, **k: None))
    monkeypatch.setitem(sys.modules, "streamd.analysis.run_analysis", dummy_module(run_rmsd_analysis=lambda *a, **k: None))
    monkeypatch.setitem(sys.modules, "streamd.preparation.complex_preparation", dummy_module(run_complex_preparation=lambda *a, **k: None))
    monkeypatch.setitem(sys.modules, "streamd.preparation.ligand_preparation", dummy_module(prepare_input_ligands=lambda *a, **k: None, check_mols=lambda *a, **k: None))
    monkeypatch.setitem(sys.modules, "streamd.utils.dask_init", dummy_module(init_dask_cluster=lambda *a, **k: None, calc_dask=lambda *a, **k: None))
    monkeypatch.setitem(
        sys.modules,
        "streamd.utils.utils",
        dummy_module(
            filepath_type=lambda *a, **k: a[0],
            run_check_subprocess=lambda *a, **k: True,
            get_protein_resid_set=lambda *a, **k: None,
            backup_prev_files=lambda *a, **k: None,
            check_to_continue_simulation_time=lambda *a, **k: True,
            merge_parts_of_simulation=lambda *a, **k: None,
        ),
    )
    monkeypatch.setitem(sys.modules, "streamd.mcpbpy_md", dummy_module(mcbpy_md=lambda *a, **k: None))

    from streamd import run_md

    captured = {}

    def fake_start(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(run_md, "start", fake_start)

    protein = tmp_path / "protein.pdb"
    protein.write_text("p")
    ligand = tmp_path / "ligand.mol"
    ligand.write_text("l")
    config = {"ligand": str(ligand), "ncpu": 4, "unused": "x"}
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.dump(config))

    monkeypatch.setattr(
        sys,
        "argv",
        ["run_md", "--config", str(config_path), "--protein", str(protein), "--ncpu", "8"],
    )
    run_md.main()

    assert captured["protein"] == str(protein)
    assert captured["ligand"] == str(ligand)
    assert captured["ncpu"] == 8
    assert "unused" not in captured


def test_run_gbsa_cli_overrides_config(monkeypatch, tmp_path):
    """Run run_gbsa.main checking CLI precedence and config usage."""

    dummy_utils = types.SimpleNamespace(
        filepath_type=lambda *a, **k: a[0],
        get_index=lambda *a, **k: ["Protein", "LIG"],
        make_group_ndx=lambda *a, **k: True,
        run_check_subprocess=lambda *a, **k: True,
        get_number_of_frames=lambda *a, **k: 1,
        temporary_directory_debug=lambda *a, **k: types.SimpleNamespace(__enter__=lambda self: tmp_path, __exit__=lambda *a: False),
    )
    dummy_dask = types.SimpleNamespace(init_dask_cluster=lambda *a, **k: None, calc_dask=lambda *a, **k: None)
    monkeypatch.setitem(sys.modules, "streamd.utils.utils", dummy_utils)
    monkeypatch.setitem(sys.modules, "streamd.utils.dask_init", dummy_dask)
    monkeypatch.setitem(sys.modules, "pandas", types.ModuleType("pandas"))

    from streamd import run_gbsa

    captured = {}

    def fake_start(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(run_gbsa, "start", fake_start)

    tpr = tmp_path / "file.tpr"
    tpr.write_text("t")
    xtc = tmp_path / "file.xtc"
    xtc.write_text("x")
    topol = tmp_path / "file.top"
    topol.write_text("top")
    index = tmp_path / "index.ndx"
    index.write_text("Protein\n")

    config = {"ligand_id": "ABC", "ncpu": 4, "unused": "x"}
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.dump(config))

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_gbsa",
            "--config",
            str(config_path),
            "--tpr",
            str(tpr),
            "--xtc",
            str(xtc),
            "--topol",
            str(topol),
            "--index",
            str(index),
            "--ncpu",
            "2",
        ],
    )
    run_gbsa.main()

    assert captured["tpr"] == str(tpr)
    assert captured["ligand_resid"] == "ABC"
    assert captured["ncpu"] == 2
    assert "unused" not in captured


def test_run_prolif_cli_overrides_config(monkeypatch, tmp_path):
    """Run run_prolif.main validating CLI/config interplay."""

    monkeypatch.setitem(sys.modules, "MDAnalysis", types.ModuleType("MDAnalysis"))
    monkeypatch.setitem(sys.modules, "pandas", types.ModuleType("pandas"))
    monkeypatch.setitem(sys.modules, "prolif", types.ModuleType("prolif"))

    barcode_mod = types.ModuleType("prolif.plotting.barcode")
    barcode_mod.Barcode = type(
        "Barcode",
        (),
        {
            "from_fingerprint": classmethod(
                lambda cls, fp: types.SimpleNamespace(
                    display=lambda *a, **k: types.SimpleNamespace(
                        figure=types.SimpleNamespace(savefig=lambda *a, **k: None)
                    )
                )
            )
        },
    )
    network_mod = types.ModuleType("prolif.plotting.network")
    network_mod.LigNetwork = type(
        "LigNetwork",
        (),
        {
            "from_fingerprint": classmethod(
                lambda cls, fp, ligand_mol=None, threshold=None: types.SimpleNamespace(save=lambda *a, **k: None)
            )
        },
    )
    monkeypatch.setitem(sys.modules, "prolif.plotting", types.ModuleType("prolif.plotting"))
    monkeypatch.setitem(sys.modules, "prolif.plotting.barcode", barcode_mod)
    monkeypatch.setitem(sys.modules, "prolif.plotting.network", network_mod)

    plt_mod = types.ModuleType("matplotlib.pyplot")
    plt_mod.ioff = lambda: None
    monkeypatch.setitem(sys.modules, "matplotlib", types.ModuleType("matplotlib"))
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", plt_mod)

    monkeypatch.setitem(
        sys.modules,
        "streamd.utils.dask_init",
        types.SimpleNamespace(init_dask_cluster=lambda *a, **k: None, calc_dask=lambda *a, **k: None),
    )
    monkeypatch.setitem(
        sys.modules, "streamd.utils.utils", types.SimpleNamespace(filepath_type=lambda *a, **k: a[0])
    )
    monkeypatch.setitem(sys.modules, "streamd.prolif.prolif2png", types.ModuleType("streamd.prolif.prolif2png"))
    monkeypatch.setitem(
        sys.modules,
        "streamd.prolif.prolif_frame_map",
        types.ModuleType("streamd.prolif.prolif_frame_map"),
    )

    from streamd.prolif import run_prolif

    captured = {}

    def fake_start(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(run_prolif, "start", fake_start)

    tpr = tmp_path / "file.tpr"
    tpr.write_text("t")
    xtc = tmp_path / "file.xtc"
    xtc.write_text("x")
    config = {"ligand": "LIG1", "ncpu": 4, "unused": "x"}
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.dump(config))

    monkeypatch.setattr(
        sys,
        "argv",
        ["run_prolif", "--config", str(config_path), "--tpr", str(tpr), "--xtc", str(xtc), "--ncpu", "2"],
    )
    run_prolif.main()

    assert captured["tpr"] == str(tpr)
    assert captured["ligand_resid"] == "LIG1"
    assert captured["ncpu"] == 2
    assert "unused" not in captured

