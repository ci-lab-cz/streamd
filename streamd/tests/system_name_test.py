from streamd.analysis.md_system_analysis import _parse_system_name


def test_parse_system_name_replica():
    protein, replica = _parse_system_name("protein_HIS_replica7")
    assert protein == "protein_HIS"
    assert replica == 7
