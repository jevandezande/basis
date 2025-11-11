"""Test the reading and writing of basis sets."""

from basis.basis import count, csv_table, difference


def test_count() -> None:
    """Test that the function counts are correct."""
    def2_svp = count("def2-svp")
    assert len(def2_svp) == 86
    assert def2_svp[1] == ([2, 1], [4, 1])


def test_diff() -> None:
    """Test that difference in counts works."""
    def2_svp = count("def2-svp")
    def2_tzvp = count("def2-tzvp")
    diff = difference(def2_svp, def2_tzvp)

    assert len(diff) == 86
    assert diff[1] == ([1, 0], [1, 0])
    assert diff[85] == ([2, 1, 1, 2], [1, 3, 2, 2])


def test_generally_contracted() -> None:
    """Test generally contracted basis sets."""
    cc_pVTZ = count("cc-pVTZ")

    assert cc_pVTZ[1] == ([3, 2, 1], [5, 2, 1])
    assert cc_pVTZ[3] == ([4, 3, 2, 1], [11, 5, 2, 1])
    assert cc_pVTZ[11] == ([5, 4, 2, 1], [16, 10, 2, 1])
    assert cc_pVTZ[21] == ([7, 6, 4, 2, 1], [20, 16, 8, 2, 1])


def test_sto() -> None:
    """Test STO-nG basis sets."""
    sto_3g = count("sto-3g")

    assert sto_3g[1] == ([1], [3])
    assert sto_3g[3] == ([2, 1], [6, 3])
    assert sto_3g[11] == ([3, 2], [9, 6])
    assert sto_3g[21] == ([4, 3, 1], [12, 9, 3])


def test_csv_table() -> None:
    """CSV output should list elements available for the requested basis."""
    basis_sets = ["sto-3g", "cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "aug-cc-pCVQZ"]
    counts = {k: count(k) for k in basis_sets}
    result = csv_table(counts, [1, 6, 9, 31])

    expected = """\
basis,element,contracted,uncontracted
sto-3g,1,"[1]","[3]"
sto-3g,6,"[2, 1]","[6, 3]"
sto-3g,9,"[2, 1]","[6, 3]"
sto-3g,31,"[4, 3, 1]","[12, 9, 3]"
cc-pVDZ,1,"[2, 1]","[4, 1]"
cc-pVDZ,6,"[3, 2, 1]","[9, 4, 1]"
cc-pVDZ,9,"[3, 2, 1]","[9, 4, 1]"
cc-pVDZ,31,"[5, 4, 2]","[14, 11, 6]"
cc-pVTZ,1,"[3, 2, 1]","[5, 2, 1]"
cc-pVTZ,6,"[4, 3, 2, 1]","[10, 5, 2, 1]"
cc-pVTZ,9,"[4, 3, 2, 1]","[10, 5, 2, 1]"
cc-pVTZ,31,"[6, 5, 3, 1]","[20, 13, 9, 1]"
cc-pVQZ,1,"[4, 3, 2, 1]","[6, 3, 2, 1]"
cc-pVQZ,6,"[5, 4, 3, 2, 1]","[12, 6, 3, 2, 1]"
cc-pVQZ,9,"[5, 4, 3, 2, 1]","[12, 6, 3, 2, 1]"
cc-pVQZ,31,"[7, 6, 4, 2, 1]","[21, 16, 12, 2, 1]"
aug-cc-pCVQZ,6,"[9, 8, 6, 4, 2]","[16, 10, 6, 4, 2]"
aug-cc-pCVQZ,9,"[9, 8, 6, 4, 2]","[16, 10, 6, 4, 2]"\
"""
    assert result == expected
