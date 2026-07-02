"""Test the reading and writing of basis sets."""

from pathlib import Path

import basis_set_exchange as bse
import pytest

from basis.basis import (
    am_letter_to_int,
    count,
    csv_table,
    difference,
    edit_basis,
    parse_elements,
    remove_angular_momentum,
    spherical_count,
)


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


def test_spherical_count() -> None:
    """Test that spherical basis function counting works correctly."""
    sto3g = count("sto-3g")
    spherical = spherical_count(sto3g)

    # H: 1 s-function → 1 x 1 = 1 spherical function
    assert spherical[1] == ([1], [])

    # C: 2 s-functions, 1 p-function → 2 x 1 = 2, 1 x 3 = 3
    assert spherical[6] == ([2, 3], [])

    # Test d-orbitals: Sc has s, p, and d functions
    assert spherical[21] == ([4, 9, 5], [])  # 4sx1, 3px3, 1dx5

    # Test that uncontracted is always empty
    assert spherical[1][1] == []
    assert spherical[6][1] == []
    assert spherical[21][1] == []


def test_spherical_count_higher_am() -> None:
    """Test spherical counting with higher angular momentum basis sets."""
    cc_pvtz = count("cc-pVTZ")
    spherical = spherical_count(cc_pvtz)

    # H: 3 s-functions, 2 p-functions, 1 d-function
    # → 3x1=3, 2x3=6, 1x5=5
    assert spherical[1] == ([3, 6, 5], [])

    # C: 4 s, 3 p, 2 d, 1 f → 4x1=4, 3x3=9, 2x5=10, 1x7=7
    assert spherical[6] == ([4, 9, 10, 7], [])


def test_am_letter_to_int() -> None:
    """Test angular momentum letter to integer conversion."""
    assert am_letter_to_int("s") == 0
    assert am_letter_to_int("p") == 1
    assert am_letter_to_int("d") == 2
    assert am_letter_to_int("f") == 3
    assert am_letter_to_int("g") == 4

    with pytest.raises(ValueError, match="Unknown angular momentum label"):
        am_letter_to_int("a")


def test_remove_angular_momentum() -> None:
    """Test that angular momentum shells are correctly removed from a basis dict."""
    bd = bse.get_basis("def2-TZVP", elements=[6])
    shells_before = bd["elements"]["6"]["electron_shells"]
    ams_before = {am for shell in shells_before for am in shell["angular_momentum"]}
    assert 3 in ams_before  # C has f functions in def2-TZVP

    result = remove_angular_momentum(bd, {3})
    shells_after = result["elements"]["6"]["electron_shells"]
    ams_after = {am for shell in shells_after for am in shell["angular_momentum"]}

    assert 3 not in ams_after
    assert {0, 1, 2}.issubset(ams_after)  # s, p, d still present


def test_remove_angular_momentum_sp_shell() -> None:
    """Test that removing p from an sp shell keeps the s component."""
    # sto-3g C has an sp shell
    bd = bse.get_basis("sto-3g", elements=[6])
    sp_shells = [
        s for s in bd["elements"]["6"]["electron_shells"] if s["angular_momentum"] == [0, 1]
    ]
    assert sp_shells, "Expected an sp shell in sto-3g C"

    result = remove_angular_momentum(bd, {1})  # remove p
    shells_after = result["elements"]["6"]["electron_shells"]
    ams_after = {am for shell in shells_after for am in shell["angular_momentum"]}

    assert 1 not in ams_after  # p removed
    assert 0 in ams_after  # s still present


def test_edit_basis_remove_f() -> None:
    """Test that edit_basis removes f functions only from specified elements."""
    # def2-TZVP: C has f, H does not; only edit C
    result = edit_basis("def2-TZVP", elements=["C"], remove=["f"], fmt="nwchem")

    # C loses its f shell
    assert "C    F" not in result
    assert "C    S" in result
    assert "C    P" in result
    assert "C    D" in result

    # H is untouched and still present in the full basis output
    assert "H    S" in result


def test_edit_basis_elements_scope() -> None:
    """Test that AM removal is scoped to the specified elements; others are unchanged."""
    # sto-3g: both C and N have p functions; only remove p from C
    result = edit_basis("sto-3g", elements=["C"], remove=["p"], fmt="nwchem")

    assert "C    P" not in result  # edited
    assert "C    S" in result  # C s-functions preserved
    assert "N    P" in result  # N is untouched


def test_edit_basis_from_file(tmp_path: Path) -> None:
    """Test round-trip: write a basis to file, read it back via --input."""
    # Write a basis set to a temp file
    bd = bse.get_basis("sto-3g", elements=[1, 6])
    basis_str = bse.write_formatted_basis_str(bd, "nwchem")
    input_file = tmp_path / "sto-3g.nw"
    input_file.write_text(basis_str)

    # Read it back and remove p functions
    result = edit_basis(basis=None, remove=["p"], fmt="nwchem", input_file=str(input_file))

    assert "H    S" in result
    assert "C    S" in result
    assert "C    P" not in result


def test_edit_basis_multi_round(tmp_path: Path) -> None:
    """Test multi-round editing: remove f then remove d from an intermediate file."""
    # Round 1: remove f from C only; full basis written to file
    round1 = edit_basis("def2-TZVP", elements=["C"], remove=["f"], fmt="nwchem")
    intermediate = tmp_path / "intermediate.nw"
    intermediate.write_text(round1)

    assert "C    F" not in round1
    assert "C    D" in round1
    assert "H    S" in round1  # other elements present

    # Round 2: remove d from C in the intermediate file, leaving other elements intact
    round2 = edit_basis(
        basis=None, elements=["C"], remove=["d"], fmt="nwchem", input_file=str(intermediate)
    )

    assert "C    D" not in round2
    assert "C    S" in round2
    assert "C    P" in round2
    assert "H    S" in round2  # H still untouched


def test_parse_elements_single() -> None:
    """Test parsing individual element tokens."""
    assert parse_elements(["H"]) == [1]
    assert parse_elements(["He"]) == [2]
    assert parse_elements(["1"]) == [1]
    assert parse_elements([6]) == [6]


def test_parse_elements_range_symbols() -> None:
    """Test expanding a symbol range."""
    assert parse_elements(["Li-Ne"]) == [3, 4, 5, 6, 7, 8, 9, 10]


def test_parse_elements_range_numbers() -> None:
    """Test expanding a numeric range."""
    assert parse_elements(["3-5"]) == [3, 4, 5]


def test_parse_elements_mixed() -> None:
    """Test a mix of singles, numeric ranges, and symbol ranges."""
    assert parse_elements(["H", "Li-Be", "11-12", "Al"]) == [1, 3, 4, 11, 12, 13]


def test_parse_elements_dedup_and_sorted() -> None:
    """Test that results are deduplicated and sorted."""
    assert parse_elements(["C", "B-N", "6"]) == [5, 6, 7]


def test_parse_elements_invalid_range() -> None:
    """Test that a reversed range raises ValueError."""
    with pytest.raises(ValueError, match="Invalid range"):
        parse_elements(["Ne-Li"])


def test_parse_elements_in_edit_basis() -> None:
    """Test that range tokens work end-to-end through edit_basis."""
    # Remove p from Li through Ne (period 2 main group); N and O should lose P shells
    result = edit_basis("sto-3g", elements=["Li-Ne"], remove=["p"], fmt="nwchem")

    assert "N    P" not in result  # N is in Li-Ne, p removed
    assert "N    S" in result  # N s-functions kept
    assert "H    S" in result  # H outside range, untouched
