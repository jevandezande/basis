from basis.basis import count, difference


def test_count() -> None:
    def2_svp = count("def2-svp")
    assert len(def2_svp) == 86
    assert def2_svp[1] == ([2, 1], [4, 1])


def test_diff() -> None:
    def2_svp = count("def2-svp")
    def2_tzvp = count("def2-tzvp")
    diff = difference(def2_svp, def2_tzvp)

    assert len(diff) == 86
    assert diff[1] == ([1, 0], [1, 0])
    assert diff[85] == ([2, 1, 1], [1, 3, 2])


def test_generally_contracted() -> None:
    cc_pVTZ = count("cc-pVTZ")

    assert cc_pVTZ[1] == ([3, 2, 1], [5, 2, 1])
    assert cc_pVTZ[3] == ([4, 3, 2, 1], [11, 5, 2, 1])
    assert cc_pVTZ[11] == ([5, 4, 2, 1], [16, 10, 2, 1])
    assert cc_pVTZ[21] == ([7, 6, 4, 2, 1], [20, 16, 8, 2, 1])


def test_sto() -> None:
    sto_3g = count("sto-3g")

    assert sto_3g[1] == ([1], [3])
    assert sto_3g[3] == ([2, 1], [6, 3])
    assert sto_3g[11] == ([3, 2], [9, 6])
    assert sto_3g[21] == ([4, 3, 1], [12, 9, 3])
