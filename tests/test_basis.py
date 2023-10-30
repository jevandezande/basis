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
