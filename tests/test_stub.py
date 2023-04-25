from basis.basis import count


def test_count() -> None:
    def2_svp = count("def2-svp")
    assert len(def2_svp) == 86
    assert def2_svp[1] == ([2, 1], [4, 1])
