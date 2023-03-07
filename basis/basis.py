from collections import defaultdict

import basis_set_exchange as bse  # type:ignore


def count(basis: str) -> dict[str, tuple[list[int], list[int]]]:
    data = bse.get_basis(basis)["elements"]
    counts = {}
    for element, values in data.items():
        contracted: dict[int, int] = defaultdict(int)
        uncontracted: dict[int, int] = defaultdict(int)
        for function in values["electron_shells"]:
            for am in function["angular_momentum"]:
                contracted[am] += 1
                uncontracted[am] += len(function["exponents"])

        con = [0] * (max(contracted) + 1)
        uncon = [0] * (max(contracted) + 1)
        for key, val in contracted.items():
            con[key] = val
        for key, val in uncontracted.items():
            uncon[key] = val

        counts[element] = (con, uncon)

    return counts


def find_max_am(counts: dict[str, dict[str, tuple[list[int], list[int]]]]) -> int:
    return max(
        len(element_data[0])
        for basis_data in counts.values()
        for element_data in basis_data.values()
    )


def table(basis_sets: list[str], elements: list[int] | None = None) -> None:
    counts = {basis: count(basis) for basis in basis_sets}
    elements = elements or list(range(1, 37))
    element_list = list(map(str, elements))
    counts = {
        basis: {
            element: element_data
            for element, element_data in basis_data.items()
            if element in element_list
        }
        for basis, basis_data in counts.items()
    }

    AM = "spdfgh"
    max_am = find_max_am(counts)
    BASIS_WIDTH = 6 * max_am + 1
    COL_WIDTH = 3 * max_am
    HLINE = "-" * (len(basis_sets) * (BASIS_WIDTH + 3) + 2) + "\n"

    out = "   | " + " | ".join(f"{basis:^{BASIS_WIDTH}s}" for basis in basis_sets) + "\n"
    out += "  " + f" |  {'  '.join(AM[:max_am])}" * 2 * len(basis_sets) + "\n"
    out += HLINE

    for element in element_list:
        if element in ["3", "11", "19", "37", "55", "87"]:  # starting a new row
            out += HLINE
        out += f"{element:>2} |"

        def count_str(basis: str) -> str:
            if element not in counts[basis]:
                return " " * BASIS_WIDTH

            con = "".join(f"{c:>3d}" for c in counts[basis][element][0])
            uncon = "".join(f"{c:>3d}" for c in counts[basis][element][1])

            return f"{uncon:<{COL_WIDTH}} →{con:<{COL_WIDTH}}"

        out += " |".join(map(count_str, basis_sets)) + "\n"

    print(out)


if __name__ == "__main__":
    # print(count("def2-SVP"))
    table(["def2-SVP", "def2-TZVP", "def2-QZVP"])
