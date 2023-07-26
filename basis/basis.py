from collections import defaultdict
from typing import Iterable, TypeVar

import basis_set_exchange as bse  # type:ignore

# fmt:off
atomic_numbers = [
    'X', 'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', # noqa: E501
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', # noqa: E501
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', # noqa: E501
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cp', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo', # noqa: E501
]
# fmt:on
atomic_dict = dict(zip(atomic_numbers, range(len(atomic_numbers))))


def count(basis: str) -> dict[int, tuple[list[int], list[int]]]:
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

        counts[int(element)] = (con, uncon)

    return counts


def find_max_am(counts: dict[str, dict[int, tuple[list[int], list[int]]]]) -> int:
    if not (
        ams := [
            len(element_data[0])
            for basis_data in counts.values()
            for element_data in basis_data.values()
        ]
    ):
        return 0
    return max(ams)


def table(basis_sets: list[str], elements: Iterable[int | str] | None = None) -> str:
    counts = {basis: count(basis) for basis in basis_sets}
    if elements is None:
        element_list = list(range(1, 37))
    else:
        element_list = sorted(elements_to_an(elements))

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

    row = 0
    rows = [0, 2, 10, 18, 36, 54, 86]
    for element in element_list:
        if element > rows[row]:
            out += HLINE
            row = searchsorted(element, rows)
        out += f"{atomic_numbers[element]:2} |"

        def count_str(basis: str) -> str:
            if element not in counts[basis]:
                return " " * BASIS_WIDTH

            con = "".join(f"{c:>3d}" for c in counts[basis][element][0])
            uncon = "".join(f"{c:>3d}" for c in counts[basis][element][1])

            return f"{uncon:<{COL_WIDTH}} →{con:<{COL_WIDTH}}"

        out += " |".join(map(count_str, basis_sets)) + "\n"

    return out


def elements_to_an(elements: Iterable[int | str]) -> list[int]:
    es: list[int] = []
    for element in elements:
        if isinstance(element, int):
            es.append(element)
        elif element.isdigit():
            es.append(int(element))
        else:
            es.append(atomic_dict[element])
    return es


T = TypeVar("T", int, float, str)


def searchsorted(value: T, target: Iterable[T], reversed: bool = False) -> int:
    """
    Find where in a sorted iterable a value would fit.

    Note: values matching existing values are placed after
    :param value: value to insert
    :param target: the iterable to examine
    :param reversed: is the target sorted in descending order
    """
    i = 0
    for i, v in enumerate(target):
        if reversed:
            if value > v:
                return i
        elif value < v:
            return i
    return i
