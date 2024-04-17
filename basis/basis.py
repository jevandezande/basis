from collections import defaultdict
from typing import Container, Iterable, TypeVar

import basis_set_exchange as bse  # type:ignore

# {element: (contracted_counts, uncontracted_counts)}
BASIS_COUNT = dict[int, tuple[list[int], list[int]]]

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


def count(basis: str) -> BASIS_COUNT:
    """
    Count the number of contracted and uncontracted basis functions for each element in a basis set

    :param basis: the basis set to count
    :return: a dictionary of element to a tuple of contracted and uncontracted counts

    >>> sto3g = count("sto-3g")
    >>> sto3g[1], sto3g[6], sto3g[9], sto3g[18]
    (([1], [3]), ([2, 1], [6, 3]), ([2, 1], [6, 3]), ([3, 2], [9, 6]))
    """
    data = bse.get_basis(basis)["elements"]
    counts = {}
    for element, values in data.items():
        contracted: dict[int, int] = defaultdict(int)
        uncontracted: dict[int, int] = defaultdict(int)
        for function in values["electron_shells"]:
            angular_momenta = function["angular_momentum"]
            coefficients = function["coefficients"]
            exponents = function["exponents"]
            for am in angular_momenta:
                # Generally contracted basis sets only specify a single angular momentum,
                # but have multiple sets of coefficients
                if len(angular_momenta) == 1 and isinstance(coefficients, list):
                    contracted[am] += len(coefficients)
                else:
                    contracted[am] += 1

                # Generally contracted basis sets share exponents
                uncontracted[am] += len(exponents)

        con = [0] * (max(contracted) + 1)
        uncon = [0] * (max(contracted) + 1)
        for key, val in contracted.items():
            con[key] = val
        for key, val in uncontracted.items():
            uncon[key] = val

        counts[int(element)] = (con, uncon)

    return counts


def find_max_am(counts: dict[str, BASIS_COUNT]) -> int:
    """
    Find the maximum angular momentum in a basis set

    :param counts: the basis sets to examine
    :return: the maximum angular momentum

    >>> find_max_am({"sto-3g": count("sto-3g")})
    3
    """
    if not (
        ams := [
            len(element_data[0])
            for basis_data in counts.values()
            for element_data in basis_data.values()
        ]
    ):
        return 0
    return max(ams)


def filter_unused_elements(counts: BASIS_COUNT, elements: Container[int]) -> BASIS_COUNT:
    """
    Filter out elements not in the list of elements

    :param counts: the basis set(s) to filter
    :param elements: the elements to keep
    :return: the filtered basis set

    >>> filter_unused_elements(count("sto-3g"), [1, 6, 9, 18])
    {1: ([1], [3]), 6: ([2, 1], [6, 3]), 9: ([2, 1], [6, 3]), 18: ([3, 2], [9, 6])}
    """
    return {element: cs for element, cs in counts.items() if element in elements}


def filter_unused_elements_multi(
    counts: dict[str, BASIS_COUNT],
    elements: Container[int],
) -> dict[str, BASIS_COUNT]:
    """
    Filter out elements not in the list of elements

    :param counts: the basis set(s) to filter
    :param elements: the elements to keep
    :return: the filtered basis sets

    >>> filter_unused_elements_multi({"sto-3g": count("sto-3g")}, [1, 6, 9, 18])
    {'sto-3g': {1: ([1], [3]), 6: ([2, 1], [6, 3]), 9: ([2, 1], [6, 3]), 18: ([3, 2], [9, 6])}}
    """
    return {
        basis: filter_unused_elements(basis_data, elements) for basis, basis_data in counts.items()
    }


def difference(basis1: BASIS_COUNT, basis2: BASIS_COUNT) -> BASIS_COUNT:
    """
    Find the difference between basis sets

    :param basis1: the first basis set
    :param basis2: the second basis set
    :return: the difference between the basis sets

    >>> sto3g = filter_unused_elements(count("sto-3g"), [1, 6, 9, 18])
    >>> sto6g = filter_unused_elements(count("sto-6g"), [1, 6, 9, 18])
    >>> difference(sto3g, sto6g)
    {1: ([0], [3]), 6: ([0, 0], [6, 3]), 9: ([0, 0], [6, 3]), 18: ([0, 0], [9, 6])}
    """

    def diff(c1: list[int], c2: list[int]) -> list[int]:
        return [b - a for a, b in zip(c1, c2)]

    b2_set = set(basis2)
    return {
        element: (
            diff(counts[0], basis2[element][0]),
            diff(counts[1], basis2[element][1]),
        )
        for element, counts in basis1.items()
        if element in b2_set
    }


def table(
    basis_sets: list[str],
    elements: Iterable[int | str] | None = None,
    diff: bool = False,
) -> str:
    """
    Generate a table of basis set counts

    :param basis_sets: the basis sets to compare
    :param elements: the elements to include
    :param diff: include a difference column
    :return: the table

    >>> print(table(["sto-3g", "sto-6g"], [1, 6, 9, 18], diff=True))
       |    sto-3g     |    sto-6g     |       Δ
       |  s  p |  s  p |  s  p |  s  p |  s  p |  s  p
    --------------------------------------------------
    H  |  3    →  1    |  6    →  1    |  3    →  0
    --------------------------------------------------
    C  |  6  3 →  2  1 | 12  6 →  2  1 |  6  3 →  0  0
    F  |  6  3 →  2  1 | 12  6 →  2  1 |  6  3 →  0  0
    --------------------------------------------------
    Ar |  9  6 →  3  2 | 18 12 →  3  2 |  9  6 →  0  0
    """
    counts: dict[str, BASIS_COUNT] = {basis: count(basis) for basis in basis_sets}
    if elements is None:
        element_list = list(range(1, 37))
    else:
        element_list = sorted(elements_to_an(elements))

    counts = filter_unused_elements_multi(counts, element_list)

    if diff:
        if len(counts) != 2:
            raise ValueError(f"Can only compare two basis sets at a time, got: {len(basis_sets)=}")

        counts["Δ"] = difference(*(counts.values()))
        basis_sets += ["Δ"]

    AM = "spdfghiklmnoqrtuvwxyz"
    max_am = find_max_am(counts)
    BASIS_WIDTH = 6 * max_am + 3
    COL_WIDTH = 3 * max_am
    HLINE = "-" * (len(basis_sets) * (BASIS_WIDTH + 1) + 2) + "\n"

    out = "   |" + "|".join(f"{basis:^{BASIS_WIDTH}s}" for basis in basis_sets) + "\n"
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

            return f"{uncon:<{COL_WIDTH}} →{con:<{COL_WIDTH}} "

        out += "|".join(map(count_str, basis_sets)).rstrip() + "\n"

    return out


def element_to_an(element: int | str) -> int:
    """
    Convert element to atomic number

    :param element: the element to convert
    :return: the atomic number

    >>> element_to_an(1)
    1
    >>> element_to_an("He")
    2
    """
    if isinstance(element, int):
        return element
    elif element.isdigit():
        return int(element)
    return atomic_dict[element]


def elements_to_an(elements: Iterable[int | str]) -> list[int]:
    """
    Convert elements to atomic number

    :param elements: the elements to convert
    :return: the atomic numbers

    >>> elements_to_an([36, "W"])
    [36, 74]
    """
    return list(map(element_to_an, elements))


T = TypeVar("T", int, float, str)


def searchsorted(value: T, target: Iterable[T], reversed: bool = False) -> int:
    """
    Find where in a sorted iterable a value would fit.

    Note: values matching existing values are placed after
    :param value: value to insert
    :param target: the iterable to examine
    :param reversed: is the target sorted in descending order

    >>> searchsorted(3, [2, 3, 4, 5])
    2
    >>> searchsorted(3, [5, 4, 3, 2], reversed=True)
    3
    """
    i = 0
    for i, v in enumerate(target):
        if reversed:
            if value > v:
                return i
        elif value < v:
            return i
    return i
