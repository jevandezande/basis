"""Module for counting basis functions in basis sets using Basis Set Exchange."""

import os
from collections import defaultdict
from itertools import zip_longest
from typing import Container, Iterable, Literal, TypeVar

import basis_set_exchange as bse  # type:ignore
from basis_set_exchange import manip, writers  # type:ignore

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
atomic_dict = dict(zip(atomic_numbers, range(len(atomic_numbers)), strict=True))

spherical_harmonics = "spdfghiklmnoqrtuvwxyz"
spherical_harmonics_counts = {am: 2 * i + 1 for i, am in enumerate(spherical_harmonics)}
cartesian_harmonics_counts = {
    am: (i + 1) * (i + 2) // 2 for i, am in enumerate(spherical_harmonics)
}


def count(basis: str) -> BASIS_COUNT:
    """Count the number of contracted and uncontracted basis functions in a basis set.

    :param basis: basis set to count
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


def count_atomic_basis_functions(contracted_counts: list[int]) -> list[int]:
    """Count the resulting number of atomic basis functions from the basis set.

    :param contracted_counts: number of contracted basis functions
    :return: number of atomic basis functions

    >>> sto3g = count("sto-3g")
    >>> H, C, F, Ar = sto3g[1][0], sto3g[6][0], sto3g[9][0], sto3g[18][0]
    >>> count_atomic_basis_functions(H)
    [1]
    >>> count_atomic_basis_functions(C)
    [2, 3]
    >>> count_atomic_basis_functions(F)
    [2, 3]
    >>> count_atomic_basis_functions(Ar)
    [3, 6]
    """
    return [
        s * c for s, c in zip(spherical_harmonics_counts.values(), contracted_counts, strict=False)
    ]


def spherical_count(basis_counts: BASIS_COUNT) -> BASIS_COUNT:
    """Convert contracted basis function counts to spherical basis function counts.

    For each contracted basis function of angular momentum l, there are 2l+1 spherical
    basis functions. For example, 3 d-type contracted functions produce 3 x 5 = 15
    spherical basis functions.

    :param basis_counts: basis set counts to convert
    :return: spherical basis function counts (contracted only, uncontracted set to [])

    >>> sto3g = count("sto-3g")
    >>> spherical = spherical_count(sto3g)
    >>> spherical[1]
    ([1], [])
    >>> spherical[6]
    ([2, 3], [])
    >>> spherical[9]
    ([2, 3], [])
    >>> spherical[18]
    ([3, 6], [])
    """
    return {
        element: (count_atomic_basis_functions(contracted), [])
        for element, (contracted, _) in basis_counts.items()
    }


def find_max_am(counts: dict[str, BASIS_COUNT]) -> int:
    """Find the maximum angular momentum in a basis set.

    :param counts: basis sets to examine
    :return: maximum angular momentum

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
    """Filter out elements not in the list of elements.

    :param counts: basis set(s) to filter
    :param elements: elements to keep
    :return: filtered basis set

    >>> filter_unused_elements(count("sto-3g"), [1, 6, 9, 18])
    {1: ([1], [3]), 6: ([2, 1], [6, 3]), 9: ([2, 1], [6, 3]), 18: ([3, 2], [9, 6])}
    """
    return {element: cs for element, cs in counts.items() if element in elements}


def filter_unused_elements_multi(
    counts: dict[str, BASIS_COUNT],
    elements: Container[int],
) -> dict[str, BASIS_COUNT]:
    """Filter out elements not in the list of elements.

    :param counts: basis set(s) to filter
    :param elements: elements to keep
    :return: filtered basis sets

    >>> filter_unused_elements_multi({"sto-3g": count("sto-3g")}, [1, 6, 9, 18])
    {'sto-3g': {1: ([1], [3]), 6: ([2, 1], [6, 3]), 9: ([2, 1], [6, 3]), 18: ([3, 2], [9, 6])}}
    """
    return {
        basis: filter_unused_elements(basis_data, elements) for basis, basis_data in counts.items()
    }


def difference(basis1: BASIS_COUNT, basis2: BASIS_COUNT) -> BASIS_COUNT:
    """Find the difference between basis sets.

    :param basis1: first basis set
    :param basis2: second basis set
    :return: difference between the basis sets

    >>> sto3g = filter_unused_elements(count("sto-3g"), [1, 6, 9, 18])
    >>> sto6g = filter_unused_elements(count("sto-6g"), [1, 6, 9, 18])
    >>> difference(sto3g, sto6g)
    {1: ([0], [3]), 6: ([0, 0], [6, 3]), 9: ([0, 0], [6, 3]), 18: ([0, 0], [9, 6])}
    """

    def diff(c1: Iterable[int], c2: Iterable[int]) -> list[int]:
        return [b - a for a, b in zip_longest(c1, c2, fillvalue=0)]

    b2_set = set(basis2)
    return {
        element: (
            diff(uncon, basis2[element][0]),
            diff(con, basis2[element][1]),
        )
        for element, (uncon, con) in basis1.items()
        if element in b2_set
    }


def table(
    basis_sets: list[str],
    elements: Iterable[int | str] | None = None,
    diff: bool = False,
    format: Literal["plain", "csv"] = "plain",
    spherical: bool = False,
) -> str:
    """Generate a table of basis set counts.

    :param basis_sets: basis sets to compare
    :param elements: elements to include
    :param diff: include a difference column
    :param format: output format
    :param spherical: show spherical basis function counts instead of contracted/uncontracted
    :return: table

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

    if spherical:
        counts = {basis: spherical_count(basis_counts) for basis, basis_counts in counts.items()}

    if elements is None:
        element_list = list(range(1, 37))
    else:
        element_list = sorted(parse_elements(elements))

    counts = filter_unused_elements_multi(counts, element_list)

    if diff:
        if len(counts) != 2:
            raise ValueError(f"Can only compare two basis sets at a time, got: {len(basis_sets)=}")

        counts["Δ"] = difference(*(counts.values()))
        basis_sets += ["Δ"]

    match format:
        case "plain":
            return plain_table(counts, element_list, spherical)
        case "csv":
            return csv_table(counts, element_list, spherical)
        case _:
            raise ValueError(f"Unsupported format: {format}")


def plain_table(
    counts: dict[str, BASIS_COUNT],
    element_list: list[int],
    spherical: bool = False,
) -> str:
    """Generate a plain text table of basis set counts."""
    max_am = find_max_am(counts)

    # Header
    if spherical:
        BASIS_WIDTH = 3 * max_am + 1
        HLINE = "-" * (len(counts) * (BASIS_WIDTH + 1) + 2) + "\n"

        out = "   |" + "|".join(f"{basis:^{BASIS_WIDTH}s}" for basis in counts) + "\n"
        out += "  " + f" |  {'  '.join(spherical_harmonics[:max_am])}" * len(counts) + "\n"
    else:
        # Normal mode: show both uncontracted and contracted with arrow
        BASIS_WIDTH = 6 * max_am + 3
        COL_WIDTH = 3 * max_am
        HLINE = "-" * (len(counts) * (BASIS_WIDTH + 1) + 2) + "\n"

        out = "   |" + "|".join(f"{basis:^{BASIS_WIDTH}s}" for basis in counts) + "\n"
        out += "  " + f" |  {'  '.join(spherical_harmonics[:max_am])}" * 2 * len(counts) + "\n"

    row = 0
    rows = [0, 2, 10, 18, 36, 54, 86]

    def count_str(element: int, basis: str) -> str:
        if element not in counts[basis]:
            return " " * BASIS_WIDTH

        contracted = counts[basis][element][0]

        if spherical:
            # Pad contracted list to max_am length and show blank spaces for zero counts
            padded = contracted + [0] * (max_am - len(contracted))
            return "".join(f"{c:>3d}" if c else "   " for c in padded)
        else:
            uncontracted = counts[basis][element][1]
            con = "".join(f"{c:>3d}" for c in contracted)
            uncon = "".join(f"{c:>3d}" for c in uncontracted)
            return f"{uncon:<{COL_WIDTH}} →{con:<{COL_WIDTH}}"

    for element in element_list:
        if element > rows[row]:
            out += HLINE
            row = searchsorted(element, rows)
        out += f"{atomic_numbers[element]:2} |"

        out += " |".join(count_str(element, basis) for basis in counts).rstrip() + "\n"

    return out.rstrip()


def csv_table(
    counts: dict[str, BASIS_COUNT],
    element_list: list[int],
    spherical: bool = False,
) -> str:
    """Generate a CSV table of basis set counts.

    >>> counts = {"sto-3g": count("sto-3g"), "def2-svp": count("def2-svp")}
    >>> print(csv_table(counts, [1, 6, 9, 18]))
    basis,element,contracted,uncontracted
    sto-3g,1,"[1]","[3]"
    sto-3g,6,"[2, 1]","[6, 3]"
    sto-3g,9,"[2, 1]","[6, 3]"
    sto-3g,18,"[3, 2]","[9, 6]"
    def2-svp,1,"[2, 1]","[4, 1]"
    def2-svp,6,"[3, 2, 1]","[7, 4, 1]"
    def2-svp,9,"[3, 2, 1]","[7, 4, 1]"
    def2-svp,18,"[4, 3, 1]","[10, 7, 1]"
    """

    def _quote(value: list[int]) -> str:
        field = str(value).replace('"', '\\"')
        return f'"{field}"'

    if spherical:
        rows = ["basis,element,spherical"]
        for basis, basis_counts in counts.items():
            for element in element_list:
                if element not in basis_counts:
                    continue

                contracted, _uncontracted = basis_counts[element]
                rows.append(f"{basis},{element},{_quote(contracted)}")
    else:
        rows = ["basis,element,contracted,uncontracted"]
        for basis, basis_counts in counts.items():
            for element in element_list:
                if element not in basis_counts:
                    continue

                contracted, uncontracted = basis_counts[element]
                rows.append(
                    f"{basis},{element},{_quote(contracted)},{_quote(uncontracted)}",
                )

    return "\n".join(rows)


def element_to_an(element: int | str) -> int:
    """Convert element to atomic number.

    :param element: element to convert
    :return: atomic number

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
    """Convert elements to atomic number.

    :param elements: elements to convert
    :return: atomic numbers

    >>> elements_to_an([36, "W"])
    [36, 74]
    """
    return list(map(element_to_an, elements))


def parse_elements(tokens: Iterable[int | str]) -> list[int]:
    """Parse element tokens into a list of atomic numbers, expanding ranges.

    Each token is either a single element (symbol or atomic number) or a range of
    the form `start-end` where both endpoints are inclusive and may be symbols or
    atomic numbers.  Tokens are whitespace- or argument-separated; the hyphen is used
    exclusively as a range separator.

    :param tokens: iterable of element tokens (e.g. `["H", "Li-Ne", "19-20"]`)
    :return: sorted, deduplicated list of atomic numbers

    >>> parse_elements(["H"])
    [1]
    >>> parse_elements(["Li-Ne"])
    [3, 4, 5, 6, 7, 8, 9, 10]
    >>> parse_elements(["H", "Li-Ne", "19"])
    [1, 3, 4, 5, 6, 7, 8, 9, 10, 19]
    >>> parse_elements(["3-5", "Ne"])
    [3, 4, 5, 10]
    """
    result: set[int] = set()
    for token in tokens:
        s = str(token)
        if "-" in s:
            left, right = s.split("-", 1)
            start = element_to_an(left)
            end = element_to_an(right)
            if start > end:
                raise ValueError(f"Invalid range {s!r}: start must be <= end")
            result.update(range(start, end + 1))
        else:
            result.add(element_to_an(s))
    return sorted(result)


def am_letter_to_int(letter: str) -> int:
    """Convert an angular momentum letter label to its integer value.

    :param letter: angular momentum label (e.g. 's', 'p', 'd', 'f')
    :return: integer angular momentum value

    >>> am_letter_to_int('s')
    0
    >>> am_letter_to_int('p')
    1
    >>> am_letter_to_int('d')
    2
    >>> am_letter_to_int('f')
    3
    """
    if letter not in spherical_harmonics:
        raise ValueError(f"Unknown angular momentum label: {letter!r}")
    return spherical_harmonics.index(letter)


def remove_angular_momentum(
    basis_dict: dict,
    am_to_remove: set[int],
    elements: set[int] | None = None,
) -> dict:
    """Remove shells with specified angular momenta from a basis set dictionary.

    Multi-angular-momentum shells (e.g. sp, spd) are split into single-AM shells
    before filtering so that only the targeted AM is removed.

    :param basis_dict: BSE basis set dictionary (from bse.get_basis or bse.read_formatted_basis_*)
    :param am_to_remove: set of integer angular momentum values to remove
    :param elements: atomic numbers of elements to edit; `None` edits all elements
    :return: new basis set dictionary with the specified shells removed from the target elements
    """
    result = manip.uncontract_spdf(basis_dict)
    for z, element_data in result["elements"].items():
        if elements is not None and int(z) not in elements:
            continue
        element_data["electron_shells"] = [
            shell
            for shell in element_data["electron_shells"]
            if not set(shell["angular_momentum"]).issubset(am_to_remove)
        ]
    return result


# BSE assigns each writer a recommended extension, but several formats share one
# (e.g. .gbs -> gaussian94/psi4/xtron).  A canonical format is preferred for those.
_EXTENSION_FORMAT_PREFERENCES = {
    ".gbs": "gaussian94",
    ".bas": "gamess_us",
    ".json": "json",
    ".molcas": "molcas",
}


def _build_extension_format_map() -> dict[str, str]:
    """Build a mapping of lowercased file extension -> BSE writer format key."""
    mapping: dict[str, str] = {}
    for fmt in bse.get_writer_formats():
        ext = writers.get_format_extension(fmt).lower()
        mapping.setdefault(ext, fmt)
    mapping.update(_EXTENSION_FORMAT_PREFERENCES)
    return mapping


EXTENSION_TO_FORMAT = _build_extension_format_map()


def guess_format(path: str) -> str | None:
    """Guess a BSE writer format key from a file path's extension.

    :param path: output file path (e.g. `'def2-TZVP.gbs'`)
    :return: BSE writer format key, or `None` if the extension is unrecognized

    >>> guess_format("def2-TZVP.gbs")
    'gaussian94'
    >>> guess_format("basis.nw")
    'nwchem'
    >>> guess_format("basis.unknown") is None
    True
    >>> guess_format("no_extension") is None
    True
    """
    ext = os.path.splitext(path)[1].lower()
    return EXTENSION_TO_FORMAT.get(ext) if ext else None


def edit_basis(
    basis: str | None,
    elements: Iterable[int | str] | None = None,
    remove: Iterable[str] | None = None,
    fmt: str = "nwchem",
    input_file: str | None = None,
) -> str:
    """Fetch or read a basis set, optionally remove angular momentum types, and format it.

    When *input_file* is given the basis set is read from that local file (enabling
    multi-round editing).  Otherwise the basis set is fetched from the Basis Set
    Exchange by name.

    When *elements* is specified, AM removal applies only to those elements; all other
    elements in the basis set are written unchanged.

    :param basis: BSE basis set name; required when *input_file* is not provided
    :param elements: elements whose shells will be edited; `None` edits all elements
    :param remove: angular momentum letter labels to remove (e.g. `['f', 'g']`)
    :param fmt: output format key accepted by BSE (default `'nwchem'`)
    :param input_file: path to a local formatted basis set file to read instead of BSE
    :return: formatted basis set string

    >>> result = edit_basis("sto-3g", elements=["H", "C"], remove=["p"], fmt="nwchem")
    >>> "H    S" in result and "C    P" not in result
    True
    """
    if input_file is not None:
        basis_dict: dict = bse.read_formatted_basis_file(input_file)
    else:
        if basis is None:
            raise ValueError("Either 'basis' or 'input_file' must be provided")
        basis_dict = bse.get_basis(basis)

    if remove:
        element_set = set(parse_elements(elements)) if elements is not None else None
        am_to_remove = {am_letter_to_int(letter) for letter in remove}
        basis_dict = remove_angular_momentum(basis_dict, am_to_remove, element_set)

    return bse.write_formatted_basis_str(basis_dict, fmt)


T = TypeVar("T", int, float, str)


def searchsorted(value: T, target: Iterable[T], reversed: bool = False) -> int:
    """Find where in a sorted iterable a value would fit.

    Note: values matching existing values are placed after
    :param value: value to insert
    :param target: iterable to examine
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
