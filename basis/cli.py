"""Command-line interface for basis set examination and editing."""

import sys
from argparse import ArgumentParser, Namespace

import basis_set_exchange as bse

from .basis import edit_basis, guess_format, table


def show_parser(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Create an ArgumentParser for the 'show' subcommand."""
    parser = parser or ArgumentParser(description="Show basis set function counts.")

    parser.add_argument("basis", nargs="+", help="Basis set(s) to examine.")
    parser.add_argument("-e", "--elements", nargs="+", help="Elements to examine.")
    parser.add_argument(
        "-d", "--diff", action="store_true", help="Find the difference between two basis sets."
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=("plain", "csv"),
        default="plain",
        help="Output format [%(default)s].",
    )
    parser.add_argument(
        "-s",
        "--spherical",
        action="store_true",
        help="Show spherical basis function counts (contracted functions x spherical harmonics).",
    )

    return parser


def edit_parser(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Create an ArgumentParser for the 'edit' subcommand."""
    parser = parser or ArgumentParser(description="Edit and export a basis set.")

    parser.add_argument(
        "basis",
        nargs="?",
        default=None,
        help="Basis set name to fetch from BSE. Required unless --input is given.",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="FILE",
        help="Read from a local formatted basis set file instead of BSE.",
    )
    parser.add_argument("-e", "--elements", nargs="+", help="Elements to include.")
    parser.add_argument(
        "-r",
        "--remove",
        nargs="+",
        metavar="AM",
        help="Angular momentum types to remove (e.g. f g).",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        default=None,
        help="Output file path. Defaults to stdout.",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=sorted(bse.get_writer_formats()),
        default=None,
        help="Output format. Guessed from the output file extension when omitted, "
        "falling back to nwchem.",
    )

    return parser


def basis_parser() -> ArgumentParser:
    """Create the top-level ArgumentParser with 'show' and 'edit' subparsers."""
    parser = ArgumentParser(description="Examine and edit basis sets from the Basis Set Exchange.")
    subparsers = parser.add_subparsers(dest="subcommand")

    show_parser(subparsers.add_parser("show", help="Show basis set function counts."))
    edit_parser(subparsers.add_parser("edit", help="Edit and export a basis set."))

    return parser


def show_cli(args: Namespace) -> None:
    """Run the 'show' subcommand."""
    print(table(args.basis, args.elements, args.diff, args.format, args.spherical))


def edit_cli(args: Namespace) -> None:
    """Run the 'edit' subcommand."""
    basis = args.basis
    input_file = args.input
    output = args.output

    if basis is None and input_file is None:
        sys.exit("error: one of 'basis' or '--input' is required")

    fmt = args.format
    if fmt is None:
        if output is not None:
            fmt = guess_format(output)
            if fmt is None:
                print(
                    f"warning: could not guess format from '{output}', using nwchem",
                    file=sys.stderr,
                )
        fmt = fmt or "nwchem"

    result = edit_basis(
        basis=basis,
        elements=args.elements,
        remove=args.remove,
        fmt=fmt,
        input_file=input_file,
    )

    if output is None:
        print(result)
    else:
        with open(output, "w") as fh:
            fh.write(result)


def basis_cli() -> None:
    """Run the basis set CLI."""
    parser = basis_parser()
    args = parser.parse_args()

    match args.subcommand:
        case "show":
            show_cli(args)
        case "edit":
            edit_cli(args)
        case _:
            parser.print_help()
            sys.exit(1)


if __name__ == "__main__":
    basis_cli()
