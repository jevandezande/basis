"""Command-line interface for basis set examination."""

from argparse import ArgumentParser

from .basis import table


def basis_parser(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Create an ArgumentParser for basis set examination."""
    parser = parser or ArgumentParser(description="Run basis")

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


def basis_cli() -> None:
    """Run the basis set examination CLI."""
    args = basis_parser().parse_args()

    print(table(args.basis, args.elements, args.diff, args.format, args.spherical))


if __name__ == "__main__":
    basis_cli()
