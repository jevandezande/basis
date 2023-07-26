from argparse import ArgumentParser

from .basis import table


def basis_parser(parser: ArgumentParser | None = None) -> ArgumentParser:
    parser = parser or ArgumentParser(description="Run basis")

    parser.add_argument(
        "basis",
        nargs="+",
    )
    parser.add_argument(
        "-e",
        "--elements",
        nargs="+",
    )
    return parser


def basis_cli() -> None:
    args = basis_parser().parse_args()

    table(args.basis, args.elements)
