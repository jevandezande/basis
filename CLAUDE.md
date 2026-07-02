# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Basis is a Python CLI tool for examining, comparing, and editing quantum chemistry basis sets using the Basis Set Exchange (BSE) library. It counts contracted and uncontracted basis functions, presents them in tabular format, and can selectively remove angular momentum types from basis sets and write the result to formatted output files.

## Development Commands

This project uses `uv` for dependency management and Python environment handling.

### Setup

```bash
uv sync                      # Install dependencies and the package (required for entry point)
```

### Testing

```bash
uv run pytest                                    # Run all tests (includes doctests)
uv run pytest tests/test_basis.py               # Run specific test file
uv run pytest -k test_count                     # Run specific test by name
uv run pytest --cov=basis --cov-report=term-missing  # Run with coverage report
```

### Linting and Formatting

```bash
uv run ruff format basis tests   # Format code
uv run ruff check .              # Lint code
uv run ty check basis            # Type check with ty
```

### Running the CLI

```bash
# Show basis set function counts (original functionality)
uv run basis show def2-SVP def2-TZVP -e H He 3 6 11
uv run basis show sto-3g def2-SVP -e 1 6 9 18 --diff
uv run basis show cc-pVDZ -f csv

# Edit a basis set fetched from BSE, write to file
# (output format is guessed from the file extension: .nw -> nwchem, .gbs -> gaussian94, ...)
uv run basis edit def2-TZVP \
    -e Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Ga Ge As Se Br Kr \
    --remove f -o def2-TZVP-f.nw

# Write Gaussian format just by naming the file .gbs (no -f needed)
uv run basis edit def2-TZVP --remove f -o def2-TZVP.gbs

# Edit a local basis set file (enables multi-round editing)
uv run basis edit --input def2-TZVP-f.nw --remove d -o def2-TZVP-fd.nw

# Print edited basis to stdout instead of a file
uv run basis edit def2-SVP -e C N O --remove d
```

### Pre-commit Hooks

```bash
pre-commit install           # Install pre-commit hooks
pre-commit run --all-files   # Run hooks manually
```

The pre-commit hooks run: YAML/TOML checks, ruff format, ruff check, ty check, and pytest.

## Code Architecture

### Core Components

**basis/basis.py** - Main module containing all basis set analysis and editing logic:

- `count(basis: str) -> BASIS_COUNT`: Queries BSE and counts contracted/uncontracted functions per element
- `table(basis_sets, elements, diff, format, spherical) -> str`: Orchestrates table generation
- `difference(basis1, basis2) -> BASIS_COUNT`: Computes differences between two basis sets
- `am_letter_to_int(letter: str) -> int`: Maps AM letter labels to integers (`s→0, p→1, d→2, f→3, ...`)
- `remove_angular_momentum(basis_dict, am_to_remove) -> dict`: Removes shells of specified AM types from a BSE basis dict; splits multi-AM shells (sp, spd) via `bse.manip.uncontract_spdf()` before filtering so only the targeted AM is removed
- `edit_basis(basis, elements, remove, fmt, input_file) -> str`: Orchestrates fetching/reading a basis set, element filtering, AM removal, and formatted string output
- `guess_format(path: str) -> str | None`: Maps an output file extension to a BSE writer format key (`.gbs → gaussian94`, `.nw → nwchem`, ...); returns `None` for unknown extensions. Uses `EXTENSION_TO_FORMAT`, a reverse of BSE's format→extension map with `_EXTENSION_FORMAT_PREFERENCES` resolving extensions shared by several writers (e.g. `.gbs`, `.bas`, `.json`)
- Output formatters: `plain_table()` and `csv_table()`

**basis/cli.py** - Command-line interface:

- Uses argparse with two subcommands: `show` and `edit`
- `show_parser()` / `show_cli()`: the original count/table functionality
- `edit_parser()` / `edit_cli()`: new edit/export functionality
- `basis_parser()`: top-level parser that wires the two subparsers together
- `basis_cli()`: entry point dispatching to `show_cli()` or `edit_cli()`

### Key Data Structures

- `BASIS_COUNT = dict[int, tuple[list[int], list[int]]]`: Maps element atomic numbers to (contracted_counts, uncontracted_counts) tuples
- Lists are indexed by angular momentum (s=0, p=1, d=2, etc.)
- Example: `{1: ([2, 1], [4, 1])}` means H has 2 contracted s, 1 contracted p, 4 uncontracted s, 1 uncontracted p

### Basis Set Exchange Integration

The package relies on the `basis_set_exchange` library. BSE ships all basis set data as bundled JSON files -- no network access is required. Key BSE API used:

- `bse.get_basis(name, elements=None)` → basis dict (or formatted string if `fmt` is given)
- `bse.read_formatted_basis_file(path)` → basis dict from a local file
- `bse.write_formatted_basis_str(basis_dict, fmt)` → formatted string
- `bse.get_writer_formats()` → dict of `{format_key: display_name}` for all supported output formats
- `bse.writers.get_format_extension(fmt)` → recommended file extension for a writer format (e.g. `'.gbs'`); used to build the extension→format map for `guess_format()`
- `bse.manip.uncontract_spdf(basis_dict)` → splits multi-AM shells into single-AM shells

The `basis_set_exchange` library has no type stubs, so all imports are marked `# type:ignore`.

### BSE Basis Dict Structure

```python
{
    "elements": {
        "6": {  # Z-number as string key
            "electron_shells": [
                {
                    "function_type": "gto",
                    "angular_momentum": [0],   # list of ints; [0,1] for sp shells
                    "exponents": ["..."],       # strings
                    "coefficients": [["..."]]   # list-of-lists; outer = general contractions
                },
                ...
            ]
        }
    }
}
```

All numeric values (exponents, coefficients) are stored as strings, not floats.

### Edit Subcommand Flow

1. Load basis dict: from BSE (`bse.get_basis`) or local file (`bse.read_formatted_basis_file`)
2. Filter elements if `-e/--elements` provided (passed to `bse.get_basis` for BSE source; post-read dict filtering for local files)
3. Remove angular momentum shells if `-r/--remove` provided via `remove_angular_momentum()`
4. Resolve the output format: use `-f/--format` if given, else `guess_format()` on the `-o/--output` extension, else fall back to `nwchem` (warning to stderr if an output extension was present but unrecognized)
5. Format and write output via `bse.write_formatted_basis_str()` to stdout or `-o/--output` file

Multi-round editing is enabled by writing to a file then reading it back with `--input`:

```bash
uv run basis edit def2-QZVP -e C --remove f -o intermediate.nw
uv run basis edit --input intermediate.nw --remove g -o final.nw
```

### Show Subcommand Flow

1. Count basis functions for each requested basis set via `count()`
2. Optionally convert to spherical counts via `spherical_count()`
3. Filter to requested elements using `filter_unused_elements_multi()`
4. Optionally compute differences with `difference()`
5. Format output via `plain_table()` or `csv_table()`

### Testing

Tests use doctests embedded in functions plus explicit test cases in `tests/test_basis.py`. The test suite covers:

- Standard basis sets (STO-nG) and correlation consistent basis sets (cc-pVXZ)
- Generally contracted basis sets
- Difference calculations and CSV output formatting
- `am_letter_to_int` including `ValueError` for invalid labels
- `remove_angular_momentum` with single-AM shells and with sp shells (verifies correct splitting)
- `edit_basis` with f-removal, element filtering, local file input, and multi-round editing

Doctests run automatically via pytest with `--doctest-modules` (configured in `pyproject.toml`).

## Linting Configuration

Ruff is configured with strict rules (line-length=100, Google-style docstrings). Several rules are disabled for complexity metrics (PLR09xx, PLR1702, PLR2004). Variable naming allows non-lowercase (`N806` disabled).

Type checking uses `ty` (Astral's type checker). Run as `uv run ty check basis`.

## Project Configuration

- **Build system**: `hatchling` (defined in `pyproject.toml`; required for `uv sync` to install the `basis` console script entry point)
- **Entry point**: `basis = "basis.cli:basis_cli"` under `[project.scripts]`
- **Main branch**: `master`
