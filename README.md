# basis

[![License](https://img.shields.io/github/license/jevandezande/basis)](https://github.com/jevandezande/basis/blob/master/LICENSE)
[![Powered by: uv](https://img.shields.io/badge/-uv-purple)](https://docs.astral.sh/uv)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Markdown style: rumdl](https://img.shields.io/badge/md%20style-rumdl-000000.svg)](https://rumdl.dev)
[![Typing: ty](https://img.shields.io/badge/typing-ty-EFC621.svg)](https://github.com/astral-sh/ty)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/jevandezande/basis/test.yml?branch=master&logo=github-actions)](https://github.com/jevandezande/basis/actions)
[![Codecov](https://img.shields.io/codecov/c/github/jevandezande/basis)](https://codecov.io/gh/jevandezande/basis)

## Usage

### Show

Compare basis set function counts across elements:

```text
❯ uv run basis show def2-SVP def2-TZVP -e H He 3 6 11
   |         def2-SVP          |         def2-TZVP
   |  s  p  d  f |  s  p  d  f |  s  p  d  f |  s  p  d  f
----------------------------------------------------------
H  |  4  1       →  2  1       |  5  1       →  3  1
He |  4  1       →  2  1       |  5  1       →  3  1
----------------------------------------------------------
Li |  7  3       →  3  2       | 11  3       →  5  3
C  |  7  4  1    →  3  2  1    | 11  6  2  1 →  5  3  2  1
----------------------------------------------------------
Na | 10  6  1    →  4  2  1    | 14  8  3    →  5  4  3
```

Element ranges are supported: `-e H-Ne Na-Ar` is equivalent to listing each element individually.

### Edit

Remove, filter, and export basis sets. The `--remove` flag accepts one or more angular
momentum labels (`s p d f g ...`). Only the specified elements are modified; all others
are written unchanged.

#### Example: def2-TZVP(-f)

Remove f functions from all main-group elements to produce def2-TZVP(-f):

```text
❯ uv run basis edit def2-TZVP \
    -e H-He Li-Ne Na-Ar K Ca Ga-Kr Rb Sr In-Xe Cs Ba Tl-Rn \
    --remove f -o def2-TZVP-f.nw
```

Multi-round editing is supported by piping through a local file with `--input`:

```text
❯ uv run basis edit def2-QZVP -e C-Ne --remove f -o intermediate.nw
❯ uv run basis edit --input intermediate.nw -e C-Ne --remove g -o final.nw
```

### Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [jevandezande/uv-cookiecutter](https://github.com/jevandezande/uv-cookiecutter) project template.
