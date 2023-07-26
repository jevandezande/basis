# Basis

[![License](https://img.shields.io/github/license/jevandezande/basis)](https://github.com/jevandezande/basis/blob/master/LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/jevandezande/basis/test.yml?branch=master)](https://github.com/jevandezande/basis/actions/)
[![Codecov](https://img.shields.io/codecov/c/github/jevandezande/basis)](https://codecov.io/gh/jevandezande/basis)

# Usage
```
❯ poetry run basis def2-SVP def2-TZVP -e H He 3 6 11
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


## Credits
This package was created with Cookiecutter and the [jevandezande/poetry-cookiecutter](https://github.com/jevandezande/poetry-cookiecutter) project template.
