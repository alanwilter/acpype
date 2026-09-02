# Instructions for setting up a development environment

## For Coding

If using `Linux` with Intel processors, you can skip the `conda` step and directly to [Setting ACPYPE](#setting-acpype). Otherwise:

### Setting `conda`

For `Linux` (Ubuntu 22.04 or newer) and `macOS`. Anyway, `CONDA` is strongly recommended.
Also recommended is GPG key, so do accordingly in [GitHub](https://docs.github.com/articles/generating-a-gpg-key/).

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh

conda create -n acpype python=3.12 ambertools openbabel ocl-icd-system ipython gromacs=2019.1 -y
# ocl-icd-system: case you have GPU

conda activate acpype
```

### Setting ACPYPE

```bash
git clone https://github.com/alanwilter/acpype.git

cd acpype

uv sync

pre-commit install

pre-commit run -a

sys=$(perl -e '`uname` eq "Darwin\n" ? print "macos" : print "linux"')
cp ./acpype/amber_${sys}/bin/charmmgen $(dirname $(which antechamber))

git config user.email _email_ # as seen in 'gpg --list-secret-keys --keyid-format=long'

git config user.signingkey _singing_key_ # as seen in 'gpg --list-secret-keys --keyid-format=long'

git config commit.gpgsign true

uv run pytest --cov=tests --cov=acpype --cov-report=term-missing:skip-covered --cov-report=xml
```

### Linting and type checking

```bash
uv run ruff check .        # lint
uv run ruff format .       # format
uv run ty check            # type check
```

`ty` runs in CI and as a `pre-commit` hook, and the tree is expected to stay free of
diagnostics. `acpype/amber_linux`, `acpype/amber_macos`, `legacy` and `recipe` are
excluded from it (see `[tool.ty.src]` in `pyproject.toml`).

### Running the tests

Every test that touches the filesystem takes the `janitor` fixture, which gives it an
empty working directory and links the test inputs in by name. Nothing is written into
`tests/`, and two `pytest` runs can happen at once without colliding -- ACPYPE and
antechamber both scatter scratch files through the working directory, so sharing one
would make runs fail unpredictably.

If you add a test that reads an input file or runs `acpype`, take `janitor` too.

Coverage settings live in `[tool.pytest]`, so a plain `uv run pytest` already produces
the same report CI does, including the `fail_under` gate in `[tool.coverage.report]`.

### Refreshing the vendored AmberTools

`acpype` ships a trimmed AmberTools so `antechamber` works out of the box.

```bash
./update_macos_bins.sh -f   # needs conda/mamba, run on macOS
./update_linux_bins.sh      # needs Docker, runs a linux/amd64 container
```

Both call `scripts/vendor_amber.py`, which derives the shared libraries from the
executables rather than a hand-written list, so a library rename between AmberTools
releases is picked up automatically. On macOS `scripts/fix_macos_rpaths.py` then
removes the duplicate `LC_RPATH` entries that conda-forge's arm64 builds carry, which
dyld refuses to load.

`charmmgen` is a special case: modern AmberTools dropped it, and conda-forge's
`ambertools` package has never contained it, so ACPYPE builds its own for legacy
CHARMM output. It is untarred *after* vendoring and is therefore **not** part of the
dependency closure, so nothing guarantees a re-vendor leaves its libraries in place.

The macOS binary is rebuilt by `scripts/build_charmmgen.sh` from the source still
maintained in [Amber-MD/AmberClassic](https://github.com/Amber-MD/AmberClassic), as a
universal arm64 + x86_64 binary so Apple Silicon no longer needs Rosetta 2. The script
patches one line, because AmberClassic renamed `AMBERHOME` to `AMBERCLASSICHOME`.

`scripts/check_amber_bundle.py` guards this: it runs every executable in the bundle
and fails if the dynamic loader rejects any of them. Both update scripts run it, and
so does CI. On Linux it runs against a stock Ubuntu image carrying only the packages
documented as host requirements, so a library that should have been bundled cannot
hide behind conda's copy.

If using `VSCode`:

- Enable `enableCommitSigning` in `settings.json` (**_Workspace_** recommended):

  ```yml
  "git.enableCommitSigning": true,
  ```

- `uv sync` creates the virtualenv at `.venv` in the repository root, which `VSCode`
  and `Kiro` pick up automatically as the `Python Interpreter` (it is also set
  explicitly in `.vscode/settings.json`). No extra configuration is needed.

  Note that `ambertools`, `gromacs` and `openbabel` are **not** installed into `.venv`:
  `acpype` ships its own `antechamber` under `acpype/amber_${sys}`, and the rest come
  from the `conda` environment. If you need those on your `PATH`, run the tests from an
  activated `conda activate acpype` shell -- `uv run` will still use `.venv` for the
  Python dependencies.

## Releasing

`./release.sh -p` publishes to PyPI and `-d` builds and pushes the Docker images.

The PyPI step calls `scripts/build_wheels.py`, which produces **one wheel per
platform** rather than a single universal one. A combined wheel carries both vendored
AmberTools trees and comes to roughly 129 MB, over PyPI's 100 MB per-file limit; split,
they are about 49 MB (Linux) and 69 MB (macOS), and each installs only where its
binaries can actually run. The script fails if either wheel creeps back over the limit,
and CI builds them on every run so that cannot go unnoticed.

## For Documenting

Using [Sphinx](https://www.sphinx-doc.org) with the Read the Docs theme; both come from
the `docs` dependency group.

Online documentation provided by [Read the Docs](http://acpype.readthedocs.io).

To test it locally:

```bash
uv run --group docs sphinx-build -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html` in a browser. The `docs/Makefile` route
(`cd docs && make html`) also works from an environment with the `docs` group synced.
