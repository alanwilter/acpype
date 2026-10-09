# Justfile for ACPYPE
#
# `just` itself comes from the `rust-just` dev dependency, so `uv sync` provides it.

# Show available commands
list:
    @just --list

# Run the formatting, linting and type checking commands
qa:
    uv run ruff format .
    uv run ruff check . --fix
    uv run ty check
    uv audit

# Check formatting, linting and type checking (no fixes, for CI)
ci:
    uv run ruff format --check .
    uv run ruff check .
    uv run ty check
    uv audit

# This is not a read-only check. Besides ruff and ty, which `qa` already runs, the
# ver_today hook stamps today's date as the version in pyproject.toml,
# src/acpype/__init__.py and uv.lock, repoints the README badges at the newest tag,
# and `git add`s all four. That is what you want immediately before a commit, and not
# what you want in the middle of a QA sweep -- hence its own recipe. To run the hooks
# without stamping a version: SKIP=ver_today just pre-commit

# Run every pre-commit hook over the whole tree (stamps the version, see above)
pre-commit:
    uv run pre-commit run -a

# Run all the tests, but allow for arguments to be passed
test *ARGS:
    uv run pytest {{ ARGS }}

# Run all the tests, but on failure, drop into the debugger
pdb *ARGS:
    uv run pytest --pdb --maxfail=10 {{ ARGS }}

# The hooks add what `qa` does not cover: trailing whitespace, missing final newlines,
# oversized files and stray debug statements. ver_today is skipped because it would
# stamp today's date as the version and stage four files, which a QA sweep must not do.

# Run the formatting, linting, type checking, hooks and tests commands
qa-all:
    just qa
    SKIP=ver_today uv run pre-commit run -a
    uv run pytest

# autoupdate moves the hook revs in .pre-commit-config.yaml, which are pinned
# independently of the dependency bounds in pyproject.toml, so bumping both here keeps
# them in step. Expect a gap in either direction anyway: autoupdate has no equivalent of
# --exclude-newer and always takes the newest tag, while the uv steps above hold the
# project a week behind, so any ruff released in that window lands in the hook first.
# Both read the same [tool.ruff] config, so a gap only bites when a release changes
# formatting or adds a rule.

# Upgrade the project libraries and pre-commit hooks, then rebuild
up *ARGS:
    uv lock --exclude-newer "7 days" -U {{ ARGS }}
    uv audit
    uv sync --exclude-newer "7 days" --all-groups
    uv run pre-commit autoupdate
    just build

# Build the per-platform wheels and check they stay under the PyPI size limit
build:
    rm -rf dist
    uv run python scripts/build_dists.py --out-dir dist

# Build the documentation
docs:
    uv run --group docs sphinx-build -b html docs docs/_build/html
    @echo "Open docs/_build/html/index.html"

# Confirm every executable in a vendored AmberTools bundle still loads
check-bundle SYS=os():
    uv run python scripts/check_amber_bundle.py src/acpype/amber_{{ if SYS == "macos" { "macos" } else { "linux" } }}

# Re-vendor AmberTools for macOS (needs conda/mamba, must run on macOS)
vendor-macos:
    ./update_macos_bins.sh -f

# Re-vendor AmberTools for Linux (needs Docker)
vendor-linux:
    ./update_linux_bins.sh

# Rebuild the bundled charmmgen from AmberClassic (macos | linux | all)
charmmgen TARGET="all":
    ./scripts/build_charmmgen.sh {{ TARGET }}

# Print the current version of the project
version:
    @uv run python -c "import acpype; print(acpype.__version__)"

# Remove build, test and coverage artefacts
clean:
    rm -rf dist build .pytest_cache .ruff_cache htmlcov docs/_build
    rm -f .coverage .coverage.* coverage.xml
    find . -name '__pycache__' -not -path './.venv/*' -not -path './acpype/amber_*' -exec rm -rf {} +
    find . -name '*.pyc' -not -path './.venv/*' -not -path './acpype/amber_*' -delete
