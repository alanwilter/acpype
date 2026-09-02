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

# Run all the tests, but allow for arguments to be passed
test *ARGS:
    uv run pytest {{ ARGS }}

# Run all the tests, but on failure, drop into the debugger
pdb *ARGS:
    uv run pytest --pdb --maxfail=10 {{ ARGS }}

# Run the formatting, linting, type checking and tests commands
qa-all:
    just qa
    uv run pytest

# Upgrade the project libraries and rebuild
up *ARGS:
    uv lock --exclude-newer "7 days" -U {{ ARGS }}
    uv audit
    uv sync --exclude-newer "7 days" --all-groups
    just build

# Build the per-platform wheels and check they stay under the PyPI size limit
build:
    rm -rf dist
    uv run python scripts/build_wheels.py --out-dir dist

# Build the documentation
docs:
    uv run --group docs sphinx-build -b html docs docs/_build/html
    @echo "Open docs/_build/html/index.html"

# Confirm every executable in a vendored AmberTools bundle still loads
check-bundle SYS=os():
    uv run python scripts/check_amber_bundle.py acpype/amber_{{ if SYS == "macos" { "macos" } else { "linux" } }}

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
