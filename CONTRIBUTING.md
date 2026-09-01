# Instructions for setting up a development environment

## For Coding

If using `Linux` with Intel processors, you can skip the `conda` step and directly to [Setting ACPYPE](#setting-acpype). Otherwise:

### Setting `conda`

For `Linux` (Ubuntu 20 recommended) and `macOS`. Anyway, `CONDA` is strongly recommended.
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

pytest --cov=tests --cov=acpype --cov-report=term-missing:skip-covered --cov-report=xml
```

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

## For Documenting

Using [Sphinx](https://www.sphinx-doc.org) with theme from `pip install sphinx-rtd-theme`.

Online documentation provided by [Read the Docs](http://acpype.readthedocs.io).

To test it locally:

```bash
cd docs/
make clean
make html
```

Then open `_build/html/index.html` in a browser.
