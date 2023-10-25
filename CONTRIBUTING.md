# Instructions for setting up a development environment

## For Coding

If using `Linux` with Intel processors, you can skip the `conda` step and directly to [Setting ACPYPE](#setting-acpype). Otherwise:

### Setting `conda`

For `Linux` (Ubuntu 20 recommended) and `macOS`. Anyway, `CONDA` is strongly recommended.
Also recommended is GPG key, so do accordingly in [GitHub](https://docs.github.com/articles/generating-a-gpg-key/).

```bash
curl -sSL https://install.python-poetry.org | python3 -

conda create -n acpype python=3.9 ambertools openbabel ocl-icd-system ipython gromacs=2019.1 -y
# ocl-icd-system: case you have GPU

conda activate acpype
```

### Setting ACPYPE

```bash
git clone https://github.com/alanwilter/acpype.git

cd acpype

poetry install

pre-commit install

pre-commit run -a

cp ./acpype/amber_linux/bin/charmmgen $(dirname $(which antechamber))

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

- Another `VSCode` nuisance, if using its graphic buttons for commit etc.: its environment won't necessarily recognise the dependencies installed via `poetry` under `conda` environment named `acpype` (unless you have started `VSCode` from folder repository with command `code .`). To avoid troubles do:

  ```bash
  conda deactivate
  poetry install # it will create its own virtualenv
  ```

  You could use this `poetry virtualenv` as the `Python Interpreter` for `VSCode`, however `ambertools`, `gromacs` and `openbabel` won't be available (unless you've had installed them system wide by other means rather than `conda`).
  To avoid further troubles, go back to `conda activate acpype` and remember to do the instructions above if you add new dependencies to the project via `poetry`.

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
