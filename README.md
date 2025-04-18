# ACPYPE

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg?style=plastic)](https://GitHub.com/alanwilter/acpype/graphs/commit-activity)<!-- [![Open Source Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103&style=plastic)](https://github.com/ellerbrock/open-source-badges/) -->
[![GitHub](https://img.shields.io/github/license/alanwilter/acpype?style=plastic)](https://github.com/alanwilter/acpype)
[![python](https://img.shields.io/badge/python-3.8%2E%2E%2E3.12-blue.svg?style=plastic&logo=python)](https://github.com/alanwilter/acpype)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/alanwilter/acpype?display_name=tag&logo=github&style=plastic)](https://github.com/alanwilter/acpype)
[![GitHub Release](https://img.shields.io/github/release-date/alanwilter/acpype?style=plastic&logo=github)](https://github.com/alanwilter/acpype)<!-- ![GitHub All Releases](https://img.shields.io/github/downloads/alanwilter/acpype/total?style=plastic) -->
[![Docker Pulls](https://img.shields.io/docker/pulls/acpype/acpype?style=plastic&logo=docker)](https://hub.docker.com/r/acpype/acpype)
[![Docker Image Size (tag)](https://img.shields.io/docker/image-size/acpype/acpype/latest?style=plastic&logo=docker)](https://hub.docker.com/r/acpype/acpype/tags)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/acpype.svg?style=plastic&logo=conda-forge)](https://anaconda.org/conda-forge/acpype)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/acpype.svg?style=plastic&logo=conda-forge)](https://anaconda.org/conda-forge/acpype/files)<!-- ![Conda](https://img.shields.io/conda/pn/conda-forge/acpype?logo=conda-forge&style=plastic) -->
[![PyPI](https://img.shields.io/pypi/v/acpype?style=plastic&logo=pypi)](https://pypi.org/project/acpype/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/acpype?style=plastic&logo=pypi)](https://pypi.org/project/acpype/#files)
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/alanwilter/acpype/.github%2Fworkflows%2Fcheck_acpype.yml?style=plastic)
[![Poetry](https://img.shields.io/endpoint?style=plastic&url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![Ruff](https://img.shields.io/endpoint?style=plastic&url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white&style=plastic)](https://github.com/pre-commit/pre-commit)
[![Commits since release](https://img.shields.io/github/commits-since/alanwilter/acpype/2023.10.27/master?style=plastic)](https://github.com/alanwilter/acpype/commits/master)
[![Codecov](https://img.shields.io/codecov/c/github/alanwilter/acpype?style=plastic)](https://app.codecov.io/gh/alanwilter/acpype)
[![Documentation Status](https://readthedocs.org/projects/acpype/badge/?version=latest&style=plastic)](https://acpype.readthedocs.io/en/latest/?badge=latest)
[![Citations Badge](https://img.shields.io/endpoint?label=citations&logo=googlescholar&style=plastic&url=https%3A%2F%2Fapi.juleskreuer.eu%2Fcitation-badge.php%3Fshield%26doi%3D10.1186%2F1756-0500-5-367)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=c68TiIUAAAAJ&citation_for_view=c68TiIUAAAAJ:UeHWp8X0CEIC)

<!-- [![Citations](https://img.shields.io/endpoint?label=citations&logo=googlescholar&style=plastic&url=https%3A%2F%2Fcitations-acpype-1g51xat28cdv.runkit.sh%2F)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=c68TiIUAAAAJ&citation_for_view=c68TiIUAAAAJ:UeHWp8X0CEIC) -->

<!-- ![Scrutinizer code quality (GitHub/Bitbucket)](https://img.shields.io/scrutinizer/quality/g/alanwilter/acpype) -->
<!-- ![Scrutinizer coverage (GitHub/BitBucket)](https://img.shields.io/scrutinizer/coverage/g/alanwilter/acpype) -->

## AnteChamber PYthon Parser interfacE

A tool based in **Python** to use **Antechamber** to generate topologies for chemical
compounds and to interface with others python applications like CCPN and ARIA.

`acpype` is pronounced as **_ace + pipe_**

Topologies files to be generated so far: CNS/XPLOR, GROMACS, CHARMM and AMBER.

**NB:** Topologies generated by `acpype/Antechamber` are based on General Amber Force
Field (GAFF) and should be used only with compatible forcefields like AMBER and
its variant.

Several flavours of AMBER FF are ported already for GROMACS (see [ffamber](http://ffamber.cnsm.csulb.edu/)) as well as to XPLOR/CNS (see [`xplor-nih`](http://ambermd.org/xplor-nih.html)) and [CHARMM](https://www.charmm.org/).

This code is released under **[GNU General Public Licence V3](https://www.gnu.org/licenses/gpl-3.0.en.html)**.

See online [documentation](https://acpype.readthedocs.io/) for more.

### **NO WARRANTY AT ALL**

It was inspired by:

- `amb2gmx.pl` (Eric Sorin, David Mobley and John Chodera)
  and depends on `Antechamber` and `OpenBabel`

- [YASARA Autosmiles](http://www.yasara.org/autosmiles.htm) (Elmar Krieger)

- `topolbuild` (Bruce Ray)

- `xplo2d` (G.J. Kleywegt)

For Non-uniform 1-4 scale factor conversion (e.g. if using **GLYCAM06**), please cite:

> BERNARDI, A., FALLER, R., REITH, D., and KIRSCHNER, K. N. ACPYPE update for
> nonuniform 1–4 scale factors: Conversion of the GLYCAM06 force field from AMBER
> to GROMACS. SoftwareX 10 (2019), 100241. Doi: [10.1016/j.softx.2019.100241](https://doi.org/10.1016/j.softx.2019.100241)

For `Antechamber`, please cite:

> 1. WANG, J., WANG, W., KOLLMAN, P. A., and CASE, D. A. Automatic atom type and
>    bond type perception in molecular mechanical calculations. Journal of Molecular
>    Graphics and Modelling 25, 2 (2006), 247–260. Doi: [10.1016/j.jmgm.2005.12.005](https://doi.org/10.1016/j.jmgm.2005.12.005)

> 2. WANG, J., WOLF, R. M., CALDWELL, J. W., KOLLMAN, P. A., and CASE, D. A.
>    Development and testing of a General Amber Force Field. Journal of Computational
>    Chemistry 25, 9 (2004), 1157–1174. Doi: [10.1002/jcc.20035](https://doi.org/10.1002/jcc.20035)

If you use this code, I am glad if you cite:

> SOUSA DA SILVA, A. W. & VRANKEN, W. F.
> ACPYPE - AnteChamber PYthon Parser interfacE.
> BMC Research Notes 5 (2012), 367 Doi: [10.1186/1756-0500-5-367](https://doi.org/10.1186/1756-0500-5-367)

and (optionally)

> BATISTA, P. R.; WILTER, A.; DURHAM, E. H. A. B. & PASCUTTI, P. G. Molecular
> Dynamics Simulations Applied to the Study of Subtypes of HIV-1 Protease.
> Cell Biochemistry and Biophysics 44 (2006), 395-404. Doi: [10.1385/CBB:44:3:395](https://doi.org/10.1385/CBB:44:3:395)

Alan Silva, DSc

alanwilter _at_ gmail _dot_ com

#### How To Use ACPYPE

##### Introduction

We now have an up-to-date _web service_ at **[Bio2Byte](http://bio2byte.be/acpype/)** (but it **does not** have the `amb2gmx` functionality).

To run `acpype`, locally, with its all functionalities, you need **ANTECHAMBER** from package
[AmberTools](http://ambermd.org/) and
[Open Babel](http://openbabel.org) if your input files are of PDB
format.

However, if one wants `acpype` just to emulate _amb2gmx.pl_, one needs nothing
at all but _[Python](http://www.python.org)_.

There are several ways of obtaining `acpype`:

1. Via **[CONDA](https://anaconda.org/search?q=acpype)**:

   _(It should be wholesome, fully functional, all batteries included)_

   ```bash
   conda install -c conda-forge acpype
   ```

2. Via **[PyPI](https://pypi.org/project/acpype/)**:

   If you're using Linux with Intel processors then

   ```bash
   pip install acpype
   ```

   is enough and you should have a complete solution. Oterwise ...

   _(Make sure you have `AmberTools` and, optionally but highly recommended, `OpenBabel`)_

   ```bash
   # You can use conda to get the needed 3rd parties for example
   conda create -n acpype --channel conda-forge ambertools openbabel

   # Or for Ubuntu 20:
   apt-get install -y openbabel python3-openbabel libarpack++2-dev libgfortran5

   pip install acpype

   # or if you feel daring

   pip install git+https://github.com/alanwilter/acpype.git
   ```

   **NB:** If using OpenBabel python module, it's really **_CRITICAL_** to have it installed in the same `Python` environment of `acpype`.

3. By downloading it via `git`:

   _(Make sure you have `AmberTools` and, optionally but highly recommended, `OpenBabel`)_

   ```bash
   # You can use conda to get the needed 3rd parties for example
   conda create -n acpype --channel conda-forge ambertools openbabel

   # Or for Ubuntu 20:
   apt-get install -y openbabel python3-openbabel libarpack++2-dev libgfortran5

   git clone https://github.com/alanwilter/acpype.git
   ```

   **NB:** Using this mode, CHARMM topology files will not be generated.

4. Via **[Docker](https://hub.docker.com/repository/docker/acpype/acpype/)**:

   _(It should be wholesome, fully functional, all batteries included)_

   If you have Docker installed, you can run `acpype_docker.sh` by:

   NOTE: first time may take some time as it pulls the `acpype` docker image.

   On Linux / macOS:

   ```bash
   ln -fsv "$PWD/acpype_docker.sh" /usr/local/bin/acpype_docker
   ```

   On Windows:
   Using Command Prompt:

   In the directory where the `acpype_docker.bat` file is found:

   ```bash
   setx /M path "%path%;%cd%"
   ```

   Commands:

   ```bash
   acpype_docker -i CCCC

   acpype_docker -i tests/DDD.pdb -c gas
   ```

**NB:**

- By installing via `conda` or using via `docker` you get `AmberTools v.21.11` and `OpenBabel v3.1.1`. Our `AmberTools v.21.11` is a stripped version from the original containing only the necessary binaries and libraries and comes with the `charmmgen` binary from `AmberTools17` in order to generate CHARMM topologies.
- By installing via `pip` you get `AmberTools` (as described above) embedded. However, the included binaries may not work in your system (library dependencies issues) and with only provide binaries for Linux (Ubuntu20) and macOS (Intel).

##### To Test, if doing via `git`

At folder `acpype/`, type:

```bash
./run_acpype.py -i tests/FFF.pdb
```

It'll create a folder called _FFF.acpype_, and inside it one may find topology
files for GROMACS and CNS/XPLOR.

Or using a molecule in [SMILES](https://archive.epa.gov/med/med_archive_03/web/html/smiles.html) notation:

```bash
./run_acpype.py -i CCCC # smiles for C4H6 1,3-Butadiene compound
```

It'll create a folder called _smiles_molecule.acpype_.

To get help and more information, type:

```bash
./run_acpype.py -h
```

##### To Install

At folder `acpype/`, type:

```bash
  ln -fsv "$PWD/run_acpype.py" /usr/local/bin/acpype
```

Then re-login or start another shell session.

If via `conda` or `pip`, `acpype` should be in your `$PATH`.

##### To Verify with GMX

GROMACS < v.5.0

```bash
cd FFF.acpype/
grompp -c FFF_GMX.gro -p FFF_GMX.top -f em.mdp -o em.tpr
mdrun -v -deffnm em
# And if you have VMD
vmd em.gro em.trr
```

GROMACS > v.5.0

```bash
cd FFF.acpype/
gmx grompp -c FFF_GMX.gro -p FFF_GMX.top -f em.mdp -o em.tpr
gmx mdrun -v -deffnm em
# And if you have VMD
vmd em.gro em.trr
```

##### For MD, do

GROMACS < v.5.0

```bash
grompp -c em.gro -p FFF_GMX.top -f md.mdp -o md.tpr
mdrun -v -deffnm md
vmd md.gro md.trr
```

GROMACS > v.5.0

```bash
gmx grompp -c em.gro -p FFF_GMX.top -f md.mdp -o md.tpr
gmx mdrun -v -deffnm md
vmd md.gro md.trr
```

#### To Emulate `amb2gmx.pl`

For any given _prmtop_ and _inpcrd_ files (outputs from AMBER LEaP), type:

```bash
acpype -p FFF_AC.prmtop -x FFF_AC.inpcrd
```

The output files `FFF_GMX.gro` and `FFF_GMX.top` will be generated inside folder _FFF_GMX.amb2gmx_

#### To Verify with CNS/XPLOR

At folder _FFF.acpype_, type:

```bash
cns < FFF_CNS.inp
```

#### To Verify with NAMD

- see [TutorialNAMD](../../wiki/Tutorial-NAMD)
