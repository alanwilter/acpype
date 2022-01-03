"""
    The Package

    Requirements:
        - ``Python 3.6`` or higher
        - ``Antechamber`` (from ``AmberTools`` preferably)
        - ``OpenBabel`` (optional, but strongly recommended)

    This code is released under **GNU General Public License V3**.

        **<<<  NO WARRANTY AT ALL!!!  >>>**

    It was inspired by:

        - ``amb2gmx.pl`` (Eric Sorin, David Mobley and John Chodera) and depends on Antechamber and Openbabel
        - `YASARA Autosmiles <http://www.yasara.org/autosmiles.htm>`_ (Elmar Krieger)
        - ``topolbuild`` (Bruce Ray)
        - ``xplo2d`` (G.J. Kleywegt)

    For Non-uniform 1-4 scale factor conversion (e.g. if using ``GLYCAM06``), please cite:

        BERNARDI, A., FALLER, R., REITH, D., and KIRSCHNER, K. N. ACPYPE update for
        nonuniform 1-4 scale factors: Conversion of the GLYCAM06 force field from AMBER
        to GROMACS. SoftwareX 10 (2019), 100241.
        doi: `10.1016/j.softx.2019.100241 <https://doi.org/10.1016/j.softx.2019.100241>`_

    For Antechamber, please cite:

        1.  WANG, J., WANG, W., KOLLMAN, P. A., and CASE, D. A. Automatic atom type and
            bond type perception in molecular mechanical calculations. Journal of Molecular
            Graphics and Modelling 25, 2 (2006), 247-260.
            doi: `10.1016/j.jmgm.2005.12.005 <https://doi.org/10.1016/j.jmgm.2005.12.005>`_
        2.  WANG, J., WOLF, R. M., CALDWELL, J. W., KOLLMAN, P. A., and CASE, D. A.
            Development and testing of a General Amber Force Field. Journal of Computational
            Chemistry 25, 9 (2004), 1157-1174.
            doi: `10.1002/jcc.20035 <https://doi.org/10.1002/jcc.20035>`_

    If you use this code, I am glad if you cite:

        SOUSA DA SILVA, A. W. & VRANKEN, W. F.
        ACPYPE - AnteChamber PYthon Parser interfacE.
        BMC Research Notes 5 (2012), 367
        doi: `10.1186/1756-0500-5-367 <http://www.biomedcentral.com/1756-0500/5/367>`_

    and (optionally)

        BATISTA, P. R.; WILTER, A.; DURHAM, E. H. A. B. & PASCUTTI, P. G. Molecular
        Dynamics Simulations Applied to the Study of Subtypes of HIV-1 Protease.
        Cell Biochemistry and Biophysics 44 (2006), 395-404.
        doi: `10.1385/CBB:44:3:395 <https://doi.org/10.1385/CBB:44:3:395>`_

    Alan Silva, D.Sc. <alanwilter _at_ gmail _dot_ com>
"""

# from https://packaging.python.org/guides/single-sourcing-package-version/
# using option 2
# updated automatically via pre-commit git-hook
__version__ = "2022.1.3"
