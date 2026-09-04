import hashlib
import io
import math
import os
import platform
import re
import subprocess as sub
import sys
import tarfile
import tempfile
from shutil import which
from urllib.request import urlopen

from acpype.errors import AcpypeError
from acpype.params import Pi

# charmmgen was dropped by AmberTools, so ACPYPE keeps an old build of its own. The
# wheels bundle it; installs without a bundle, conda's among them, can fetch it with
# `acpype --fetch-charmmgen`. Pinned to a tag and checksummed so the download is
# reproducible: bump the ref and the digests together whenever the binaries change.
CHARMMGEN_REF = "2026.9.3"
CHARMMGEN_URL = "https://raw.githubusercontent.com/alanwilter/acpype/{ref}/charmmgen_{key}.tgz"
CHARMMGEN_SHA256 = {
    "linux": "851408fa0e42fe53e283acda6bd1675e47835284a4134b1250e01a4c1299c93a",
    "macos": "67057e988d7b620f9460fa816000fa5a560d2bc8e5607b2ce38b3d012c442c01",
}


def find_bin(abin):
    return which(abin) or ""


def checkOpenBabelVersion():
    """Check openbabel version"""
    import warnings

    from openbabel import openbabel as obl

    warnings.filterwarnings("ignore")
    return int(obl.OBReleaseVersion().replace(".", ""))


def dotproduct(aa, bb):
    """Scalar product"""
    return sum((a * b) for a, b in zip(aa, bb, strict=False))


def cross_product(a, b):
    """Cross product"""
    c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
    return c


def length(v):
    """Distance between 2 vectors"""
    return math.sqrt(dotproduct(v, v))


def vec_sub(aa, bb):
    """Vector A - B"""
    return [a - b for a, b in zip(aa, bb, strict=False)]


def imprDihAngle(a, b, c, d):
    """Calculate improper dihedral angle"""
    ba = vec_sub(a, b)
    bc = vec_sub(c, b)
    cb = vec_sub(b, c)
    cd = vec_sub(d, c)
    n1 = cross_product(ba, bc)
    n2 = cross_product(cb, cd)
    angle = math.acos(dotproduct(n1, n2) / (length(n1) * length(n2))) * 180 / Pi
    cp = cross_product(n1, n2)
    if dotproduct(cp, bc) < 0:
        angle = -1 * angle
    return angle


def distanceAA(c1, c2):
    """Distance between two atoms"""
    # print c1, c2
    dist2 = (c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[2] - c2[2]) ** 2
    # dist2 = math.sqrt(dist2)
    return dist2


def elapsedTime(seconds, add_s=False, separator=" "):
    """Takes an amount of seconds and turns it into a human-readable amount of time."""
    suffixes = ["y", "w", "d", "h", "m", "s"]
    # the formatted time string to be returned
    atime = []

    # the pieces of time to iterate over (days, hours, minutes, etc)
    # - the first piece in each tuple is the suffix (d, h, w)
    # - the second piece is the length in seconds (a day is 60s * 60m * 24h)
    parts = [
        (suffixes[0], 60 * 60 * 24 * 7 * 52),
        (suffixes[1], 60 * 60 * 24 * 7),
        (suffixes[2], 60 * 60 * 24),
        (suffixes[3], 60 * 60),
        (suffixes[4], 60),
        (suffixes[5], 1),
    ]

    # for each time piece, grab the value and remaining seconds, and add it to
    # the time string
    for suffix, alength in parts:
        value = seconds // alength
        if value > 0:
            seconds = seconds % alength
            atime.append(f"{value!s}{(suffix, (suffix, suffix + 's')[value > 1])[add_s]}")
        if seconds < 1:
            break

    return separator.join(atime)


def splitBlock(dat):
    """Split a amber parm dat file in blocks
    0 = mass, 1 = extra + bond, 2 = angle, 3 = dihedral, 4 = improp, 5 = hbond
    6 = equiv nbon, 7 = nbon, 8 = END, 9 = etc.
    """
    dict_ = {}
    count = 0
    for line in dat:
        line = line.rstrip()
        if count in dict_:
            dict_[count].append(line)
        else:
            dict_[count] = [line]
        if not line:
            count += 1
    return dict_


def parseFrcmod(lista):
    """Parse FRCMOD file"""
    heads = ["MASS", "BOND", "ANGL", "DIHE", "IMPR", "HBON", "NONB"]
    dict_ = {}
    dd = {}
    ahead = None
    for line in lista[1:]:
        line = line.strip()
        if line[:4] == "CMAP":
            # ff19SB's frcmod ends with a CMAP block of %FLAG records and grids that has
            # no place in a parm .dat; it used to be swept into NONB as junk keys.
            ahead = None
            continue
        if line[:4] in heads:
            ahead = line[:4]
            dict_[ahead] = []
            dd = {}
            continue
        if line and ahead:
            key = line.replace(" -", "-").replace("- ", "-").split()[0]
            if key in dd:
                if not dd[key].count(line):
                    dd[key].append(line)
            else:
                dd[key] = [line]
            dict_[ahead] = dd
    # Built as a comprehension rather than popping in place: mutating the dict while
    # iterating it raises RuntimeError, which frcmod.ff14SB never triggered because it
    # has no empty section, but frcmod.ff99SB does.
    return {kk: vv for kk, vv in dict_.items() if vv}


def readAmberLibTypes(libFiles):
    """Read ``(residue, atom name) -> atom type`` from AMBER OFF library files.

    Args:
        libFiles: paths to ``.lib`` files such as ``amino12.lib``.

    Returns:
        dict: the atom type each residue template assigns to each atom name.
    """
    types = {}
    for path in libFiles:
        residue = None
        with open(path) as fh:
            for line in fh:
                if line.startswith("!"):
                    match = re.match(r"!entry\.(\S+)\.unit\.atoms table", line)
                    residue = match.group(1) if match else None
                    continue
                if residue is None:
                    continue
                tokens = line.split()
                if len(tokens) >= 2:
                    types[(residue, tokens[0].strip('"'))] = tokens[1].strip('"')
    return types


RESIDUE_ALIASES = {"HIS": "HIE", "CYS2": "CYX", "ASPH": "ASH", "GLUH": "GLH", "LYSH": "LYS"}


def retypeMol2Atoms(mol2Lines, templates, allowedTypes, aliases=RESIDUE_ALIASES):
    """Apply residue-template atom types to the matching atoms of a mol2.

    An atom is retyped when its residue and atom names match a template and the
    template type is in ``allowedTypes``. The first and last residues try their
    N- and C-terminal templates first. Everything else keeps its current type, so
    non-protein atoms are untouched and hydrogen naming conventions never matter.

    Args:
        mol2Lines: the lines of a Tripos mol2 file.
        templates: ``(residue, atom name) -> type`` as returned by
            :func:`readAmberLibTypes`.
        allowedTypes: the only types this function may assign.
        aliases: residue name synonyms tried when a name is not in the templates.

    Returns:
        tuple: the new lines, and a list of ``(residue, atom, old type, new type)``.
    """
    lines = list(mol2Lines)
    try:
        start = next(i for i, line in enumerate(lines) if line.startswith("@<TRIPOS>ATOM")) + 1
    except StopIteration:
        return lines, []
    end = next((i for i in range(start, len(lines)) if lines[i].startswith("@<TRIPOS>")), len(lines))

    atoms = []
    for i in range(start, end):
        tokens = lines[i].split()
        if len(tokens) >= 8:
            atoms.append((i, tokens))
    if not atoms:
        return lines, []
    substIds = [int(tokens[6]) for _, tokens in atoms]
    first, last = min(substIds), max(substIds)

    changes = []
    for i, tokens in atoms:
        name, current, substId, residue = tokens[1], tokens[5], int(tokens[6]), tokens[7]
        candidates = [residue, aliases.get(residue, residue)]
        if substId == first:
            candidates = ["N" + c for c in candidates] + candidates
        if substId == last:
            candidates = ["C" + c for c in candidates] + candidates
        newType = next((templates[(c, name)] for c in candidates if (c, name) in templates), None)
        if newType is None or newType not in allowedTypes or newType == current:
            continue
        # Replace the type token in place, keeping the line's own column layout.
        parts = re.split(r"(\s+)", lines[i].rstrip("\n"))
        fields = [k for k, part in enumerate(parts) if part and not part.isspace()]
        parts[fields[5]] = newType.ljust(len(current))
        lines[i] = "".join(parts) + "\n"
        changes.append((residue, name, current, newType))
    return lines, changes


_FRCMOD_SECTIONS = (
    ("MASS", "MASS"),
    ("BOND", "BOND"),
    ("ANGL", "ANGLE"),
    ("DIHE", "DIHE"),
    ("IMPR", "IMPROPER"),
    ("NONB", "NONBON"),
)


def _reverseKey(key):
    return "-".join(reversed(key.split("-")))


def mergeFrcmodGaps(amberOnly, withGaff):
    """Keep GAFF-analogue parameters only where the AMBER set has a genuine gap.

    Both arguments are the lines of a parmchk2 frcmod for the same molecule. The first
    comes from a run against the AMBER parameter set alone, so its entries are the
    terms AMBER cannot supply. The second comes from a run against AMBER merged with
    GAFF; it also contains GAFF analogues that parmchk2 prefers over AMBER wildcards
    it already had, which is what must not reach tleap. The result carries the first
    run's keys with the second run's values.

    Args:
        amberOnly: frcmod lines from the AMBER-only parmchk2 run.
        withGaff: frcmod lines from the AMBER + GAFF parmchk2 run.

    Returns:
        list: the lines of the merged frcmod.
    """
    gaps = parseFrcmod(amberOnly)
    filled = parseFrcmod(withGaff)
    out = ["Remark: terms AMBER lacks, filled by parmchk2 from GAFF where it could\n"]
    for section, header in _FRCMOD_SECTIONS:
        out.append(header + "\n")
        available = filled.get(section, {})
        for key, lines in gaps.get(section, {}).items():
            chosen = available.get(key) or available.get(_reverseKey(key)) or lines
            out.extend(line.rstrip("\n") + "\n" for line in chosen)
        out.append("\n")
    return out


def readParmAtomTypes(parmFiles):
    """Collect the atom types declared in the MASS block of parm ``.dat`` and frcmod files.

    A ``.dat`` opens with a title line and then its MASS entries up to the first blank
    line; a frcmod opens with a remark line and a ``MASS`` header before the same block.

    Args:
        parmFiles: paths to the files to read.

    Returns:
        set: every atom type name declared.
    """
    types = set()
    for path in parmFiles:
        with open(path) as fh:
            lines = fh.read().splitlines()
        start = 1
        for index, line in enumerate(lines[:3]):
            if line.strip().upper().startswith("MASS"):
                start = index + 1
                break
        for line in lines[start:]:
            if not line.strip():
                break
            types.add(line.split()[0])
    return types


def unknownMol2Types(mol2Lines, knownTypes):
    """Find the atoms of a mol2 whose type no parameter set defines.

    ``DU`` is antechamber's dummy type for an atom it could not type at all; anything
    else outside ``knownTypes`` is a type the force field being used has no entry for.

    Args:
        mol2Lines: the lines of a Tripos mol2 file.
        knownTypes: atom types the parameter files declare.

    Returns:
        list: ``(atom name, atom type)`` for every offending atom, in file order.
    """
    found = []
    inAtoms = False
    for line in mol2Lines:
        if line.startswith("@<TRIPOS>ATOM"):
            inAtoms = True
            continue
        if line.startswith("@<TRIPOS>"):
            inAtoms = False
            continue
        tokens = line.split()
        if inAtoms and len(tokens) >= 6 and (tokens[5] == "DU" or tokens[5] not in knownTypes):
            found.append((tokens[1], tokens[5]))
    return found


def parmMerge(fdat1, fdat2, frcmod=False):
    """Merge two amber parm dat/frcmod files, caching the result in the temp directory.

    The cache name includes a digest of both inputs. Keying it on the file names alone
    meant an AmberTools upgrade silently reused a merge built from the previous
    release's parameters, since the names never change.

    Args:
        fdat1: first parm dat file.
        fdat2: second parm dat file, or a frcmod file when ``frcmod`` is set.
        frcmod: whether ``fdat2`` is a frcmod rather than a dat file.

    Returns:
        str: path of the merged file.
    """
    name1 = os.path.basename(fdat1).split(".dat")[0]
    if frcmod:
        name2 = os.path.basename(fdat2).split(".")[1]
    else:
        name2 = os.path.basename(fdat2).split(".dat")[0]
    digest = hashlib.sha256()
    for path in (fdat1, fdat2):
        with open(path, "rb") as fh:
            digest.update(fh.read())
    mname = os.path.join(tempfile.gettempdir(), f"{name1}{name2}-{digest.hexdigest()[:12]}.dat")
    if os.path.exists(mname) and os.path.getsize(mname):
        return mname
    mdatFile = open(mname, "w")
    mdat = [f"merged {name1} {name2}"]

    dat1 = splitBlock(open(fdat1).readlines())

    if frcmod:
        dHeads = {
            "MASS": 0,
            "BOND": 1,
            "ANGL": 2,
            "DIHE": 3,
            "IMPR": 4,
            "HBON": 5,
            "NONB": 7,
        }
        dat2 = parseFrcmod(open(fdat2).readlines())  # dict
        for kk in dat2:
            for parEntry in dat2[kk]:
                idFirst = None
                for line in dat1[dHeads[kk]][:]:
                    if line:
                        key = line.replace(" -", "-").replace("- ", "-").split()[0]
                        if key == parEntry:
                            if not idFirst:
                                idFirst = dat1[dHeads[kk]].index(line)
                            dat1[dHeads[kk]].remove(line)
                rev = dat2[kk][parEntry][:]
                rev.reverse()
                if idFirst is None:
                    idFirst = 0
                for ll in rev:
                    if dHeads[kk] in [
                        0,
                        1,
                        7,
                    ]:  # MASS has title in index 0 and so BOND, NONB
                        dat1[dHeads[kk]].insert(idFirst + 1, ll)
                    else:
                        dat1[dHeads[kk]].insert(idFirst, ll)
        dat1[0][0] = mdat[0]
        for kk in dat1:
            for line in dat1[kk]:
                mdatFile.write(line + "\n")
        return mname

    dat2 = splitBlock(open(fdat2).readlines())
    id1 = 0
    id2 = 0
    for kk in list(dat1)[:8]:
        if kk == 0:
            lines = dat1[kk][1:-1] + dat2[kk][1:-1] + [""]
            for line in lines:
                mdat.append(line)
        if kk == 1:
            for i in dat1[kk]:
                if "-" in i:
                    id1 = dat1[kk].index(i)
                    break
            for j in dat2[kk]:
                if "-" in j:
                    id2 = dat2[kk].index(j)
                    break
            l1 = dat1[kk][:id1]
            l2 = dat2[kk][:id2]
            line = ""
            for item in l1 + l2:
                line += item.strip() + " "
            mdat.append(line)
            lines = dat1[kk][id1:-1] + dat2[kk][id2:-1] + [""]
            for line in lines:
                mdat.append(line)
        if kk in [2, 3, 4, 5, 6]:  # angles, p dih, imp dih
            lines = dat1[kk][:-1] + dat2[kk][:-1] + [""]
            for line in lines:
                mdat.append(line)
        if kk == 7:
            lines = dat1[kk][:-1] + dat2[kk][1:-1] + [""]
            for line in lines:
                mdat.append(line)
    for kk in list(dat1)[8:]:
        for line in dat1[kk]:
            mdat.append(line)
    for kk in list(dat2)[9:]:
        for line in dat2[kk]:
            mdat.append(line)
    for line in mdat:
        mdatFile.write(line + "\n")
    mdatFile.close()

    return mname


def job_pids_family(jpid):
    """INTERNAL: Return all job processes (PIDs)"""
    apid = repr(jpid)
    dict_pids = {}
    pids = [apid]
    cmd = f"ps -A -o uid,pid,ppid|grep {os.getuid()}"
    out = _getoutput(cmd).split("\n")  # getoutput("ps -A -o uid,pid,ppid|grep %i" % os.getuid()).split('\n')
    for item in out:
        vec = item.split()
        dict_pids[vec[2]] = vec[1]
    while True:
        try:
            apid = dict_pids[apid]
            pids.append(apid)
        except KeyError:
            break
    return " ".join(pids)


def _getoutput(cmd):
    """To simulate commands.getoutput
    shell=True is necessary despite security issues
    """
    out = sub.Popen(cmd, shell=True, stderr=sub.STDOUT, stdout=sub.PIPE).communicate()[0][:-1]
    return out.decode()


def while_replace(string):
    while "  " in string:
        string = string.replace("  ", " ")
    return string


#: Vendored AmberTools directory name per platform, for the wheels that carry one.
BUNDLE_BY_PLATFORM = {"linux": "amber_linux", "darwin": "amber_macos"}


def bundled_amber_dir():
    """Locate the AmberTools bundled inside this installation, if there is one.

    Only the Linux x86_64 and macOS arm64 wheels carry AmberTools. Installs from the
    source distribution, and platforms with no wheel at all, have none.

    Returns:
        str | None: path of the bundle, or ``None`` when this installation has none.
    """
    name = BUNDLE_BY_PLATFORM.get(sys.platform)
    if name is None:
        return None
    path = os.path.join(os.path.dirname(__file__), name)
    return path if os.path.isdir(os.path.join(path, "bin")) else None


def charmmgen_platform():
    """Name the charmmgen build this machine can run, if there is one.

    Returns:
        str | None: ``"macos"`` (a universal binary, so it serves Apple Silicon and
        Intel alike), ``"linux"`` for x86_64, or ``None`` where no build exists, such
        as Linux on aarch64.
    """
    if sys.platform == "darwin":
        return "macos"
    if sys.platform.startswith("linux") and platform.machine().lower() in ("x86_64", "amd64"):
        return "linux"
    return None


def amber_bin_dir():
    """The AmberTools ``bin`` the antechamber that will actually run belongs to.

    PATH comes first, and the bundle only as a fallback, mirroring
    :func:`set_for_pip`, which wires the bundle in only when nothing is on PATH
    already. Getting that order wrong would check one AmberTools while antechamber
    ran from another, and antechamber resolves charmmgen against its own AMBERHOME.

    Returns:
        str | None: the directory, or ``None`` when no AmberTools can be found.
    """
    antechamber = which("antechamber")
    if antechamber:
        return os.path.dirname(antechamber)
    bundle = bundled_amber_dir()
    return os.path.join(bundle, "bin") if bundle else None


def charmmgen_path():
    """Path of the charmmgen antechamber would actually run, or ``None``.

    antechamber invokes ``$AMBERHOME/bin/charmmgen`` by absolute path, so a copy merely
    on PATH is never found. This looks where antechamber looks.
    """
    directory = amber_bin_dir()
    if directory is None:
        return None
    path = os.path.join(directory, "charmmgen")
    return path if os.access(path, os.X_OK) else None


def download_charmmgen(destDir=None):
    """Download this platform's charmmgen and install it beside antechamber.

    Args:
        destDir: where to write it; the AmberTools ``bin`` in use by default.

    Returns:
        str: path of the installed executable.

    Raises:
        AcpypeError: when no build exists for this platform, no AmberTools can be
            found, the destination is not writable, the download fails, or the archive
            does not match its recorded checksum.
    """
    key = charmmgen_platform()
    if key is None:
        raise AcpypeError(
            f"no charmmgen build for {sys.platform} on {platform.machine()}; "
            "the Docker image ships a full AmberTools and can write CHARMM files"
        )
    directory = destDir or amber_bin_dir()
    if directory is None:
        raise AcpypeError("no AmberTools found to install charmmgen into; is antechamber on PATH?")

    url = CHARMMGEN_URL.format(ref=CHARMMGEN_REF, key=key)
    if not os.access(directory, os.W_OK):
        raise AcpypeError(
            f"'{directory}' is not writable. Install it by hand with:\n"
            f"    curl -sL {url} | tar xz --strip-components=4 -C '{directory}'"
        )
    try:
        with urlopen(url, timeout=60) as response:
            blob = response.read()
    except OSError as exc:
        raise AcpypeError(f"could not download '{url}': {exc}") from exc

    digest = hashlib.sha256(blob).hexdigest()
    if digest != CHARMMGEN_SHA256[key]:
        raise AcpypeError(f"'{url}' does not match its recorded checksum ({digest}); refusing to install it")

    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as archive:
        member = next((m for m in archive.getmembers() if m.name.endswith("bin/charmmgen")), None)
        extracted = archive.extractfile(member) if member else None
        if extracted is None:
            raise AcpypeError(f"'{url}' holds no charmmgen binary")
        binary = extracted.read()

    path = os.path.join(directory, "charmmgen")
    with open(path, "wb") as fh:
        fh.write(binary)
    os.chmod(path, 0o755)  # noqa: S103 - it sits beside antechamber, which is world-executable too
    return path


def set_for_pip(binaries):
    """Point the environment at the bundled AmberTools when one is installed.

    Args:
        binaries: mapping of tool name to executable, as in ``acpype.params.binaries``.
    """
    # For pip package
    if which(binaries["ac_bin"]) is None:
        bundle = bundled_amber_dir()
        if bundle is None:
            # Nothing to wire up: either an sdist install or an unsupported platform.
            # The caller reports this properly once it finds no antechamber.
            return
        os.environ["PATH"] += os.pathsep + os.path.join(bundle, "bin")
        os.environ["AMBERHOME"] = bundle + os.sep
        os.environ["LD_LIBRARY_PATH"] = os.path.join(bundle, "lib") + os.sep
        if sys.platform == "darwin":
            os.environ["DYLD_LIBRARY_PATH"] = os.path.join(bundle, "lib") + os.sep
