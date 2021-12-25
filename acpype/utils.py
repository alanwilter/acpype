import math
import os
import subprocess as sub
import sys
from shutil import which

from acpype.params import Pi


def find_bin(abin):
    return which(abin) or ""


def checkOpenBabelVersion():
    "check openbabel version"
    import warnings

    import openbabel as obl

    warnings.filterwarnings("ignore")
    return int(obl.OBReleaseVersion().replace(".", ""))


def dotproduct(aa, bb):
    """scalar product"""
    return sum((a * b) for a, b in zip(aa, bb))


def cross_product(a, b):
    """cross product"""
    c = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    return c


def length(v):
    """distance between 2 vectors"""
    return math.sqrt(dotproduct(v, v))


def vec_sub(aa, bb):
    """vector A - B"""
    return [a - b for a, b in zip(aa, bb)]


def imprDihAngle(a, b, c, d):
    """calculate improper dihedral angle"""
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
    dist2 = (c1[0] - c2[0]) ** 2 + (c1[1] - c2[1]) ** 2 + (c1[0] - c2[0]) ** 2 + (c1[2] - c2[2]) ** 2
    # dist2 = math.sqrt(dist2)
    return dist2


def elapsedTime(seconds, add_s=False, separator=" "):
    """
    Takes an amount of seconds and turns it into a human-readable amount of time.
    """
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
            atime.append(f'{str(value)}{(suffix, (suffix, suffix + "s")[value > 1])[add_s]}')
        if seconds < 1:
            break

    return separator.join(atime)


def splitBlock(dat):
    """split a amber parm dat file in blocks
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
    for line in lista[1:]:
        line = line.strip()
        if line[:4] in heads:
            ahead = line[:4]
            dict_[ahead] = []
            dd = {}
            continue
        if line:
            key = line.replace(" -", "-").replace("- ", "-").split()[0]
            if key in dd:
                if not dd[key].count(line):
                    dd[key].append(line)
            else:
                dd[key] = [line]
            dict_[ahead] = dd
    for kk in dict_:
        if not dict_[kk]:
            dict_.pop(kk)
    return dict_


def parmMerge(fdat1, fdat2, frcmod=False):
    """merge two amber parm dat/frcmod files and save in /tmp"""
    name1 = os.path.basename(fdat1).split(".dat")[0]
    if frcmod:
        name2 = os.path.basename(fdat2).split(".")[1]
    else:
        name2 = os.path.basename(fdat2).split(".dat")[0]
    mname = "/tmp/" + name1 + name2 + ".dat"
    mdatFile = open(mname, "w")
    mdat = [f"merged {name1} {name2}"]

    dat1 = splitBlock(open(fdat1).readlines())

    if frcmod:
        dHeads = {"MASS": 0, "BOND": 1, "ANGL": 2, "DIHE": 3, "IMPR": 4, "HBON": 5, "NONB": 7}
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
                    if dHeads[kk] in [0, 1, 7]:  # MASS has title in index 0 and so BOND, NONB
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
    """
    To simulate commands.getoutput
    shell=True is necessary despite security issues
    """
    out = sub.Popen(cmd, shell=True, stderr=sub.STDOUT, stdout=sub.PIPE).communicate()[0][:-1]
    return out.decode()


def while_replace(string):
    while "  " in string:
        string = string.replace("  ", " ")
    return string


def set_for_pip(binaries):
    # For pip package
    if which(binaries["ac_bin"]) is None:
        LOCAL_PATH = os.path.dirname(__file__)
        if sys.platform == "linux":
            os.environ["PATH"] += os.pathsep + LOCAL_PATH + "/amber21-11_linux/bin"
            os.environ["AMBERHOME"] = LOCAL_PATH + "/amber21-11_linux/"
            os.environ["LD_LIBRARY_PATH"] = LOCAL_PATH + "/amber21-11_linux/lib/"
        elif sys.platform == "darwin":
            os.environ["PATH"] += os.pathsep + LOCAL_PATH + "/amber21-11_os/bin"
            os.environ["AMBERHOME"] = LOCAL_PATH + "/amber21-11_os/"
            os.environ["LD_LIBRARY_PATH"] = LOCAL_PATH + "/amber21-11_os/lib/"
            os.environ["DYLD_LIBRARY_PATH"] = LOCAL_PATH + "/amber21-11_os/lib/"
