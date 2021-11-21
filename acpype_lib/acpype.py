#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Requirements: Python 3.6 or higher
                  Antechamber (from AmberTools preferably)
                  OpenBabel (optional, but strongly recommended)

    This code is released under GNU General Public License V3.

          <<<  NO WARRANTY AT ALL!!!  >>>

    It was inspired by:

    - amb2gmx.pl (Eric Sorin, David Mobley and John Chodera)
      and depends on Antechamber and Openbabel

    - YASARA Autosmiles:
      http://www.yasara.org/autosmiles.htm (Elmar Krieger)

    - topolbuild (Bruce Ray)

    - xplo2d (G.J. Kleywegt)

    For Non-uniform 1-4 scale factor conversion (e.g. if using GLYCAM06), please cite:

    BERNARDI, A., FALLER, R., REITH, D., and KIRSCHNER, K. N. ACPYPE update for
    nonuniform 1-4 scale factors: Conversion of the GLYCAM06 force field from AMBER
    to GROMACS. SoftwareX 10 (2019), 100241. doi: 10.1016/j.softx.2019.100241

    For Antechamber, please cite:
    1. WANG, J., WANG, W., KOLLMAN, P. A., and CASE, D. A. Automatic atom type and
       bond type perception in molecular mechanical calculations. Journal of Molecular
       Graphics and Modelling 25, 2 (2006), 247-260. doi: 10.1016/j.jmgm.2005.12.005
    2. WANG, J., WOLF, R. M., CALDWELL, J. W., KOLLMAN, P. A., and CASE, D. A.
       Development and testing of a General Amber Force Field. Journal of Computational
       Chemistry 25, 9 (2004), 1157-1174. doi: 10.1002/jcc.20035

    If you use this code, I am glad if you cite:

    SOUSA DA SILVA, A. W. & VRANKEN, W. F.
    ACPYPE - AnteChamber PYthon Parser interfacE.
    BMC Research Notes 5 (2012), 367 doi: 10.1186/1756-0500-5-367
    http://www.biomedcentral.com/1756-0500/5/367

    BATISTA, P. R.; WILTER, A.; DURHAM, E. H. A. B. & PASCUTTI, P. G. Molecular
    Dynamics Simulations Applied to the Study of Subtypes of HIV-1 Protease.
    Cell Biochemistry and Biophysics 44 (2006), 395-404. doi: 10.1385/CBB:44:3:395

    Alan Wilter Sousa da Silva, D.Sc.
    Bioinformatician, UniProt, EMBL-EBI
    Hinxton, Cambridge CB10 1SD, UK.
    >>http://www.ebi.ac.uk/~awilter<<

    alanwilter _at_ gmail _dot_ com
"""

import traceback
import signal
import time
import argparse
import math
import os
import pickle
import sys
import subprocess as sub
import abc
import array  # to pacify PyLint
import sysconfig
from datetime import datetime
from shutil import copy2, rmtree, which
from acpype_lib import version

MAXTIME = 3 * 3600
# For pip package
if which("antechamber") is None:
    LOCAL_PATH = sysconfig.get_paths()["purelib"]
    if sys.platform == "linux":
        os.environ["PATH"] += os.pathsep + LOCAL_PATH + "/amber21-11_linux/bin:" + LOCAL_PATH + "/amber21-11_linux/dat/"
        os.environ["AMBERHOME"] = LOCAL_PATH + "/amber21-11_linux/"
        os.environ["ACHOME"] = LOCAL_PATH + "/amber21-11_linux/bin/"
        os.environ["LD_LIBRARY_PATH"] = LOCAL_PATH + "/amber21-11_linux/lib/"
    elif sys.platform == "darwin":
        os.environ["PATH"] += os.pathsep + LOCAL_PATH + "/amber21-11_os/bin:" + LOCAL_PATH + "/amber21-11_os/dat/"
        os.environ["AMBERHOME"] = LOCAL_PATH + "/amber21-11_os/"
        os.environ["ACHOME"] = LOCAL_PATH + "/amber21-11_os/bin/"
        os.environ["LD_LIBRARY_PATH"] = LOCAL_PATH + "/amber21-11_os/lib/"

if sys.version_info < (3, 6):
    raise Exception("Sorry, you need python 3.6 or higher")

year = datetime.today().year
# tag = "2019-09-26T19:44:00UTC"
tag = version

lineHeader = f"""
| ACPYPE: AnteChamber PYthon Parser interfacE v. {tag} (c) {year} AWSdS |
"""
frameLine = (len(lineHeader) - 2) * "="
header = f"{frameLine}{lineHeader}{frameLine}"

#    TODO:
#        Howto Charmm and Amber with NAMD
#        Howto build topology for a modified amino acid
#        CYANA topology files

# List of Topology Formats created by acpype so far:
outTopols = ["gmx", "cns", "charmm"]
qDict = {"mopac": 0, "divcon": 1, "sqm": 2}

# Residues that are not solute, to be avoided when balancing charges in
# amb2gmx mode
ionOrSolResNameList = ["Cl-", "Na+", "K+", "CIO", "Cs+", "IB", "Li+", "MG2", "Rb+", "WAT", "MOH", "NMA"]

# leapAmberFile = 'leaprc.ff99SB'  # 'leaprc.ff10' and 'leaprc.ff99bsc0' has extra Atom Types not in parm99.dat
leapAmberFile = "leaprc.protein.ff14SB"  # 'leaprc.ff14SB'

# "qm_theory='AM1', grms_tol=0.0002, maxcyc=999, tight_p_conv=1, scfconv=1.d-10,"
# "AM1 ANALYT MMOK GEO-OK PRECISE"

cal = 4.184
Pi = 3.141593
qConv = 18.2223
radPi = 57.295780  # 180/Pi
maxDist = 3.0
minDist = 0.5
maxDist2 = maxDist ** 2  # squared Ang.
minDist2 = minDist ** 2  # squared Ang.
diffTol = 0.01

dictAmbAtomType2AmbGmxCode = {
    "BR": "1",
    "C": "2",
    "CA": "3",
    "CB": "4",
    "CC": "5",
    "CK": "6",
    "CM": "7",
    "CN": "8",
    "CQ": "9",
    "CR": "10",
    "CT": "11",
    "CV": "12",
    "CW": "13",
    "C*": "14",
    "Ca": "15",
    "F": "16",
    "H": "17",
    "HC": "18",
    "H1": "19",
    "H2": "20",
    "H3": "21",
    "HA": "22",
    "H4": "23",
    "H5": "24",
    "HO": "25",
    "HS": "26",
    "HW": "27",
    "HP": "28",
    "I": "29",
    "Cl": "30",
    "Na": "31",
    "IB": "32",
    "Mg": "33",
    "N": "34",
    "NA": "35",
    "NB": "36",
    "NC": "37",
    "N2": "38",
    "N3": "39",
    "N*": "40",
    "O": "41",
    "OW": "42",
    "OH": "43",
    "OS": "44",
    "O2": "45",
    "P": "46",
    "S": "47",
    "SH": "48",
    "CU": "49",
    "FE": "50",
    "K": "51",
    "Rb": "52",
    "Cs": "53",
    "Li": "56",
    "Zn": "57",
    "Sr": "58",
    "Ba": "59",
    "MCH3A": "MCH3A",
    "MCH3B": "MCH3B",
    "MNH2": "MNH2",
    "MNH3": "MNH3",
    "MW": "MW",
}

dictOplsAtomType2OplsGmxCode = {
    "Ac3+": ["697"],
    "Am3+": ["699"],
    "Ar": ["Ar", "097"],
    "Ba2+": ["414"],
    "Br": ["722", "730"],
    "Br-": ["402"],
    "CT": [
        "064",
        "076",
        "122",
        "135",
        "136",
        "137",
        "138",
        "139",
        "148",
        "149",
        "152",
        "157",
        "158",
        "159",
        "161",
        "173",
        "174",
        "175",
        "181",
        "182",
        "183",
        "184",
        "206",
        "207",
        "208",
        "209",
        "210",
        "211",
        "212",
        "213",
        "214",
        "215",
        "216",
        "217",
        "218",
        "219",
        "220",
        "223",
        "224",
        "225",
        "229",
        "230",
        "242",
        "243",
        "244",
        "256",
        "257",
        "258",
        "259",
        "273",
        "274",
        "275",
        "276",
        "291",
        "292",
        "293",
        "294",
        "297",
        "305",
        "306",
        "307",
        "308",
        "331",
        "371",
        "373",
        "375",
        "391",
        "396",
        "421",
        "431",
        "443",
        "448",
        "453",
        "455",
        "458",
        "461",
        "468",
        "476",
        "482",
        "484",
        "486",
        "490",
        "491",
        "492",
        "498",
        "499",
        "505",
        "515",
        "516",
        "645",
        "670",
        "671",
        "672",
        "673",
        "674",
        "675",
        "676",
        "677",
        "678",
        "679",
        "680",
        "681",
        "701",
        "725",
        "747",
        "748",
        "755",
        "756",
        "757",
        "758",
        "762",
        "764",
        "765",
        "766",
        "774",
        "775",
        "776",
        "782",
        "783",
        "903",
        "904",
        "905",
        "906",
        "907",
        "908",
        "912",
        "913",
        "914",
        "915",
        "942",
        "943",
        "944",
        "945",
        "951",
        "957",
        "959",
        "960",
        "961",
        "962",
        "963",
        "964",
    ],
    "CA": [
        "053",
        "145",
        "147",
        "166",
        "199",
        "221",
        "228",
        "260",
        "263",
        "266",
        "302",
        "312",
        "315",
        "317",
        "336",
        "351",
        "362",
        "380",
        "457",
        "460",
        "463",
        "472",
        "488",
        "521",
        "522",
        "523",
        "528",
        "532",
        "533",
        "538",
        "539",
        "551",
        "582",
        "590",
        "591",
        "592",
        "593",
        "604",
        "605",
        "606",
        "607",
        "608",
        "609",
        "610",
        "611",
        "612",
        "625",
        "644",
        "647",
        "648",
        "649",
        "650",
        "651",
        "652",
        "714",
        "716",
        "718",
        "720",
        "724",
        "727",
        "729",
        "731",
        "735",
        "736",
        "737",
        "738",
        "739",
        "742",
        "752",
        "768",
        "916",
        "917",
        "918",
    ],
    "C3": [
        "007",
        "010",
        "036",
        "039",
        "063",
        "065",
        "067",
        "068",
        "069",
        "070",
        "080",
        "088",
        "090",
        "092",
        "096",
        "106",
        "107",
        "109",
        "126",
        "132",
        "415",
        "418",
        "425",
        "429",
    ],
    "C": [
        "001",
        "017",
        "026",
        "058",
        "095",
        "131",
        "231",
        "234",
        "235",
        "247",
        "252",
        "267",
        "320",
        "322",
        "334",
        "366",
        "378",
        "470",
        "471",
        "772",
        "952",
    ],
    "C2": [
        "005",
        "009",
        "015",
        "016",
        "019",
        "022",
        "027",
        "028",
        "031",
        "034",
        "037",
        "056",
        "057",
        "061",
        "071",
        "081",
        "089",
        "091",
        "093",
        "110",
    ],
    "CT_2": ["223B", "224B", "225B", "246", "283", "284", "285", "292B", "293B", "295", "298", "299", "906B", "912B"],
    "CM": ["141", "142", "143", "227", "323", "324", "337", "338", "381", "382", "517", "518", "708"],
    "CW": ["508", "514", "543", "552", "561", "567", "575", "583", "588", "637"],
    "CB": ["050", "349", "350", "364", "365", "501", "595", "623", "624"],
    "CH": ["006", "008", "014", "025", "029", "030", "060", "073"],
    "CZ": ["261", "423", "754", "925", "927", "928", "929", "931"],
    "CO": ["189", "191", "193", "195", "197", "198"],
    "C_2": ["232", "233", "277", "280", "465"],
    "CR": ["506", "509", "558", "572", "634"],
    "CQ": ["347", "531", "621", "642"],
    "CV": ["507", "560", "574", "636"],
    "CY": ["711", "712", "713", "733"],
    "CS": ["544", "568", "589"],
    "CK": ["353", "627"],
    "CN": ["502", "594"],
    "CP": ["043", "048"],
    "CU": ["550", "581"],
    "CT_3": ["245", "296"],
    "C=": ["150", "178"],
    "CD": ["011", "075"],
    "C4": ["066"],
    "C7": ["077"],
    "C8": ["074"],
    "C9": ["072"],
    "CX": ["510"],
    "C!": ["145B"],
    "C*": ["500"],
    "C+": ["700"],
    "C_3": ["271"],
    "CC": ["045"],
    "CF": ["044"],
    "CG": ["049"],
    "CT_4": ["160"],
    "Ca2+": ["412"],
    "Cl": ["123", "151", "226", "264"],
    "Cl-": ["401", "709"],
    "Cs+": ["410"],
    "Cu2+": ["Cu2+"],
    "Eu3+": ["705"],
    "F": ["164", "719", "721", "726", "728", "786", "956", "965"],
    "F-": ["400"],
    "Fe2+": ["Fe2+"],
    "Gd3+": ["706"],
    "HA": [
        "146",
        "316",
        "318",
        "389",
        "524",
        "525",
        "526",
        "529",
        "534",
        "535",
        "536",
        "540",
        "541",
        "546",
        "547",
        "554",
        "555",
        "556",
        "563",
        "564",
        "565",
        "569",
        "570",
        "576",
        "577",
        "578",
        "584",
        "585",
        "586",
        "597",
        "598",
        "599",
        "600",
        "601",
        "602",
        "613",
        "614",
        "615",
        "616",
        "617",
        "618",
        "619",
        "629",
        "630",
        "631",
        "638",
        "639",
        "640",
        "643",
        "653",
        "654",
        "655",
        "656",
        "715",
        "717",
        "740",
        "741",
        "746",
    ],
    "HC": [
        "140",
        "144",
        "153",
        "156",
        "165",
        "176",
        "185",
        "190",
        "192",
        "194",
        "196",
        "255",
        "279",
        "282",
        "329",
        "330",
        "332",
        "344",
        "372",
        "374",
        "376",
        "392",
        "416",
        "419",
        "422",
        "426",
        "430",
        "432",
        "444",
        "449",
        "454",
        "456",
        "459",
        "462",
        "469",
        "477",
        "483",
        "485",
        "487",
        "702",
        "710",
        "759",
        "763",
        "777",
        "778",
        "779",
        "784",
        "911",
        "926",
        "930",
        "950",
        "958",
    ],
    "H": [
        "004",
        "013",
        "041",
        "047",
        "128",
        "240",
        "241",
        "250",
        "254",
        "314",
        "325",
        "327",
        "339",
        "342",
        "343",
        "357",
        "358",
        "360",
        "367",
        "369",
        "383",
        "385",
        "387",
        "388",
        "428",
        "479",
        "481",
        "504",
        "513",
        "545",
        "553",
        "562",
        "596",
        "632",
        "744",
        "745",
        "909",
        "910",
    ],
    "H3": ["021", "052", "055", "104", "105", "289", "290", "301", "304", "310", "941", "955"],
    "HO": ["024", "079", "155", "163", "168", "170", "172", "188", "270", "435"],
    "HS": ["033", "086", "087", "204", "205"],
    "HW": ["112", "114", "117", "119", "796"],
    "H4": ["345", "390"],
    "H5": ["355", "359"],
    "He": ["130"],
    "I": ["732"],
    "I-": ["403"],
    "K+": ["408"],
    "Kr": ["098"],
    "LP": ["433", "797"],
    "La3+": ["703"],
    "Li+": ["404", "406"],
    "MCH3A": ["MCH3A"],
    "MCH3B": ["MCH3B"],
    "MNH2": ["MNH2"],
    "MNH3": ["MNH3"],
    "MW": ["MW", "115"],
    "Mg2+": ["411"],
    "NA": [
        "040",
        "046",
        "319",
        "321",
        "333",
        "354",
        "361",
        "377",
        "379",
        "503",
        "512",
        "542",
        "548",
        "557",
        "587",
        "628",
    ],
    "NC": ["311", "335", "346", "348", "363", "520", "527", "530", "537", "603", "620", "622", "641", "646"],
    "N": ["003", "012", "094", "237", "238", "239", "249", "251", "265", "478", "480", "787"],
    "N3": ["020", "101", "102", "103", "286", "287", "288", "309", "427", "940", "953"],
    "N2": ["051", "054", "300", "303", "313", "341", "356", "368", "386", "743"],
    "NB": ["042", "352", "511", "549", "559", "573", "580", "626", "635"],
    "N*": ["319B", "333B", "354B", "377B"],
    "NT": ["127", "900", "901", "902"],
    "NZ": ["262", "424", "750", "753"],
    "NO": ["760", "767"],
    "NY": ["749", "751"],
    "Na+": ["405", "407"],
    "Nd3+": ["704"],
    "Ne": ["129"],
    "OS": ["062", "108", "179", "180", "186", "395", "442", "447", "452", "467", "473", "566", "571", "579", "773"],
    "O": ["002", "059", "236", "248", "253", "326", "328", "340", "370", "384", "771", "788"],
    "OH": ["023", "078", "154", "162", "167", "169", "171", "187", "268", "420", "434"],
    "O2": ["018", "125", "272", "394", "441", "446", "451", "954"],
    "OW": ["111", "113", "116", "118", "795"],
    "O_2": ["278", "281", "466"],
    "OY": ["475", "494", "497"],
    "OL": ["120"],
    "ON": ["761"],
    "OU": ["437"],
    "O_3": ["269"],
    "P": ["393", "440", "445", "450", "785"],
    "P+": ["781"],
    "Rb+": ["409"],
    "S": ["035", "038", "084", "085", "124", "202", "203", "222", "633"],
    "SH": ["032", "082", "083", "200", "201", "417", "734"],
    "SI": ["SI"],
    "SY": ["474"],
    "SY2": ["493"],
    "SZ": ["496"],
    "Sr2+": ["413"],
    "Th4+": ["698"],
    "U": ["436"],
    "Xe": ["099"],
    "Yb3+": ["707"],
    "Zn2+": ["Zn2+"],
}

# reverse dictOplsAtomType2OplsGmxCode
oplsCode2AtomTypeDict = {}
for k, vv in list(dictOplsAtomType2OplsGmxCode.items()):
    for code in vv:
        oplsCode2AtomTypeDict[code] = k

# Cross dictAmbAtomType2AmbGmxCode with dictOplsAtomType2OplsGmxCode & add H1,HP,H2
dictAtomTypeAmb2OplsGmxCode = {"H1": ["140", "1.00800"], "HP": ["140", "1.00800"], "H2": ["140", "1.00800"]}
dictOplsMass = {
    "SY2": ["32.06000"],
    "Zn2+": ["65.37000"],
    "CQ": ["12.01100"],
    "CP": ["12.01100"],
    "Nd3+": ["144.24000"],
    "Br-": ["79.90400"],
    "Cu2+": ["63.54600"],
    "Br": ["79.90400"],
    "H": ["1.00800"],
    "P": ["30.97376"],
    "Sr2+": ["87.62000"],
    "ON": ["15.99940"],
    "OL": ["0.00000"],
    "OH": ["15.99940"],
    "OY": ["15.99940"],
    "OW": ["15.99940"],
    "OU": ["15.99940"],
    "OS": ["15.99940"],
    "Am3+": ["243.06000"],
    "HS": ["1.00800"],
    "HW": ["1.00800"],
    "HO": ["1.00800"],
    "HC": ["1.00800"],
    "HA": ["1.00800"],
    "O2": ["15.99940"],
    "Ca2+": ["40.08000"],
    "Th4+": ["232.04000"],
    "He": ["4.00260"],
    "C": ["12.01100"],
    "Cs+": ["132.90540"],
    "O": ["15.99940"],
    "Gd3+": ["157.25000"],
    "S": ["32.06000"],
    "P+": ["30.97376"],
    "La3+": ["138.91000"],
    "H3": ["1.00800"],
    "H4": ["1.00800"],
    "MNH2": ["0.00000"],
    "MW": ["0.00000"],
    "NB": ["14.00670"],
    "K+": ["39.09830"],
    "Ne": ["20.17970"],
    "Rb+": ["85.46780"],
    "C+": ["12.01100"],
    "C*": ["12.01100"],
    "NO": ["14.00670"],
    "CT_4": ["12.01100"],
    "NA": ["14.00670"],
    "C!": ["12.01100"],
    "NC": ["14.00670"],
    "NZ": ["14.00670"],
    "CT_2": ["12.01100"],
    "CT_3": ["12.01100"],
    "NY": ["14.00670"],
    "C9": ["14.02700"],
    "C8": ["13.01900"],
    "C=": ["12.01100"],
    "Yb3+": ["173.04000"],
    "C3": ["15.03500", "12.01100"],
    "C2": ["14.02700"],
    "C7": ["12.01100"],
    "C4": ["16.04300"],
    "CK": ["12.01100"],
    "Cl-": ["35.45300"],
    "N*": ["14.00670"],
    "CH": ["13.01900"],
    "CO": ["12.01100"],
    "CN": ["12.01100"],
    "CM": ["12.01100"],
    "F": ["18.99840"],
    "CC": ["12.01100"],
    "CB": ["12.01100"],
    "CA": ["12.01100"],
    "CG": ["12.01100"],
    "CF": ["12.01100"],
    "N": ["14.00670"],
    "CZ": ["12.01100"],
    "CY": ["12.01100"],
    "CX": ["12.01100"],
    "Ac3+": ["227.03000"],
    "CS": ["12.01100"],
    "CR": ["12.01100"],
    "N2": ["14.00670"],
    "N3": ["14.00670"],
    "CW": ["12.01100"],
    "CV": ["12.01100"],
    "CU": ["12.01100"],
    "CT": ["12.01100"],
    "SZ": ["32.06000"],
    "SY": ["32.06000"],
    "Cl": ["35.45300"],
    "NT": ["14.00670"],
    "O_2": ["15.99940"],
    "Xe": ["131.29300"],
    "SI": ["28.08000"],
    "SH": ["32.06000"],
    "Eu3+": ["151.96000"],
    "F-": ["18.99840"],
    "MNH3": ["0.00000"],
    "H5": ["1.00800"],
    "C_3": ["12.01100"],
    "C_2": ["12.01100"],
    "I-": ["126.90450"],
    "LP": ["0.00000"],
    "I": ["126.90450"],
    "Na+": ["22.98977"],
    "Li+": ["6.94100"],
    "U": ["0.00000"],
    "MCH3A": ["0.00000"],
    "MCH3B": ["0.00000"],
    "CD": ["13.01900", "12.01100"],
    "O_3": ["15.99940"],
    "Kr": ["83.79800"],
    "Fe2+": ["55.84700"],
    "Ar": ["39.94800"],
    "Mg2+": ["24.30500"],
    "Ba2+": ["137.33000"],
}
for ambKey in dictAmbAtomType2AmbGmxCode:
    if ambKey in dictOplsAtomType2OplsGmxCode:
        dictAtomTypeAmb2OplsGmxCode[ambKey] = dictOplsAtomType2OplsGmxCode[ambKey] + list(dictOplsMass[ambKey])

# learnt from 22 residues test.
dictAtomTypeAmb2OplsGmxCode = {
    "HS": ["204", "1.008"],
    "HP": ["140", "1.008"],
    "HO": ["155", "168", "1.008"],
    "HC": ["140", "1.008"],
    "HA": ["146", "1.008"],
    "O2": ["272", "15.9994"],
    "C*": ["500", "12.011"],
    "NA": ["503", "512", "14.0067"],
    "NB": ["511", "14.0067"],
    "CB": ["501", "12.011"],
    "C": ["235", "271", "12.011"],
    "CN": ["502", "12.011"],
    "CM": ["302", "12.011"],
    "CC": ["507", "508", "510", "12.011"],
    "H": ["240", "241", "290", "301", "304", "310", "504", "513", "1.008"],
    "CA": ["145", "166", "12.011"],
    "O": ["236", "15.9994"],
    "N": ["237", "238", "239", "14.0067"],
    "S": ["202", "32.06"],
    "CR": ["506", "509", "12.011"],
    "N2": ["300", "303", "14.0067"],
    "N3": ["287", "309", "14.0067"],
    "CW": ["508", "510", "514", "12.011"],
    "CV": ["507", "12.011"],
    "CT": [
        "135",
        "136",
        "137",
        "149",
        "157",
        "158",
        "206",
        "209",
        "210",
        "223B",
        "224B",
        "245",
        "246",
        "274",
        "283",
        "284",
        "285",
        "292",
        "292B",
        "293B",
        "296",
        "307",
        "308",
        "505",
        "12.011",
    ],
    "OH": ["154", "167", "15.9994"],
    "H1": ["140", "1.008"],
    "H4": ["146", "1.008"],
    "H5": ["146", "1.008"],
    "SH": ["200", "32.06"],
}

# learnt from 22 residues test.
dictAtomTypeGaff2OplsGmxCode = {
    "cc": ["500", "506", "507", "508", "514", "12.011"],
    "ca": ["145", "166", "501", "502", "12.011"],
    "h1": ["140", "1.008"],
    "h4": ["146", "1.008"],
    "h5": ["146", "1.008"],
    "cz": ["302", "12.011"],
    "c2": ["509", "510", "12.011"],
    "nh": ["300", "303", "14.0067"],
    "ha": ["146", "1.008"],
    "na": ["503", "512", "14.0067"],
    "nc": ["511", "14.0067"],
    "nd": ["511", "14.0067"],
    "hx": ["140", "1.008"],
    "hs": ["204", "1.008"],
    "hn": ["240", "241", "290", "301", "304", "310", "504", "513", "1.008"],
    "ho": ["155", "168", "1.008"],
    "c3": [
        "135",
        "136",
        "137",
        "149",
        "157",
        "158",
        "206",
        "209",
        "210",
        "223B",
        "224B",
        "245",
        "246",
        "274",
        "283",
        "284",
        "285",
        "292",
        "292B",
        "293B",
        "296",
        "307",
        "308",
        "505",
        "12.011",
    ],
    "hc": ["140", "1.008"],
    "cd": ["500", "506", "507", "508", "514", "12.011"],
    "c": ["235", "271", "12.011"],
    "oh": ["154", "167", "15.9994"],
    "ss": ["202", "32.06"],
    "o": ["236", "272", "15.9994"],
    "n": ["237", "238", "239", "14.0067"],
    "sh": ["200", "32.06"],
    "n4": ["287", "309", "14.0067"],
}

# draft
atomTypeAmber2oplsDict = {
    "HS": ["HS"],
    "HP": ["HC"],
    "HO": ["HO"],
    "HC": ["HC"],
    "HA": ["HA"],
    "O2": ["O2"],
    "C*": ["C*"],
    "NA": ["NA"],
    "NB": ["NB"],
    "CB": ["CB"],
    "CN": ["CN"],
    "CV": ["CV"],
    "CM": ["CA"],
    "CA": ["CA"],
    "CR": ["CR"],
    "OH": ["OH"],
    "H1": ["HC"],
    "H4": ["HA"],
    "N2": ["N2"],
    "N3": ["N3"],
    "H5": ["HA"],
    "SH": ["SH"],
    "N": ["N"],
    "S": ["S"],
    "O": ["O"],
    "C": ["C", "C_3"],
    "CW": ["CW", "CX"],
    "H": ["H", "H3"],
    "CC": ["CX", "CW", "CV"],
    "CT": ["CT", "CT_2", "CT_3"],
}

# draft
a2oD = {
    "amber99_2": ["opls_235", "opls_271"],
    "amber99_3": ["opls_302", "opls_145"],
    "amber99_5": ["opls_507", "opls_508", "opls_510"],
    "amber99_11": [
        "opls_209",
        "opls_158",
        "opls_283",
        "opls_223B",
        "opls_293B",
        "opls_284",
        "opls_292B",
        "opls_274",
        "opls_136",
        "opls_135",
        "opls_292",
        "opls_157",
        "opls_206",
        "opls_137",
        "opls_505",
        "opls_224B",
        "opls_307",
        "opls_308",
        "opls_210",
        "opls_149",
    ],
    "amber99_13": ["opls_514"],
    "amber99_14": ["opls_500"],
    "amber99_17": ["opls_504", "opls_241", "opls_240", "opls_290", "opls_301", "opls_310", "opls_304", "opls_513"],
    "amber99_18": ["opls_140"],
    "amber99_19": ["opls_140"],
    "amber99_22": ["opls_146"],
    "amber99_23": ["opls_146"],
    "amber99_25": ["opls_155"],
    "amber99_26": ["opls_204"],
    "amber99_28": ["opls_140"],
    "amber99_34": ["opls_238", "opls_239", "opls_237"],
    "amber99_35": ["opls_512", "opls_503"],
    "amber99_36": ["opls_511"],
    "amber99_38": ["opls_300", "opls_303"],
    "amber99_39": ["opls_309", "opls_287"],
    "amber99_41": ["opls_236"],
    "amber99_43": ["opls_154"],
    "amber99_45": ["opls_272"],
    "amber99_47": ["opls_202"],
    "amber99_48": ["opls_200"],
}

# pid = 0 # not sure if I still need it.

head = "%s created by acpype (v: " + tag + ") on %s\n"

date = datetime.now().ctime()

usage = """
    acpype -i _file_ | _SMILES_string_ [-c _string_] [-n _int_] [-m _int_] [-a _string_] [-f] etc. or
    acpype -p _prmtop_ -x _inpcrd_ [-d | -w]"""

epilog = """

    output: assuming 'root' is the basename of either the top input file,
            the 3-letter residue name or user defined (-b option)
    root_bcc_gaff.mol2:  final mol2 file with 'bcc' charges and 'gaff' atom type
    root_AC.inpcrd    :  coord file for AMBER
    root_AC.prmtop    :  topology and parameter file for AMBER
    root_AC.lib       :  residue library file for AMBER
    root_AC.frcmod    :  modified force field parameters
    root_GMX.gro      :  coord file for GROMACS
    root_GMX.top      :  topology file for GROMACS
    root_GMX.itp      :  molecule unit topology and parameter file for GROMACS
    root_GMX_OPLS.itp :  OPLS/AA mol unit topol & par file for GROMACS (experimental!)
    em.mdp, md.mdp    :  run parameters file for GROMACS
    root_NEW.pdb      :  final pdb file generated by ACPYPE
    root_CNS.top      :  topology file for CNS/XPLOR
    root_CNS.par      :  parameter file for CNS/XPLOR
    root_CNS.inp      :  run parameters file for CNS/XPLOR
    root_CHARMM.rtf   :  topology file for CHARMM
    root_CHARMM.prm   :  parameter file for CHARMM
    root_CHARMM.inp   :  run parameters file for CHARMM"""

TLEAP_TEMPLATE = """
verbosity 1
source %(leapAmberFile)s
source %(leapGaffFile)s
mods = loadamberparams %(acBase)s.frcmod
%(res)s = loadmol2 %(acMol2FileName)s
check %(res)s
saveamberparm %(res)s %(acBase)s.prmtop %(acBase)s.inpcrd
saveoff %(res)s %(acBase)s.lib
quit
"""


def checkOpenBabelVersion():
    "check openbabel version"
    import openbabel as obl
    import warnings

    warnings.filterwarnings("ignore")
    return int(obl.OBReleaseVersion().replace(".", ""))


def checkSmiles(smiles):

    if checkOpenBabelVersion() >= 300:
        from openbabel import openbabel as ob
        from openbabel import pybel

        ob.cvar.obErrorLog.StopLogging()
    elif checkOpenBabelVersion() >= 200 and checkOpenBabelVersion() < 300:
        import openbabel as ob
        import pybel  # type: ignore

        ob.cvar.obErrorLog.StopLogging()

    " Check if input is a smiles string "

    try:
        ob.obErrorLog.SetOutputLevel(0)
        pybel.readstring("smi", smiles)
        return True
    except Exception:
        ob.obErrorLog.SetOutputLevel(0)
        return False


def dotproduct(aa, bb):
    """scalar product"""
    return sum((a * b) for a, b in zip(aa, bb))


def crosproduct(a, b):
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
    n1 = crosproduct(ba, bc)
    n2 = crosproduct(cb, cd)
    angle = math.acos(dotproduct(n1, n2) / (length(n1) * length(n2))) * 180 / Pi
    cp = crosproduct(n1, n2)
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
    """to simulate commands.getoutput in order to work with python 2.6 up to 3.x"""
    out = sub.Popen(cmd, shell=True, stderr=sub.STDOUT, stdout=sub.PIPE).communicate()[0][:-1]
    return out.decode()


def while_replace(string):
    while "  " in string:
        string = string.replace("  ", " ")
    return string


class Topology_14:

    """
    Amber topology abstraction for non-uniform 1-4 scale factors

    """

    def __init__(self) -> None:
        self.pointers = array.array("d")
        self.charge = array.array("d")
        self.atom_type_index = array.array("d")
        self.nonbonded_parm_index = array.array("d")
        self.scee_scale_factor = array.array("d")
        self.scnb_scale_factor = array.array("d")
        self.dihedral_force_constants = array.array("d")
        self.dihedral_periodicity = array.array("d")
        self.dihedral_phase = array.array("d")
        self.dihedral_yes_H = array.array("d")
        self.dihedral_no_H = array.array("d")
        self.lennard_jones_acoef = array.array("d")
        self.lennard_jones_bcoef = array.array("d")

    def read_amber_topology(self, buff):
        """Read AMBER topology file"""
        flag_strings = [
            "%FLAG POINTERS",
            "%FLAG CHARGE",
            "%FLAG ATOM_TYPE_INDEX",
            "%FLAG NONBONDED_PARM_INDEX",
            "%FLAG SCEE_SCALE_FACTOR",
            "%FLAG SCNB_SCALE_FACTOR",
            "%FLAG DIHEDRAL_FORCE_CONSTANT",
            "%FLAG DIHEDRAL_PERIODICITY",
            "%FLAG DIHEDRAL_PHASE",
            "%FLAG DIHEDRALS_INC_HYDROGEN",
            "%FLAG DIHEDRALS_WITHOUT_HYDROGEN",
            "%FLAG LENNARD_JONES_ACOEF",
            "%FLAG LENNARD_JONES_BCOEF",
        ]
        attributes = [
            "pointers",
            "charge",
            "atom_type_index",
            "nonbonded_parm_index",
            "scee_scale_factor",
            "scnb_scale_factor",
            "dihedral_force_constants",
            "dihedral_periodicity",
            "dihedral_phase",
            "dihedral_yes_H",
            "dihedral_no_H",
            "lennard_jones_acoef",
            "lennard_jones_bcoef",
        ]
        for i, _item in enumerate(attributes):
            try:
                setattr(self, attributes[i], self.p7_array_read(buff, flag_strings[i]))
            except Exception:
                print("Skipping non-existent attributes", attributes[i], flag_strings[i])

    @staticmethod
    def skipline(buff, index):
        """skip line"""
        while buff[index] != "\n":
            index += 1
        index += 1
        return index

    def p7_array_read(self, buff, flag_string):
        """Convert AMBER topology data to python array"""
        myarray = array.array("d")
        i = buff.index(flag_string)
        i = self.skipline(buff, i)
        i = self.skipline(buff, i)
        while 1:
            while buff[i] == " " or buff[i] == "\t" or buff[i] == "\n":
                i += 1
            j = i
            if buff[i] == "%":
                break
            while buff[i] != " " and buff[i] != "\t" and buff[i] != "\n":
                i += 1
            myarray.append(float(buff[j:i]))
        return myarray

    def print_gmx_pairs(self):
        """Generate non-bonded pairs list"""
        pair_list = []
        pair_buff = "[ pairs_nb ]\n;       ai         aj  funct         qi         qj           sigma         epsilon\n"
        pair_list.append(pair_buff)
        dihedrals = self.dihedral_yes_H + self.dihedral_no_H
        dih_number = len(dihedrals)
        j = int(0)
        while j < dih_number:
            if dihedrals[j + 2] > 0:
                parm_idx = int(dihedrals[j + 4]) - 1
                scee_scale_factor = self.scee_scale_factor[parm_idx]
                if scee_scale_factor == 0:
                    scee_scale_factor = 1.2
                ai = int(abs(dihedrals[j]) / 3)
                al = int(abs(dihedrals[j + 3]) / 3)
                qi = self.charge[ai] / qConv
                ql = self.charge[al] / qConv / scee_scale_factor
                ntypes = int(self.pointers[1])
                ai_index = int(self.atom_type_index[ai])
                al_index = int(self.atom_type_index[al])
                nb_parm_index = int(self.nonbonded_parm_index[ntypes * (ai_index - 1) + al_index - 1]) - 1
                scnb_scale_factor = self.scnb_scale_factor[parm_idx]
                if scnb_scale_factor == 0:
                    scnb_scale_factor = 2
                lj_acoeff = self.lennard_jones_acoef[nb_parm_index] / scnb_scale_factor
                lj_bcoeff = self.lennard_jones_bcoef[nb_parm_index] / scnb_scale_factor
                if lj_bcoeff != 0:
                    sigma6 = lj_acoeff / lj_bcoeff
                else:
                    sigma6 = 1  # arbitrary and doesnt matter
                epsilon = lj_bcoeff / 4 / sigma6 * 4.184
                sigma = sigma6 ** (1 / 6) / 10
                pair_buff = (
                    "{:>10.0f} {:>10.0f} {:>6.0f} ".format(ai + 1, al + 1, 1)
                    + "{:>10.6f} {:>10.6f} ".format(qi, ql)
                    + "{:>15.5e} {:>15.5e}\n".format(sigma, epsilon)
                )
                pair_list.append(pair_buff)
            j += 5
        return "".join(pair_list)

    def hasNondefault14(self):
        """Check non-uniform 1-4 scale factor"""
        for val in self.scee_scale_factor:
            if val not in (0, 1.2):
                return True
        for val in self.scnb_scale_factor:
            if val not in (0, 2):
                return True
        return False

    def patch_gmx_topol14(self, gmx_init_top):
        """Patch GMX topology file for non-uniform 1-4 scale factor"""
        pair_buff = self.print_gmx_pairs()
        jdefault = gmx_init_top.index("\n[ atomtypes ]")
        ipair = gmx_init_top.index("[ pairs ]")
        jpair = gmx_init_top.index("\n[ angles ]")
        init_buff = (
            "\n\n[ defaults ]\n"
            + "; nbfunc        comb-rule       gen-pairs       \n"
            + "1               2               no              \n"
        )
        return (
            gmx_init_top.splitlines()[0]
            + init_buff
            + gmx_init_top[jdefault:ipair]
            + pair_buff
            + gmx_init_top[jpair : len(gmx_init_top)]
        )


class AbstractTopol:

    """
    Super class to build topologies
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self):
        if self.__class__ is AbstractTopol:
            raise TypeError("Attempt to create istance of abstract class AbstractTopol")
        self.debug = None
        self.verbose = None
        self.chargeVal = None
        self.tmpDir = None
        self.absInputFile = None
        self.chargeType = None
        self.babelExe = None
        self.baseName = None
        self.acExe = None
        self.force = None
        self.acBaseName = None
        self.atomType = None
        self.acMol2FileName = None
        self.multiplicity = None
        self.qFlag = None
        self.ekFlag = None
        self.timeTol = None
        self.acXyzFileName = None
        self.acTopFileName = None
        self.acParDict = None
        self.tleapExe = None
        self.parmchkExe = None
        self.acFrcmodFileName = None
        self.gaffDatfile = None
        self.homeDir = None
        self.rootDir = None
        self.extOld = None
        self.direct = None
        self.disam = None
        self.gmx4 = None
        self.sorted = None
        self.chiral = None
        self.outTopols = None
        self.ext = None
        self.xyzFileData = None
        self.charmmBase = None
        self.allhdg = None
        self.topo14Data = None
        self.atomPairs = None
        self.properDihedralsGmx45 = None
        self.properDihedralsAlphaGamma = None
        self.properDihedralsCoefRB = None
        self.resName = None
        self.acLog = None
        self.tleapLog = None
        self.parmchkLog = None
        self.inputFile = None
        self.babelLog = None
        self.absHomeDir = None
        self.molTopol = None
        self.topFileData = None
        self.residueLabel = None
        self._atomTypeNameList = None
        self.atomTypeSystem = None
        self.totalCharge = None
        self.atoms = None
        self.atomTypes = None
        self.pbc = None
        self.bonds = None
        self.angles = None
        self.properDihedrals = None
        self.improperDihedrals = None
        self.condensedProperDihedrals = None
        self.chiralGroups = None
        self.excludedAtoms = None
        self.atomsGromacs = None
        self.atomTypesGromacs = None
        self.CnsTopFileName = None
        self.CnsInpFileName = None
        self.CnsParFileName = None
        self.CnsPdbFileName = None
        self.is_smiles = None
        self.smiles = None
        self.amb2gmx = None

    def printDebug(self, text=""):
        """Debug log level"""
        if self.debug:
            print(f"DEBUG: {while_replace(text)}")

    def printWarn(self, text=""):
        """Warn log level"""
        if self.verbose:
            print(f"WARNING: {while_replace(text)}")

    def printError(self, text=""):
        """Error log level"""
        if self.verbose:
            print(f"ERROR: {while_replace(text)}")

    def printMess(self, text=""):
        """Info log level"""
        if self.verbose:
            print(f"==> {while_replace(text)}")

    def printQuoted(self, text=""):
        """Print quoted messages"""
        if self.verbose:
            print(10 * "+" + "start_quote" + 59 * "+")
            print(while_replace(text))
            print(10 * "+" + "end_quote" + 61 * "+")

    def search(self, name=None, alist=False):
        """
        returns a list with all atomName matching 'name'
        or just the first case
        """
        ll = [x for x in self.atoms if x.atomName == name.upper()]
        if ll and not alist:
            ll = ll[0]
        return ll

    def guessCharge(self):
        """
        Guess the charge of a system based on antechamber
        Returns None in case of error
        """
        done = False
        error = False
        charge = self.chargeVal
        localDir = os.path.abspath(".")
        if not os.path.exists(self.tmpDir):
            os.mkdir(self.tmpDir)
        if not os.path.exists(os.path.join(self.tmpDir, self.inputFile)):
            copy2(self.absInputFile, self.tmpDir)
        os.chdir(self.tmpDir)

        if self.chargeType == "user":
            if self.ext == ".mol2":
                self.printMess("Reading user's charges from mol2 file...")
                charge = self.readMol2TotalCharge(self.inputFile)
                done = True
            else:
                self.chargeType = "bcc"
                self.printWarn("cannot read charges from a PDB file")
                self.printWarn("using now 'bcc' method for charge")
        if self.chargeVal is None and not done:
            self.printWarn("no charge value given, trying to guess one...")
            mol2FileForGuessCharge = self.inputFile
            if self.ext == ".pdb":
                cmd = f"{self.babelExe} -ipdb {self.inputFile} -omol2 -O {self.baseName}.mol2"
                self.printDebug(f"guessCharge: {cmd}")
                out = _getoutput(cmd)
                self.printDebug(out)
                mol2FileForGuessCharge = os.path.abspath(f"{self.baseName}.mol2")
                in_mol = "mol2"
            else:
                in_mol = self.ext[1:]
                if in_mol == "mol":
                    in_mol = "mdl"

            cmd = f"{self.acExe} -dr no -i {mol2FileForGuessCharge} -fi {in_mol} -o tmp -fo mol2 -c gas -pf y"

            if self.debug:
                self.printMess("Debugging...")
                cmd = cmd.replace("-pf y", "-pf n")
                print(while_replace(cmd))

            log = _getoutput(cmd).strip()

            if os.path.exists("tmp"):
                charge = self.readMol2TotalCharge("tmp")
            else:
                try:
                    charge = float(
                        log.strip()
                        .split("equal to the total charge (")[-1]
                        .split(") based on Gasteiger atom type, exit")[0]
                    )
                except Exception:
                    error = True

            if error:
                self.printError("guessCharge failed")
                os.chdir(localDir)
                self.printQuoted(log)
                self.printMess("Trying with net charge = 0 ...")
        charge = float(charge)
        charge2 = int(round(charge))
        drift = abs(charge2 - charge)
        self.printDebug("Net charge drift '%3.6f'" % drift)
        if drift > diffTol:
            self.printError("Net charge drift '%3.5f' bigger than tolerance '%3.5f'" % (drift, diffTol))
            if not self.force:
                raise Exception("Error with calculated charge")
        self.chargeVal = str(charge2)
        self.printMess(f"... charge set to {charge2}")
        os.chdir(localDir)

    def setResNameCheckCoords(self):
        """Set a 3 letter residue name
        and check coords duplication
        """
        exit_ = False
        localDir = os.path.abspath(".")
        if not os.path.exists(self.tmpDir):
            os.mkdir(self.tmpDir)

        copy2(self.absInputFile, self.tmpDir)
        os.chdir(self.tmpDir)

        exten = self.ext[1:]
        if self.ext == ".pdb":
            tmpFile = open(self.inputFile, "r")
        else:
            if exten == "mol":
                exten = "mdl"
            cmd = f"{self.acExe} -dr no -i {self.inputFile} -fi {exten} -o tmp -fo ac -pf y"
            self.printDebug(cmd)
            out = _getoutput(cmd)
            if not out.isspace():
                self.printDebug(out)
            try:
                tmpFile = open("tmp", "r")
            except Exception:
                rmtree(self.tmpDir)
                raise

        tmpData = tmpFile.readlines()
        residues = set()
        coords = {}
        for line in tmpData:
            if "ATOM  " in line or "HETATM" in line:
                residues.add(line[17:20])
                at = line[0:17]
                cs = line[30:54]
                if cs in coords:
                    coords[cs].append(at)
                else:
                    coords[cs] = [at]
        # self.printDebug(coords)

        if len(residues) > 1:
            self.printError(f"more than one residue detected '{str(residues)}'")
            self.printError(f"verify your input file '{self.inputFile}'. Aborting ...")
            raise Exception("Only ONE Residue is allowed for ACPYPE to work")

        dups = ""
        shortd = ""
        longd = ""
        longSet = set()
        id_ = 0
        items = list(coords.items())
        ll = len(items)
        for item in items:
            id_ += 1
            if len(item[1]) > 1:  # if True means atoms with same coordinates
                for i in item[1]:
                    dups += f"{i} {item[0]}\n"
            for id2 in range(id_, ll):
                item2 = items[id2]
                c1 = list(map(float, [item[0][i : i + 8] for i in range(0, 24, 8)]))
                c2 = list(map(float, [item2[0][i : i + 8] for i in range(0, 24, 8)]))
                dist2 = distanceAA(c1, c2)
                if dist2 < minDist2:
                    dist = math.sqrt(dist2)
                    shortd += f"{dist:8.5f}       {item[1]} {item2[1]}\n"
                if dist2 < maxDist2:  # and not longOK:
                    longSet.add(str(item[1]))
                    longSet.add(str(item2[1]))
            if str(item[1]) not in longSet and ll > 1:
                longd += f"{item[1]}\n"

        if dups:
            self.printError(f"Atoms with same coordinates in '{self.inputFile}'!")
            self.printQuoted(dups[:-1])
            exit_ = True

        if shortd:
            self.printError(f"Atoms TOO close (< {minDist} Ang.)")
            self.printQuoted(f"Dist (Ang.)    Atoms\n{shortd[:-1]}")
            exit_ = True

        if longd:
            self.printError(f"Atoms TOO alone (> {maxDist} Ang.)")
            self.printQuoted(longd[:-1])
            exit_ = True

        if exit_:
            if self.force:
                self.printWarn("You chose to proceed anyway with '-f' option. GOOD LUCK!")
            else:
                self.printError("Use '-f' option if you want to proceed anyway. Aborting ...")
                rmtree(self.tmpDir)
                raise Exception("Coordinates issues with your system")
        try:  # scape resname list index out of range
            resname = list(residues)[0].strip()
            newresname = resname
        except Exception:
            self.printError("resname list index out of range, using default resname: 'LIG'")
            resname = "LIG"
            newresname = resname

        # To avoid resname likes: 001 (all numbers), 1e2 (sci number), ADD : reserved terms for leap
        leapWords = [
            "_cmd_options_",
            "_types_",
            "add",
            "addAtomTypes",
            "addIons",
            "addIons2",
            "addPath",
            "addPdbAtomMap",
            "addPdbResMap",
            "alias",
            "alignAxes",
            "bond",
            "bondByDistance",
            "center",
            "charge",
            "check",
            "clearPdbAtomMap",
            "clearPdbResMap",
            "clearVariables",
            "combine",
            "copy",
            "createAtom",
            "createParmset",
            "createResidue",
            "createUnit",
            "crossLink",
            "debugOff",
            "debugOn",
            "debugStatus",
            "deleteBond",
            "deleteOffLibEntry",
            "deleteRestraint",
            "desc",
            "deSelect",
            "displayPdbAtomMap",
            "displayPdbResMap",
            "edit",
            "flip",
            "groupSelectedAtoms",
            "help",
            "impose",
            "list",
            "listOff",
            "loadAmberParams",
            "loadAmberPrep",
            "loadMol2",
            "loadOff",
            "loadPdb",
            "loadPdbUsingSeq",
            "logFile",
            "matchVariables",
            "measureGeom",
            "quit",
            "relax",
            "remove",
            "restrainAngle",
            "restrainBond",
            "restrainTorsion",
            "saveAmberParm",
            "saveAmberParmPert",
            "saveAmberParmPol",
            "saveAmberParmPolPert",
            "saveAmberPrep",
            "saveMol2",
            "saveOff",
            "saveOffParm",
            "savePdb",
            "scaleCharges",
            "select",
            "sequence",
            "set",
            "setBox",
            "solvateBox",
            "solvateCap",
            "solvateDontClip",
            "solvateOct",
            "solvateShell",
            "source",
            "transform",
            "translate",
            "verbosity",
            "zMatrix",
        ]
        isLeapWord = False
        for word in leapWords:
            if resname.upper().startswith(word.upper()):
                self.printDebug(f"Residue name is a reserved word: '{word.upper()}'")
                isLeapWord = True
        try:
            float(resname)
            self.printDebug(f"Residue name is a 'number': '{resname}'")
            isNumber = True
        except ValueError:
            isNumber = False

        if resname[0].isdigit() or isNumber or isLeapWord:
            newresname = "R" + resname
        if not resname.isalnum():
            newresname = "MOL"
        if newresname != resname:
            self.printWarn(
                f"In {self.acBaseName}.lib, residue name will be '{newresname}' instead of '{resname}' elsewhere"
            )

        self.resName = newresname

        os.chdir(localDir)
        self.printDebug("setResNameCheckCoords done")

    def readMol2TotalCharge(self, mol2File):
        """Reads the charges in given mol2 file and returns the total"""
        charge = 0.0
        ll = []
        cmd = f"{self.acExe} -dr no -i {mol2File} -fi mol2 -o tmp -fo mol2 -c wc -cf tmp.crg -pf y"
        if self.debug:
            self.printMess("Debugging...")
            cmd = cmd.replace("-pf y", "-pf n")

        self.printDebug(cmd)

        log = _getoutput(cmd)

        if os.path.exists("tmp.crg"):
            tmpFile = open("tmp.crg", "r")
            tmpData = tmpFile.readlines()
            for line in tmpData:
                ll += line.split()
            charge = sum(map(float, ll))
        if not log.isspace() and self.debug:
            self.printQuoted(log)

        self.printDebug("readMol2TotalCharge: " + str(charge))

        return charge

    def execAntechamber(self, chargeType=None, atomType=None):
        """To call Antechamber and execute it

        Welcome to antechamber 17.3: molecular input file processor.

        Usage: antechamber -i   input file name
                           -fi  input file format
                           -o   output file name
                           -fo  output file format
                           -c   charge method
                           -cf  charge file name
                           -nc  net molecular charge (int)
                           -a   additional file name
                           -fa  additional file format
                           -ao  additional file operation
                                crd   : only read in coordinate
                                crg   : only read in charge
                                radius: only read in radius
                                name  : only read in atom name
                                type  : only read in atom type
                                bond  : only read in bond type
                           -m   multiplicity (2S+1), default is 1
                           -rn  residue name, overrides input file, default is MOL
                           -rf  residue toplogy file name in prep input file,
                                default is molecule.res
                           -ch  check file name for gaussian, default is 'molecule'
                           -ek  mopac or sqm keyword, inside quotes; overwrites previous ones
                           -gk  gaussian job keyword, inside quotes
                           -gm  gaussian memory keyword, inside quotes, such as "%mem=1000MB"
                           -gn  gaussian number of processors keyword, inside quotes, such as "%nproc=8"
                           -gv  add keyword to generate gesp file (for Gaussian 09 only)
                                1    : yes
                                0    : no, the default
                           -ge  gaussian esp file generated by iop(6/50=1), default is g09.gesp
                           -df  am1-bcc precharge flag, 2 - use sqm(default); 0 - use mopac
                           -at  atom type
                                gaff : the default
                                gaff2: for gaff2 (beta-version)
                                amber: for PARM94/99/99SB
                                bcc  : bcc
                                sybyl: sybyl
                           -du  fix duplicate atom names: yes(y)[default] or no(n)
                           -bk  component/block Id, for ccif
                           -an  adjust atom names: yes(y) or no(n)
                                the default is 'y' for 'mol2' and 'ac' and 'n' for the other formats
                           -j   atom type and bond type prediction index, default is 4
                                0    : no assignment
                                1    : atom type
                                2    : full  bond types
                                3    : part  bond types
                                4    : atom and full bond type
                                5    : atom and part bond type
                           -s   status information: 0(brief), 1(default) or 2(verbose)
                           -eq  equalizing atomic charge, default is 1 for '-c resp' and '-c bcc' and 0
                                for the other charge methods
                                0    : no use
                                1    : by atomic paths
                                2    : by atomic paths and structural information, i.e. E/Z configurations
                           -pf  remove intermediate files: yes(y) or no(n)[default]
                           -pl  maximum path length to determin equivalence of atomic charges for resp and bcc,
                                the smaller the value, the faster the algorithm, default is -1 (use full length),
                                set this parameter to 10 to 30 if your molecule is big (# atoms >= 100)
                           -dr  acdoctor mode: yes(y)[default] or no(n)
                           -i -o -fi and -fo must appear; others are optional
                           Use 'antechamber -L' to list the supported file formats and charge methods

                              List of the File Formats

                 file format type  abbre. index | file format type abbre. index
                --------------------------------------------------------------
                Antechamber        ac       1  | Sybyl Mol2         mol2    2
                PDB                pdb      3  | Modified PDB       mpdb    4
                AMBER PREP (int)   prepi    5  | AMBER PREP (car)   prepc   6
                Gaussian Z-Matrix  gzmat    7  | Gaussian Cartesian gcrt    8
                Mopac Internal     mopint   9  | Mopac Cartesian    mopcrt 10
                Gaussian Output    gout    11  | Mopac Output       mopout 12
                Alchemy            alc     13  | CSD                csd    14
                MDL                mdl     15  | Hyper              hin    16
                AMBER Restart      rst     17  | Jaguar Cartesian   jcrt   18
                Jaguar Z-Matrix    jzmat   19  | Jaguar Output      jout   20
                Divcon Input       divcrt  21  | Divcon Output      divout 22
                SQM Input          sqmcrt  23  | SQM Output         sqmout 24
                Charmm             charmm  25  | Gaussian ESP       gesp   26
                Component cif      ccif    27  |
                --------------------------------------------------------------
                        AMBER restart file can only be read in as additional file.

                              List of the Charge Methods

                charge method     abbre. index | charge method    abbre. index
                --------------------------------------------------------------
                RESP               resp     1  |  AM1-BCC           bcc     2
                CM1                cm1      3  |  CM2               cm2     4
                ESP (Kollman)      esp      5  |  Mulliken          mul     6
                Gasteiger          gas      7  |  Read in charge    rc      8
                Write out charge   wc       9  |  Delete Charge     dc     10
                --------------------------------------------------------------
        """
        global pid

        self.printMess("Executing Antechamber...")

        self.makeDir()

        ct = chargeType or self.chargeType
        at = atomType or self.atomType
        if "amber2" in at:
            at = "amber"

        if ct == "user":
            ct = ""
        else:
            ct = f"-c {ct}"

        exten = self.ext[1:]
        if exten == "mol":
            exten = "mdl"

        cmd = "%s -dr no -i %s -fi %s -o %s -fo mol2 %s -nc %s -m %s -s 2 -df %s -at %s -pf y %s" % (
            self.acExe,
            self.inputFile,
            exten,
            self.acMol2FileName,
            ct,
            self.chargeVal,
            self.multiplicity,
            self.qFlag,
            at,
            self.ekFlag,
        )

        if self.debug:
            self.printMess("Debugging...")
            cmd = cmd.replace("-pf y", "-pf n")

        self.printDebug(cmd)

        if os.path.exists(self.acMol2FileName) and not self.force:
            self.printMess("AC output file present... doing nothing")
        else:
            try:
                os.remove(self.acMol2FileName)
            except Exception:
                pass
            signal.signal(signal.SIGALRM, self.signal_handler)
            signal.alarm(self.timeTol)
            p = sub.Popen(cmd, shell=True, stderr=sub.STDOUT, stdout=sub.PIPE)
            pid = p.pid

            out = str(p.communicate()[0].decode())  # p.stdout.read()
            self.acLog = out

        if os.path.exists(self.acMol2FileName):
            self.printMess("* Antechamber OK *")
        else:
            self.printQuoted(self.acLog)
            return True
        return False

    def signal_handler(self, _signum, _frame):  # , pid = 0):
        """Signal handler"""
        global pid
        pids = job_pids_family(pid)
        self.printDebug("PID: %s, PIDS: %s" % (pid, pids))
        self.printMess("Timed out! Process %s killed, max exec time (%ss) exceeded" % (pids, self.timeTol))
        # os.system('kill -15 %s' % pids)
        for i in pids.split():
            os.kill(int(i), 15)
        raise Exception("Semi-QM taking too long to finish... aborting!")

    def delOutputFiles(self):
        """Delete temporary output files"""
        delFiles = [
            "mopac.in",
            "tleap.in",
            "fixbo.log",
            "addhs.log",
            "ac_tmp_ot.mol2",
            "frcmod.ac_tmp",
            "fragment.mol2",
            self.tmpDir,
        ]  # , 'divcon.pdb', 'mopac.pdb', 'mopac.out'] #'leap.log'
        self.printMess("Removing temporary files...")
        for file_ in delFiles:
            file_ = os.path.join(self.absHomeDir, file_)
            if os.path.exists(file_):
                if os.path.isdir(file_):
                    rmtree(file_)
                else:
                    os.remove(file_)

    def checkXyzAndTopFiles(self):
        """Cehck XYZ and TOP files"""
        fileXyz = self.acXyzFileName
        fileTop = self.acTopFileName
        if os.path.exists(fileXyz) and os.path.exists(fileTop):
            return True
        return False

    def execTleap(self):
        """Execute tleap"""
        fail = False

        self.makeDir()

        if self.ext == ".pdb":
            self.printMess("... converting pdb input file to mol2 input file")
            if self.convertPdbToMol2():
                self.printError("convertPdbToMol2 failed")

        if self.execAntechamber():
            self.printError("Antechamber failed")
            fail = True

        if self.execParmchk():
            self.printError("Parmchk failed")
            fail = True

        if fail:
            return True

        tleapScpt = TLEAP_TEMPLATE % self.acParDict

        fp = open("tleap.in", "w")
        fp.write(tleapScpt)
        fp.close()

        cmd = "%s -f tleap.in" % self.tleapExe

        if self.checkXyzAndTopFiles() and not self.force:
            self.printMess("Topologies files already present... doing nothing")
        else:
            try:
                os.remove(self.acTopFileName)
                os.remove(self.acXyzFileName)
            except Exception:
                pass
            self.printMess("Executing Tleap...")
            self.printDebug(cmd)
            self.tleapLog = _getoutput(cmd)
            self.checkLeapLog(self.tleapLog)

        if self.checkXyzAndTopFiles():
            self.printMess("* Tleap OK *")
        else:
            self.printQuoted(self.tleapLog)
            return True
        return False

    def checkLeapLog(self, log):
        """Check Leap log"""
        log = log.splitlines(True)
        check = ""
        block = False
        for line in log:
            # print "*"+line+"*"
            if "Checking '" in line:
                # check += line
                block = True
            if "Checking Unit." in line:
                block = False
            if block:
                check += line
        self.printQuoted(check[:-1])

    def locateDat(self, aFile):
        """locate a file pertinent to $AMBERHOME/dat/leap/parm/"""
        amberhome = os.environ.get("AMBERHOME")
        if amberhome:
            aFileF = os.path.join(amberhome, "dat/leap/parm", aFile)
            if os.path.exists(aFileF):
                return aFileF
        aFileF = os.path.join(os.path.dirname(self.acExe), "../dat/leap/parm", aFile)
        if os.path.exists(aFileF):
            return aFileF
        return None

    def execParmchk(self):
        """Execute parmchk"""
        self.makeDir()
        cmd = f"{self.parmchkExe} -i {self.acMol2FileName} -f mol2 -o {self.acFrcmodFileName}"

        if "amber" in self.atomType:
            gaffFile = self.locateDat(self.gaffDatfile)
            parmfile = self.locateDat("parm10.dat")
            frcmodffxxSB = self.locateDat("frcmod.ff14SB")
            # frcmodparmbsc0 = self.locateDat('frcmod.parmbsc0')
            parmGaffFile = parmMerge(parmfile, gaffFile)
            parmGaffffxxSBFile = parmMerge(parmGaffFile, frcmodffxxSB, frcmod=True)
            # parm99gaffff99SBparmbsc0File = parmMerge(parm99gaffff99SBFile, frcmodparmbsc0, frcmod = True)
            # parm10file = self.locateDat('parm10.dat') # PARM99 + frcmod.ff99SB + frcmod.parmbsc0 in AmberTools 1.4

            cmd += f" -p {parmGaffffxxSBFile}"  # Ignoring BSC0
        elif "gaff2" in self.atomType:
            cmd += " -s 2"

        self.parmchkLog = _getoutput(cmd)

        self.printDebug(cmd)

        if os.path.exists(self.acFrcmodFileName):
            check = self.checkFrcmod()
            if check:
                self.printWarn("Couldn't determine all parameters:")
                self.printMess("From file '%s'\n" % self.acFrcmodFileName + check)
            else:
                self.printMess("* Parmchk OK *")
        else:
            self.printQuoted(self.parmchkLog)
            return True
        return False

    def checkFrcmod(self):
        """Check FRCMOD file"""
        check = ""
        frcmodContent = open(self.acFrcmodFileName, "r").readlines()
        for line in frcmodContent:
            if "ATTN, need revision" in line:
                check += line
        return check

    def convertPdbToMol2(self):
        """Convert PDB to MOL2 by using babel"""
        if self.ext == ".pdb":
            if self.execBabel():
                self.printError("convert pdb to mol2 via babel failed")
                return True
        return False

    def convertSmilesToMol2(self):
        if not self.babelExe:
            raise Exception("SMILES needs openbabel python module")
        if checkOpenBabelVersion() >= 300:
            from openbabel import pybel

        elif checkOpenBabelVersion() >= 200 and checkOpenBabelVersion() < 300:
            import pybel  # type: ignore

        """Convert Smiles to MOL2 by using babel"""
        try:
            mymol = pybel.readstring("smi", str(self.smiles))
            mymol.addh()
            mymol.make3D()
            mymol.write(self.ext.replace(".", ""), self.absInputFile, overwrite=True)
            return True
        except Exception:
            return False

    def execBabel(self):
        """Execute babel"""
        self.makeDir()

        cmd = f"{self.babelExe} -ipdb {self.inputFile} -omol2 -O {self.baseName}.mol2"
        self.printDebug(cmd)
        self.babelLog = _getoutput(cmd)
        self.ext = ".mol2"
        self.inputFile = self.baseName + self.ext
        self.acParDict["ext"] = "mol2"
        if os.path.exists(self.inputFile):
            self.printMess("* Babel OK *")
        else:
            self.printQuoted(self.babelLog)
            return True
        return False

    def makeDir(self):
        """Make Dir"""
        os.chdir(self.rootDir)
        self.absHomeDir = os.path.abspath(self.homeDir)
        if not os.path.exists(self.homeDir):
            os.mkdir(self.homeDir)
        os.chdir(self.homeDir)
        if self.absInputFile:
            copy2(self.absInputFile, ".")
        return True

    def createACTopol(self):
        """
        If successful, Amber Top and Xyz files will be generated
        """
        if self.execTleap():
            self.printError("Tleap failed")
        if not self.debug:
            self.delOutputFiles()

    def createMolTopol(self):
        """
        Create molTop obj
        """
        self.topFileData = open(self.acTopFileName, "r").readlines()
        self.molTopol = MolTopol(
            self,
            verbose=self.verbose,
            debug=self.debug,
            gmx4=self.gmx4,
            disam=self.disam,
            direct=self.direct,
            is_sorted=self.sorted,
            chiral=self.chiral,
        )
        if self.outTopols:
            if "cns" in self.outTopols:
                self.molTopol.writeCnsTopolFiles()
            if "gmx" in self.outTopols:
                self.molTopol.writeGromacsTopolFiles()
            if "charmm" in self.outTopols:
                self.writeCharmmTopolFiles()
        try:  # scape the pickle save error
            self.pickleSave()
        except Exception:
            self.printError("pickleSave failed")
            pass

        self.delOutputFiles()  # required to use on Jupyter Notebook

    def pickleSave(self):
        """
        To restore:
            from acpype import *
            #import cPickle as pickle
            import pickle
            o = pickle.load(open('DDD.pkl','rb'))
            NB: It fails to restore with ipython in Mac (Linux OK)
        """
        pklFile = self.baseName + ".pkl"
        dumpFlag = False
        if not os.path.exists(pklFile):
            mess = "Writing pickle file %s" % pklFile
            dumpFlag = True
        elif self.force:
            mess = "Overwriting pickle file %s" % pklFile
            dumpFlag = True
        else:
            mess = "Pickle file %s already present... doing nothing" % pklFile
        self.printMess(mess)
        if dumpFlag:
            with open(pklFile, "wb") as f:  # for python 3.3 or higher
                pickle.dump(self, f)

    def getFlagData(self, flag):
        """
        For a given acFileTop flag, return a list of the data related
        """
        block = False
        tFlag = "%FLAG " + flag
        data = ""

        if not self.topFileData:
            raise Exception("PRMTOP file empty?")

        for rawLine in self.topFileData:
            if "%COMMENT" in rawLine:
                continue
            line = rawLine.replace("\r", "").replace("\n", "")
            if tFlag in line:
                block = True
                continue
            if block and "%FLAG " in line:
                break
            if block:
                if "%FORMAT" in line:
                    line = line.strip().strip("%FORMAT()").split(".")[0]
                    for c in line:
                        if c.isalpha():
                            f = int(line.split(c)[1])
                            break
                    continue
                data += line
        # data need format
        sdata = [data[i : i + f].strip() for i in range(0, len(data), f)]
        if "+" and "." in data and flag != "RESIDUE_LABEL":  # it's a float
            ndata = list(map(float, sdata))
        elif flag != "RESIDUE_LABEL":
            try:  # try if it's integer
                ndata = list(map(int, sdata))
            except Exception:
                ndata = sdata
        else:
            ndata = sdata
        if flag == "AMBER_ATOM_TYPE":
            nn = []
            ll = set()
            prefixed = False
            for ii in ndata:
                prefixed = True
                if ii[0].isdigit():
                    ll.add(ii)
                    ii = "A" + ii
                nn.append(ii)
            if prefixed and ll:
                self.printDebug("GMX does not like atomtype starting with Digit")
                self.printDebug("prefixing AtomType %s with 'A'." % list(ll))
            ndata = nn
        return ndata  # a list

    def getResidueLabel(self):
        """
        Get a 3 capital letters code from acFileTop
        Returns a list.
        """
        residueLabel = self.getFlagData("RESIDUE_LABEL")
        residueLabel = list(map(str, residueLabel))
        if residueLabel[0] != residueLabel[0].upper():
            self.printWarn("residue label '%s' in '%s' is not all UPPERCASE" % (residueLabel[0], self.inputFile))
            self.printWarn("this may raise problem with some applications like CNS")
        self.residueLabel = residueLabel

    def getCoords(self):
        """
        For a given acFileXyz file, return a list of coords as:
        [[x1,y1,z1],[x2,y2,z2], etc.]
        """
        if not self.xyzFileData:
            raise Exception("INPCRD file empty?")
        data = ""
        for rawLine in self.xyzFileData[2:]:
            line = rawLine.replace("\r", "").replace("\n", "")
            data += line
        ll = len(data)
        ndata = list(map(float, [data[i : i + 12] for i in range(0, ll, 12)]))

        gdata = []
        for i in range(0, len(ndata), 3):
            gdata.append([ndata[i], ndata[i + 1], ndata[i + 2]])

        self.printDebug("getCoords done")

        return gdata

    def getAtoms(self):
        """
        Set a list with all atoms objects build from dat in acFileTop
        Set also if molTopol atom type system is gaff or amber
        Set also list atomTypes
        Set also resid
        Set also molTopol total charge
        """
        atomNameList = self.getFlagData("ATOM_NAME")
        atomTypeNameList = self.getFlagData("AMBER_ATOM_TYPE")
        self._atomTypeNameList = atomTypeNameList
        massList = self.getFlagData("MASS")
        chargeList = self.getFlagData("CHARGE")
        resIds = self.getFlagData("RESIDUE_POINTER") + [0]
        coords = self.getCoords()
        ACOEFs, BCOEFs = self.getABCOEFs()

        atoms = []
        atomTypes = []
        tmpList = []  # a list with unique atom types
        totalCharge = 0.0
        countRes = 0
        id_ = 0
        FirstNonSoluteId = None
        for atomName in atomNameList:
            if atomName != atomName.upper():
                self.printDebug("atom name '%s' HAS to be all UPPERCASE... Applying this here." % atomName)
                atomName = atomName.upper()
            atomTypeName = atomTypeNameList[id_]
            if id_ + 1 == resIds[countRes]:
                resid = countRes
                countRes += 1
            resName = self.residueLabel[resid]
            if resName in ionOrSolResNameList and not FirstNonSoluteId:
                FirstNonSoluteId = id_
            mass = massList[id_]
            charge = chargeList[id_]
            chargeConverted = charge / qConv
            totalCharge += charge
            coord = coords[id_]
            ACOEF = ACOEFs[id_]
            BCOEF = BCOEFs[id_]
            atomType = AtomType(atomTypeName, mass, ACOEF, BCOEF)
            if atomTypeName not in tmpList:
                tmpList.append(atomTypeName)
                atomTypes.append(atomType)
            atom = Atom(atomName, atomType, id_ + 1, resid, mass, chargeConverted, coord)
            atoms.append(atom)
            id_ += 1

        balanceChargeList, balanceValue, balanceIds = self.balanceCharges(chargeList, FirstNonSoluteId)

        for id_ in balanceIds:
            atoms[id_].charge = balanceValue / qConv

        if atomTypeName[0].islower():
            self.atomTypeSystem = "gaff"
        else:
            self.atomTypeSystem = "amber"
        self.printDebug("Balanced TotalCharge %13.10f" % float(sum(balanceChargeList) / qConv))

        self.totalCharge = int(round(totalCharge / qConv))

        self.atoms = atoms
        self.atomTypes = atomTypes

        self.pbc = None
        if len(coords) == len(atoms) + 2 or len(coords) == len(atoms) * 2 + 2:
            self.pbc = [coords[-2], coords[-1]]
        self.printDebug("PBC = %s" % self.pbc)
        self.printDebug("getAtoms done")

    def getBonds(self):
        """Get Bonds"""
        uniqKbList = self.getFlagData("BOND_FORCE_CONSTANT")
        uniqReqList = self.getFlagData("BOND_EQUIL_VALUE")
        bondCodeHList = self.getFlagData("BONDS_INC_HYDROGEN")
        bondCodeNonHList = self.getFlagData("BONDS_WITHOUT_HYDROGEN")
        bondCodeList = bondCodeHList + bondCodeNonHList
        bonds = []
        for i in range(0, len(bondCodeList), 3):
            idAtom1 = bondCodeList[i] // 3  # remember python starts with id 0
            idAtom2 = bondCodeList[i + 1] // 3
            bondTypeId = bondCodeList[i + 2] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            kb = uniqKbList[bondTypeId]
            req = uniqReqList[bondTypeId]
            atoms = [atom1, atom2]
            bond = Bond(atoms, kb, req)
            bonds.append(bond)
        self.bonds = bonds
        self.printDebug("getBonds done")

    def getAngles(self):
        """Get Angles"""
        uniqKtList = self.getFlagData("ANGLE_FORCE_CONSTANT")
        uniqTeqList = self.getFlagData("ANGLE_EQUIL_VALUE")
        # for list below, true atom number = index/3 + 1
        angleCodeHList = self.getFlagData("ANGLES_INC_HYDROGEN")
        angleCodeNonHList = self.getFlagData("ANGLES_WITHOUT_HYDROGEN")
        angleCodeList = angleCodeHList + angleCodeNonHList
        angles = []
        for i in range(0, len(angleCodeList), 4):
            idAtom1 = angleCodeList[i] // 3  # remember python starts with id 0
            idAtom2 = angleCodeList[i + 1] // 3
            idAtom3 = angleCodeList[i + 2] // 3
            angleTypeId = angleCodeList[i + 3] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            kt = uniqKtList[angleTypeId]
            teq = uniqTeqList[angleTypeId]  # angle given in rad in prmtop
            atoms = [atom1, atom2, atom3]
            angle = Angle(atoms, kt, teq)
            angles.append(angle)
        self.angles = angles
        self.printDebug("getAngles done")

    def getDihedrals(self):
        """
        Get dihedrals (proper and imp), condensed list of prop dih and
        atomPairs
        """
        uniqKpList = self.getFlagData("DIHEDRAL_FORCE_CONSTANT")
        uniqPeriodList = self.getFlagData("DIHEDRAL_PERIODICITY")
        uniqPhaseList = self.getFlagData("DIHEDRAL_PHASE")
        # for list below, true atom number = abs(index)/3 + 1
        dihCodeHList = self.getFlagData("DIHEDRALS_INC_HYDROGEN")
        dihCodeNonHList = self.getFlagData("DIHEDRALS_WITHOUT_HYDROGEN")
        dihCodeList = dihCodeHList + dihCodeNonHList
        properDih = []
        improperDih = []
        condProperDih = []  # list of dihedrals condensed by the same quartet
        # atomPairs = []
        atomPairs = set()
        for i in range(0, len(dihCodeList), 5):
            idAtom1 = dihCodeList[i] // 3  # remember python starts with id 0
            idAtom2 = dihCodeList[i + 1] // 3
            # 3 and 4 indexes can be negative: if id3 < 0, end group interations
            # in amber are to be ignored; if id4 < 0, dihedral is improper
            idAtom3raw = dihCodeList[i + 2] // 3  # can be negative -> exclude from 1-4vdw
            idAtom4raw = dihCodeList[i + 3] // 3  # can be negative -> Improper
            idAtom3 = abs(idAtom3raw)
            idAtom4 = abs(idAtom4raw)
            dihTypeId = dihCodeList[i + 4] - 1
            atom1 = self.atoms[idAtom1]
            atom2 = self.atoms[idAtom2]
            atom3 = self.atoms[idAtom3]
            atom4 = self.atoms[idAtom4]
            kPhi = uniqKpList[dihTypeId]  # already divided by IDIVF
            period = int(uniqPeriodList[dihTypeId])  # integer
            phase = uniqPhaseList[dihTypeId]  # angle given in rad in prmtop
            if phase == kPhi == 0:
                period = 0  # period is set to 0
            atoms = [atom1, atom2, atom3, atom4]
            dihedral = Dihedral(atoms, kPhi, period, phase)
            if idAtom4raw > 0:
                try:
                    atomsPrev = properDih[-1].atoms
                except Exception:
                    atomsPrev = []
                properDih.append(dihedral)
                if idAtom3raw < 0 and atomsPrev == atoms:
                    condProperDih[-1].append(dihedral)
                else:
                    condProperDih.append([dihedral])
                pair = (atom1, atom4)
                # if atomPairs.count(pair) == 0 and idAtom3raw > 0:
                if idAtom3raw > 0:
                    atomPairs.add(pair)
            else:
                improperDih.append(dihedral)
        try:
            atomPairs = sorted(atomPairs)
        except Exception:
            pass
        self.properDihedrals = properDih
        self.improperDihedrals = improperDih
        self.condensedProperDihedrals = condProperDih  # [[],[],...]
        self.atomPairs = atomPairs  # set((atom1, atom2), ...)
        self.printDebug("getDihedrals done")

    def getChirals(self):
        """
        Get chiral atoms, its 4 neighbours and improper dihedral angle
        to store non-planar improper dihedrals for CNS (and CNS only!)
        """
        self.chiralGroups = []
        if not self.babelExe:
            self.printWarn("No Openbabel python module, no chiral groups")

        if checkOpenBabelVersion() >= 300:
            from openbabel import openbabel as ob
            from openbabel import pybel

        elif checkOpenBabelVersion() >= 200 and checkOpenBabelVersion() < 300:
            import openbabel as ob
            import pybel  # type: ignore

        self.printMess("Using OpenBabel v." + ob.OBReleaseVersion() + "\n")

        " obchiral script - replace the obchiral executable"
        out = []
        _filename, file_extension = os.path.splitext(self.inputFile)
        mol = pybel.readfile(file_extension.replace(".", ""), self.inputFile)
        for ml in mol:
            for at in ml:
                if ob.OBStereoFacade(ml.OBMol).HasTetrahedralStereo(at.idx):
                    out.append(at.idx)
        " end of obchiral script "

        chiralGroups = []
        for id_ in out:
            atChi = self.atoms[id_]
            quad = []
            for bb in self.bonds:
                bAts = bb.atoms[:]
                if atChi in bAts:
                    bAts.remove(atChi)
                    quad.append(bAts[0])
            if len(quad) != 4:
                if self.chiral:
                    self.printWarn("Atom %s has less than 4 connections to 4 different atoms. It's NOT Chiral!" % atChi)
                continue
            v1, v2, v3, v4 = [x.coords for x in quad]
            chiralGroups.append((atChi, quad, imprDihAngle(v1, v2, v3, v4)))
        self.chiralGroups = chiralGroups

    def sortAtomsForGromacs(self):
        """
        Re-sort atoms for gromacs, which expects all hydrogens to immediately
        follow the heavy atom they are bonded to and belong to the same charge
        group.

        Currently, atom mass < 1.2 is taken to denote a proton.  This behavior
        may be changed by modifying the 'is_hydrogen' function within.

        JDC 2011-02-03
        """

        # Build dictionary of bonded atoms.
        bonded_atoms = dict()
        for atom in self.atoms:
            bonded_atoms[atom] = list()
        for bond in self.bonds:
            [atom1, atom2] = bond.atoms
            bonded_atoms[atom1].append(atom2)
            bonded_atoms[atom2].append(atom1)

        # Define hydrogen and heavy atom classes.
        def is_hydrogen(atom):
            """Check for H"""
            return atom.mass < 1.2

        def is_heavy(atom):
            """Check for non H"""
            return not is_hydrogen(atom)

        # Build list of sorted atoms, assigning charge groups by heavy atom.
        sorted_atoms = list()
        cgnr = 1  # charge group number: each heavy atoms is assigned its own charge group
        # First pass: add heavy atoms, followed by the hydrogens bonded to them.
        for atom in self.atoms:
            if is_heavy(atom):
                # Append heavy atom.
                atom.cgnr = cgnr
                sorted_atoms.append(atom)
                # Append all hydrogens.
                for bonded_atom in bonded_atoms[atom]:
                    if is_hydrogen(bonded_atom) and bonded_atom not in sorted_atoms:
                        # Append bonded hydrogen.
                        bonded_atom.cgnr = cgnr
                        sorted_atoms.append(bonded_atom)
                cgnr += 1

        # Second pass: Add any remaining atoms.
        if len(sorted_atoms) < len(self.atoms):
            for atom in self.atoms:
                if atom not in sorted_atoms:
                    atom.cgnr = cgnr
                    sorted_atoms.append(atom)
                    cgnr += 1

        # Replace current list of atoms with sorted list.
        self.atoms = sorted_atoms

        # Renumber atoms in sorted list, starting from 1.
        for (index, atom) in enumerate(self.atoms):
            atom.id = index + 1

    def balanceCharges(self, chargeList, FirstNonSoluteId=None):
        """
        Note that python is very annoying about floating points.
        Even after balance, there will always be some residue of order e-12
        to e-16, which is believed to vanished once one writes a topology
        file, say, for CNS or GMX, where floats are represented with 4 or 5
        maximum decimals.
        """
        limIds = []
        total = sum(chargeList)
        totalConverted = total / qConv
        self.printDebug("charge to be balanced: total %13.10f" % (totalConverted))
        maxVal = max(chargeList[:FirstNonSoluteId])
        minVal = min(chargeList[:FirstNonSoluteId])
        if abs(maxVal) >= abs(minVal):
            lim = maxVal
        else:
            lim = minVal
        nLims = chargeList.count(lim)
        # limId = chargeList.index(lim)
        diff = totalConverted - round(totalConverted)
        fix = lim - diff * qConv / nLims
        id_ = 0
        for c in chargeList:
            if c == lim:
                limIds.append(id_)
                chargeList[id_] = fix
            id_ += 1
        # self.printDebug(chargeList)
        self.printDebug("balanceCharges done")
        return chargeList, fix, limIds

    def getABCOEFs(self):
        """Get non-bonded coeficients"""
        uniqAtomTypeIdList = self.getFlagData("ATOM_TYPE_INDEX")
        nonBonIdList = self.getFlagData("NONBONDED_PARM_INDEX")
        rawACOEFs = self.getFlagData("LENNARD_JONES_ACOEF")
        rawBCOEFs = self.getFlagData("LENNARD_JONES_BCOEF")
        # print nonBonIdList, len(nonBonIdList), rawACOEFs, len(rawACOEFs)
        ACOEFs = []
        BCOEFs = []
        ntypes = max(uniqAtomTypeIdList)
        # id_ = 0
        # for atName in self._atomTypeNameList:
        for id_ in range(len(self._atomTypeNameList)):
            # id_ = self._atomTypeNameList.index(atName)
            atomTypeId = uniqAtomTypeIdList[id_]
            index = ntypes * (atomTypeId - 1) + atomTypeId
            nonBondId = nonBonIdList[index - 1]
            # print "*****", index, ntypes, atName, id_, atomTypeId, nonBondId
            ACOEFs.append(rawACOEFs[nonBondId - 1])
            BCOEFs.append(rawBCOEFs[nonBondId - 1])
            # id_ += 1
        # print ACOEFs
        self.printDebug("getABCOEFs done")
        return ACOEFs, BCOEFs

    def setProperDihedralsCoef(self):
        """
        It takes self.condensedProperDihedrals and returns
        self.properDihedralsCoefRB, a reduced list of quartet atoms + RB.
        Coeficients ready for GMX (multiplied by 4.184)

        self.properDihedralsCoefRB = [ [atom1,..., atom4], C[0:5] ]

        For proper dihedrals: a quartet of atoms may appear with more than
        one set of parameters and to convert to GMX they are treated as RBs.

        The resulting coefs calculated here may look slighted different from
        the ones calculated by amb2gmx.pl because python is taken full float
        number from prmtop and not rounded numbers from rdparm.out as
        amb2gmx.pl does.
        """
        properDihedralsCoefRB = []
        properDihedralsAlphaGamma = []
        properDihedralsGmx45 = []
        for item in self.condensedProperDihedrals:
            V = 6 * [0.0]
            C = 6 * [0.0]
            for dih in item:
                period = dih.period  # Pn
                kPhi = dih.kPhi  # in rad
                phaseRaw = dih.phase * radPi  # in degree
                phase = int(phaseRaw)  # in degree
                if period > 4 and self.gmx4:
                    raise Exception("Likely trying to convert ILDN to RB, DO NOT use option '-z'")
                if phase in [0, 180]:
                    properDihedralsGmx45.append([item[0].atoms, phaseRaw, kPhi, period])
                    if self.gmx4:
                        if kPhi != 0:
                            V[period] = 2 * kPhi * cal
                        if period == 1:
                            C[0] += 0.5 * V[period]
                            if phase == 0:
                                C[1] -= 0.5 * V[period]
                            else:
                                C[1] += 0.5 * V[period]
                        elif period == 2:
                            if phase == 180:
                                C[0] += V[period]
                                C[2] -= V[period]
                            else:
                                C[2] += V[period]
                        elif period == 3:
                            C[0] += 0.5 * V[period]
                            if phase == 0:
                                C[1] += 1.5 * V[period]
                                C[3] -= 2 * V[period]
                            else:
                                C[1] -= 1.5 * V[period]
                                C[3] += 2 * V[period]
                        elif period == 4:
                            if phase == 180:
                                C[2] += 4 * V[period]
                                C[4] -= 4 * V[period]
                            else:
                                C[0] += V[period]
                                C[2] -= 4 * V[period]
                                C[4] += 4 * V[period]
                else:
                    properDihedralsAlphaGamma.append([item[0].atoms, phaseRaw, kPhi, period])
                    # print phaseRaw, kPhi, period
            if phase in [0, 180]:
                properDihedralsCoefRB.append([item[0].atoms, C])

        self.printDebug("setProperDihedralsCoef done")

        self.properDihedralsCoefRB = properDihedralsCoefRB
        self.properDihedralsAlphaGamma = properDihedralsAlphaGamma
        self.properDihedralsGmx45 = properDihedralsGmx45

    def writeCharmmTopolFiles(self):
        """Write CHARMM topology files"""
        self.printMess("Writing CHARMM files\n")

        at = self.atomType
        self.getResidueLabel()
        res = self.resName  # self.residueLabel[0]
        # print res, self.residueLabel, type(self.residueLabel)

        cmd = (
            "%s -dr no -i %s -fi mol2 -o %s -fo charmm -s 2 -at %s \
        -pf y -rn %s"
            % (self.acExe, self.acMol2FileName, self.charmmBase, at, res)
        )

        if self.debug:
            cmd = cmd.replace("-pf y", "-pf n")
            self.printDebug(cmd)

        log = _getoutput(cmd)
        if self.debug:
            self.printQuoted(log)

    def writePdb(self, file_):
        """
        Write a new PDB file_ with the atom names defined by Antechamber
        Input: file_ path string
        The format generated here use is slightly different from
        http://www.wwpdb.org/documentation/format23/sect9.html respected to
        atom name
        """
        # TODO: assuming only one residue ('1')
        pdbFile = open(file_, "w")
        fbase = os.path.basename(file_)
        pdbFile.write("REMARK " + head % (fbase, date))
        id_ = 1
        for atom in self.atoms:
            # id_ = self.atoms.index(atom) + 1
            aName = atom.atomName
            if len(aName) == 2 or len(aName) == 3:
                aName = "%s " % aName.capitalize()
                s = aName

            elif len(aName) == 1:
                aName = " %s  " % aName
                if "C" or "H" in aName:
                    s = "Xx"

            #        for ll in aName:
            #            if ll.isalpha():
            #                s = ll
            #                break

            rName = self.residueLabel[0]
            x = atom.coords[0]
            y = atom.coords[1]
            z = atom.coords[2]
            line = "%-6s%5d %4s %3s Z%4d%s%8.3f%8.3f%8.3f%6.2f%6.2f%s%2s\n" % (
                "ATOM",
                id_,
                aName,
                rName,
                1,
                4 * " ",
                x,
                y,
                z,
                1.0,
                0.0,
                10 * " ",
                "".join(filter(str.isalpha, s)),
            )
            pdbFile.write(line)
            id_ += 1
        pdbFile.write("END\n")

    def writeGromacsTopolFiles(self):
        """
        # from ~/Programmes/amber10/dat/leap/parm/gaff.dat
        #atom type        atomic mass        atomic polarizability        comments
        ca                12.01                 0.360                    Sp2 C in pure aromatic systems
        ha                1.008                 0.135                    H bonded to aromatic carbon

        #bonded atoms        harmonic force kcal/mol/A^2       eq. dist. Ang.  comments
        ca-ha                  344.3*                           1.087**         SOURCE3  1496    0.0024    0.0045
        * for gmx: 344.3 * 4.184 * 100 * 2 = 288110 kJ/mol/nm^2 (why factor 2?)
        ** convert Ang to nm ( div by 10) for gmx: 1.087 A = 0.1087 nm
        # CA HA      1    0.10800   307105.6 ; ged from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from ffamber99bon.itp)
        # CA-HA  367.0    1.080       changed from 340. bsd on C6H6 nmodes; PHE,TRP,TYR (from parm99.dat)

        # angle        HF kcal/mol/rad^2    eq angle degrees     comments
        ca-ca-ha        48.5*             120.01                SOURCE3 2980   0.1509   0.2511
        * to convert to gmx: 48.5 * 4.184 * 2 = 405.848 kJ/mol/rad^2 (why factor 2?)
        # CA  CA  HA           1   120.000    418.400 ; new99 (from ffamber99bon.itp)
        # CA-CA-HA    50.0      120.00 (from parm99.dat)

        # dihedral    idivf        barrier hight/2 kcal/mol  phase degrees       periodicity     comments
        X -ca-ca-X    4           14.500*                    180.000             2.000           intrpol.bsd.on C6H6
        *convert 2 gmx: 14.5/4 * 4.184 * 2 (?) (yes in amb2gmx, not in topolbuild, why?) = 30.334 or 15.167 kJ/mol
        # X -CA-CA-X    4   14.50        180.0     2.         intrpol.bsd.on C6H6 (from parm99.dat)
        # X CA CA  X    3   30.334       0.000   -30.33400     0.000     0.000     0.000   ; intrpol.bsd.on C6H6
        ;propers treated as RBs in GMX to use combine multiple AMBER torsions per quartet (from ffamber99bon.itp)

        # impr. dihedral        barrier hight/2      phase degrees       periodicity     comments
        X -X -ca-ha             1.1*                  180.                      2.       bsd.on C6H6 nmodes
        * to convert to gmx: 1.1 * 4.184 = 4.6024 kJ/mol/rad^2
        # X -X -CA-HA         1.1          180.          2.           bsd.on C6H6 nmodes (from parm99.dat)
        # X   X   CA  HA       1      180.00     4.60240     2      ; bsd.on C6H6 nmodes
        ;impropers treated as propers in GROMACS to use correct AMBER analytical function (from ffamber99bon.itp)

        # 6-12 parms     sigma = 2 * r * 2^(-1/6)    epsilon
        # atomtype        radius Ang.                    pot. well depth kcal/mol      comments
          ha                  1.4590*                      0.0150**                         Spellmeyer
          ca                  1.9080                    0.0860                            OPLS
        * to convert to gmx:
            sigma = 1.4590 * 2^(-1/6) * 2 = 2 * 1.29982 Ang. = 2 * 0.129982 nm  = 1.4590 * 2^(5/6)/10 =  0.259964 nm
        ** to convert to gmx: 0.0150 * 4.184 = 0.06276 kJ/mol
        # amber99_3    CA     0.0000  0.0000  A   3.39967e-01  3.59824e-01 (from ffamber99nb.itp)
        # amber99_22   HA     0.0000  0.0000  A   2.59964e-01  6.27600e-02 (from ffamber99nb.itp)
        # C*          1.9080  0.0860             Spellmeyer
        # HA          1.4590  0.0150             Spellmeyer (from parm99.dat)
        # to convert r and epsilon to ACOEF and BCOEF
        # ACOEF = sqrt(e1*e2) * (r1 + r2)^12 ; BCOEF = 2 * sqrt(e1*e2) * (r1 + r2)^6 = 2 * ACOEF/(r1+r2)^6
        # to convert ACOEF and BCOEF to r and epsilon
        # r = 0.5 * (2*ACOEF/BCOEF)^(1/6); ep = BCOEF^2/(4*ACOEF)
        # to convert ACOEF and BCOEF to sigma and epsilon (GMX)
        # sigma = (ACOEF/BCOEF)^(1/6) * 0.1 ; epsilon = 4.184 * BCOEF^2/(4*ACOEF)
        #   ca   ca       819971.66        531.10
        #   ca   ha        76245.15        104.66
        #   ha   ha         5716.30         18.52

        For proper dihedrals: a quartet of atoms may appear with more than
        one set of parameters and to convert to GMX they are treated as RBs;
        use the algorithm:
          for(my $j=$i;$j<=$lines;$j++){
            my $period = $pn{$j};
            if($pk{$j}>0) {
              $V[$period] = 2*$pk{$j}*$cal;
            }
            # assign V values to C values as predefined #
            if($period==1){
              $C[0]+=0.5*$V[$period];
              if($phase{$j}==0){
                $C[1]-=0.5*$V[$period];
              }else{
                $C[1]+=0.5*$V[$period];
              }
            }elsif($period==2){
              if(($phase{$j}==180)||($phase{$j}==3.14)){
                $C[0]+=$V[$period];
                $C[2]-=$V[$period];
              }else{
                $C[2]+=$V[$period];
              }
            }elsif($period==3){
              $C[0]+=0.5*$V[$period];
              if($phase{$j}==0){
                $C[1]+=1.5*$V[$period];
                $C[3]-=2*$V[$period];
              }else{
                $C[1]-=1.5*$V[$period];
                $C[3]+=2*$V[$period];
              }
            }elsif($period==4){
              if(($phase{$j}==180)||($phase{$j}==3.14)){
                $C[2]+=4*$V[$period];
                $C[4]-=4*$V[$period];
              }else{
                $C[0]+=$V[$period];
                $C[2]-=4*$V[$period];
                $C[4]+=4*$V[$period];
              }
            }
          }
        """

        self.printMess("Writing GROMACS files\n")

        self.setAtomType4Gromacs()

        self.writeGroFile()

        self.writeGromacsTop()

        self.writeMdpFiles()

    def setAtomType4Gromacs(self):
        """Atom types names in Gromacs TOP file are not case sensitive;
        this routine will append a '_' to lower case atom type.
        E.g.: CA and ca -> CA and ca_
        """
        if self.disam:
            self.printMess("Disambiguating lower and uppercase atomtypes in GMX top file.\n")
            self.atomTypesGromacs = self.atomTypes
            self.atomsGromacs = self.atoms
            return

        atNames = [at.atomTypeName for at in self.atomTypes]
        delAtomTypes = []
        modAtomTypes = []
        atomTypesGromacs = []
        dictAtomTypes = {}
        for at in self.atomTypes:
            atName = at.atomTypeName
            dictAtomTypes[atName] = at
            if atName.islower() and atName.upper() in atNames:
                atUpper = self.atomTypes[atNames.index(atName.upper())]
                if at.ACOEF is atUpper.ACOEF and at.BCOEF is at.BCOEF:
                    delAtomTypes.append(atName)
                else:
                    newAtName = atName + "_"
                    modAtomTypes.append(atName)
                    atomType = AtomType(newAtName, at.mass, at.ACOEF, at.BCOEF)
                    atomTypesGromacs.append(atomType)
                    dictAtomTypes[newAtName] = atomType
            else:
                atomTypesGromacs.append(at)

        atomsGromacs = []
        for a in self.atoms:
            atName = a.atomType.atomTypeName
            if atName in delAtomTypes:
                atom = Atom(a.atomName, dictAtomTypes[atName.upper()], a.id, a.resid, a.mass, a.charge, a.coords)
                atom.cgnr = a.cgnr
                atomsGromacs.append(atom)
            elif atName in modAtomTypes:
                atom = Atom(a.atomName, dictAtomTypes[atName + "_"], a.id, a.resid, a.mass, a.charge, a.coords)
                atom.cgnr = a.cgnr
                atomsGromacs.append(atom)
            else:
                atomsGromacs.append(a)

        self.atomTypesGromacs = atomTypesGromacs
        self.atomsGromacs = atomsGromacs

    def writeGromacsTop(self):
        """Write GMX topology file"""
        if self.atomTypeSystem == "amber":
            d2opls = dictAtomTypeAmb2OplsGmxCode
        else:
            d2opls = dictAtomTypeGaff2OplsGmxCode

        topText = []
        itpText = []
        oitpText = []
        otopText = []
        top = self.baseName + "_GMX.top"
        itp = self.baseName + "_GMX.itp"
        otop = self.baseName + "_GMX_OPLS.top"
        oitp = self.baseName + "_GMX_OPLS.itp"

        headDefault = """
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333
"""
        headItp = """
; Include %s topology
#include "%s"
"""
        headOpls = """
; Include forcefield parameters
#include "ffoplsaa.itp"
"""
        headSystem = """
[ system ]
 %s
"""
        headMols = """
[ molecules ]
; Compound        nmols
"""
        headAtomtypes = """
[ atomtypes ]
;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb
"""
        headAtomtypesOpls = """
; For OPLS atomtypes manual fine tuning
; AC_at:OPLS_at:OPLScode: Possible_Aternatives (see ffoplsaa.atp and ffoplsaanb.itp)
"""
        headMoleculetype = """
[ moleculetype ]
;name            nrexcl
 %-16s 3
"""
        headAtoms = """
[ atoms ]
;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
"""
        headBonds = """
[ bonds ]
;   ai     aj funct   r             k
"""
        headPairs = """
[ pairs ]
;   ai     aj    funct
"""
        headAngles = """
[ angles ]
;   ai     aj     ak    funct   theta         cth
"""
        headProDih = """
[ dihedrals ] ; propers
; treated as RBs in GROMACS to use combine multiple AMBER torsions per quartet
;    i      j      k      l   func    C0         C1         C2         C3         C4         C5
"""

        headProDihAlphaGamma = """; treated as usual propers in GROMACS since Phase angle diff from 0 or 180 degrees
;    i      j      k      l   func   phase     kd      pn
"""

        headProDihGmx45 = """
[ dihedrals ] ; propers
; for gromacs 4.5 or higher, using funct 9
;    i      j      k      l   func   phase     kd      pn
"""

        headImpDih = """
[ dihedrals ] ; impropers
; treated as propers in GROMACS to use correct AMBER analytical function
;    i      j      k      l   func   phase     kd      pn
"""

        # NOTE: headTopWaterTip3p and headTopWaterSpce actually do NOTHING
        # ==============================================================================================================
        #         _headTopWaterTip3p = """
        # [ bondtypes ]
        #   ; i    j      func       b0          kb
        #   OW    HW         1    0.09572   462750.4 ; TIP3P water
        #   HW    HW         1    0.15139   462750.4 ; TIP3P water
        #
        # [ angletypes ]
        #   ;  i    j    k  func       th0       cth
        #   HW  OW  HW           1   104.520    836.800 ; TIP3P water
        #   HW  HW  OW           1   127.740      0.000 ; (found in crystallographic water with 3 bonds)
        # """
        #
        #         _headTopWaterSpce = """
        # [ bondtypes ]
        #   ; i    j      func       b0          kb
        #   OW    HW         1    0.1       462750.4 ; SPCE water
        #   HW    HW         1    0.1633    462750.4 ; SPCE water
        #
        # [ angletypes ]
        #   ;  i    j    k  func       th0       cth
        #   HW  OW  HW           1   109.47      836.800 ; SPCE water
        #   HW  HW  OW           1   125.265     0.000 ; SPCE water
        # """
        # ==============================================================================================================

        headNa = """
[ moleculetype ]
  ; molname       nrexcl
  NA+             1

[ atoms ]
  ; id_    at type res nr  residu name     at name  cg nr  charge   mass
    1       %s      1          NA+         NA+       1      1     22.9898
"""
        headCl = """
[ moleculetype ]
  ; molname       nrexcl
  CL-             1

[ atoms ]
  ; id_    at type res nr  residu name     at name  cg nr  charge   mass
    1       %s      1         CL-           CL-      1     -1     35.45300
"""
        headK = """
[ moleculetype ]
  ; molname       nrexcl
  K+             1

[ atoms ]
  ; id_    at type res nr  residu name     at name  cg nr  charge   mass
    1       %s       1          K+         K+       1      1     39.100
"""
        headWaterTip3p = """
[ moleculetype ]
; molname       nrexcl ; TIP3P model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.834   16.00000
     2     HW      1     WAT    H1      1      0.417    1.00800
     3     HW      1     WAT    H2      1      0.417    1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.09572   462750.4 0.09572   462750.4
1   3   1   0.09572   462750.4 0.09572   462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   104.520    836.800  104.520    836.800
#else
[ settles ]
; i j   funct   length
1   1   0.09572 0.15139

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
"""

        headWaterSpce = """
[ moleculetype ]
; molname       nrexcl ; SPCE model
  WAT             2

[ atoms ]
;   nr   type  resnr residue  atom   cgnr     charge       mass
     1     OW      1     WAT     O      1     -0.8476  15.99940
     2     HW      1     WAT    H1      1      0.4238   1.00800
     3     HW      1     WAT    H2      1      0.4238   1.00800

#ifdef FLEXIBLE
[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.1 462750.4  0.1     462750.4
1   3   1   0.1 462750.4  0.1     462750.4

[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   109.47  836.800 109.47  836.800
#else
[ settles ]
; OW    funct   doh dhh
1   1   0.1 0.16330

[ exclusions ]
1   2   3
2   1   3
3   1   2
#endif
"""
        if self.direct and self.amb2gmx:
            self.printMess("Converting directly from AMBER to GROMACS (EXPERIMENTAL).\n")

        # Dict of ions dealt by acpype emulating amb2gmx
        ionsDict = {"Na+": headNa, "Cl-": headCl, "K+": headK}
        ionsSorted = []
        # NOTE: headWaterTip3p and headWaterSpce actually do the real thing
        #      so, skipping headTopWaterTip3p and headWaterTip3p
        # headTopWater = headTopWaterTip3p
        headWater = headWaterTip3p

        nWat = 0
        # topFile.write("; " + head % (top, date))
        topText.append("; " + head % (top, date))
        otopText.append("; " + head % (otop, date))
        # topFile.write(headDefault)
        topText.append(headDefault)

        nSolute = 0
        if not self.amb2gmx:
            topText.append(headItp % (itp, itp))
            otopText.append(headOpls)
            otopText.append(headItp % (itp, itp))
            itpText.append("; " + head % (itp, date))
            oitpText.append("; " + head % (oitp, date))

        self.printDebug("atomTypes %i" % len(self.atomTypesGromacs))
        temp = []
        otemp = []
        for aType in self.atomTypesGromacs:
            aTypeName = aType.atomTypeName
            oaCode = d2opls.get(aTypeName, ["x", "0"])[:-1]
            aTypeNameOpls = oplsCode2AtomTypeDict.get(oaCode[0], "x")
            A = aType.ACOEF
            B = aType.BCOEF
            # one cannot infer sigma or epsilon for B = 0, assuming 0 for them
            if B == 0.0:
                sigma, epsilon, r0, epAmber = 0, 0, 0, 0
            else:
                r0 = 0.5 * math.pow((2 * A / B), (1.0 / 6))
                epAmber = 0.25 * B * B / A
                sigma = 0.1 * math.pow((A / B), (1.0 / 6))
                epsilon = cal * epAmber
            if aTypeName == "OW":
                if A == 629362.166 and B == 625.267765:
                    # headTopWater = headTopWaterSpce
                    headWater = headWaterSpce
            # OW 629362.166 625.267765 spce
            # OW 581935.564 594.825035 tip3p
            #       print aTypeName, A, B
            line = " %-8s %-11s %3.5f  %3.5f   A   %13.5e %13.5e" % (
                aTypeName,
                aTypeName,
                0.0,
                0.0,
                sigma,
                epsilon,
            ) + " ; %4.2f  %1.4f\n" % (r0, epAmber)
            oline = "; %s:%s:opls_%s: %s\n" % (aTypeName, aTypeNameOpls, oaCode[0], repr(oaCode[1:]))
            # tmpFile.write(line)
            temp.append(line)
            otemp.append(oline)
        if self.amb2gmx:
            topText.append(headAtomtypes)
            topText += temp
            nWat = self.residueLabel.count("WAT")
            for ion in ionsDict:
                nIon = self.residueLabel.count(ion)
                if nIon > 0:
                    idIon = self.residueLabel.index(ion)
                    ionType = self.search(name=ion).atomType.atomTypeName
                    ionsSorted.append((idIon, nIon, ion, ionType))
            ionsSorted.sort()
        else:
            itpText.append(headAtomtypes)
            itpText += temp
            oitpText.append(headAtomtypesOpls)
            oitpText += otemp
        self.printDebug("GMX atomtypes done")

        if len(self.atoms) > 3 * nWat + sum([x[1] for x in ionsSorted]):
            nSolute = 1

        if nWat:
            # topText.append(headTopWater)
            self.printDebug("type of water '%s'" % headWater[43:48].strip())

        if nSolute:
            if self.amb2gmx:
                topText.append(headMoleculetype % self.baseName)
            else:
                itpText.append(headMoleculetype % self.baseName)
                oitpText.append(headMoleculetype % self.baseName)

        self.printDebug("atoms %i" % len(self.atoms))
        qtot = 0.0
        count = 1
        temp = []
        otemp = []
        id2oplsATDict = {}
        for atom in self.atomsGromacs:
            resid = atom.resid
            resname = self.residueLabel[resid]
            if not self.direct:
                if resname in list(ionsDict) + ["WAT"]:
                    break
            aName = atom.atomName
            aType = atom.atomType.atomTypeName
            oItem = d2opls.get(aType, ["x", 0])
            oplsAtName = oplsCode2AtomTypeDict.get(oItem[0], "x")
            id_ = atom.id
            id2oplsATDict[id_] = oplsAtName
            oaCode = "opls_" + oItem[0]
            cgnr = id_
            if self.sorted:
                cgnr = atom.cgnr  # JDC
            charge = atom.charge
            mass = atom.mass
            omass = float(oItem[-1])
            qtot += charge
            resnr = resid + 1
            line = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot %1.3f\n" % (
                id_,
                aType,
                resnr,
                resname,
                aName,
                cgnr,
                charge,
                mass,
                qtot,
            )  # JDC
            oline = "%6d %4s %5d %5s %5s %4d %12.6f %12.5f ; qtot % 3.3f  %-4s\n" % (
                id_,
                oaCode,
                resnr,
                resname,
                aName,
                cgnr,
                charge,
                omass,
                qtot,
                oplsAtName,
            )  # JDC
            count += 1
            temp.append(line)
            otemp.append(oline)
        if temp:
            if self.amb2gmx:
                topText.append(headAtoms)
                topText += temp
            else:
                itpText.append(headAtoms)
                itpText += temp
                oitpText.append(headAtoms)
                oitpText += otemp
        self.printDebug("GMX atoms done")

        # remove bond of water
        self.printDebug("bonds %i" % len(self.bonds))
        temp = []
        otemp = []
        for bond in self.bonds:
            res1 = self.residueLabel[bond.atoms[0].resid]
            res2 = self.residueLabel[bond.atoms[0].resid]
            if "WAT" in [res1, res2]:
                continue
            a1Name = bond.atoms[0].atomName
            a2Name = bond.atoms[1].atomName
            id1 = bond.atoms[0].id
            id2 = bond.atoms[1].id
            oat1 = id2oplsATDict.get(id1)
            oat2 = id2oplsATDict.get(id2)
            line = "%6i %6i %3i %13.4e %13.4e ; %6s - %-6s\n" % (
                id1,
                id2,
                1,
                bond.rEq * 0.1,
                bond.kBond * 200 * cal,
                a1Name,
                a2Name,
            )
            oline = "%6i %6i %3i ; %13.4e %13.4e ; %6s - %-6s %6s - %-6s\n" % (
                id1,
                id2,
                1,
                bond.rEq * 0.1,
                bond.kBond * 200 * cal,
                a1Name,
                a2Name,
                oat1,
                oat2,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                topText.append(headBonds)
                topText += temp
            else:
                itpText.append(headBonds)
                itpText += temp
                oitpText.append(headBonds)
                oitpText += otemp
        self.printDebug("GMX bonds done")

        self.printDebug("atomPairs %i" % len(self.atomPairs))
        temp = []
        for pair in self.atomPairs:
            # if not printed:
            #    tmpFile.write(headPairs)
            #    printed = True
            a1Name = pair[0].atomName
            a2Name = pair[1].atomName
            id1 = pair[0].id
            id2 = pair[1].id
            # id1 = self.atoms.index(pair[0]) + 1
            # id2 = self.atoms.index(pair[1]) + 1
            line = "%6i %6i %6i ; %6s - %-6s\n" % (id1, id2, 1, a1Name, a2Name)
            temp.append(line)
        temp.sort()
        if temp:
            if self.amb2gmx:
                topText.append(headPairs)
                topText += temp
            else:
                itpText.append(headPairs)
                itpText += temp
                oitpText.append(headPairs)
                oitpText += temp
        self.printDebug("GMX pairs done")

        self.printDebug("angles %i" % len(self.angles))
        temp = []
        otemp = []
        for angle in self.angles:
            a1 = angle.atoms[0].atomName
            a2 = angle.atoms[1].atomName
            a3 = angle.atoms[2].atomName
            id1 = angle.atoms[0].id
            id2 = angle.atoms[1].id
            id3 = angle.atoms[2].id
            oat1 = id2oplsATDict.get(id1)
            oat2 = id2oplsATDict.get(id2)
            oat3 = id2oplsATDict.get(id3)
            line = "%6i %6i %6i %6i %13.4e %13.4e ; %6s - %-6s - %-6s\n" % (
                id1,
                id2,
                id3,
                1,
                angle.thetaEq * radPi,
                2 * cal * angle.kTheta,
                a1,
                a2,
                a3,
            )
            oline = "%6i %6i %6i %6i ; %13.4e %13.4e ; %6s - %-4s - %-6s %4s - %+4s - %-4s\n" % (
                id1,
                id2,
                id3,
                1,
                angle.thetaEq * radPi,
                2 * cal * angle.kTheta,
                a1,
                a2,
                a3,
                oat1,
                oat2,
                oat3,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                topText.append(headAngles)
                topText += temp
            else:
                itpText.append(headAngles)
                itpText += temp
                oitpText.append(headAngles)
                oitpText += otemp
        self.printDebug("GMX angles done")

        self.setProperDihedralsCoef()
        self.printDebug("properDihedralsCoefRB %i" % len(self.properDihedralsCoefRB))
        self.printDebug("properDihedralsAlphaGamma %i" % len(self.properDihedralsAlphaGamma))
        self.printDebug("properDihedralsGmx45 %i" % len(self.properDihedralsGmx45))
        temp = []
        otemp = []
        if self.gmx4:
            self.printMess("Writing RB dihedrals for old GMX 4.\n")
            for dih in self.properDihedralsCoefRB:
                a1 = dih[0][0].atomName
                a2 = dih[0][1].atomName
                a3 = dih[0][2].atomName
                a4 = dih[0][3].atomName
                id1 = dih[0][0].id
                id2 = dih[0][1].id
                id3 = dih[0][2].id
                id4 = dih[0][3].id
                oat1 = id2oplsATDict.get(id1)
                oat2 = id2oplsATDict.get(id2)
                oat3 = id2oplsATDict.get(id3)
                oat4 = id2oplsATDict.get(id4)
                c0, c1, c2, c3, c4, c5 = dih[1]
                line = "%6i %6i %6i %6i %6i %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    3,
                    c0,
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                ) + " ; %6s-%6s-%6s-%6s\n" % (a1, a2, a3, a4)
                oline = "%6i %6i %6i %6i %6i ; %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    3,
                    c0,
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                ) + " ; %6s-%6s-%6s-%6s    %4s-%4s-%4s-%4s\n" % (a1, a2, a3, a4, oat1, oat2, oat3, oat4)
                temp.append(line)
                otemp.append(oline)
            temp.sort()
            otemp.sort()
            if temp:
                if self.amb2gmx:
                    topText.append(headProDih)
                    topText += temp
                else:
                    itpText.append(headProDih)
                    itpText += temp
                    oitpText.append(headProDih)
                    oitpText += otemp
            self.printDebug("GMX proper dihedrals done")
        else:
            self.printMess("Writing GMX dihedrals for GMX 4.5 and higher.\n")
            funct = 9  # 9
            for dih in self.properDihedralsGmx45:
                a1 = dih[0][0].atomName
                a2 = dih[0][1].atomName
                a3 = dih[0][2].atomName
                a4 = dih[0][3].atomName
                id1 = dih[0][0].id
                id2 = dih[0][1].id
                id3 = dih[0][2].id
                id4 = dih[0][3].id
                ph = dih[1]  # phase already in degree
                kd = dih[2] * cal  # kPhi PK
                pn = dih[3]  # .period
                line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    funct,
                    ph,
                    kd,
                    pn,
                    a1,
                    a2,
                    a3,
                    a4,
                )
                oline = "%6i %6i %6i %6i %6i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                    id1,
                    id2,
                    id3,
                    id4,
                    funct,
                    ph,
                    kd,
                    pn,
                    a1,
                    a2,
                    a3,
                    a4,
                )
                temp.append(line)
                otemp.append(oline)
            temp.sort()
            otemp.sort()
            if temp:
                if self.amb2gmx:
                    topText.append(headProDihGmx45)
                    topText += temp
                else:
                    itpText.append(headProDihGmx45)
                    itpText += temp
                    oitpText.append(headProDihGmx45)
                    oitpText += otemp

        # for properDihedralsAlphaGamma
        if not self.gmx4:
            funct = 4  # 4
        else:
            funct = 1
        temp = []
        otemp = []
        for dih in self.properDihedralsAlphaGamma:
            a1 = dih[0][0].atomName
            a2 = dih[0][1].atomName
            a3 = dih[0][2].atomName
            a4 = dih[0][3].atomName
            id1 = dih[0][0].id
            id2 = dih[0][1].id
            id3 = dih[0][2].id
            id4 = dih[0][3].id
            ph = dih[1]  # phase already in degree
            kd = dih[2] * cal  # kPhi PK
            pn = dih[3]  # .period
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            oline = "%6i %6i %6i %6i %6i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                topText.append(headProDihAlphaGamma)
                topText += temp
            else:
                itpText.append(headProDihAlphaGamma)
                itpText += temp
                oitpText.append(headProDihAlphaGamma)
                oitpText += otemp
        self.printDebug("GMX special proper dihedrals done")

        self.printDebug("improperDihedrals %i" % len(self.improperDihedrals))
        temp = []
        otemp = []
        for dih in self.improperDihedrals:
            a1 = dih.atoms[0].atomName
            a2 = dih.atoms[1].atomName
            a3 = dih.atoms[2].atomName
            a4 = dih.atoms[3].atomName
            id1 = dih.atoms[0].id
            id2 = dih.atoms[1].id
            id3 = dih.atoms[2].id
            id4 = dih.atoms[3].id
            kd = dih.kPhi * cal
            pn = dih.period
            ph = dih.phase * radPi
            line = "%6i %6i %6i %6i %6i %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            oline = "%6i %6i %6i %6i %6i ; %8.2f %9.5f %3i ; %6s-%6s-%6s-%6s\n" % (
                id1,
                id2,
                id3,
                id4,
                funct,
                ph,
                kd,
                pn,
                a1,
                a2,
                a3,
                a4,
            )
            temp.append(line)
            otemp.append(oline)
        temp.sort()
        otemp.sort()
        if temp:
            if self.amb2gmx:
                topText.append(headImpDih)
                topText += temp
            else:
                itpText.append(headImpDih)
                itpText += temp
                oitpText.append(headImpDih)
                oitpText += otemp
        self.printDebug("GMX improper dihedrals done")

        if not self.direct:
            for ion in ionsSorted:
                topText.append(ionsDict[ion[2]] % ion[3])

            if nWat:
                topText.append(headWater)

        topText.append(headSystem % (self.baseName))
        topText.append(headMols)
        otopText.append(headSystem % (self.baseName))
        otopText.append(headMols)

        if nSolute > 0:
            topText.append(" %-16s %-6i\n" % (self.baseName, nSolute))
            otopText.append(" %-16s %-6i\n" % (self.baseName, nSolute))

        if not self.direct:
            for ion in ionsSorted:
                topText.append(" %-16s %-6i\n" % (ion[2].upper(), ion[1]))

            if nWat:
                topText.append(" %-16s %-6i\n" % ("WAT", nWat))

        if self.topo14Data.hasNondefault14():
            citation = (
                "     BERNARDI, A., FALLER, R., REITH, D., and KIRSCHNER, K. N. ACPYPE update for\n"
                + "     nonuniform 1-4 scale factors: Conversion of the GLYCAM06 force field from AMBER\n"
                + '     to GROMACS. SoftwareX 10 (2019), 100241. doi: 10.1016/j.softx.2019.100241"\n'
            )

            msg = "Non-default 1-4 scale parameters detected.  Converting individually. Please cite:\n\n" + citation

            self.printMess(msg)
            topText = self.topo14Data.patch_gmx_topol14("".join(topText))

        gmxDir = os.path.abspath(".")
        topFileName = os.path.join(gmxDir, top)
        topFile = open(topFileName, "w")
        topFile.writelines(topText)
        self.topText = topText

        if not self.amb2gmx:
            itpFileName = os.path.join(gmxDir, itp)
            itpFile = open(itpFileName, "w")
            itpFile.writelines(itpText)
            oitpFileName = os.path.join(gmxDir, oitp)
            oitpFile = open(oitpFileName, "w")
            oitpFile.writelines(oitpText)
            otopFileName = os.path.join(gmxDir, otop)
            otopFile = open(otopFileName, "w")
            otopFile.writelines(otopText)

    def writeGroFile(self):
        """Write GRO files"""
        # print "Writing GROMACS GRO file\n"
        self.printDebug("writing GRO file")
        gro = self.baseName + "_GMX.gro"
        gmxDir = os.path.abspath(".")
        groFileName = os.path.join(gmxDir, gro)
        groFile = open(groFileName, "w")
        groFile.write(head % (gro, date))
        groFile.write(" %i\n" % len(self.atoms))
        count = 1
        for atom in self.atoms:
            coords = [c * 0.1 for c in atom.coords]
            resid = atom.resid
            line = "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" % (
                resid + 1,
                self.residueLabel[resid],
                atom.atomName,
                count,
                coords[0],
                coords[1],
                coords[2],
            )
            count += 1
            if count == 100000:
                count = 0
            groFile.write(line)
        if self.pbc:
            boxX = self.pbc[0][0] * 0.1
            boxY = self.pbc[0][1] * 0.1
            boxZ = self.pbc[0][2] * 0.1
            vX = self.pbc[1][0]
            # vY = self.pbc[1][1]
            # vZ = self.pbc[1][2]
            if vX == 90.0:
                self.printDebug("PBC triclinic")
                text = "%11.5f %11.5f %11.5f\n" % (boxX, boxY, boxZ)
            elif round(vX, 2) == 109.47:
                self.printDebug("PBC octahedron")
                f1 = 0.471405  # 1/3 * sqrt(2)
                f2 = 0.333333 * boxX
                v22 = boxY * 2 * f1
                v33 = boxZ * f1 * 1.73205  # f1 * sqrt(3)
                v21 = v31 = v32 = 0.0
                v12 = f2
                v13 = -f2
                v23 = f1 * boxX
                text = "%11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n" % (
                    boxX,
                    v22,
                    v33,
                    v21,
                    v31,
                    v12,
                    v32,
                    v13,
                    v23,
                )
        else:
            self.printDebug("Box size estimated")
            X = [a.coords[0] * 0.1 for a in self.atoms]
            Y = [a.coords[1] * 0.1 for a in self.atoms]
            Z = [a.coords[2] * 0.1 for a in self.atoms]
            boxX = max(X) - min(X)  # + 2.0 # 2.0 is double of rlist
            boxY = max(Y) - min(Y)  # + 2.0
            boxZ = max(Z) - min(Z)  # + 2.0
            text = "%11.5f %11.5f %11.5f\n" % (boxX * 20.0, boxY * 20.0, boxZ * 20.0)
        groFile.write(text)

    def writeMdpFiles(self):
        """Write MDP for test with GROMACS"""
        emMdp = f"""; to test
; gmx grompp -f em.mdp -c {self.baseName}_GMX.gro -p {self.baseName}_GMX.top -o em.tpr -v
; gmx mdrun -ntmpi 1 -v -deffnm em
integrator               = steep
nsteps                   = 500
"""
        mdMdp = f"""; to test
; gmx grompp -f md.mdp -c em.gro -p {self.baseName}_GMX.top -o md.tpr
; gmx mdrun -ntmpi 1 -v -deffnm md
integrator               = md
nsteps                   = 10000
nstxout                  = 10
; vmd md.gro md.trr
"""
        rungmx = f"""
gmx grompp -f em.mdp -c {self.baseName}_GMX.gro -p {self.baseName}_GMX.top -o em.tpr -v
gmx mdrun -ntmpi 1 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p {self.baseName}_GMX.top -o md.tpr
gmx mdrun -ntmpi 1 -v -deffnm md
"""
        emMdpFile = open("em.mdp", "w")
        mdMdpFile = open("md.mdp", "w")
        runGmxFile = open("rungmx.sh", "w")
        emMdpFile.write(emMdp)
        mdMdpFile.write(mdMdp)
        runGmxFile.write(rungmx)
        os.chmod("rungmx.sh", 0o744)

    def writeCnsTopolFiles(self):
        """Write CNS topology files"""
        autoAngleFlag = True
        autoDihFlag = True
        cnsDir = os.path.abspath(".")

        pdb = self.baseName + "_NEW.pdb"
        par = self.baseName + "_CNS.par"
        top = self.baseName + "_CNS.top"
        inp = self.baseName + "_CNS.inp"

        pdbFileName = os.path.join(cnsDir, pdb)
        parFileName = os.path.join(cnsDir, par)
        topFileName = os.path.join(cnsDir, top)
        inpFileName = os.path.join(cnsDir, inp)

        self.CnsTopFileName = topFileName
        self.CnsInpFileName = inpFileName
        self.CnsParFileName = parFileName
        self.CnsPdbFileName = pdbFileName

        parFile = open(parFileName, "w")
        topFile = open(topFileName, "w")
        inpFile = open(inpFileName, "w")

        self.printMess("Writing NEW PDB file\n")
        self.writePdb(pdbFileName)

        self.printMess("Writing CNS/XPLOR files\n")

        # print "Writing CNS PAR file\n"
        parFile.write("Remarks " + head % (par, date))
        parFile.write("\nset echo=false end\n")

        parFile.write("\n{ Bonds: atomType1 atomType2 kb r0 }\n")
        lineSet = []
        for bond in self.bonds:
            a1Type = bond.atoms[0].atomType.atomTypeName + "_"
            a2Type = bond.atoms[1].atomType.atomTypeName + "_"
            kb = 1000.0
            if not self.allhdg:
                kb = bond.kBond
            r0 = bond.rEq
            line = "BOND %5s %5s %8.1f %8.4f\n" % (a1Type, a2Type, kb, r0)
            lineRev = "BOND %5s %5s %8.1f %8.4f\n" % (a2Type, a1Type, kb, r0)
            if line not in lineSet:
                if lineRev not in lineSet:
                    lineSet.append(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Angles: aType1 aType2 aType3 kt t0 }\n")
        lineSet = []
        for angle in self.angles:
            a1 = angle.atoms[0].atomType.atomTypeName + "_"
            a2 = angle.atoms[1].atomType.atomTypeName + "_"
            a3 = angle.atoms[2].atomType.atomTypeName + "_"
            kt = 500.0
            if not self.allhdg:
                kt = angle.kTheta
            t0 = angle.thetaEq * radPi
            line = "ANGLe %5s %5s %5s %8.1f %8.2f\n" % (a1, a2, a3, kt, t0)
            lineRev = "ANGLe %5s %5s %5s %8.1f %8.2f\n" % (a3, a2, a1, kt, t0)
            if line not in lineSet:
                if lineRev not in lineSet:
                    lineSet.append(line)
        for item in lineSet:
            parFile.write(item)

        parFile.write(
            "\n{ Proper Dihedrals: aType1 aType2 aType3 aType4 kt per\
iod phase }\n"
        )
        lineSet = set()
        for item in self.condensedProperDihedrals:
            seq = ""
            id_ = 0
            for dih in item:
                # id_ = item.index(dih)
                ll = len(item)
                a1 = dih.atoms[0].atomType.atomTypeName + "_"
                a2 = dih.atoms[1].atomType.atomTypeName + "_"
                a3 = dih.atoms[2].atomType.atomTypeName + "_"
                a4 = dih.atoms[3].atomType.atomTypeName + "_"
                kp = 750.0
                if not self.allhdg:
                    kp = dih.kPhi
                p = dih.period
                ph = dih.phase * radPi
                if ll > 1:
                    if id_ == 0:
                        line = (
                            "DIHEdral %5s %5s %5s %5s  MULT %1i %7.3f %4i %8\
.2f\n"
                            % (a1, a2, a3, a4, ll, kp, p, ph)
                        )
                    else:
                        line = "%s %7.3f %4i %8.2f\n" % (40 * " ", kp, p, ph)
                else:
                    line = "DIHEdral %5s %5s %5s %5s %15.3f %4i %8.2f\n" % (a1, a2, a3, a4, kp, p, ph)
                seq += line
                id_ += 1
            lineSet.add(seq)
        for item in lineSet:
            parFile.write(item)

        parFile.write(
            "\n{ Improper Dihedrals: aType1 aType2 aType3 aType4 kt p\
eriod phase }\n"
        )
        lineSet = set()
        for idh in self.improperDihedrals:
            a1 = idh.atoms[0].atomType.atomTypeName + "_"
            a2 = idh.atoms[1].atomType.atomTypeName + "_"
            a3 = idh.atoms[2].atomType.atomTypeName + "_"
            a4 = idh.atoms[3].atomType.atomTypeName + "_"
            kp = 750.0
            if not self.allhdg:
                kp = idh.kPhi
            p = idh.period
            ph = idh.phase * radPi
            line = "IMPRoper %5s %5s %5s %5s %13.1f %4i %8.2f\n" % (a1, a2, a3, a4, kp, p, ph)
            lineSet.add(line)
        if self.chiral:
            for idhc in self.chiralGroups:
                _atc, neig, angle = idhc
                a1 = neig[0].atomType.atomTypeName + "_"
                a2 = neig[1].atomType.atomTypeName + "_"
                a3 = neig[2].atomType.atomTypeName + "_"
                a4 = neig[3].atomType.atomTypeName + "_"
                kp = 11000.0
                p = 0
                ph = angle
                line = "IMPRoper %5s %5s %5s %5s %13.1f %4i %8.2f\n" % (a1, a2, a3, a4, kp, p, ph)
                lineSet.add(line)

        for item in lineSet:
            parFile.write(item)

        parFile.write("\n{ Nonbonded: Type Emin sigma; (1-4): Emin/2 sigma }\n")
        for at in self.atomTypes:
            A = at.ACOEF
            B = at.BCOEF
            atName = at.atomTypeName + "_"
            if B == 0.0:
                sigma = epAmber = ep2 = sig2 = 0.0
            else:
                epAmber = 0.25 * B * B / A
                ep2 = epAmber / 2.0
                sigma = math.pow((A / B), (1.0 / 6))
                sig2 = sigma
            line = "NONBonded %5s %11.6f %11.6f %11.6f %11.6f\n" % (atName, epAmber, sigma, ep2, sig2)
            parFile.write(line)
        parFile.write("\nset echo=true end\n")

        # print "Writing CNS TOP file\n"
        topFile.write("Remarks " + head % (top, date))
        topFile.write("\nset echo=false end\n")
        topFile.write("\nautogenerate angles=%s dihedrals=%s end\n" % (autoAngleFlag, autoDihFlag))
        topFile.write("\n{ atomType  mass }\n")

        for at in self.atomTypes:
            atType = at.atomTypeName + "_"
            mass = at.mass
            line = "MASS %-5s %8.3f\n" % (atType, mass)
            topFile.write(line)

        topFile.write("\nRESIdue %s\n" % self.residueLabel[0])
        topFile.write("\nGROUP\n")

        topFile.write("\n{ atomName  atomType  Charge }\n")

        for at in self.atoms:
            atName = at.atomName
            atType = at.atomType.atomTypeName + "_"
            charge = at.charge
            line = "ATOM %-5s TYPE= %-5s CHARGE= %8.4f END\n" % (atName, atType, charge)
            topFile.write(line)

        topFile.write("\n{ Bonds: atomName1  atomName2 }\n")
        for bond in self.bonds:
            a1Name = bond.atoms[0].atomName
            a2Name = bond.atoms[1].atomName
            line = "BOND %-5s %-5s\n" % (a1Name, a2Name)
            topFile.write(line)

        if not autoAngleFlag or 1:  # generating angles anyway
            topFile.write("\n{ Angles: atomName1 atomName2 atomName3}\n")
            for angle in self.angles:
                a1Name = angle.atoms[0].atomName
                a2Name = angle.atoms[1].atomName
                a3Name = angle.atoms[2].atomName
                line = "ANGLe %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name)
                topFile.write(line)

        if not autoDihFlag or 1:  # generating angles anyway
            topFile.write("\n{ Proper Dihedrals: name1 name2 name3 name4 }\n")
            for item in self.condensedProperDihedrals:
                for dih in item:
                    a1Name = dih.atoms[0].atomName
                    a2Name = dih.atoms[1].atomName
                    a3Name = dih.atoms[2].atomName
                    a4Name = dih.atoms[3].atomName
                    line = "DIHEdral %-5s %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name, a4Name)
                    break
                topFile.write(line)

        topFile.write("\n{ Improper Dihedrals: aName1 aName2 aName3 aName4 }\n")
        for dih in self.improperDihedrals:
            a1Name = dih.atoms[0].atomName
            a2Name = dih.atoms[1].atomName
            a3Name = dih.atoms[2].atomName
            a4Name = dih.atoms[3].atomName
            line = "IMPRoper %-5s %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name, a4Name)
            topFile.write(line)

        if self.chiral:
            for idhc in self.chiralGroups:
                _atc, neig, angle = idhc
                a1Name = neig[0].atomName
                a2Name = neig[1].atomName
                a3Name = neig[2].atomName
                a4Name = neig[3].atomName
                line = "IMPRoper %-5s %-5s %-5s %-5s\n" % (a1Name, a2Name, a3Name, a4Name)
                topFile.write(line)

        topFile.write("\nEND {RESIdue %s}\n" % self.residueLabel[0])

        topFile.write("\nset echo=true end\n")

        inpFile.write("Remarks " + head % (inp, date))
        inpData = """
topology
  @%(CNS_top)s
end

parameters
  @%(CNS_par)s
  nbonds
      atom cdie shift eps=1.0  e14fac=0.4   tolerance=0.5
      cutnb=9.0 ctonnb=7.5 ctofnb=8.0
      nbxmod=5 vswitch wmin 1.0
  end
  remark dielectric constant eps set to 1.0
end

flags exclude elec ? end

segment name="    "
  chain
   coordinates @%(NEW_pdb)s
  end
end
coordinates @%(NEW_pdb)s
coord copy end

! Remarks If you want to shake up the coordinates a bit ...
 vector do (x=x+6*(rand()-0.5)) (all)
 vector do (y=y+6*(rand()-0.5)) (all)
 vector do (z=z+6*(rand()-0.5)) (all)
 write coordinates output=%(CNS_ran)s end

! Remarks RMS diff after randomisation and before minimisation
coord rms sele=(known and not hydrogen) end

print threshold=0.02 bonds
print threshold=3.0 angles
print threshold=3.0 dihedrals
print threshold=3.0 impropers

! Remarks Do Powell energy minimisation
minimise powell
  nstep=250 drop=40.0
end

write coordinates output=%(CNS_min)s end
write structure   output=%(CNS_psf)s end

! constraints interaction (not hydro) (not hydro) end

print threshold=0.02 bonds
print threshold=3.0 angles
print threshold=3.0 dihedrals
print threshold=3.0 impropers

flags exclude * include vdw end energy end
distance from=(not hydro) to=(not hydro) cutoff=2.6 end

! Remarks RMS fit after minimisation
coord fit sele=(known and not hydrogen) end

stop
"""
        dictInp = {}
        dictInp["CNS_top"] = top
        dictInp["CNS_par"] = par
        dictInp["NEW_pdb"] = pdb
        dictInp["CNS_min"] = self.baseName + "_NEW_min.pdb"
        dictInp["CNS_psf"] = self.baseName + "_CNS.psf"
        dictInp["CNS_ran"] = self.baseName + "_rand.pdb"
        line = inpData % dictInp
        inpFile.write(line)
        self.printDebug("chiralGroups %i" % len(self.chiralGroups))


class ACTopol(AbstractTopol):

    """
    Class to build the AC topologies (Antechamber AmberTools)
    """

    def __init__(
        self,
        inputFile,
        chargeType="bcc",
        chargeVal=None,
        multiplicity="1",
        atomType="gaff",
        force=False,
        basename=None,
        debug=False,
        outTopol="all",
        allhdg=False,
        timeTol=MAXTIME,
        qprog="sqm",
        ekFlag=None,
        verbose=True,
        gmx4=False,
        disam=False,
        direct=False,
        is_sorted=False,
        chiral=False,
        amb2gmx=False,
    ):
        super().__init__()
        self.amb2gmx = amb2gmx
        self.debug = debug
        self.verbose = verbose
        self.gmx4 = gmx4
        self.disam = disam
        self.direct = direct
        self.sorted = is_sorted
        self.chiral = chiral
        self.babelExe = which("obabel") or which("babel") or ""
        if not os.path.exists(self.babelExe):
            if self.ext != ".mol2" and self.ext != ".mdl":
                self.printError("no 'babel' executable; you need it if input is PDB or SMILES")
                self.printError("otherwise use only MOL2 or MDL file as input ... aborting!")
                raise Exception("Missing BABEL")
            else:
                self.printWarn("no 'babel' executable, no PDB file as input can be used!")
        self.inputFile = os.path.basename(inputFile)
        self.rootDir = os.path.abspath(".")
        self.absInputFile = os.path.abspath(inputFile)
        if not os.path.exists(self.absInputFile) and checkSmiles(inputFile):
            self.is_smiles = True
            self.smiles = inputFile
            if not basename:
                self.inputFile = "smiles_molecule.mol2"
            else:
                self.inputFile = f"{basename}.mol2"
            self.absInputFile = os.path.abspath(self.inputFile)
        elif not os.path.exists(self.absInputFile):
            raise Exception(f"Input file {inputFile} DOES NOT EXIST")
        baseOriginal, ext = os.path.splitext(self.inputFile)
        base = basename or baseOriginal
        self.baseOriginal = baseOriginal
        self.ext = ext
        self.baseName = base  # name of the input file without ext.
        if self.is_smiles:
            self.convertSmilesToMol2()
        self.timeTol = timeTol
        self.printDebug("Max execution time tolerance is %s" % elapsedTime(self.timeTol))
        if ekFlag == '"None"' or ekFlag is None:
            self.ekFlag = ""
        else:
            self.ekFlag = "-ek %s" % ekFlag
        self.extOld = ext
        self.homeDir = self.baseName + ".acpype"
        self.chargeType = chargeType
        self.chargeVal = chargeVal
        self.multiplicity = multiplicity
        self.atomType = atomType
        self.gaffDatfile = "gaff.dat"
        leapGaffFile = "leaprc.gaff"
        if "2" in self.atomType:
            leapGaffFile = "leaprc.gaff2"
            self.gaffDatfile = "gaff2.dat"
        self.force = force
        self.allhdg = allhdg
        self.acExe = ""
        dirAmber = os.getenv("AMBERHOME", os.getenv("ACHOME"))
        if dirAmber:
            for ac_bin in ["bin", "exe"]:
                ac_path = os.path.join(dirAmber, ac_bin, "antechamber")
                if os.path.exists(ac_path):
                    self.acExe = ac_path
                    break
        if not self.acExe:
            self.acExe = which("antechamber") or ""
        if not os.path.exists(self.acExe):
            self.printError("no 'antechamber' executable... aborting ! ")
            hint1 = "HINT1: is 'AMBERHOME' or 'ACHOME' environment variable set?"
            hint2 = (
                "HINT2: is 'antechamber' in your $PATH?"
                + "\n    What 'which antechamber' in your terminal says?"
                + "\n    'alias' doesn't work for ACPYPE."
            )
            self.printMess(hint1)
            self.printMess(hint2)
            raise Exception("Missing ANTECHAMBER")
        self.tleapExe = which("tleap") or ""
        self.parmchkExe = which("parmchk2") or ""
        acBase = base + "_AC"
        self.acBaseName = acBase
        self.acXyzFileName = acBase + ".inpcrd"
        self.acTopFileName = acBase + ".prmtop"
        self.acFrcmodFileName = acBase + ".frcmod"
        self.tmpDir = os.path.join(self.rootDir, ".acpype_tmp_%s" % os.path.basename(base))
        self.setResNameCheckCoords()
        self.guessCharge()
        acMol2FileName = "%s_%s_%s.mol2" % (base, chargeType, atomType)
        self.acMol2FileName = acMol2FileName
        self.charmmBase = "%s_CHARMM" % base
        self.qFlag = qDict[qprog]
        self.outTopols = [outTopol]
        if outTopol == "all":
            self.outTopols = outTopols
        self.acParDict = {
            "base": base,
            "ext": ext[1:],
            "acBase": acBase,
            "acMol2FileName": acMol2FileName,
            "res": self.resName,
            "leapAmberFile": leapAmberFile,
            "baseOrg": self.baseOriginal,
            "leapGaffFile": leapGaffFile,
        }


class MolTopol(AbstractTopol):

    """ "
    Class to write topologies and parameters files for several applications

    http://amber.scripps.edu/formats.html (not updated to amber 10 yet)

    Parser, take information in AC xyz and top files and convert to objects

    INPUTS: acFileXyz and acFileTop
    RETURN: molTopol obj or None
    """

    def __init__(
        self,
        acTopolObj=None,
        acFileXyz=None,
        acFileTop=None,
        debug=False,
        basename=None,
        verbose=True,
        gmx4=False,
        disam=False,
        direct=False,
        is_sorted=False,
        chiral=False,
        amb2gmx=False,
    ):
        super().__init__()
        self.amb2gmx = amb2gmx
        self.chiral = chiral
        self.allhdg = False
        self.debug = debug
        self.gmx4 = gmx4
        self.disam = disam
        self.direct = direct
        self.sorted = is_sorted
        self.verbose = verbose
        self.inputFile = acFileTop

        if acTopolObj:
            if not acFileXyz:
                acFileXyz = acTopolObj.acXyzFileName
            if not acFileTop:
                acFileTop = acTopolObj.acTopFileName
            self._parent = acTopolObj
            self.allhdg = self._parent.allhdg
            self.debug = self._parent.debug
            self.inputFile = self._parent.inputFile
        if not os.path.exists(acFileXyz) and not os.path.exists(acFileTop):
            self.printError("Files '%s' and/or '%s' don't exist" % (acFileXyz, acFileTop))
            self.printError("molTopol object won't be created")

        self.xyzFileData = open(acFileXyz, "r").readlines()
        self.topFileData = [x for x in open(acFileTop, "r").readlines() if not x.startswith("%COMMENT")]
        self.topo14Data = Topology_14()
        self.topo14Data.read_amber_topology("".join(self.topFileData))
        self.printDebug("prmtop and inpcrd files loaded")

        self.getResidueLabel()
        if len(self.residueLabel) > 1:
            self.baseName = basename or os.path.splitext(os.path.basename(acFileTop))[0]  # 'solute'
        else:
            self.baseName = basename or self.residueLabel[0]  # 3 caps letters
        if acTopolObj:
            self.baseName = basename or acTopolObj.baseName
        self.printDebug("basename defined = '%s'" % self.baseName)

        self.getAtoms()

        self.getBonds()

        self.getAngles()

        self.getDihedrals()

        if self.amb2gmx:
            self.rootDir = os.path.abspath(".")
            self.homeDir = f"{self.baseName}.amb2gmx"
            self.makeDir()
        else:
            self.getChirals()

        # Sort atoms for gromacs output. # JDC
        if self.sorted:
            self.printMess("Sorting atoms for gromacs ordering.\n")
            self.sortAtomsForGromacs()


class Atom:

    """
    Charges in prmtop file has to be divide by 18.2223 to convert to charge
    in units of the electron charge.
    To convert ACOEF and BCOEF to r0 (Ang.) and epsilon (kcal/mol), as seen
    in gaff.dat for example; same atom type (i = j):
        r0 = 1/2 * (2 * ACOEF/BCOEF)^(1/6)
        epsilon = 1/(4 * A) * BCOEF^2
    To convert r0 and epsilon to ACOEF and BCOEF
        ACOEF = sqrt(ep_i * ep_j) * (r0_i + r0_j)^12
        BCOEF = 2 * sqrt(ep_i * ep_j) * (r0_i + r0_j)^6
              = 2 * ACOEF/(r0_i + r0_j)^6
    where index i and j for atom types.
    Coord is given in Ang. and mass in Atomic Mass Unit.
    """

    def __init__(self, atomName, atomType, id_, resid, mass, charge, coord):
        self.atomName = atomName
        self.atomType = atomType
        self.id = id_
        self.cgnr = id_
        self.resid = resid
        self.mass = mass
        self.charge = charge  # / qConv
        self.coords = coord

    def __str__(self):
        return "<Atom id=%s, name=%s, %s>" % (self.id, self.atomName, self.atomType)

    def __repr__(self):
        return "<Atom id=%s, name=%s, %s>" % (self.id, self.atomName, self.atomType)


class AtomType:

    """
    AtomType per atom in gaff or amber.
    """

    def __init__(self, atomTypeName, mass, ACOEF, BCOEF):
        self.atomTypeName = atomTypeName
        self.mass = mass
        self.ACOEF = ACOEF
        self.BCOEF = BCOEF

    def __str__(self):
        return "<AtomType=%s>" % self.atomTypeName

    def __repr__(self):
        return "<AtomType=%s>" % self.atomTypeName


class Bond:

    """
    attributes: pair of Atoms, spring constant (kcal/mol), dist. eq. (Ang)
    """

    def __init__(self, atoms, kBond, rEq):
        self.atoms = atoms
        self.kBond = kBond
        self.rEq = rEq

    def __str__(self):
        return "<%s, r=%s>" % (self.atoms, self.rEq)

    def __repr__(self):
        return "<%s, r=%s>" % (self.atoms, self.rEq)


class Angle:

    """
    attributes: 3 Atoms, spring constant (kcal/mol/rad^2), angle eq. (rad)
    """

    def __init__(self, atoms, kTheta, thetaEq):
        self.atoms = atoms
        self.kTheta = kTheta
        self.thetaEq = thetaEq  # rad, to convert to degree: thetaEq * 180/Pi

    def __str__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.thetaEq * 180 / Pi)

    def __repr__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.thetaEq * 180 / Pi)


class Dihedral:

    """
    attributes: 4 Atoms, spring constant (kcal/mol), periodicity,
    phase (rad)
    """

    def __init__(self, atoms, kPhi, period, phase):
        self.atoms = atoms
        self.kPhi = kPhi
        self.period = period
        self.phase = phase  # rad, to convert to degree: kPhi * 180/Pi

    def __str__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.phase * 180 / Pi)

    def __repr__(self):
        return "<%s, ang=%.2f>" % (self.atoms, self.phase * 180 / Pi)


def get_option_parser():
    parser = argparse.ArgumentParser(usage=usage + epilog)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="input file type like '.pdb', '.mdl', '.mol2' or SMILES string (mandatory if -p and -x not set)",
    )
    parser.add_argument(
        "-b",
        "--basename",
        action="store",
        dest="basename",
        help="a basename for the project (folder and output files)",
    )
    parser.add_argument(
        "-x", "--inpcrd", action="store", dest="inpcrd", help="amber inpcrd file name (always used with -p)",
    )
    parser.add_argument(
        "-p", "--prmtop", action="store", dest="prmtop", help="amber prmtop file name (always used with -x)",
    )
    parser.add_argument(
        "-c",
        "--charge_method",
        choices=["gas", "bcc", "user"],
        action="store",
        default="bcc",
        dest="charge_method",
        help="charge method: gas, bcc (default), user (user's charges in mol2 file)",
    )
    parser.add_argument(
        "-n",
        "--net_charge",
        action="store",
        type=int,
        default=None,
        dest="net_charge",
        help="net molecular charge (int), for gas default is 0",
    )
    parser.add_argument(
        "-m",
        "--multiplicity",
        action="store",
        type=int,
        default=1,
        dest="multiplicity",
        help="multiplicity (2S+1), default is 1",
    )
    parser.add_argument(
        "-a",
        "--atom_type",
        choices=["gaff", "amber", "gaff2", "amber2"],
        action="store",
        default="gaff",
        dest="atom_type",
        help="atom type, can be gaff, gaff2, amber (AMBER14SB) or amber2 (AMBER14SB + GAFF2), default is gaff",
    )
    parser.add_argument(
        "-q",
        "--qprog",
        choices=["mopac", "sqm", "divcon"],
        action="store",
        default="sqm",
        dest="qprog",
        help="am1-bcc flag, sqm (default), divcon, mopac",
    )
    parser.add_argument(
        "-k", "--keyword", action="store", dest="keyword", help="mopac or sqm keyword, inside quotes",
    )
    parser.add_argument(
        "-f", "--force", action="store_true", dest="force", help="force topologies recalculation anew",
    )
    group.add_argument(
        "-d",
        "--debug",
        action="store_true",
        dest="debug",
        help="for debugging purposes, keep any temporary file created (not allowed with arg -w)",
    )
    parser.add_argument(
        "-o",
        "--outtop",
        choices=["all"] + outTopols,
        action="store",
        default="all",
        dest="outtop",
        help="output topologies: all (default), gmx, cns or charmm",
    )
    parser.add_argument(
        "-z", "--gmx4", action="store_true", dest="gmx4", help="write RB dihedrals old GMX 4.0",
    )
    parser.add_argument(
        "-t",
        "--cnstop",
        action="store_true",
        dest="cnstop",
        help="write CNS topology with allhdg-like parameters (experimental)",
    )
    parser.add_argument(
        "-s",
        "--max_time",
        type=int,
        action="store",
        default=MAXTIME,
        dest="max_time",
        help="max time (in sec) tolerance for sqm/mopac, default is %i hours" % (MAXTIME // 3600),
    )
    parser.add_argument(
        "-y", "--ipython", action="store_true", dest="ipython", help="start iPython interpreter",
    )
    group.add_argument(
        "-w",
        "--verboseless",
        action="store_false",
        default=True,
        dest="verboseless",
        help="print nothing (not allowed with arg -d)",
    )
    parser.add_argument(
        "-g",
        "--disambiguate",
        action="store_true",
        dest="disambiguate",
        help="disambiguate lower and uppercase atomtypes in GMX top file",
    )
    parser.add_argument(
        "-u",
        "--direct",
        action="store_true",
        dest="direct",
        help="for 'amb2gmx' mode, does a direct conversion, for any solvent (EXPERIMENTAL)",
        # NOTE: when solvent is present, gmx mdrun is not working, lack solvent topology
    )
    parser.add_argument(
        "-l", "--sorted", action="store_true", dest="sorted", help="sort atoms for GMX ordering",
    )
    parser.add_argument(
        "-j",
        "--chiral",
        action="store_true",
        dest="chiral",
        help="create improper dihedral parameters for chiral atoms in CNS",
    )
    parser.add_argument(
        "-v", "--version", action="store_true", dest="version", help="Show the Acpype version and exit",
    )
    return parser


def init_main(argv=None):

    """
    Main function, to satisfy Conda
    """
    if argv is None:
        argv = sys.argv[1:]

    parser = get_option_parser()
    args = parser.parse_args(argv)

    at0 = time.time()

    amb2gmxF = False

    if args.version:
        print(header)
        sys.exit(0)

    print(header)

    if not args.input:
        amb2gmxF = True
        if not args.inpcrd or not args.prmtop:
            parser.error("missing input files")
    elif args.inpcrd or args.prmtop:
        parser.error("either '-i' or ('-p', '-x'), but not both")

    if args.debug:
        texta = "Python Version %s" % sys.version
        print("DEBUG: %s" % while_replace(texta))

    if args.direct and not amb2gmxF:
        parser.error("option -u is only meaningful in 'amb2gmx' mode (args '-p' and '-x')")

    try:
        if amb2gmxF:
            print("Converting Amber input files to Gromacs ...")
            system = MolTopol(
                acFileXyz=args.inpcrd,
                acFileTop=args.prmtop,
                amb2gmx=True,
                debug=args.debug,
                basename=args.basename,
                verbose=args.verboseless,
                gmx4=args.gmx4,
                disam=args.disambiguate,
                direct=args.direct,
                is_sorted=args.sorted,
                chiral=args.chiral,
            )

            system.printDebug("prmtop and inpcrd files parsed")
            system.writeGromacsTopolFiles()
        else:
            molecule = ACTopol(
                args.input,
                chargeType=args.charge_method,
                chargeVal=args.net_charge,
                debug=args.debug,
                multiplicity=args.multiplicity,
                atomType=args.atom_type,
                force=args.force,
                outTopol=args.outtop,
                allhdg=args.cnstop,
                basename=args.basename,
                timeTol=args.max_time,
                qprog=args.qprog,
                ekFlag='''"%s"''' % args.keyword,
                verbose=args.verboseless,
                gmx4=args.gmx4,
                disam=args.disambiguate,
                direct=args.direct,
                is_sorted=args.sorted,
                chiral=args.chiral,
                amb2gmx=False,
            )

            molecule.createACTopol()
            molecule.createMolTopol()

        acpypeFailed = False
    except Exception:
        _exceptionType, exceptionValue, exceptionTraceback = sys.exc_info()
        print("ACPYPE FAILED: %s" % exceptionValue)
        if args.debug:
            traceback.print_tb(exceptionTraceback, file=sys.stdout)
        acpypeFailed = True

    execTime = int(round(time.time() - at0))
    if execTime == 0:
        amsg = "less than a second"
    else:
        amsg = elapsedTime(execTime)
    print("Total time of execution: %s" % amsg)

    if args.ipython:
        import IPython  # pylint: disable=import-outside-toplevel

        IPython.embed()

    try:
        rmtree(molecule.tmpDir)
    except Exception:
        pass
    if acpypeFailed:
        sys.exit(19)
    try:
        os.chdir(molecule.rootDir)
    except Exception:
        pass

    if checkSmiles(args.input):
        afile = "smiles_molecule.mol2"
        if os.path.exists(afile):
            os.remove(afile)


if __name__ == "__main__":
    init_main()  # necessary for to call in anaconda package;
