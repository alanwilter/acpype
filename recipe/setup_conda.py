"""
AnteChamber PYthon Parser interfacE

A tool based in Python to use Antechamber to generate topologies for chemical compounds
and to interface with others python applications like CCPN or ARIA.
"""
from typing import Dict
from setuptools import setup

dver: Dict = {}
with open("./acpype_lib/__init__.py") as fp:
    exec(fp.read(), dver)

setup(
    name="acpype",
    version=dver["__version__"],
    description="ACPYPE - AnteChamber PYthon Parser interfacE",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    url="https://github.com/alanwilter/acpype",
    author="Alan Silva",
    author_email="alanwilter@gmail.com",
    license="GPL-3.0-or-later",
    packages=["acpype_lib"],
    keywords=["acpype"],
    include_package_data=True,
    entry_points={"console_scripts": ["acpype = acpype_lib.acpype:init_main"]},
    zip_safe=False,
)
