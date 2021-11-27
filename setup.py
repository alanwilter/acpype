#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File needed for PyPI
"""
from typing import Dict
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

dver: Dict = {}
with open("./acpype_lib/__init__.py") as fp:
    exec(fp.read(), dver)

setup(
    name="acpype",
    version=dver["__version__"],
    scripts=["acpype_lib/acpype.py"],
    author="Alan Silva",
    author_email="alanwilter@gmail.com",
    description="ACPYPE - AnteChamber PYthon Parser interfacE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alanwilter/acpype/tree/new_pip_package",
    project_urls={
        "Bug Tracker": "https://github.com/alanwilter/acpype/issues",
        "Wiki": "https://github.com/alanwilter/acpype/wiki",
    },
    packages=["acpype_lib", "amber21-11_linux", "amber21-11_os"],
    package_dir={"acpype": "acpype"},
    include_package_data=True,
    keywords="acpype",
    entry_points={"console_scripts": ["acpype = acpype_lib.acpype:init_main"]},
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
    install_requires=[],
    zip_safe=False,
)
