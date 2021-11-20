#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File needed for PyPI
"""
import setuptools
from acpype_lib._version import version

with open("README.md", "r") as fh:

    long_description = fh.read()

# updated__ = "2021.11.14"
# version = updated__[:19].replace("-", "").replace("T", "").replace(":", "")

setuptools.setup(
    name="acpype",
    version=version,
    scripts=["acpype_lib/acpype.py"],
    author="Alan Wilter Sousa da Silva",
    author_email="alanwilter@gmail.com",
    description="ACPYPE - AnteChamber PYthon Parser interfacE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alanwilter/acpype",
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
    install_requires=['pbr>=5.7.0'],
    zip_safe=False,
)
