#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:

    long_description = fh.read()

__updated__ = "2020-03-14T11:44:00CET"
version = __updated__[:19].replace('-', '').replace('T', '').replace(':', '')

setuptools.setup(
    name='acpype',
    version=version,
    scripts=['scripts/acpype.py', ],
    author='Alan Wilter Sousa da Silva',
    author_email='alanwilter@gmail.com',
    description='ACPYPE - AnteChamber PYthon Parser interfacE',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alanwilter/acpype",
    packages=['scripts', 'amber17-6_linux', 'amber17-6_os'],
    package_dir={'acpype':
                 'acpype'},
    include_package_data=True,
    keywords='acpype',
    entry_points={'console_scripts': ['acpype = scripts.acpype:init_main']},
    classifiers=[
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[],
    zip_safe=False
)
