# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:19:59 2023

@author: mwodring
"""

from setuptools import setup, find_packages

DESCRIPTION = 'Angua, a Bioinformatics pipeline using Blast, and its extension, Luggage.'
LONG_DESCRIPTION = 'Angua, a Bioinformatics pipeline using Blast, and its extension, Luggage.'

setup(
        package_data={"data": ["/data/"]},
        include_package_data=True,
        author="Sam McGreig and Morgan Wodring",
        author_email="morgan.wodring@gmail.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION
)