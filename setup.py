# -*- coding: utf-8 -*-
# CUWP: Chemical Upcycling of Waste Plastics Process Models
# Copyright (C) 2025-2027, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the MIT open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/CUWP/blob/master/LICENSE.txt
# for license details.
from setuptools import setup

setup(
    name='cuwp',
    packages=['cuwp'],
    license='MIT',
    version='0.0.4',
    description="Chemical Upcycling of Waste Plastics Process Models",
    long_description=open('README.rst', encoding='utf-8').read(),
    author='Yoel Cortes-Pena',
    install_requires=['biosteam>=2.53.0',
                      'biorefineries>=2.35.0'],
    python_requires=">=3.12",
    package_data={
        'plastics': [
            'cytiva',
            'strap/*',
            'strap/data/*',
            'pyrolysis',
            'pyrolysis/data/*',
        ]
    },
    platforms=['Windows', 'Mac', 'Linux'],
    author_email='yoelcortes@gmail.com',
    url='https://github.com/BioSTEAMDevelopmentGroup/cuwp',
    download_url='https://github.com/BioSTEAMDevelopmentGroup/cuwp',
    classifiers=['Development Status :: 3 - Alpha',
                 'Environment :: Console',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.12',
                 'Programming Language :: Python :: 3.13',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Mathematics'],
    keywords='chemical process simulation plastic bioprocess engineering STRAP solvent targeted dissolution precipitation',
)
