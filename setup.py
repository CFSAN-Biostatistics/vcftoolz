#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
    'biopython',
    'numpy',
    'pandas',
    'PyVCF',
    'matplotlib',
    'matplotlib_venn',
]

test_requirements = [
    "pytest",
]

setup(
    name='vcftoolz',
    version='1.2.0',
    description="Tools for working with Variant Call Format files.",
    long_description=readme + '\n\n' + history,
    author="Steve Davis",
    author_email='steven.davis@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/vcftoolz',
    packages=[
        'pyvenn',
        'vcftoolz',
    ],
    package_dir={'vcftoolz':
                 'vcftoolz'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['bioinformatics', 'NGS', 'vcftoolz'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: Freely Distributable',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={'console_scripts': ['vcftoolz = vcftoolz.cli:main']},
    setup_requires=["pytest-runner"],
    tests_require=test_requirements
)
