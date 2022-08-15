"""
---
    Copyright (c) 2018 Baskar Ganapathysubramanian, Balaji Sesha Sarath Pokuri
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
---
"""

## --- end license text --- ##
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name="PARyOpt",
    version="1.0.2.1",
    license='LICENSE.txt',
    packages=['PARyOpt'],
    package_data={
        'PARyOpt': ['PARyOpt/*.py', 'PARyOpt/evaluators/*.py', 'PARyOpt/kernel/*.py']
    },
    include_package_data=True,
    description="A flexible Bayesian optimization framework.",
    long_description=open('README.rst').read(),
    url="https://bitbucket.org/baskargroup/paryopt",
    python_requires='>=3.5',  # need parameter typing
    install_requires=[
        'numpy>=1.12.1',
        'paramiko>=2.1.2,<3',
        'scipy>=1.0.1',
        'sphinx==1.7.0'       # hard requirement for apidoc autobuild
    ],
    extras_require={
        'allow executing parse_result functions directly on remote hosts': ['dill>=0.2.5'],
        'for visualization': ['matplotlib>=3.0.0']
    },
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
        ],
    project_urls={
        'Documentation': 'https://paryopt.rtfd.io',
        'Source': 'https://bitbucket.org/baskargroup/paryopt',
        'Tracker': 'https://bitbucket.org/baskargroup/paryopt/issues',
        },
    author="Balaji Pokuri",
    author_email="balajip@iastate.edu",
    maintainer="Balaji Pokuri",
    maintainer_email="balajip@iastate.edu"
)
