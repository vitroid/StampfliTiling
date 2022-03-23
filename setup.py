#!/usr/bin/env python3

# from distutils.core import setup, Extension
from setuptools import setup, find_packages
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(os.path.dirname(__file__), 'genice2_stampfli', '__init__.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))



long_desc = "".join(open("README.md").readlines())


setup(
    name='genice2-stampfli',
    version=metadata['version'],
    description='Lattice generator for GenIce2.',
    long_description=long_desc,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3.6",
    ],
    author='Masakazu Matsumoto',
    author_email='vitroid@gmail.com',
    url='https://github.com/vitroid/stampfli/',
    keywords=['genice', 'stampfli', 'quasicrystal'],

    # package_dir={"": "genice2_stampfli"},
    packages=['genice2_stampfli',
              'genice2_stampfli.lattices',
    ],

    entry_points = {
        'genice2_lattice': [
            'stampfli = genice2_stampfli.lattices.stampfli',
        ],
    },
    install_requires=['GenIce2',],

    license='MIT',
)
