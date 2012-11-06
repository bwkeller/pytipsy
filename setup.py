#!/usr/bin/env python

from distutils.core import setup

setup(
		name="PyTipsy",
		version='1.0',
		author='BW Keller',
		author_email='malzraa@gmail.com',
		url='https://github.com/bwkeller/pytipsy',
		py_modules=['pytipsy'],
		licence='LICENCE.txt',
		long_description=open('README.md').read(),
		install_requires=[
			"Numpy >= 1.5.1"
			],
)
