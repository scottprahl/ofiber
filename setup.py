"""
Copyright 2018 Scott Prahl

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from setuptools import setup

# use README as the long description
# make sure to use the syntax that works in both ReST and markdown
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
	name='ofiber',
	packages=['ofiber'],
	version='0.4.2',
	description='Light Propagation in Optical Fibers',
	url='https://github.com/scottprahl/ofiber.git',  
	author='Scott Prahl',
	author_email='scott.prahl@oit.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'License :: OSI Approved :: MIT License',
		'Intended Audience :: Science/Research',
		'Programming Language :: Python',
		'Topic :: Scientific/Engineering :: Physics',
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 3.7',
	],
	keywords=['refraction', 'Sellmeier', 'modes', 'dispersion'],
	install_requires=['numpy','scipy'],
	test_suite='ofiber.test.test',
    long_description=long_description,
    long_description_content_type='text/markdown'
)