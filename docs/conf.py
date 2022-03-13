# pylint: disable=invalid-name
# pylint: disable=consider-using-f-string
"""
Configuration file for building documentation.

Sphinx builds the docs using couple of external modules: napoleon and nbsphinx.

The overall format is controlled by `.rst` files. The top level file is `index.rst`

`napoleon` builds the API in HTML assuming that the code is documented with
docstrings that follow the Google docstring format.

`nbsphinx` convert the Jupyter notebooks to html with nbsphinx, will
"""

import re
import os.path

project = 'ofiber'
<<<<<<< HEAD
copyright = '2018-2022 Scott Prahl'
author = 'Scott Prahl'
=======
>>>>>>> bf58f3a (put __version__ in __init__.py)

def get_init_property(prop):
    """Return property from __init__.py."""
    here = os.path.abspath(os.path.dirname(__file__))
    file_name = os.path.join(here, '..', project, '__init__.py')
    regex = r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop)
    with open(file_name, 'r', encoding='utf-8') as file:
        result = re.search(regex, file.read())
    return result.group(1)

release = get_init_property("__version__")
author = get_init_property("__author__")

master_doc = 'index'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx_automodapi.automodapi',
    'nbsphinx',
]
numpydoc_show_class_members = False
napoleon_use_rtype = False

<<<<<<< HEAD
# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.tox', '**.ipynb_checkpoints',
                    'sources/*', 'first_version/*', '8-Optical-Fiber-Amplifiers.ipynb']
=======
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '.DS_Store', '.tox', '**.ipynb_checkpoints']
>>>>>>> bf58f3a (put __version__ in __init__.py)

# I execute the notebooks manually in advance. If notebooks test the code,
# they should be run at build time.
nbsphinx_execute = 'never'
nbsphinx_allow_errors = True

# Add type of source files
source_suffix = ['.rst', '.ipynb']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

<<<<<<< HEAD
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

=======
>>>>>>> bf58f3a (put __version__ in __init__.py)
html_scaled_image_link = False

html_sourcelink_suffix = ''