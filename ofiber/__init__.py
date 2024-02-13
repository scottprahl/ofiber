"""
Calculate light propagation through optical fibers.

For full documentation see <https://ofiber.readthedocs.io>

Info about calculating simple optical fibers parameters::

    help(ofiber.basics)

Info about modes and other characteristics of cylindrical fibers::

    help(ofiber.cylinder_step)

Info about material and waveguide dispersion::

    help(ofiber.dispersion)

Info about noise calculations relevant to optical communications::

    help(ofiber.dispersion)

Info about modes in planar waveguides with step index profiles::

    help(ofiber.planar_step)

Info about modes in planar waveguides with parabolic index profiles::

    help(ofiber.planar_parabolic)

Info about refractive index of glasses::

    help(ofiber.refraction)
"""
__version__ = '0.8.0'
__author__ = 'Scott Prahl'
__email__ = 'scott.prahl@oit.edu'
__copyright__ = '2018-23, Scott Prahl'
__license__ = 'MIT'
__url__ = 'https://github.com/scottprahl/ofiber'

from .basics import *
from .cylinder_step import *
from .dispersion import *
from .graded_index import *
from .noise import *
from .planar_parabolic import *
from .planar_step import *
from .refraction import *
