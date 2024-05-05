.. |pypi-badge| image:: https://img.shields.io/pypi/v/ofiber?color=68CA66
   :target: https://pypi.org/project/ofiber/
   :alt: pypi
.. |github-badge| image:: https://img.shields.io/github/v/tag/scottprahl/ofiber?label=github&color=68CA66
   :target: https://github.com/scottprahl/ofiber
   :alt: github
.. |conda-badge| image:: https://img.shields.io/conda/vn/conda-forge/ofiber?label=conda&color=68CA66
   :target: https://github.com/conda-forge/ofiber-feedstock
   :alt: conda
.. |doi-badge| image:: https://zenodo.org/badge/122556263.svg
   :target: https://zenodo.org/doi/10.5281/zenodo.8368598
   :alt: doi  

.. |license-badge| image:: https://img.shields.io/github/license/scottprahl/ofiber?color=68CA66
   :target: https://github.com/scottprahl/ofiber/blob/master/LICENSE.txt
   :alt: License
.. |test-badge| image:: https://github.com/scottprahl/ofiber/actions/workflows/test.yaml/badge.svg
   :target: https://github.com/scottprahl/ofiber/actions/workflows/test.yaml
   :alt: Testing
.. |readthedocs-badge| image:: https://readthedocs.org/projects/ofiber/badge?color=68CA66
   :target: https://ofiber.readthedocs.io
   :alt: Docs
.. |downloads-badge| image:: https://img.shields.io/pypi/dm/ofiber?color=68CA66
   :target: https://pypi.org/project/ofiber/
   :alt: Downloads

ofiber
=======

by Scott Prahl

|pypi-badge| |github-badge| |conda-badge| |doi-badge|

|license-badge| |test-badge| |readthedocs-badge| |downloads-badge|

Python code to calculate light propagation through optical fibers following
the approach presented in `Ghatak and Thyagarajan, An Introduction to Fiber Optics <https://doi.org/10.1017/CBO9781139174770>`_.  Far-field fiber calculations are based on `Chen, Foundations for 
Guided-Wave Optics <https://doi.org/10.1002/0470042222>`_.


Installation
-------------

Use ``pip``::

    pip install ofiber

or ``conda``::

    conda install -c conda-forge ofiber

Usage
-----

A few examples are shown below. For all examples, see `ofiber documentation <https://ofiber.readthedocs.io>`_

Symmetric planar waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/planarwaveguide.svg
   :target: https://ofiber.readthedocs.io/en/latest/3-Planar-Waveguide-Modes.html
   :align: center
   :alt: Planar Waveguide
   
Cylindrical fibers with step index profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/stepindexmodes.svg
   :target: https://ofiber.readthedocs.io/en/latest/4-Circular-Step-Index-Fiber.html
   :align: center
   :alt: Modes in Step Index Fiber

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/modeirradiance.svg
   :target: https://ofiber.readthedocs.io/en/latest/4-Circular-Step-Index-Fiber.html
   :align: center
   :alt: Mode Irradiance

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/internalmodes.svg
   :target: https://ofiber.readthedocs.io/en/latest/4-Circular-Step-Index-Fiber.html
   :align: center
   :alt: Internal Modes

Far-field emission for step index fibers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/farfieldirradiance.svg
   :target: https://ofiber.readthedocs.io/en/latest/9-Far-field-irradiance.html
   :align: center
   :alt: Far-field Irradiance

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/theta01.svg
   :target: https://ofiber.readthedocs.io/en/latest/9-Far-field-irradiance.html
   :align: center
   :alt: polar angle of the minimum of the central irradiance lobe

Fiber design
^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/fiberdesign.svg
   :target: https://ofiber.readthedocs.io/en/latest/6-Zero-Dispersion.html
   :align: center
   :alt: Fiber Design

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/master/docs/dispersion.svg
   :target: https://ofiber.readthedocs.io/en/latest/6-Zero-Dispersion.html
   :align: center
   :alt: Dispersion

Google Colaboratory
^^^^^^^^^^^^^^^^^^^^

Use a Jupyter notebook immediately by clicking the Google Colaboratory button below

.. image:: https://colab.research.google.com/assets/colab-badge.svg
  :target: https://colab.research.google.com/github/scottprahl/ofiber/blob/master
  :alt: Colab


License
-------

``ofiber`` is licensed under the terms of the MIT license.