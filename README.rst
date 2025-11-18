.. |pypi-badge| image:: https://img.shields.io/pypi/v/ofiber?color=68CA66
   :target: https://pypi.org/project/ofiber/
   :alt: PyPI

.. |github-badge| image:: https://img.shields.io/github/v/tag/scottprahl/ofiber?label=github&color=68CA66
   :target: https://github.com/scottprahl/ofiber
   :alt: GitHub

.. |conda-badge| image:: https://img.shields.io/conda/vn/conda-forge/ofiber?label=conda&color=68CA66
   :target: https://github.com/conda-forge/ofiber-feedstock
   :alt: Conda

.. |doi-badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.8368598.svg
   :target: https://doi.org/10.5281/zenodo.8368598
   :alt: DOI

.. |license-badge| image:: https://img.shields.io/github/license/scottprahl/ofiber?color=68CA66
   :target: https://github.com/scottprahl/ofiber/blob/main/LICENSE.txt
   :alt: License

.. |test-badge| image:: https://github.com/scottprahl/ofiber/actions/workflows/test.yaml/badge.svg
   :target: https://github.com/scottprahl/ofiber/actions/workflows/test.yaml
   :alt: Testing

.. |readthedocs-badge| image:: https://readthedocs.org/projects/ofiber/badge?color=68CA66
   :target: https://ofiber.readthedocs.io
   :alt: Documentation

.. |downloads-badge| image:: https://img.shields.io/pypi/dm/ofiber?color=68CA66
   :target: https://pypi.org/project/ofiber/
   :alt: Downloads

.. |lite| image:: https://img.shields.io/badge/try-JupyterLite-68CA66.svg
   :target: https://scottprahl.github.io/ofiber/
   :alt: Try Online

ofiber
=======

|pypi-badge| |github-badge| |conda-badge| |doi-badge|

|license-badge| |test-badge| |readthedocs-badge| |downloads-badge|

|lite|

``ofiber`` is a Python library for analyzing guided-wave propagation in optical fibers and related dielectric waveguiding structures.
It provides analytical and numerical tools for mode analysis, dispersion engineering, and far-field radiation modeling.

The theoretical framework closely follows well-established treatments in:

- **Ghatak & Thyagarajan**, *An Introduction to Fiber Optics* (Cambridge University Press)
  `<https://doi.org/10.1017/CBO9781139174770>`_

- **Chen**, *Foundations for Guided-Wave Optics* (Wiley)
  `<https://doi.org/10.1002/0470042222>`_

The accompanying examples and notebooks are intended to support instruction, research reproducibility, and optical design workflows.

Installation
------------

With ``pip``::

    pip install ofiber

Or with ``conda`` from conda-forge::

    conda install -c conda-forge ofiber


Documentation and Examples
---------------------------

Comprehensive user documentation, theory notes, and executable Jupyter examples are available at:

📄 https://ofiber.readthedocs.io

or use immediately in your browser via the JupyterLite button below

    |lite|

Examples
--------

Symmetric planar waveguides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/main/docs/planarwaveguide.svg
   :target: https://ofiber.readthedocs.io/en/latest/3-Planar-Waveguide-Modes.html
   :align: center

Cylindrical step-index fibers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/main/docs/stepindexmodes.svg
   :target: https://ofiber.readthedocs.io/en/latest/4-Circular-Step-Index-Fiber.html
   :align: center

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/main/docs/modeirradiance.svg
   :target: https://ofiber.readthedocs.io/en/latest/4-Circular-Step-Index-Fiber.html
   :align: center

Far-field radiation patterns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/main/docs/farfieldirradiance.svg
   :target: https://ofiber.readthedocs.io/en/latest/9-Far-field-irradiance.html
   :align: center

Fiber design and dispersion control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: https://raw.githubusercontent.com/scottprahl/ofiber/main/docs/fiberdesign.svg
   :target: https://ofiber.readthedocs.io/en/latest/6-Zero-Dispersion.html
   :align: center


Citation
--------

If you use ``ofiber`` in academic, instructional, or applied technical work, please cite:

Prahl, S. (2025). *ofiber: A Python module for modeling guided-wave light propagation in optical fibers* (Version 0.9.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.8368598

BibTeX
^^^^^^

.. code-block:: bibtex

   @software{ofiber_prahl_2025,
     author    = {Scott Prahl},
     title     = {ofiber: A Python module for modeling guided-wave light propagation in optical fibers},
     year      = {2025},
     version   = {0.9.0},
     doi       = {10.5281/zenodo.8368598},
     url       = {https://github.com/scottprahl/ofiber},
     publisher = {Zenodo}
   }

---

License
-------

``ofiber`` is released under the MIT License.
