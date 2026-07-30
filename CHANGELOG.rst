Changelog
==========

Unreleased ()
-------------
* migrate Makefile to a uv-first workflow with shared RUN wrappers and no .venv/.ready bootstrap prerequisite
* add RM/RMR cleanup variables and replace hardcoded rm cleanup commands
* raise package Python support metadata to 3.10 through 3.14 and update docs/lite dependency ranges
* update CI/runtime Python defaults to 3.14 and the test matrix to 3.10 and 3.14
* add release dates to each changelog heading using release or tag dates
* fix missing += in _d2_sellmeier so d2n() sums all three Sellmeier terms again
* restore the full Schott names of 44 glasses that were truncated to six characters
* find_glass() prefers an exact name and raises ValueError instead of quietly returning SiO2
* add unit tests for refraction.py
* make test runs the unit tests instead of doing nothing
* add make coverage to report unit test coverage of the package
* add make update-notebooks to re-execute the docs notebooks in place
* share NB_TIMEOUT between note-test and update-notebooks
* fix first_derivative() appending its padding value in the wrong place and dividing by a negative step
* add unit tests for graded_index.py
* fix power_law_profile() so a scalar radius outside the core also returns the cladding index
* export ray_delay() from the package, previously called velocity() and missing from __all__
* fix the material dispersion of ray_delay(), which did not reduce to the plain form at zero dispersion
* add a profile_dispersion argument to ray_delay() so the delay-equalising exponent comes out at q = 2 - 2P
* fix __all__ in cylinder_step.py, misspelled _all_, which leaked np, plt, scipy and special into ofiber
* add unit tests for cylinder_step.py
* fix integer arrays of V, ell or em returning zeros from LP_mode_value(), PetermannW(), V_d2bV_by_V(), bending_loss_db() and FF_node_polar_angle()
* fix the em check in LP_mode_value() and FF_node_polar_angle() raising on an array, which made the em branch unreachable
* correct the accuracy range in the V_d2bV_by_V_Approx() docstring, 1% holds for 1.35<V<2.08 rather than 1.4<V<2.4

0.9.1 (2026-01-14)
------------------
* remove requirements*.txt, all deps in pyproject.toml
* use np.trapezoid()
* version info only in __init__.py
* update readthedocs configuration
* update docs/conf.py
* update github action to publish to pypi
* fix pylint warnings in update_citation.py
* move jupyter_lite_config.json to ofiber folder
* update zenodo url

0.9.0 (2025-11-18)
-------------------
* support for jupyterlite
* use pyproject.toml and requirements-dev.txt
* better citation
* normalizing code formatting with black

0.8.1 (2024-05-05)
-------------------
* Implemented FF_node_polar_angle (@matt8s)
* Added FF_node_polar_angle() example to 9-Far-field-irradiance.ipynb (@matt8s)
* Modified _FF_polar_x to take argument kasin (@matt8s)
* Improved configuration for ruff linting tool
* Edited Jupyter notebooks
* Added images and improved README

0.8.0 (2024-02-13)
-------------------
* add functions for far-field irradiance (thanks @matt8s)
* improve calc of b for cylindrical step fibers for large V
* add new notebook for far-field irradiance
* enable use of arrays in `cylinder_step.LP_mode_value()`

0.7.1 (2023-09-21)
-------------------
* fix long-description for pypi

0.7.0 (2023-09-21)
-------------------
* better handling of propagation factor b
* add import boilerplate needed for Jupyterlite
* Makefile packaging for Jupyterlite file
* use __version__ and other dunders
* flake8 now passes
* fix various docstrings
* readthedocs fixes
* ruff support
* add erbium fiber doc
* clean up notebooks
* add copyright stuff
* add .github actions
* update badges
* conda support

v0.6.2 (2021-08-06)
-------------------
* create pure python packages
* include wheel file
* package as python3 only

v0.6.1 (2021-03-26)
-------------------
* add colab and binder badges
* improve docs
* rename dispersion api
* start using more descriptive names in notebooks

v0.6.0 (2021-01-05)
-------------------
* improve help()
* improve sphinx docs
* new names

v0.5.0 (2020-05-25)
-------------------
* sphinx docs
* reorganize and improve docs
* add lint target to Makefile
* add rcheck target to Makefile

v0.4.0 (2020-01-22)
--------------------
* improve packaging
* add some documentation
* more glasses
* more functions

v0.3.0 (2018-03-12)
-------------------
* add basics.py
* add documentation for functions

v0.2.1 (2018-03-05)
-------------------
* fix typo

v0.2.0 (2018-03-05)
-------------------
* Add noise.py
* Improve packaging

v0.1.0 (2018-02-26)
-------------------
* Initial release
