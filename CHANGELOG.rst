Changelog
==========

Unreleased
----------
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
* add make coverage to report unit test coverage of the package in the terminal
* add make update-notebooks to format the docs notebooks with black and re-execute them in place
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
* fix TE_planar_parabolic_field() calling np.math.factorial, removed in numpy 2.0
* fix parabolic_propagation_constants() counting modes with an unrelated formula and passing a float to np.empty
* add unit tests for planar_parabolic.py
* add unit tests for planar_step.py
* fix TE_field() and TM_field() returning a flat field for a mode the waveguide cannot guide, nan is now returned
* fix TE_propagation_constant() and TM_propagation_constant() raising TypeError on a scalar V
* fix integer arrays of V returning zeros from TE_propagation_constant() and TM_propagation_constant()
* correct the TM signatures in the planar_step docstring and label them magnetic rather than electric
* fix the symmetric branch of the mode plots being drawn every pi/2 instead of every pi
* fix TM_mode_plot() scaling its x axis by (n1/n2)**2, xi never exceeds V/2
* V_d2bV_by_V() returns nan rather than 0 for an unguided mode, since 0 is a legal value
* document that FF_node_polar_angle() returns k*a*sin(theta), not the polar angle itself
* correct the Dispersion() docstring, it returns material and waveguide dispersion but not their sum
* correct the glass() docstring, which referenced a nonexistent ofiber.glass_index and inverted find_glass and glass_name
* correct LP_mode_value() pointing at itself instead of _LP_mode_value for the calculation details
* correct the cutoff_wavelength() docstring, the default is the LP11 cutoff and LP01 has no cutoff at all
* add V_parameter() to the basics docstring, drop the duplicated critical_angle entry and stray colons
* describe m in R_par(), R_per() and R_unpolarized() as the relative index of the two media
* fix R_par(), R_per() and R_unpolarized() returning nan past the critical angle instead of total reflection
* add unit tests for basics.py
* note in the noise docstrings that I0 is the post-gain current for shot_noise but the primary one for best_APD_gain
* correct the thermal_min_power() summary, it reaches the requested signal-to-noise ratio rather than one
* list the noise functions in the module docstring and explain Np in quantum_min_power()
* use scipy.constants in noise.py instead of four-digit values for q, k, h and c, matching the rest of the package
* replace the 4.34 and 4.343 dB conversions in cylinder_step.py with an exact 10/ln(10), they disagreed with each other
* use scipy.constants.speed_of_light rather than 3e8 in 1-Refractive-Index.ipynb
* explain in Waveguide_Dispersion() why q defaults to 1e20 rather than np.inf, which makes esi_Delta nan
* correct the accuracy claim for approx=True in Waveguide_Dispersion(), the error peaks near 6% at V=2.44
* list the functions in the dispersion module docstring and note where the silica material zero falls
* add unit tests for dispersion.py, completing 100% statement and branch coverage of the package
* build Delta and V in Waveguide_Dispersion() from the basics routines instead of repeating their formulas
* run the unit tests in CI rather than only pip check, and execute every notebook in a second job
* require the tests to pass before the pypi workflow publishes a release
* drop the package-wide automodapi from docs/ofiber.rst, which documented every function a second time
* accept q=np.inf in esi_Delta(), esi_radius() and esi_V_parameter() so dispersion and basics share one step-index sentinel
* raise RuntimeError rather than StopIteration when the far-field node search finds no sign change
* fix an RST substitution error from writing ``|x|`` in the power_law_profile() docstring
* add a black-check target that reports formatting without rewriting files
* add a lint target that runs black, ruff, pylint, rstcheck and yamllint, and have rcheck defer to it
* run make lint in CI as its own job, using the official astral-sh/setup-uv action with caching
* update the CI actions, checkout v4 to v7 and setup-python v5 to v7
* reformat with black, which had drifted in six files
* add unit tests for noise.py

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
