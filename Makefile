SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docs
BUILDDIR      = docs/_build

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

clean:
	rm -rf docs/api 
	rm -rf docs/_build 
	rm -rf docs/_build/.buildinfo
	rm -rf docs/_build/.doctrees
	rm -rf dist
	rm -rf ofiber.egg-info
	rm -rf ofiber/__pycache__
	rm -rf ofiber/*.pyc
	rm -rf .tox

rstcheck:
	-rstcheck README.rst
	-rstcheck CHANGELOG.rst
	-rstcheck docs/index.rst
	-rstcheck docs/changelog.rst
	-rstcheck --ignore-directives automodule docs/ofiber.rst

pylint:
	-pylint ofiber/basics.py
	-pydocstyle ofiber/basics.py
	-pylint ofiber/cylinder_step.py
	-pydocstyle ofiber/cylinder_step.py
	-pylint ofiber/dispersion.py
	-pydocstyle ofiber/dispersion.py
	-pylint ofiber/graded_index.py
	-pydocstyle ofiber/graded_index.py
	-pylint ofiber/noise.py
	-pydocstyle ofiber/noise.py
	-pylint ofiber/planar_parabolic.py
	-pydocstyle ofiber/planar_parabolic.py
	-pylint ofiber/planar_step.py
	-pydocstyle ofiber/planar_step.py
	-pylint ofiber/refraction.py
	-pydocstyle ofiber/refraction.py
	-pylint ofiber/__init__.py
	-pydocstyle ofiber/__init__.py

rcheck:
	make clean
	touch docs/*ipynb
	touch docs/*rst
	make pylint
	make rstcheck
	make html
	check-manifest
	pyroma -d .
#	tox

.PHONY: clean rcheck html lint