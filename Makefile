# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docs
BUILDDIR      = docs/_build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

clean:
	rm -rf docs/_build 
	rm -rf docs/api 
	rm -rf docs/_build/.buildinfo
	rm -rf docs/_build/.doctrees
	rm -rf dist
	rm -rf ofiber.egg-info
	rm -rf ofiber/__pycache__
	rm -rf ofiber/*.pyc
	rm -rf .tox

lint:
	-pylint ofiber/basics.py
	-pep257 ofiber/basics.py
	-pylint ofiber/cylinder_step.py
	-pep257 ofiber/cylinder_step.py
	-pylint ofiber/dispersion.py
	-pep257 ofiber/dispersion.py
	-pylint ofiber/graded_index.py
	-pep257 ofiber/graded_index.py
	-pylint ofiber/noise.py
	-pep257 ofiber/noise.py
	-pylint ofiber/parabolic.py
	-pep257 ofiber/parabolic.py
	-pylint ofiber/planar_step.py
	-pep257 ofiber/planar_step.py
	-pylint ofiber/refraction.py
	-pep257 ofiber/refraction.py

rcheck:
	make clean
	touch docs/*ipynb
	touch docs/*rst
	make lint
	make html
	check-manifest
	pyroma -d .
#	tox

.PHONY: clean check rcheck html