SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docs
BUILDDIR      = docs/_build

html:
	$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

rstcheck:
	-rstcheck README.rst
	-rstcheck CHANGELOG.rst
	-rstcheck docs/index.rst
	-rstcheck docs/changelog.rst
	-rstcheck --ignore-directives automodule docs/ofiber.rst

lintcheck:
	-pylint ofiber/basics.py
	-pylint ofiber/cylinder_step.py
	-pylint ofiber/dispersion.py
	-pylint ofiber/graded_index.py
	-pylint ofiber/noise.py
	-pylint ofiber/planar_parabolic.py
	-pylint ofiber/planar_step.py
	-pylint ofiber/refraction.py
	-pylint ofiber/__init__.py

doccheck:
	-pydocstyle ofiber/basics.py
	-pydocstyle ofiber/cylinder_step.py
	-pydocstyle ofiber/dispersion.py
	-pydocstyle ofiber/graded_index.py
	-pydocstyle ofiber/noise.py
	-pydocstyle ofiber/planar_parabolic.py
	-pydocstyle ofiber/planar_step.py
	-pydocstyle ofiber/refraction.py
	-pydocstyle ofiber/__init__.py

notecheck:
	make clean
	pytest --verbose -n 4 test_all_notebooks.py

rcheck:
	make notecheck
	make lintcheck
	make doccheck
	make rstcheck
	make html
	check-manifest
	pyroma -d .
#	tox

clean:
	rm -rf docs/api 
	rm -rf docs/_build 
	rm -rf dist
	rm -rf ofiber.egg-info
	rm -rf ofiber/__pycache__
	rm -rf ofiber/*.pyc
	rm -rf .tox
	rm -rf __pycache__
	rm -rf build


.PHONY: clean rcheck html rstcheck lintcheck doccheck rcheck