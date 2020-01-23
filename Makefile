# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf _build 
	rm -rf api 
	rm -rf build 
	rm -rf dist
	rm -rf ofiber.egg-info
	rm -rf ofiber/__pycache__
	rm -rf ofiber/*.pyc
	
check:
	-check-manifest
	-pyroma -d .
	-pylint ofiber/basics.py
	-pep257 ofiber/basics.py
	-pylint ofiber/cylinder_step.py
	-pep257 ofiber/cylinder_step.py
	-pylint ofiber/graded_index.py
	-pep257 ofiber/graded_index.py
	-pylint ofiber/noise.py
	-pep257 ofiber/noise.py
	-pylint ofiber/parabolic.py
	-pep257 ofiber/parabolic.py
	
.PHONEY: clean check