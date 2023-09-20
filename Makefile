html:
	cd docs && python -m sphinx -T -E -b html -d _build/doctrees -D language=en . _build
	open docs/_build/index.html

rstcheck:
	-rstcheck README.rst
	-rstcheck CHANGELOG.rst
	-rstcheck docs/index.rst
	-rstcheck docs/changelog.rst
	-rstcheck --ignore-directives automodule docs/ofiber.rst

lint:
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
	pytest --verbose tests/test_all_notebooks.py

rcheck:
	make lintcheck
	make doccheck
	make rstcheck
	make html
	check-manifest
	pyroma -d .
	make notecheck

lite:
	mkdir files
	cp docs/0-Basics.ipynb files
	cp docs/1-Refractive-Index.ipynb files
	cp docs/2-Materials.ipynb files
	cp docs/3-Planar-Waveguide-Modes.ipynb files
	cp docs/4-Circular-Step-Index-Fiber.ipynb files
	cp docs/5-Dispersion.ipynb files
	cp docs/6-Zero-Dispersion.ipynb files
	cp docs/7-Detectors.ipynb files
	jupyter lite init
	jupyter lite build
	jupyter lite serve
#	cp docs/8-Optical-Fiber-Amplifiers.ipynb

clean:
	rm -rf docs/api 
	rm -rf docs/_build 
	rm -rf dist
	rm -rf ofiber.egg-info
	rm -rf ofiber/__pycache__
	rm -rf ofiber/*.pyc
	rm -rf ofiber/.ipynb_checkpoints
	rm -rf .tox
	rm -rf __pycache__
	rm -rf build
	rm -rf files
	rm -rf _output
	rm -rf .jupyterlite.doit.db
	rm -rf .pytest_cache


.PHONY: clean rcheck html rstcheck lintcheck doccheck rcheck