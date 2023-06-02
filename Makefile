.PHONY: test clean-pyc clean-build docs clean

# https://blog.ganssle.io/articles/2021/10/setup-py-deprecated.html

help:
	@echo "init - install submodules"
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "test-all - run tests on every Python version with tox"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "release - package and upload a release"
	@echo "dist - package"
	@echo "install - install the package to the active Python's site-packages"

init:
	git submodule update --init --recursive

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/

lint:
	flake8 vcftoolz tests

test:
	python setup.py test

test-all:
	tox

coverage:
	coverage run --source vcftoolz setup.py test
	coverage report -m
	coverage html
	open htmlcov/index.html

docs:
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	xdg-open docs/_build/html/index.html

release: clean dist
	twine upload dist/*

dist: clean
	python -m build
	ls -l dist

install: clean
	python -m pip install -e .
