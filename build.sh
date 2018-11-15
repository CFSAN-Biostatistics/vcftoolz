#!/bin/bash

rm -r vcftoolz.egg-info/
python setup.py clean --all
python setup.py sdist --verbose

