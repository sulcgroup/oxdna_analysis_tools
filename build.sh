#!/bin/bash

# Autoincrement the version, remove old files, rebuild, and reinstall
awk -F'.' '{printf("%s.%s.%s'\''",$1,$2,$3+1)}' ./src/oxDNA_analysis_tools/__version__.py > ./src/oxDNA_analysis_tools/__version__.py.tmp
mv ./src/oxDNA_analysis_tools/__version__.py.tmp ./src/oxDNA_analysis_tools/__version__.py
rm dist/*
python -m build
python -m pip install .