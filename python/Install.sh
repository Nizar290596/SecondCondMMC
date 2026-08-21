#!/bin/bash

MMCFOAMREADER_VERSION="1.0.0"

python3 setup.py bdist_wheel
pip3 install --force-reinstall dist/mmcFoamReader-${MMCFOAMREADER_VERSION}-py3-none-any.whl
