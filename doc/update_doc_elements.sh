#!/bin/sh

echo "Update elements for documentation"

echo "Convert IPython notebooks to rst"

python -m IPython nbconvert --to rst ../notebooks/*ipynb

echo "Convert README.md to rst"

pandoc -s -t rst --toc ../README.md -o readme.rst
