#!/bin/bash

conda
brew install imagemagick
conda install -c conda-forge img2pdf
webshot::install_phantomjs()
Rscript scripts/install_dep.R