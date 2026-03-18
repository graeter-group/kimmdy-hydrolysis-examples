#!/bin/bash
# combine data for publishing into a zip file to upload with a release

# This release contains results.zip, which is meant to be unzipped here, as an artifact.
# in such a way that the paths are maintained
#
#     results/rates.csv
#     results/sample_rates.csv
#     examples/gly-hydrolysis/dist.xvg

zip results.zip \
    results/rates.csv \
    results/sample_rates.csv \
    examples/gly-hydrolysis/dist.xvg
