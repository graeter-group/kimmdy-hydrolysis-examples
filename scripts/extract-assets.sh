#! /bin/env bash

# those files are too large to be included in the repository
# so we zip them first and then extract them

# create individual zip files, but only of the large files, if they do not exist yet
files='examples/collagen-hydrolysis/assets/collagen.ndx examples/collagen-hydrolysis/assets/collagen.top examples/collagen-hydrolysis/assets/npt.gro examples/collagen-hydrolysis/assets/nvt.gro'
for file in $files; do
  if [ ! -f "${file}.zip" ]; then
    echo "Creating zip file for $file"
    zip "${file}.zip" "$file"
  fi
done

# extract the files from the zip files if they do not exist yet
for file in $files; do
  if [ ! -f "$file" ]; then
    echo "Extracting $file from ${file}.zip"
    unzip "${file}.zip" $file
  fi
done

