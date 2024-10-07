## Summary

Code to reproduce figures from Williams, Oliphant and Au et al. "Luminal breast epithelial cells of BRCA1 or BRCA2 mutation carriers and non-carriers harbor common breast cancer copy number alterations".

## Steps to run code

1. Download the processed data from [zenodo](https://zenodo.org/records/13645602) and call the folder zenododirectory.
2. Change the path of `basedir` in `config.yaml` to point to the directory above where the zenodo folder has been placed.
3. Download processed cancer data from [zenodo](https://zenodo.org/records/13898350) and set `cancerdir` in `config.yaml`.
4. Run code (in `scripts/`) for each of the figures.

## Environment

The [docker](https://hub.docker.com/repository/docker/marcjwilliams1/signals) file associated with [signals](https://github.com/shahcompbio/signals) should have everything you need. Code was tested with v0.11.0. There is also session info at the end of each Rmd.