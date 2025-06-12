#/bin/bash
snakemake --configfile config.yaml --configfile config.local.yaml --use-conda --cores 4 all

# NB delete --configfile config.local.yaml if you don't have a local config file