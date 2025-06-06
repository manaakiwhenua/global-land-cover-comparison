```bash
conda env create -f workflow/envs/base.yaml
conda activate snakemake9
./run.sh
```

NB: you'll also need to add a (git-ignored) file called `secrets.env` which defines the following environment variables for secrets, in the form `KEY=value`:
- `LRIS_KEY` (API Key for [LRIS](https://lris.scinfo.org.nz/))