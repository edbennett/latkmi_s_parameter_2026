# The Peskin-Takeuchi $S$ parameter and vector meson decay constants in $N_f = 8$ QCD

## Analysis workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21297277.svg)](https://doi.org/10.5281/zenodo.21297277)

The workflow in this repository performs
the analyses presented in the paper
[The Peskin-Takeuchi $S$ parameter and vector meson decay constants
in $N_f = 8$ QCD][paper].

## Requirements

- Conda, for example, installed from [Miniforge][miniforge]
- [Snakemake][snakemake], which may be installed using Conda
- LaTeX, for example, from [TeX Live][texlive]

## Setup

1. Install the dependencies above.
2. Clone this repository including submodules
   (or download its Zenodo release and `unzip` it)
   and `cd` into it:

   ```shellsession
   git clone --recurse-submodules https://github.com/edbennett/latkmi_s_parameter_2026
   cd latkmi_s_parameter_2026
   ```

3. From the [data release][datarelease],
   download the files `data.zip` and `ensembles.csv`.
   Unzip the former in the repository root to create the `data/` directory,
   and place the latter in the `metadata` directory.

## Running the workflow

The workflow is run using Snakemake:

``` shellsession
snakemake --cores 1 --use-conda
```

where the number `1`
may be replaced by
the number of CPU cores you wish to allocate to the computation.

Snakemake will automatically download and install
all required Python packages.
This requires an Internet connection;
if you are running in an HPC environment where you would need
to run the workflow without Internet access,
you can use the commands:

``` shellsession
snakemake --conda-create-envs-only
snakemake previous_data/spectrum_2505.08658.csv previous_data/fit_results_2505.08658.csv
```

The former command prepares the Conda environments,
while the latter downloads data from our previous work that this analysis relies upon.

The workflow takes around three hours to run
using six cores of an Apple M1 CPU.

## Output

Output plots, tables, and definitions
are placed in the `assets/plots`, `assets/tables`, and `assets/definitions` directories.

Output data assets are placed into the `data_assets` directory.

Intermediary data are placed in the `intermediary_data` directory.

## Reusability

This workflow is relatively tailored to the data
which it was originally written to analyse.
Additional ensembles may be added to the analysis
by adding relevant files to the `raw_data` directory,
and adding corresponding entries to the files in the `metadata` directory.
However,
extending the analysis in this way
has not been as fully tested as the rest of the workflow,
and is not guaranteed to be trivial for someone not already familiar with the code.

## References

This release contains data from the following sources:

- Phys.Rev.Lett. 101 (2008) 242001 [[0806.4222](https://arxiv.org/abs/0806.4222)]
  (`external_data/jlqcd_prl08_l10_r.csv`)
- Phys.Rev.D 81 (2010) 014504 [[0909.4931](https://arxiv.org/abs/0909.4931)]
  (`external_data/rbc_ukqcd_prd10_l10_r.csv`)
- Phys.Rev.D 90 (2014) 11, 114502 [[1405.4752](https://arxiv.org/abs/1405.2752)]
  (`external_data/lsd_prd14_spectra_sparameter_table_6.csv`, `lsd_prd14_spectra_sparameter_table_6.csv`)
- Phys.Rev.D 99 (2019) 1, 014509 [[1807.08411](https://arxiv.org/abs/1807.08411)]
  (`external_data/lsd_prd19_spectra_nf08_table_1_3_4.csv`)
- [https://doi.org/10.5281/zenodo.17037868](https://doi.org/10.5281/zenodo.17037868),
  updated with additional data from previous work as noted in the file
  (`previous_data/spectrum_2505.08658_updated.csv`)

[datarelease]: https://doi.org/10.5281/zenodo.21693826
[miniforge]: https://github.com/conda-forge/miniforge
[paper]: https://doi.org/10.48550/arXiv.2609.35371
[snakemake]: https://snakemake.github.io
[snakemake-conda]: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html
[texlive]: https://tug.org/texlive/
