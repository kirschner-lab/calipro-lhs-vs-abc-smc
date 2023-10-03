# CaliPro-LHS vs ABC-SMC scripts for simulations and plots

Python and R scripts used for the simulations and plots published in
Nanda and Kirschner 2023 Front Appl Math Stat.

To re-run the analysis, install the dependencies:

```shell
Rscript -e 'install.packages(c("tidyverse", "deSolve", "lhs", "plyr", "extraDistr", "RNetCDF", "Rtsne", "cowplot"))'
python3 -m pip install --user pymc seaborn ninja
```

Then in the root directory with the `build.ninja` file,
run the `ninja` command to run the model simulations and plot their outputs.
Output results are also available in a separate 1.1 GB release.
Note the the `build.ninja` file was added later
to try and automate running the entire analysis,
therefore please instead refer to the output results release
as the primary results reference.

Informal discussion of Bayesian pairs plots and parameter summaries
is documented in `doc/suppl-bayesian.pdf`
