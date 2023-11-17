Bayesian hierarchical modeling of size spectra
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

This repository contains code to recreate the analyses in Wesner et
al. *Bayesian hierarchical modeling of size spectra*. The purpose of the
paper is to demonstrate how build Bayesian generalized models using the
truncated Pareto distribution. For a detailed derivation of the
distribution, see A. M. Edwards et al. (2017) and A. Edwards et al.
(2020). All analyses in the paper depend on the `isdbayes` package,
which can be found here: <https://github.com/jswesner/isdbayes>.

## Recreate analyses

All figures and tables in the main text can be recreated by running the
R scripts in the `\code` folder. In that folder, each figure has two
associated scripts: “figure_x.R” and “figure_x_model_and_code.R”. The
**“figure_x.R”** scripts recreate the .rds and .jpg files of each
figure. They plot the model results obtained from the
**“figure_x_model_code.R”** scripts.

For example, `code/figure_4.R` creates the figure below, showing how
three modeling approaches (a-c) deal with an outlier in the size
spectrum parameter $\lambda$.

<img src="ms/fig4hierarchical_regression_plot.jpg" width="100%" />

The code to simulate data and fit these models is in
`code/figure_4_model_code.R`. That code requires the installation of
`rstan`, `brms`, and `isdbayes`, which in turn require installation of
an external C++ compiler. However, all models are already compiled and
saved in the `\models` folder. Unless you’d like to refit each model,
there is no need to run `code/figure_4_model_code.R`.

## Needing a cluster

For most of the scripts, running the modeling exercise will not work on
a personal computer. It requires a cluster, because the models are
repeatedly fit 1000 times, resulting in large files. In the code, this
is achieved using the `replicate()` function, with n = 1000. To avoid
automatically running large models, the replicates in the scripts are
set to a lower value, like n = 2.

Some of the fitted model simulations were too big to upload to GitHub.
However, the fitted models are summarized (e.g., quantiles of the
posteriors) before plotting and all of those summaries are available
either in the `\models` folder or the `\posteriors` folder.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-edwards2017" class="csl-entry">

Edwards, A. M., J. P. W. Robinson, M. J. Plank, J. K. Baum, and J. L.
Blanchard. 2017. “Testing and Recommending Methods for Fitting Size
Spectra to Data.” *Methods in Ecology and Evolution* 8 (1): 57–67.
<http://dx.doi.org/10.1111/2041-210X.12641>.

</div>

<div id="ref-edwards2020" class="csl-entry">

Edwards, Am, Jpw Robinson, Jl Blanchard, Jk Baum, and Mj Plank. 2020.
“Accounting for the Bin Structure of Data Removes Bias When Fitting Size
Spectra.” *Marine Ecology Progress Series* 636 (February): 19–33.
<https://doi.org/10.3354/meps13230>.

</div>

</div>
