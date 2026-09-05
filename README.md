# PMLi: Statistical Tests for Partially Matched Samples

`PMLi` is an educational R package implementing five hypothesis-testing procedures for two samples that contain both paired observations and unmatched observations.

The package was developed by Kai Li for AMS 597: Statistical Computing in Spring 2021.

## Implemented methods

- `weighted.z()` — Liptak weighted-Z p-value pooling
- `modified.t()` — Kim et al. modified *t* statistic
- `corrected.z()` — Looney and Jones corrected-Z test
- `mle.homo()` — Ekbohm maximum-likelihood test under homoscedasticity
- `mle.hetero()` — Lin and Stivers maximum-likelihood test under heteroscedasticity

Each function accepts a two-column matrix or data frame with missing values marking unmatched observations. If the input is not partially matched, the implementation delegates to an appropriate standard one- or two-sample procedure.

## Installation

Install the packaged source archive:

```r
install.packages("PMLi_1.0.tar.gz", repos = NULL, type = "source")
```

Or install the repository with `remotes`:

```r
install.packages("remotes")
remotes::install_github("Garylikai/PMLi-1.0-R-Package")
```

The package declares R 2.10 or later and imports base `stats`. Building the vignette also requires `knitr` and `rmarkdown`.

## Example

```r
library(PMLi)

# The package includes a partially matched sample dataset named pm.
weighted.z(pm)
modified.t(pm, alternative = "less", conf.level = 0.99)
corrected.z(pm)
mle.homo(pm)
mle.hetero(pm)
```

Returned objects use the standard `htest` structure, including a test statistic, p-value, confidence interval, estimate, and method description.

## Repository structure

- `R/` — package source
- `man/` — generated function and dataset documentation
- `data/` and `data-raw/` — sample data and its construction script
- `vignettes/vignette.Rmd` — methodological overview and examples

## Scope

This package was created as a statistical-computing project and is not on CRAN. Validate method assumptions and results independently before using it for consequential research, clinical, or regulatory decisions.

## License

MIT. See `LICENSE`.
