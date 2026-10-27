
<!-- README.md is generated from README.Rmd. Please edit that file -->

# realm

realm (Relational Analysis of Mobility) is an R package for the statistical
analysis of weighted mobility networks and distribution networks. It is based
on the [MoNAn](https://github.com/stocnet/MoNAn) package (version 2.1.1.9000).

This package is in an early stage of development.

Programming MoNAn version 2, on which realm is based, extensively used Claude
Code, which produced the code for the majority of the C++ implementation of the
MCMC. The C++ implementation is checked in all steps and effects that it exactly
reproduces the human-coded R-based implementation. Additionally, Claude Code
was used for optimization; similarly, output from the machine was checked
against previous human-coded implementation.

# Installation

realm contains C++ code, so installing it from GitHub requires a C++ compiler
(on Windows Rtools, on macOS the Xcode command line tools).

```r
# install.packages("remotes")
remotes::install_github("perblock/realm")
```

An example of how to use the package is in `realmExampleScript.R`.
