
<!-- README.md is generated from README.Rmd. Please edit that file -->

# realm

realm (relational alignment models) is an R package for the statistical
analysis of relational alignment models. Its code is based
on the [MoNAn](https://github.com/stocnet/MoNAn) package (version 2.1.2).

This package is in an early stage of development.

Development of MoNAn version 2, on which realm is based, used Claude
Code for the majority of the C++ implementation of the
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
